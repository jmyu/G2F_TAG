library(data.table)
pre_jgra <- function(tester, trait){
  env_para <- subset(EnvInd_table,
                     Tester == tester & Trait == trait)
  
  pheno <- read.table(paste("./Pheno_LbE/",tester,"_",trait,"_LbE_table",sep = ''),
                      stringsAsFactors = F, header = T)
  
  env <- read.table(paste("./Env_MeanPara/",tester,"_",
                          env_para$Env_para,"_",
                          trait,"_envMeanPara_",
                          env_para$Win_s,"_",env_para$Win_e,
                          sep = ''),
                    stringsAsFactors = F, header = T, sep='\t') %>%
    filter(env_code %in% colnames(pheno))
  
  env.list <- intersect(colnames(pheno)[-1],env$env_code)
  
  env_ind <- env[,c('env_code',env_para$Env_para)]
  env_ind <- env_ind[match(env.list,env_ind$env_code),]
  pheno <- pheno[,c('line_code',env.list)]
  env_ind$env_code == colnames(pheno)[-1]
  colnames(env_ind)[2] <- 'env_ind'
  
  geno <- as.matrix(G[match(pheno$line_code, G$X.Marker.),-1]*2 - 1)
  sum(rownames(geno) != pheno$line_code)
  
  return(list(pheno = pheno,
              geno = geno,
              env_ind = env_ind))
}


JGRA.1to2 <- function(pheno,geno,envir,
                 fold = 10,reshuffle = 50)
{
  envir.name <- colnames(pheno)[-1]

  n.line <- nrow(pheno)
  n.envir <- ncol(pheno)-1


    pheno.y <- pheno %>% gather(-line_code,key = 'env_code',value = 'Yobs', na.rm = T)
    

    lines <- unique(pheno.y$line_code)
    envirs <- unique(pheno.y$env_code)
    
    pred.ab <- data.frame()
    pred.y <- data.frame()
    
    for(k in 1:length(envirs))
    {
      envir.training <- subset(envir,env_code != envirs[k])
      for(j in 1:length(lines))
      {
        pheno.j <- subset(pheno.y,line_code == lines[j]) %>%
          filter(env_code != envirs[k]) %>%
          merge(envir,by = "env_code")
        
        fit <- lm(Yobs ~ env_ind, data = pheno.j)
        
        pred.b <- as.vector(round(fit$coefficient[2],4))
        pred.a <- as.vector(round(predict(fit, 
                                          data.frame(env_ind = mean(envir.training$env_ind))), 4)) ## adjusted by the population mean
        
        pred.ab.i <- data.frame(line_code = lines[j],
                                pred.a = pred.a,
                                pred.b = pred.b,
                                env_removed = envirs[k])
        pred.y.i <- data.frame(env_code = envirs[k],
                               line_code = lines[j],
                               pred = pred.b*(envir$env_ind[match(envirs[k],envir$env_code)] - 
                                                mean(envir$env_ind)) + pred.a)
        pred.ab <- rbind(pred.ab,pred.ab.i)
        pred.y <- rbind(pred.y,pred.y.i)
      }
    }
    
    pred.y <- merge(pheno.y,pred.y,by = c("line_code","env_code"),all.x = T)
    colnames(pred.y)[3:4] <- c('obs','pred')
    
    pred.accuracy.y <- pred.y %>%
      group_by(env_code) %>%
      summarize(r = cor(pred,obs,use = 'complete.obs')) %>% ungroup() %>%
      mutate(r_overall = cor(pred.y$pred,pred.y$obs,use = 'complete.obs'))
    
    
    return(list(pred.ab = pred.ab,
                pred.y = pred.y,
                pred.accuracy.y = pred.accuracy.y))
    
    
}

JGRA.1to3 <- function(pheno,geno,envir,fold=10,reshuffle=2){
  # colnames(pheno)=as.character(unlist(pheno[1,]));
  envir.name=colnames(pheno)[-1];
  
  # pheno=pheno[-1,];
  m=pheno[,-1];
  m=as.numeric(as.character(unlist(m)));m <- matrix(data=m, ncol=dim(pheno)[2]-1, nrow=dim(pheno)[1]);
  colnames(m)=colnames(pheno)[-1];
  pheno_=data.frame(line_code=pheno$line_code,m);
  colnames(pheno_)=c("line_code",envir.name);
  pheno=pheno_;
  
  n.line=dim(pheno)[1];
  n.envir=dim(pheno)[2]-1;
  
  intercept=numeric();
  slope=numeric();
  for(j in 1:n.line)
  {
    x1=envir[,2];
    y1=as.vector(t(pheno[j,-c(1)]));
    
    coe=lm(y~x,data=data.frame(x=x1,y=y1));
    inter=summary(coe)$coefficients[1,1]
    slop=summary(coe)$coefficients[2,1];
    intercept=c(intercept,inter);
    slope=c(slope,slop);
  }
  
  Marker=geno;
  
  intercept.hat=numeric();slope.hat=numeric();cor.within=matrix(999,reshuffle,n.envir);cor.all=numeric();
  for(i in 1:reshuffle)
  { 
    
    #      cross=sample(rep(1:fold,each=round(n.line/fold,0)),n.line);
    cross = sample(rep(1:fold,each=ceiling(n.line/fold)),n.line)
    yhat.whole.cross=numeric();yobs.whole.cross=numeric();
    for(f in 1:fold)
    {
      id.T=c(1:n.line)[cross!=f]; id.V=c(1:n.line)[cross==f];
      ##Intercept###
      y0=intercept; 
      ans<-mixed.solve(y=y0[id.T],Z=Marker[id.T,])
      e=as.matrix(ans$u)
      G.pred=Marker[id.V,]
      y_pred=as.matrix(G.pred) %*% e
      GEBV.inter=c(y_pred[,1])+c(ans$beta);
      ##Slope###
      y0=slope; 
      ans<-mixed.solve(y=y0[id.T],Z=Marker[id.T,])
      e=as.matrix(ans$u)
      G.pred=Marker[id.V,]
      y_pred=as.matrix(G.pred) %*% e
      GEBV.slope=c(y_pred[,1])+c(ans$beta);
      ###All the predicted slope and intercept
      yhat.envir=matrix(999,length(id.V),n.envir);yobs.envir=matrix(999,length(id.V),n.envir)
      for(j in 1:n.envir)
      {
        yhat=GEBV.inter+GEBV.slope*envir[j,2];
        yobs=pheno[id.V,j+1];
        yhat.envir[,j]=yhat;yobs.envir[,j]=yobs;
      }
      yhat.whole.cross=rbind(yhat.whole.cross,yhat.envir);
      yobs.whole.cross=rbind(yobs.whole.cross,yobs.envir);
    }
    for(j in 1:n.envir)
    {cor.within[i,j]=cor(yhat.whole.cross[,j],yobs.whole.cross[,j],use = "complete.obs");}
    
    cor.all=c(cor.all,cor(as.vector(yhat.whole.cross),as.vector(yobs.whole.cross),use = "complete.obs"));
  }
  #Correlation within environment 50 times
  r_within=apply(cor.within,2,mean);names(r_within)=colnames(pheno)[-1];
  r_within=data.frame(cor_within=r_within,envir=colnames(pheno)[-1]);
  #Correlation across environment 50 times
  r_across=mean(cor.all);
  #Observation and prediction last time
  outforfigure=data.frame(obs=as.vector(yobs.whole.cross),
                          pre=as.vector(yhat.whole.cross),
                          envir=rep(colnames(pheno)[-1],times=rep(n.line,n.envir)))
  colnames(cor.within)=colnames(pheno)[-1];
  r_within=cor.within;
  r_across=cor.all;
  return(list(outforfigure,r_within,r_across));
}   


JGRA.1to4 <- function(pheno,geno,envir,fold=10,reshuffle=2)
{
  
  # colnames(pheno)=as.character(unlist(pheno[1,]));
  envir.name=colnames(pheno)[-1];
  
  # pheno=pheno[-1,];
  m=pheno[,-1];
  m=as.numeric(as.character(unlist(m)));m <- matrix(data=m, ncol=dim(pheno)[2]-1, nrow=dim(pheno)[1]);
  colnames(m)=colnames(pheno)[-1];
  pheno_=data.frame(line_code=pheno$line_code,m);
  colnames(pheno_)=c("line_code",envir.name);
  pheno=pheno_;
  
  n.line=dim(pheno)[1];
  n.envir=dim(pheno)[2]-1;


    Marker=geno;
    
    cor.within=matrix(999,reshuffle,n.envir);cor.all=numeric();
    for(i in 1:reshuffle)
    {
      obs_matrix=matrix(999,n.line,n.envir);pre_matrix=matrix(999,n.line,n.envir);
      for(k in 1:n.envir)
      {
        intercept=numeric();
        slope=numeric();
        for(j in 1:n.line)
        {
          x1=envir[-k,2];
          y1=as.vector(t(pheno[j,-c(1,1+k)]));
          
          coe=lm(y~x,data=data.frame(x=x1,y=y1));
          inter=summary(coe)$coefficients[1,1]
          slop=summary(coe)$coefficients[2,1];
          intercept=c(intercept,inter);
          slope=c(slope,slop);
        }
        
        cross = sample(rep(1:fold,each=ceiling(n.line/fold)),n.line)
        yhat.whole=numeric();yobs.whole=numeric();
        
        for(f in 1:fold)
        {
          id.T=c(1:n.line)[cross!=f]; id.V=c(1:n.line)[cross==f];
          ##Intercept###
          y0=intercept; 
          ans<-mixed.solve(y=y0[id.T],Z=Marker[id.T,])
          e=as.matrix(ans$u)
          G.pred=Marker[id.V,]
          y_pred=as.matrix(G.pred) %*% e
          GEBV.inter=c(y_pred[,1])+c(ans$beta);
          ##Slope###
          y0=slope; 
          ans<-mixed.solve(y=y0[id.T],Z=Marker[id.T,])
          e=as.matrix(ans$u)
          G.pred=Marker[id.V,]
          y_pred=as.matrix(G.pred) %*% e
          GEBV.slope=c(y_pred[,1])+c(ans$beta);
          ###All the predicted slope and intercept
          yhat=GEBV.inter+GEBV.slope*envir[k,2];
          yobs=pheno[id.V,k+1];
          
          yhat.whole=c(yhat.whole,yhat);
          yobs.whole=c(yobs.whole,yobs);
        }
        cor.within[i,k]=cor(yhat.whole,yobs.whole,use = "complete.obs");
        obs_matrix[,k]=yobs.whole;
        pre_matrix[,k]=yhat.whole;
      }
      cor.shuffle=cor(as.vector(obs_matrix),as.vector(pre_matrix),use = "complete.obs")
      cor.all=c(cor.all,cor.shuffle);
    }
    
    yhat.whole.cross=pre_matrix;
    yobs.whole.cross=obs_matrix;
    #Correlation within environment 50 times
    r_within=apply(cor.within,2,mean);names(r_within)=colnames(pheno)[-1];
    #Correlation across environment 50 times
    r_across=mean(cor.all);
    #Observation and prediction last time
    outforfigure=data.frame(obs=as.vector(yobs.whole.cross),
                            pre=as.vector(yhat.whole.cross),
                            # col=rep(coloo,times=rep(n.line,n.envir)),
                            envir=rep(colnames(pheno)[-1],
                                      times=rep(n.line,n.envir)))
    colnames(cor.within)=colnames(pheno)[-1];
    r_within=cor.within;
    r_across=cor.all;
    
  
  return(list(outforfigure,r_within,r_across));
}


  