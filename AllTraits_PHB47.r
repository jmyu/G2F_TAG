library(dplyr)
library(tidyr)
library(rrBLUP)
source("./JGRA.function.r")

EnvInd_table <- read.table("./EnvInd_table.txt",header = T)
G <- data.frame(fread('./McNellie_Final_Genotype_g2f_geno_hybrid_imp_v3_numeric_reduced_snpEff.csv'))
rownames(G) <- G$X.Marker.

tester <- 'PHB47'

traits <- c("DTA","DTS","PH","EH","GM","GY")

result.table <- data.frame(Tester = character(),
                           Trait = character(),
                           TGUE_across = numeric(),
                           UGTE_across = numeric(),
                           UGUE_across = numeric(),
                           TGUE_within = numeric(),
                           UGTE_within = numeric(),
                           UGUE_within = numeric())

for (trait in traits) {
  dt_input <- pre_jgra(tester,trait)
  result1to2 <- JGRA.1to2(pheno = dt_input$pheno,
                          envir = dt_input$env_ind,
                          geno = dt_input$geno,
                          fold = 10,reshuffle = 10)
  write.table(result1to2$pred.y,
              paste("./result/",tester,"_",trait,"_1to2_out4figure",sep = ''),
              row.names = F,quote = F)
  r_1to2_within <- mean(result1to2$pred.accuracy$r)
  r_1to2_overall <- mean(result1to2$pred.accuracy$r_overall)
  
  result1to3 <- JGRA.1to3(pheno = dt_input$pheno,
                          envir = dt_input$env_ind,
                          geno = dt_input$geno,
                          fold = 10,reshuffle = 10)
  write.table(result1to3[[1]],
              paste("./result/",tester,"_",trait,"_1to3_out4figure",sep = ''),
              row.names = F,quote = F)
  r_1to3_within <- mean(result1to3[[2]])
  r_1to3_overall <- mean(result1to3[[3]])
  
  result1to4 <- JGRA.1to4(pheno = dt_input$pheno,
                          envir = dt_input$env_ind,
                          geno = dt_input$geno,
                          fold = 10,reshuffle = 10)
  write.table(result1to4[[1]],
              paste("./result/",tester,"_",trait,"_1to4_out4figure",sep = ''),
              row.names = F,quote = F)
  r_1to4_within <- mean(result1to4[[2]])
  r_1to4_overall <- mean(result1to4[[3]])
  
  result <- data.frame(Tester = tester,
                       Trait = trait,
                       TGUE_across = r_1to2_overall,
                       UGTE_across = r_1to3_overall,
                       UGUE_across = r_1to4_overall,
                       TGUE_within = r_1to2_within,
                       UGTE_within = r_1to2_within,
                       UGUE_within = r_1to2_within)
  result.table <- rbind(result.table,result)
}
write.table(result.table,
            paste("./result/table/PredictionAbility",
                  tester,sep = '_'),
            row.names = F,quote = F)
