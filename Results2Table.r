library(dplyr)
library(readr)
results <- list.files("./result/table/", 
                 full.names = TRUE) %>% 
  lapply(read.table,header = T) %>% 
  bind_rows 

results[,3:8] <- round(results[,3:8],2)

result.across <- results %>%
  select(1:5)
result.within <- results %>%
  select(1,2,6,7,8)

write.csv(result.across,
            file = "./PredictiveAbility_across.csv",
          row.names = F,quote = F)
write.csv(result.within,
          file = "./PredictiveAbility_within.csv",
          row.names = F,quote = F)
