rm(list = ls())

x<-read.csv("~/work/GEA/Convergence/dev_scripts/random_WZA_set_summaries/random_WZA_set_summary.multi_species.csv")


par( mfrow = c(2,2))

hist( x$TruePositives, 
      xlab = "True Number of Convergent Genes",
      main = "")

hist( x$NumberOfHits, 
      xlab = "Number of Convergent Genes Identified",
      main = "")


hist( x$NumberOfGoodGuesses, 
      xlab = "Number of Times Correct Configurations" ,
      main = "")


hist( x$FalsePositives, 
      xlab = "Number of False Positives" ,
      main = "")


library(ggplot2)

multi <- read.csv("~/work/GEA/Convergence/dev_scripts/random_WZA_sets/random_set_10.multi_species.csv")

multi_strim <- multi[ !is.na( multi$converge_pVal ) , ]

multi_strim$converge_pVal_adjust <- p.adjust( multi_strim$converge_pVal, 
                                              method = "BH")

sum(multi_strim$converge_pVal < 0.05)

ggplot(data = multi_strim, 
       aes( x = orthogroup, 
            y = -log10(converge_pVal_adjust),
            col = as.factor(converge)))+
  geom_point()+
  geom_hline(aes( yintercept = -log10(0.05)))


