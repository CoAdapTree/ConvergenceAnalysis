rm(list = ls ())

# in this script, we are going to apply the p_jth stats to the WZA simulations

library(ggplot2)
library(reshape2)
library(ggpubr)

formatWZAresult <- function(wza_df, reps, id){
  output = list()
  for (i in seq(length(reps))){
    set_rep <- wza_df[wza_df$rep == reps[i],]
    temp_df = data.frame( emp_Z = set_rep$kendall_Z ,
                          LA = set_rep$LA,
                          gene= set_rep$gene,
                          position = set_rep$position + (i-1)*10e6)
    names(temp_df) <- c("Z",
                        paste("LA_",reps[1],sep = ""),
                        paste("gene_",reps[1],sep = ""),
                        paste("position_",reps[1],sep = ""))
    output[[i]] =  temp_df 
  }
  boundDF <- do.call(rbind, output)
  
  boundDF$Z <-  1- ((rank( boundDF$Z ) -1)/ length(boundDF$Z))
  
  names(boundDF) <- c(paste("empZ_",id,sep = ""),
                      paste("LA_",id,sep = ""),
                      paste("gene_",id,sep = ""),
                      paste("position_",id,sep = ""))
  return(boundDF)
}


p_jth <- function(plist,j) {
  # get the list of the p-value vector
  k=length(plist)
  # return an NA if the j parameter is bigger than the length of the p vector
  if(j>k) {return(NA)}
  # sort the p-values
  pjth=sort(plist)[j]
  # scale the probabilities 
  prob = if(j<k) (pjth /sort(plist)[j+1]) else max(plist)
  # perform the test
  dbinom(x = j, size = j, prob=prob )
}

p_value_adjustment <- function(p_vec){
  p_adj = 1 - (1 - min(p_vec))^length(p_vec)
  return(p_adj)
}



LA_indicator_function_Mult <- function(LA_vec){
  
  if (length(LA_vec) !=5){
    return(NA)
  }
  if (-99 %in% LA_vec){
    return("Linked")
  }
  else if (sum(LA_vec<=0.001) == length(LA_vec)){
    return("Neutral")
  }
  else if ( sum(LA_vec>0.001) == 1){
    return("One Hit")
  }
  else if ( sum(LA_vec>0.001) == 2){
    return("Two Hits")
  }
  else if ( sum(LA_vec>0.001) == 3){
    return("Three Hits")
  }
  else if ( sum(LA_vec>0.001) == 4){
    return("Four Hits")
  }
  else if ( sum(LA_vec>0.001) == 5){
    return("Five Hits")
  }
  
}

LA_indicator_function_Quad <- function(LA_vec){
  
  if (length(LA_vec) !=5){
    return(NA)
  }
  if (-99 %in% LA_vec){
    return("Linked")
  }
  else if (sum(LA_vec<=0.005) == length(LA_vec)){
    return("Neutral")
  }
  else if ( sum(LA_vec>0.005) == 1){
    return("Single")
  }
  else if ( sum(LA_vec>0.005) >= 2){
    return("Convergent")
  }
}



FisherCombProb <- function( pVec ){
  pSort=sort(pVec)[2:length(pVec)]
  testStat = -2*sum(log(pSort))
  pVal = 1 - pchisq( testStat, 2*length(pSort) )
  return( pVal )
}


plot_pMax <- function(wzaDF1, wzaDF2, wzaDF3, shuffle = F){
  ## There are 20 simulation replicates
  
  datasets <- sample( seq( 1,20 )) 
  shuffle = T
  wzaDF1<-BC
  wzaDF2<- cline
  wzaDF3<- trunc
  set_1 <- formatWZAresult(wzaDF1, c(1:20)[c(11:20,1:10)] , 1)
  set_2 <- formatWZAresult(wzaDF1, c(1:20)[c(1:20)], 2)
  set_3 <- formatWZAresult(wzaDF3, c(1:20)[c(11:20,1:10)] , 3)
  set_4 <- formatWZAresult(wzaDF2, (1:20)[c(1:20)], 4)
  set_5 <- formatWZAresult(wzaDF2, c(1:20)[c(11:20,1:10)] , 5)

  #set_1 <- formatWZAresult(wzaDF1, datasets[sample( seq( 1,20 )) ] , 1)
  if (shuffle == T){
    set.seed(1)
    
      sampler = sample(nrow(set_1), replace = FALSE)
      set_1$empZ_1 <- set_1$empZ_1[sampler]
      set_1$LA_1 <- set_1$LA_1[sampler]
    }
  
  set_2 <- formatWZAresult(wzaDF1, datasets[sample( seq( 1,20 )) ], 2)
  if (shuffle == T){
    set.seed(13)
  
    sampler = sample(nrow(set_2), replace = FALSE)
    set_2$empZ_2 <- set_2$empZ_2[sampler]
    set_2$LA_2 <- set_2$LA_2[sampler]
  }
  
  set_3 <- formatWZAresult(wzaDF2, datasets[sample( seq( 1,20 )) ], 3)
  if (shuffle == T){
    set.seed(1151)
    
    sampler = sample(nrow(set_3), replace = FALSE)
    set_3$empZ_3 <- set_3$empZ_3[sampler]
    set_3$LA_3 <- set_3$LA_3[sampler]
  }
  
  set_4 <- formatWZAresult(wzaDF2, datasets[sample( seq( 1,20 )) ], 4)
  if (shuffle == T){
    set.seed(222)
    
    sampler = sample(nrow(set_4), replace = FALSE)
    set_4$empZ_4 <- set_4$empZ_4[sampler]
    set_4$LA_4 <- set_4$LA_4[sampler]
  }
  
  set_5 <- formatWZAresult(wzaDF3, datasets[sample( seq( 1,20 )) ], 5)
  if (shuffle == T){
    set.seed(123)
    sampler = sample(nrow(set_5), replace = FALSE )
    set_5$empZ_5 <- set_5$empZ_5[sampler]
    set_5$LA_5 <- set_5$LA_5[sampler]
  }
  
    test_set <- cbind( set_1,
                     set_2,
                     set_3,
                     set_4,
                     set_5)
  
  pDF <- test_set[c(1,5,9,13,17)]

  melted_test_set_p1 <- melt(test_set, 
                             id = c("LA_1", "gene_1","position_1", "LA_2","gene_2","position_2",   
                                    "LA_3","gene_3","position_3","LA_4","gene_4","position_4","LA_5","gene_5","position_5"),
                             value.name	= "empZ" , variable.name = "replicate")
  melted_test_set <- melt(melted_test_set_p1, 
                          id = c("gene_1","position_1", "gene_2","position_2",   
                                 "gene_3","position_3","gene_4","position_4","gene_5","position_5", "empZ", "replicate"),
                          value.name	= "LA" , variable.name = "replicate_rep")

  test_set$p_1th = apply(pDF, 1, function(x) p_jth( x[!is.na(x)], 1) )
  test_set$p_2th = apply(pDF, 1, function(x) p_jth( x[!is.na(x)], 2) )
  test_set$p_3th = apply(pDF, 1, function(x) p_jth( x[!is.na(x)], 3) )
  test_set$p_4th = apply(pDF, 1, function(x) p_jth( x[!is.na(x)], 4) )
  test_set$p_5th = apply(pDF, 1, function(x) p_jth( x[!is.na(x)], 5) )
  
  
  test_set$pMax <- apply(test_set[,c(21,22,23,24,25)], 1, p_value_adjustment)
  
  test_set$Fisher = apply(pDF, 1, FisherCombProb )
  
  # now just plop the results in the list
  
  test_set$LA_indicator_mult <- apply(test_set[,c(2,6,10,14, 18)], 1, LA_indicator_function_Mult)
  test_set$LA_indicator_quad <- apply(test_set[,c(2,6,10,14, 18)], 1, LA_indicator_function_Quad)
  
  test_set$LA_indicator_mult <- factor( test_set$LA_indicator_mult, 
                                        levels = c("Neutral", "Linked", "One Hit", "Two Hits", "Three Hits", "Four Hits", "Five Hits"))
  
  test_set$LA_indicator_quad <- factor( test_set$LA_indicator_quad, 
                                        levels = c("Neutral", "Linked", "Single", "Convergent"))
  test_set$FDR_p <- -log10(p.adjust( test_set$pMax, method = "fdr"))

  test_set$FDR_p_fisher <- -log10(p.adjust( test_set$Fisher, method = "fdr"))
  
  test_set$position <- test_set[,4]
  
  test_set$label <- "pMax"
  
  berlin <- ggplot(data = test_set, aes( x = position/1e6, y = -log10(pMax) , col = LA_indicator_quad))+
    geom_point()+
    geom_hline(aes(yintercept = -log10(0.05)), lty = 2, col = "grey")+
    scale_color_brewer(palette = "Dark2")+
    facet_grid(label~.)+
    theme_bw()+
    ggtitle("pMax Test")+
    scale_x_continuous("Position (Mbp)")+
    scale_y_continuous("-log10( p-values )",limits = c(0,7))+
    theme(
      legend.title = element_blank(),
      legend.position = "bottom"
    )
  
  test_set$marker <- 1
  bar <- ggplot(data = test_set[-log10(test_set$pMax)> -log10(0.05), ], aes( x =  marker, fill = LA_indicator_mult))+
    geom_bar(position = "stack")+
        scale_y_continuous("Count")+
    #    scale_y_continuous("Count", limits = c(0,40), breaks = 1:10*4)+
    scale_fill_brewer(palette = "Paired")+
    theme_bw()+
    theme(
      legend.title = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  melted_test_set$replicate <- factor( melted_test_set$replicate,
                                       levels = c("empZ_1","empZ_2","empZ_3","empZ_4","empZ_5"),
                                       labels = c("Species 1","Species 2","Species 3","Species 4","Species 5"))
  manhattan <- ggplot( data = melted_test_set, aes(x= position_1/1e6, y = -log10(empZ), col = LA>0.001, group = replicate))+
    geom_point()+
    facet_grid( replicate~.)+
    scale_y_continuous("-log10( empirical p-value )", limits = c(0,6) )+
    scale_color_manual("Local Adaptation Gene", values = c("black","red"))+
    scale_x_continuous("Position (Mbp)")+
    theme_bw()

  manhattan2 <- ggplot( data = melted_test_set, aes(x= position_1/1e6, y = LA, col = LA>0.001, group = replicate_rep))+
    geom_point()+
    facet_grid( replicate_rep~.)+
#    scale_y_continuous("-log10( empirical p-value )", limits = c(0,6) )+
    scale_color_manual("Local Adaptation Gene", values = c("black","red"))+
    scale_x_continuous("Position (Mbp)", limits = c(0,2))+
    theme_bw()
  print( ggarrange( ggarrange(berlin, bar, widths = c(5,1), ncol = 2), ggarrange( manhattan+guides(colour=F), get_legend(manhattan), widths = c(5,1),ncol = 2), ncol = 1, nrow = 2, heights = c(2, 5) ), align = "v" )
  
  
}

# 
trunc <- read.csv("~/UBC/GEA/Convergence/WZA_results_for_testing/analysisFiles_s0.003/trunc_sampled.WZA.csv")

## What is the average proporiton of locally adapted genes? 
mean(tapply(trunc$LA, trunc$rep, function(x) sum(x > 0.005))/12)

# 
# pdf("trunc_convergence_tests.pdf",width = 12, height = 6)
# for (i in 1:10){
#   plot_pMax(trunc)
# }
# dev.off()
# 
# 
cline <- read.csv("~/UBC/GEA/Convergence/WZA_results_for_testing/analysisFiles_s0.003/cline_sampled.WZA.csv")

mean(tapply(cline$LA, cline$rep, function(x) sum(x > 0.005))/12)
# 
# pdf("cline_convergence_tests.pdf",width = 12, height = 6)
# for (i in 1:10){
#   plot_pMax(cline)
# }
# dev.off()
# 
# 
BC <- read.csv("~/UBC/GEA/Convergence/WZA_results_for_testing/analysisFiles_s0.003/BC_Map_sampled.WZA.csv")
#
mean(tapply(BC$LA, BC$rep, function(x) sum(x > 0.005))/12)

# pdf("BC_Map_convergence_tests.pdf",width = 12, height = 6)
# for (i in 1:10){
#   plot_pMax(BC)
# }
# dev.off()
# 
# for (i in 1:2){
# png(paste("~/UBC/GEA/Convergence/WZA_results_for_testing/",i,"_s0.003_convergence_tests_pMax_NoFDR.png", sep = ""),width = 12, height = 10, units = "in", res= 300)
#   plot_pMax(BC, cline, trunc)
#   dev.off()
# }
# 
# for (i in 1:2){
#   png(paste("~/UBC/GEA/Convergence/WZA_results_for_testing/shuffle_",i,"_s0.003_convergence_tests_pMax_NoFDR.png", sep = ""),width = 12, height = 10, units = "in", res= 300)
#   plot_pMax(BC, cline, trunc, shuffle = T)
#   dev.off()
# }


