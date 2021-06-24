rm(list = ls () )


p_value_adjustment <- function(p_vec){
  if (length(na.omit(p_vec)) == 1){
    return(na.omit(p_vec)[1])
  }
  else if (length(na.omit(p_vec)) == 0){
    return(NA)
  }
  else{
    p_adj = 1 - (1 - min(na.omit(p_vec)))^length(na.omit(p_vec))
  }
  return(p_adj)
}

get_orthogroup_p_values <- function( species_DF , speciesLabel){
  gene_Score <- tapply(species_DF$Z_p, species_DF$OG, p_value_adjustment)
  gene_Score_Names <- names(gene_Score)
  orthoDF <- data.frame(orthogroup = gene_Score_Names, temp = as.numeric(gene_Score)  )
  names(orthoDF) <- c("orthogroup",speciesLabel)
  return(orthoDF)
}

cut_df_for_pMax <- function(focal_df, name){
  focal_df$Z_p <- 1 - rank(focal_df$Z )/nrow(focal_df) ## This converts Z to empirical ps. 
##### NOTE. If you want to use some variable other than the one named Z, you need to replace
  ### "focal_df$Z" in the above line with whatever your replacement it
  df_cut <- focal_df[,names(focal_df)%in%c("OG", "Z_p")] ## grab out the orthogroup and Z_p columns 
  df_ortho <- get_orthogroup_p_values( df_cut , name ) 
  return( df_ortho )
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



## Here's a wrapper function to run pMax on the simulated ortholog datasets
run_pMax <- function(species_DF, n_species, max_n){
  # Make a matrix of the p-values from the species dataframe
  species_M <- as.matrix(species_DF[,2:(n_species+1)])
  
  # create and empty list to put results in
  pMax_results = list()
  # loop over the different convergence configurations
  for (n in 2:max_n){
    # This ugly line applies the p_jth function to each line of the matrix using the current configuration 
    pMax_results_temp = apply(species_M, 1, function(x) p_jth( x[!is.na(x)], n) )
    # now just plop the results in the list
    pMax_results[[n-1]] = pMax_results_temp
  }
  # Convert the results into a nice data frame
  pMax_DF = as.data.frame( do.call(cbind, pMax_results) )
  # Give the columns infomative names
  names(pMax_DF) = paste("pMax_",2:max_n, sep = '')
  # done
  return(cbind( species_DF, pMax_DF ))
}


wl_raw <- read.csv("~/UBC/GEA/Convergence/dev_data/WL_all_wza_tc_OG.txt", sep = "\t")
wl <- wl_raw[wl_raw$env == "MSP",]

jp_raw <- read.csv("~/UBC/GEA/Convergence/dev_data/JP_all_wza_tc_OG.txt",sep = "\t")
jp <- jp_raw[jp_raw$env == "MSP",]

dfi_raw <- read.csv("~/UBC/GEA/Convergence/dev_data/DFI_all_wza_tc_OG.txt", sep = "\t")
dfi <- dfi_raw[dfi_raw$env == "MSP",]

dfc_raw <- read.csv("~/UBC/GEA/Convergence/dev_data/DFC_all_wza_tc_OG.txt", sep = "\t")
dfc <- dfc_raw[dfc_raw$env == "MSP",]

wl_cut <- cut_df_for_pMax( wl, "wl")
jp_cut <- cut_df_for_pMax( jp, "jp")
dfi_cut <- cut_df_for_pMax( dfi, "dfi")
dfc_cut <- cut_df_for_pMax( dfc, "dfc")

n_taxa = 3

hist(jp_raw[jp_raw$env == "MSP",]$Z, breaks = 100)
hist(dfi_raw[dfi_raw$env == "MSP",]$Z, breaks = 100)

cline <- read.csv("~/UBC/GEA/Convergence/WZA_results_for_testing/analysisFiles_Vs192/cline_sampled.WZA.csv")
hist(cline[cline$rep == 1,]$empR_Z, breaks = 100)

SpeciesComparison <- Reduce(function(x, y) merge(x, y, by = "orthogroup", all = T), 
                            list(dfi_cut, jp_cut, wl_cut))

SpeciesComparison_pMax <- run_pMax(SpeciesComparison, n_taxa, n_taxa)

SpeciesComparison_pMax$pMax <-  apply(as.matrix(SpeciesComparison_pMax[,5:6]), 1,  p_value_adjustment )
SpeciesComparison_pMax$FDR_p <- p.adjust( SpeciesComparison_pMax$pMax, method = "fdr")



hist(SpeciesComparison_pMax[!is.na(SpeciesComparison_pMax$pMax),]$pMax, main = "pMax results for AHM - between DF_I, JP and WL", xlab = "pMax", breaks = 100)
min(na.omit(SpeciesComparison_pMax$pMax))
abline(  h = length(na.omit(SpeciesComparison_pMax$pMax))/100 )

SpeciesComparison_pMax_results <- SpeciesComparison_pMax[!is.na(SpeciesComparison_pMax$pMax),]
SpeciesComparison_pMax$FDR_p <- p.adjust( SpeciesComparison_pMax$pMax, method = "fdr")


# Here's the list of pMax hits
SpeciesComparison_pMax_results[p.adjust( SpeciesComparison_pMax_results$pMax, method = "fdr") < 0.05,]$orthogroup

