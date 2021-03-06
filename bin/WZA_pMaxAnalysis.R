rm(list = ls () )


p_value_adjustment <- function(p_vec){
  
  p_adj = 1 - (1 - min(p_vec))^length(p_vec)
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


wl <- read.csv("~/UBC/GEA/Convergence/dev_data/WL_all_wza_tc_OG.txt", sep = "\t")
wl <- wl[wl$env == "AHM",]

jp <- read.csv("~/UBC/GEA/Convergence/dev_data/JP_all_wza_tc_OG.txt",sep = "\t")
jp <- jp[jp$env == "AHM",]

dfi <- read.csv("~/UBC/GEA/Convergence/dev_data/DFI_all_wza_tc_OG.txt", sep = "\t")
dfi <- dfi[dfi$env == "AHM",]

dfc <- read.csv("~/UBC/GEA/Convergence/dev_data/DFC_all_wza_tc_OG.txt", sep = "\t")
dfc <- dfc[dfc$env == "AHM",]

wl_cut <- cut_df_for_pMax( wl, "wl")
jp_cut <- cut_df_for_pMax( jp, "jp")
dfi_cut <- cut_df_for_pMax( dfi, "dfi")
dfc_cut <- cut_df_for_pMax( dfc, "dfc")

n_taxa = 3

SpeciesComparison <- Reduce(function(x, y) merge(x, y, by = "orthogroup", all = T), 
                            list(dfi_cut, jp_cut, wl_cut))

SpeciesComparison_pMax <- run_pMax(SpeciesComparison, n_taxa, n_taxa)

SpeciesComparison_pMax$pMax <-  apply(as.matrix(SpeciesComparison_pMax[,5:6]), 1,  p_value_adjustment )

length(na.omit( SpeciesComparison_pMax$pMax ) )
head(SpeciesComparison_pMax[order( SpeciesComparison_pMax$pMax),] )

hist(SpeciesComparison_pMax[!is.na(SpeciesComparison_pMax$pMax),]$pMax, main = "pMax results for AHM - between DF_I, DF_C, JP and WL", xlab = "pMax", breaks = 10)
abline(  h = length(na.omit(SpeciesComparison_pMax$pMax))/10 )

SpeciesComparison_pMax_results <- SpeciesComparison_pMax[!is.na(SpeciesComparison_pMax$pMax),]

# Here's the list of pMax hits
SpeciesComparison_pMax_results[SpeciesComparison_pMax_results$pMax < 0.05,]$orthogroup

