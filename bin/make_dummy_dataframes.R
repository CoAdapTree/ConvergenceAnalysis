rm(list = ls())

## Simulate data to analyse as pretend convergence data

## Number of orthogroups 30000

# A poisson dist. with mean 2.5 gives a number of orthogroups that matches well


make_species_df <- function(numOrtho, pois_lambda = 2.5){ 
  if (pois_lambda == 1){
    species_genes <- rep(1,numOrtho)
  }
  else {
    species_genes <- rpois(numOrtho, pois_lambda )
    
  }
  gene_names <- paste("gene",1:sum(species_genes), sep = "")

  orthogroups <- paste("ortho", rep( 1:numOrtho, times = species_genes), sep = "")

  gene_df = data.frame( gene = gene_names, orthogroup = orthogroups, Z = rnorm(sum(species_genes), mean = 0, sd = 2))

  return(gene_df)
}

nOrthos = 30000
species_1 <- make_species_df(nOrthos)
species_2 <- make_species_df(nOrthos)
species_3 <- make_species_df(nOrthos)
species_4 <- make_species_df(nOrthos)
species_5 <- make_species_df(nOrthos)

write.csv(species_1, "~/work/GEA/Convergence/dev_scripts/species_1.csv", row.names = F, quote = F)
write.csv(species_2, "~/work/GEA/Convergence/dev_scripts/species_2.csv", row.names = F, quote = F)
write.csv(species_3, "~/work/GEA/Convergence/dev_scripts/species_3.csv", row.names = F, quote = F)
write.csv(species_4, "~/work/GEA/Convergence/dev_scripts/species_4.csv", row.names = F, quote = F)
write.csv(species_5, "~/work/GEA/Convergence/dev_scripts/species_5.csv", row.names = F, quote = F)
