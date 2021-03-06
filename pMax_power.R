rm (list = ls())

## Let's test the power of the pMax test

## Here's a function for simulating correlated data

simulateCorrelation <- function(r, n = 100){
# r is the population correlation coefficient 
# x1 is the sample data
    x1 = rnorm(n)
# x2 is the data that will be transformed so as to be correlated with x1
    x2 = rnorm(n)

    Y = r*x1+sqrt(1-r*r)*x2     
    return( cor.test(Y, x1)$p.value )
}

## Now make the function run in a vectorised fashion 
simulateCorrelationVectorized <- Vectorize(simulateCorrelation)

# Wrap the function up so that it takes a vector of correlations for each "species"
#and returns the data in a nice way
simulateMultipleGEAs <- function(corVector, times){
  t(replicate(times,  simulateCorrelationVectorized(corVector, n = 20) )  )
}

## Note that the 'n = 40' in the above specifies that a sample of 40 is being taken.
#This is a relatively small number to make things nice and noisy

## Here's the familiar p_jth function

p_jth <- function(plist,j) {
  k=length(plist)
  if(j>k) {return(NA)}
  pjth=sort(plist)[j]
  prob = if(j<k) (pjth /sort(plist)[j+1]) else max(plist)
  dbinom(x = j, size = j, prob=prob )
}

## ...and the p-value adjustment

p_value_adjustment <- function(p_vec){
  p_adj = 1 - (1 - min(p_vec))^length(p_vec)
  return(p_adj)
}

## Here's a wrapper function that runs the pMax on the 2, 3, 4 and 5 cases

p_jth_wrapper <- function( correlationSet ){
  p2th <- apply(correlationSet, 1, function(x) p_jth(x,2))
  p3th <- apply(correlationSet, 1, function(x) p_jth(x,3))
  p4th <- apply(correlationSet, 1, function(x) p_jth(x,4))
  p5th <- apply(correlationSet, 1, function(x) p_jth(x,5))

  return( data.frame(p2 = p2th, p3 = p3th, p4 = p4th, p5 = p5th) )
}

## Let's start with a case with no correlation (100 replicates)
noCorrelationSet <- simulateMultipleGEAs( c(0.0, 0.0, 0.0, 0.0, 0.0), 10000)

noCor_result  <- p_jth_wrapper( noCorrelationSet )
  
## Hey look, no correlation among values under the p2, p3, p4 , p5 cases under 
#the null - that's what we want to see
cor(noCor_result)

## Let's now look at the top hits
noCor_combined <- apply(as.matrix(noCor_result), 1, p_value_adjustment)

## A nice unifrom distribution - that's good
hist(noCor_combined)




## Now let's look at what happens when we have a single true positive

singleHitSet <- simulateMultipleGEAs( c(0.3, 0.3, 0.3, 0.0, 0.0), 10000)

singleHit_result  <- p_jth_wrapper( singleHitSet )
## Hey look, no correlation among values under the p2, p3, p4 , p5 cases under 
#the null - that's what we want to see
cor(singleHit_result)

## Let's now look at the top hits
singleHit_combined <- apply(as.matrix(singleHit_result), 1, p_value_adjustment)

head(singleHit_result)

## A nice unifrom distribution - that's good
hist(singleHit_combined)

sum(singleHit_combined<0.05)/length(singleHit_combined)

nSpecies=4


calculatePower <- function( nSpecies, alpha ){
  outLines = list()
  correlationsToTest <-  c(0.05,1:9/10)
  
  for (i in seq(length(correlationsToTest) )){
    print(c(rep(correlationsToTest[i],nSpecies),rep(0,5 - nSpecies)) )
    HitSet <- simulateMultipleGEAs( c(rep(correlationsToTest[i],nSpecies),rep(0,5 - nSpecies)), 1000)
    Hit_result  <- p_jth_wrapper( HitSet )
    Hit_combined <- apply(as.matrix(Hit_result), 1, p_value_adjustment)
    
    Hit_minP <- apply(as.matrix(Hit_result), 1, function(x) c(2:5)[which.min(x)])

    outLines[[i]] <- c(nSpecies, 
                       correlationsToTest[i], 
                       sum(Hit_combined<alpha)/length(Hit_combined),
                       sum(Hit_minP == nSpecies)/length(Hit_combined)
    )
  }
  
  temp_df <-data.frame( do.call(rbind, outLines) )
  names( temp_df ) <- c("nSpecies", "Cor", "pMaxPower","nSpeciesPower")
  return(temp_df)
}


AllSpeciesComparison_alpha05 <- rbind( calculatePower(1, 0.05),
                               calculatePower(2, 0.05),
                               calculatePower(3, 0.05),
                               calculatePower(4, 0.05),
                               calculatePower(5, 0.05) )
AllSpeciesComparison_alpha05$alpha = 0.05

AllSpeciesComparison_alpha01 <- rbind( calculatePower(1, 0.01),
                                       calculatePower(2, 0.01),
                                       calculatePower(3, 0.01),
                                       calculatePower(4, 0.01),
                                       calculatePower(5, 0.01) )
AllSpeciesComparison_alpha01$alpha = 0.01


AllSpeciesComparison_alpha001 <- rbind( calculatePower(1, 0.001),
                                       calculatePower(2, 0.001),
                                       calculatePower(3, 0.001),
                                       calculatePower(4, 0.001),
                                       calculatePower(5, 0.001) )
AllSpeciesComparison_alpha001$alpha = 0.001

AllSpeciesComparison <- rbind( AllSpeciesComparison_alpha05, AllSpeciesComparison_alpha01, AllSpeciesComparison_alpha001)

AllSpeciesComparison$alpha <- factor(AllSpeciesComparison$alpha,
                                     levels = c(0.05, 0.01, 0.001),
                                     labels = c(expression(alpha*" = 0.05"),
                                                expression(alpha*" = 0.01"),
                                                expression(alpha*" = 0.001")))
library ( ggplot2 ) 


ggplot( data = AllSpeciesComparison, aes( x = Cor, y = pMaxPower, col = as.factor(nSpecies)))+
  geom_line(lwd = 1)+
  scale_x_continuous("Population Correlation Coefficient for True Positives", limits = c(0,1))+
  scale_y_continuous(expression("Power to to detect gene"), limits = c(0,1))+
  scale_color_brewer("Number of \nspecies\nwith a true\npositive", palette = "Dark2")+
  facet_grid(~alpha, labeller = label_parsed)+
  theme_bw()
## The plateau at 0.05 is the false positives! The plataeu at ~0.2 for a singleton is something to be discussed


ggplot( data = AllSpeciesComparison[AllSpeciesComparison$nSpecies!=1,], aes( x = Cor, y = nSpeciesPower, col = as.factor(nSpecies)))+
  geom_line(lwd = 1)+
  scale_x_continuous("Population Correlation Coefficient for True Positives", limits = c(0,1))+
  scale_y_continuous(expression("Proportion of cases where the min p-value is for the right number of species"), limits = c(0,1))+
  scale_color_brewer("Number of \nspecies\nwith a true\npositive", palette = "Dark2")+
  geom_hline(aes(yintercept = 0.25), lty = 2)+
  facet_grid(~alpha, labeller = label_parsed)+
  theme_bw()
## The plataeu at 0.25 is the false positives - we will be right 1/4 of the time by chance



## Let's now make a plot of the p-value distributions across this range, so we can see wha kind of signal isrequired to get these results


p_value_Set <- simulateMultipleGEAs( c(0.05, 0.1 ,0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), 1000)
p_value_DF <- data.frame(p_value_Set)
names(p_value_DF) <- c("rho005", "rho01" ,"rho02", "rho03", "rho04", "rho05", "rho06", "rho07", "rho08", "rho09")
library(reshape2)
p_value_melted <- melt(p_value_DF)
p_value_melted$variable <- factor(p_value_melted$variable, 
                                  levels = c("rho005", "rho01" ,"rho02", "rho03", "rho04", "rho05", "rho06", "rho07", "rho08", "rho09"),
                                  labels = c(expression(italic("r")*" = 0.05"),
                                             expression(italic("r")*" = 0.1"),
                                             expression(italic("r")*" = 0.2"),
                                             expression(italic("r")*" = 0.3"),
                                             expression(italic("r")*" = 0.4"),
                                             expression(italic("r")*" = 0.5"),
                                             expression(italic("r")*" = 0.6"),
                                             expression(italic("r")*" = 0.7"),
                                             expression(italic("r")*" = 0.8"),
                                             expression(italic("r")*" = 0.9")))

ggplot(data = p_value_melted, aes(x = value))+
  geom_histogram(alpha= 0.75, bins = 50)+
  facet_wrap(~variable, scales = "free", ncol = 5, nrow= 2, labeller = label_parsed)+
  scale_fill_brewer("Population \ncorrelation\n coefficient", palette = "Dark2")+
  xlab(expression(italic("p")*"-value"))+
  ylab("Count")+
  guides(color = F, fill = F)+
  theme_bw()

