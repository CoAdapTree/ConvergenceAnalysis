rm(list = ls())


library(reshape2)        
library(ggplot2)
par(mfrow = c(2,4))
power_vec = c(0,0,0.1,0.4,0.6,0.7,0.8,0.9)
nSpecies = 7
p_lik = 0.012
outputDF_0.2 <- as.data.frame( t( replicate( 100, summary(factor(rbinom(20000*p_lik, nSpecies, 0.2), levels = 0:nSpecies)) ) )  )
plot(melt(outputDF_0.2), xlab = "# of Species With Local Adaptation Gene", ylab = "Number of true positives", main = expression(italic(p[Lik])*" = 0.012; "* italic(p[Adap])*" = 0.2"))
outputDF_0.4 <- as.data.frame( t( replicate( 100, summary(factor(rbinom(20000*p_lik, nSpecies, 0.4), levels = 0:nSpecies)) ) )  )
plot(melt(outputDF_0.4), xlab = "# of Species With Local Adaptation Gene", ylab = "Number of true positives", main = expression(italic(p[Lik])*" = 0.012; "* italic(p[Adap])*" = 0.4"))
outputDF_0.6 <- as.data.frame( t( replicate( 100, summary(factor(rbinom(20000*p_lik, nSpecies, 0.6), levels = 0:nSpecies)) ) )  )
plot(melt(outputDF_0.6), xlab = "# of Species With Local Adaptation Gene", ylab = "Number of true positives", main = expression(italic(p[Lik])*" = 0.012; "* italic(p[Adap])*" = 0.6"))
outputDF_0.8 <- as.data.frame( t( replicate( 100, summary(factor(rbinom(20000*p_lik, nSpecies, 0.8), levels = 0:nSpecies)) ) )  )
plot(melt(outputDF_0.8), xlab = "# of Species With Local Adaptation Gene", ylab = "Number of true positives", main = expression(italic(p[Lik])*" = 0.012; "* italic(p[Adap])*" = 0.8"))

outputDF_0.2 <- as.data.frame( t( replicate( 100, summary(factor(rbinom(20000*p_lik,nSpecies, 0.2), levels = 0:nSpecies))*power_vec ) )  )
plot(melt(outputDF_0.2), xlab = "# of Species With Local Adaptation Gene", ylab = "Number detected", main = expression(italic(p[Lik])*" = 0.012; "* italic(p[Adap])*" = 0.2"))
outputDF_0.4 <- as.data.frame( t( replicate( 100, summary(factor(rbinom(20000*p_lik,nSpecies, 0.4), levels = 0:nSpecies))*power_vec ) )  )
plot(melt(outputDF_0.4), xlab = "# of Species With Local Adaptation Gene", ylab = "Number detected", main = expression(italic(p[Lik])*" = 0.012; "* italic(p[Adap])*" = 0.4"))
outputDF_0.6 <- as.data.frame( t( replicate( 100, summary(factor(rbinom(20000*p_lik,nSpecies, 0.6), levels = 0:nSpecies))*power_vec ) )  )
plot(melt(outputDF_0.6), xlab = "# of Species With Local Adaptation Gene", ylab = "Number detected", main = expression(italic(p[Lik])*" = 0.012; "* italic(p[Adap])*" = 0.6"))
outputDF_0.8 <- as.data.frame( t( replicate( 100, summary(factor(rbinom(20000*p_lik,nSpecies, 0.8), levels = 0:nSpecies))*power_vec ) )  )
plot(melt(outputDF_0.8), xlab = "# of Species With Local Adaptation Gene", ylab = "Number detected", main = expression(italic(p[Lik])*" = 0.012; "* italic(p[Adap])*" = 0.8"))
library(extraDistr)        

