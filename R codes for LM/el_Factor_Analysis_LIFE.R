# (C) Elena Cesnaite, email: e.cesnaite@gmail.com, page: https://www.researchgate.net/profile/Elena-Cesnaite

# This code was created to analyze data described in 'Alterations in rhythmic and non-rhythmic resting-state
# EEG activity and their link to cognition in older age' paper.
# The code runs Exploratory factor analysis on cognitive score data
# Last updated 22.06.2021

## ---------------------------------------------------------------------------------------------------------- ##
rm(list=ls())

setwd("")
#load cog score data
cognition = read.csv('')

#Factor Analysis

install.packages("psych")
library(psych)
install.packages("nFactors")
library(nFactors)

cog_final <- cognition[, c(2:9)] #exclude ID information

# decide for the number of latent factors based on the scree plots
ev <- eigen(cor(cog_final))
ap <- parallel(subject = nrow(cog_final), var = ncol(cog_final), rep = 100, cent = .05)
nS <- nScree(x=ev$values, aparallel = ap$eigen$qevpea)
plotnScree(nS)

library(semTools)
attr(efa.ekc(cog_final), 'nfactors') # Empirical Kaiser Criterion (EKC)

Factor_Analysis2<- factanal(cog_final,factors = 3,score ="regression", rotation = "promax") # factor analysis estimating 3 latent factors
print(Factor_Analysis2, digits=3, cutoff=0.30, sort=TRUE)
loadings = loadings(Factor_Analysis2)
plot(loadings)

factors<-Factor_Analysis2$scores
write.csv(factors, file = "Three factors promax.csv")
