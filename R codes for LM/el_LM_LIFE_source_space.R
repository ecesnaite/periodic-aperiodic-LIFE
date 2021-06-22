# (C) Elena Cesnaite, email: e.cesnaite@gmail.com, page: https://www.researchgate.net/profile/Elena-Cesnaite

# This code was created to analyze data described in 'Alterations in rhythmic and non-rhythmic resting-state
# EEG activity and their link to cognition in older age' paper.
# The code runs linear models where factors derived from the Factor analysis are dependent variables and rsEEG 
# parameters together with age, sex, and education, are independednt variables. We also include interaction terms 
# between age and rsEEG parameters.
# Last updated 22.06.2021

## ---------------------------------------------------------------------------------------------------------- ##

rm(list=ls())

# load rsEEG parameters for 10 ROIs at source space
source_all = read.csv('')

# load 3 factors that form cognition
factors = read.csv('')

#inspect data for outliers
layout(matrix(c(1,2,3,4,5,6),2,3)) 

for (i in 1:length(source_all)){
  hist(source_all[,i], main = colnames(source_all[i]))
}

sex_final = vector()
for (g in 1:length(source_all$sex)){
 if (source_all$sex[g] == 'male') # male - 1, female - 2
   sex_final[g] = 1
else
  sex_final[g] = 2
}

# scale predictors
source_sc <- lapply(source_all, scale) 

regions_source = c(1:10)

#LM for the 1st factor
coef_freq_f1 = numeric()
p_freq_f1 = numeric()
p_all = numeric()
res_f1 = list()

for (g in 1:length(regions_source)){
  print(paste(regions_source[g]))
  res_f1[g] = lm(paste0('scale(factors$Factor1) ~ (source_sc$freq_', regions_source[g],'+ source_sc$power_',regions_source[g],
                     '+ source_sc$theta_',regions_source[g],' + source_sc$slope_', regions_source[g],') * source_sc$age + source_sc$edu + sex_final'))

  print(summary(res_f1))
  coef_freq_f1[g] = summary(res_f1)$coefficients[2,1]
  p_freq_f1[g] = summary(res_f1)$coefficients[2,4]
  p_all = rbind(p_all,summary(res_f1)$coefficients[,4])
  
}

#LM for the 2nd factor
coef_freq_f2 = numeric()
p_freq_f2 = numeric()
res_f2 = list()

for (g in 1:length(regions_source)){
  print(paste(regions_source[g]))
  res_f2[g] = lm(paste0('scale(factors$Factor2) ~ (source_sc$freq_', regions_source[g],'+ source_sc$power_',regions_source[g],
                     '+ source_sc$theta_',regions_source[g],' + source_sc$slope_', regions_source[g],') * source_sc$age + source_sc$edu + sex_final'))
  print(summary(res_f2))
  coef_freq_f2[g] = summary(res_f2)$coefficients[8,1]
  p_freq_f2[g] = summary(res_f2)$coefficients[8,4]
  p_all = rbind(p_all,summary(res_f2)$coefficients[,4])
  
}

#LM for the 3rd factor
coef_freq_f3 = numeric()
p_freq_f3 = numeric()
res_f3= list()

for (g in 1:length(regions_source)){
  print(paste(regions_source[g]))
  res_f3[g] = lm(paste0('scale(factors$Factor3) ~ (source_sc$freq_', regions_source[g],'+ source_sc$power_',regions_source[g],
                     '+ source_sc$theta_',regions_source[g],' + source_sc$slope_', regions_source[g],') * source_sc$age + source_sc$edu + sex_final'))
  print(summary(res_f3))
  coef_freq_f3[g] = summary(res_f3)$coefficients[2,1]
  p_freq_f3[g] = summary(res_f3)$coefficients[2,4]
  p_all = rbind(p_all,summary(res_f3)$coefficients[,4])
  
}



