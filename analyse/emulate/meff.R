# meff.m
#
# Produce plots of the emulator main effects
#
# inputs: nominal_values, model (emulator), 
# outputs: none
# plots produced: for each variable, a plot of the model error as a function of parameter

np  <- 1000 #number of x-axis points
nvar <- length(nominal_values) #number of variables

normalized_xx      <- matrix(NaN, nrow = np, ncol = nvar) #initialize array which will store along variable values
dimensional_xx     <- matrix(NaN, nrow = np, ncol = nvar) #dimensional version of above
normalized_yy_mean <- matrix(NaN, nrow = np, ncol = nvar) #will store emulator predictions in dimensionless form
normalized_yy_u95  <- matrix(NaN, nrow = np, ncol = nvar) #will store dimensionless emulator u95 
normalized_yy_l95  <- matrix(NaN, nrow = np, ncol = nvar) #will store dimensionless emulator l95

dimensionless_prior_mean     <- (prior_mean - input_normalization_mean)/input_normalization_sd
dimensionless_prior_sd       <- prior_sd / input_normalization_sd
dimensionless_nominal_values <- (nominal_values - input_normalization_mean)/input_normalization_sd

for (i in 1:nvar){ #loop over each variable
  
  #get the x-axis variables as +- 2sd from the prior mean
  normalized_xx[,i]  <- seq(from = (dimensionless_prior_mean[i] - 3*dimensionless_prior_sd[i]), to = (dimensionless_prior_mean[i] + 3*dimensionless_prior_sd[i]) , length.out = np)
  dimensional_xx[,i] <- seq(from = (prior_mean[i]               - 3*prior_sd[i]),               to = (prior_mean[i]               + 3*prior_sd[i]) , length.out = np)
  
  #create the test values first as a repeat of the nominal values repeated, and then replace the appropriate column
  test_values <- matrix(rep(dimensionless_nominal_values, times = np), nrow = np, ncol = 7, byrow = TRUE)
  test_values[,i] <- normalized_xx[,i] #replace the ith column with normalized xx values
  
  #test the model on these values
  tmp <- predict(model,test_values)
  
  #store the mean and upper and lower quartiles
  normalized_yy_mean[,i] <- tmp$mean
  normalized_yy_u95[,i] <- tmp$upper95
  normalized_yy_l95[,i] <- tmp$lower95
}

col = hcl.colors(nvar, palette = "Dark 3")

for (i in 1:nvar){
  
  plot(0,0,xlim = c(min(dimensional_xx[,i]),max(dimensional_xx[,i])), ylim = c(-1,1),type = "n",xlab = input_colnames[i], ylab = "GrV error (/1e12)", cex = 1.1, cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
  abline(v = 0, lwd = 0.5)
  abline(h = 0, lwd = 0.5)
  lines(dimensional_xx[,i], normalized_yy_mean[,i], col = col[1], type = 'l', lwd = 2)
  polygon(c(dimensional_xx[,i], rev(dimensional_xx[,i])), c(normalized_yy_l95[,i], rev(normalized_yy_u95[,i])), col = adjustcolor(col[1], alpha.f=0.5),  border = NA)
  
}

