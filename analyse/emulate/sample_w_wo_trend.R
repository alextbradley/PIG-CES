# sample_w_wo_trend.R
#
# Draw samples from the posterior distributions with and without a trend
#

nvar <- dim(posterior_samples_noindex)[2]
nbins <- 100 #number of bins for the posterior histogram
n_samples <- 10000

samples <- matrix(NaN, nrow = n_samples, ncol = nvar)  #initialize matrix

for (i in 1:nvar){
  data <- posterior_samples_noindex[,i]
  hist_result <- hist(data, plot = FALSE, breaks = nbins) # plot=FALSE to avoid plotting
  
  breaks <- hist_result$breaks
  counts <- hist_result$counts
  probabilities <- counts / sum(counts)
  
  samples[,i] <- sample(x = breaks[-length(breaks)], size = n_samples, replace = TRUE, prob = probabilities)

  
  #hh <- hist(sampled_values, freq = FALSE, col = rgb(0.5, 0, 0, 0.5), add =TRUE, breaks = nbins) #uncomment if you want to plot samples
}

## loop over the valued and samples
ice_loss_with_trend <- c()
ice_loss_no_trend <- c()

for (i in 1:n_samples){
  test_values <- samples[i,] 
  
  #normalize these values
  normalized_test_values <-  (test_values - input_normalization_mean)/input_normalization_sd
  
  #get the model prediction of the ice loss with trend
  tmp <- predict(model, array(unname(normalized_test_values), dim = c(1, length(normalized_test_values))))
  ice_loss_with_trend <- c(ice_loss_with_trend,tmp$mean)
  
  # replace the trend with zero
  test_values[4] <- 0
  
  # replace the bump with zero
  #test_values[5] <- 0
  
  normalized_test_values <-  (test_values - input_normalization_mean)/input_normalization_sd
  
  
  #model prediciton with no trend
  tmp <- predict(model, array(unname(normalized_test_values), dim = c(1, length(normalized_test_values))))
  ice_loss_no_trend <- c(ice_loss_no_trend,tmp$mean)
  
  
}

# make histograms of them both 
hist(ice_loss_with_trend, plot = TRUE, freq = FALSE, breaks = 50) # plot=FALSE to avoid plotting
hist(ice_loss_no_trend, plot = TRUE, freq = FALSE, breaks = 50, add = TRUE, col = rgb(1, 0, 0, 0.5)) # plot=FALSE to avoid plotting

