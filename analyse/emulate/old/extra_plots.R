#make extra plots, need to run emulate_sample first

# plot the emulated model output as a function of input parameter for each of the parameters, at the naive posterior mean
np <- 1000 #number of x points
dimensionless_prior_mean <- normalized_prior_mean
dimensionless_prior_sd <- normalized_prior_sd
normalized_xx <- matrix(NaN, nrow = np, ncol= length(dimensional_prior_mean))
dimensional_xx <- matrix(NaN, nrow = np, ncol= length(dimensional_prior_mean))
normalized_yy_mean <- matrix(NaN, nrow = np, ncol= length(dimensional_prior_mean))
normalized_yy_u95 <- matrix(NaN, nrow = np, ncol= length(dimensional_prior_mean))
normalized_yy_l95 <- matrix(NaN, nrow = np, ncol= length(dimensional_prior_mean))


for (i in 1:7){
  normalized_xx[,i] <- seq(from = (dimensionless_prior_mean[i] - 3*dimensionless_prior_sd[i]), to = (dimensionless_prior_mean[i] + 3*dimensionless_prior_sd[i]) , length.out = np)
  dimensional_xx[,i] <- (normalized_xx[,i]*dimensional_model_input_sd[i] + dimensional_model_input_mean[i])
  test_values <- matrix(rep(normalized_naive_posterior_mean, times = np), nrow = np, ncol = 7, byrow = TRUE) #test values at naive posterior mean everywhere
  test_values[,i] <- normalized_xx[,i] #replace the ith column with normalized xx values
  
  #test the model on these values
  tmp <- predict(model,test_values)
  
  normalized_yy_mean[,i] <- tmp$mean
  normalized_yy_u95[,i] <- tmp$upper95
  normalized_yy_l95[,i] <- tmp$lower95
}

normalized_xx <- data.frame(normalized_xx)
dimensional_xx <- data.frame(dimensional_xx)
normalized_yy_mean <- data.frame(normalized_yy_mean)
normalized_yy_l95 <- data.frame(normalized_yy_l95)
normalized_yy_u95 <- data.frame(normalized_yy_u95)

plot_list <- list()

for (i in 1:7) {
  # Create a data frame for the current column
  df <- data.frame(
    x = dimensional_xx[[i]],
    y = normalized_yy_mean[[i]] - normalized_actual,
    ymin = normalized_yy_l95[[i]] - normalized_actual,
    ymax = normalized_yy_u95[[i]] - normalized_actual
  )
  
  
  # Create the plot
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_line() +
    geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = dimensional_naive_posterior_mean[i], linetype = "dashed") +
    
    ylim(-.1, .1) +
    xlim(xlims[i,]) + 
    labs(x = input_headers[i], y = "emu error") +
    theme_minimal()
  
  # Add the plot to the list
  plot_list[[i]] <- p
}

grid.arrange(grobs = plot_list, ncol = 3, nrow = 3)

  

