extract_features <- function(df){
  measurement <- df$measurement
  time <- df$time
  
  min <- min(measurement, na.rm = T)
  max <- max(measurement, na.rm = T)
  max_min <- max - min
  median <- median(measurement, na.rm = T)
  mean <- mean(measurement, na.rm = T)
  quantile_5th <- quantile(measurement, probs = .05)
  quantile_95th <- quantile(measurement, probs = .95)
  variance <- var(measurement, na.rm = T)
  amplitude <- (max + min) / 2
  rms <- sqrt(mean(measurement ^ 2, na.rm = T))
  skewness <- skewness(measurement, na.rm = T)
  kurtosis <- kurtosis(measurement, na.rm = T)
  slope <- lm(measurement ~ seq_along(measurement))$coefficients[2]
  entropy <- -sum((measurement / sum(measurement)) * log(measurement / sum(measurement)))
  autocorrelation <- cor(measurement[-length(measurement)], measurement[-1])
  
  features <- data.frame(min, max, max_min, median, mean, quantile_5th, quantile_95th,
                         variance, amplitude, rms, skewness, kurtosis, slope,
                         entropy, autocorrelation)
  
  return(features)
}
