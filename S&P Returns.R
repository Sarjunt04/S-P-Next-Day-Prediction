#ST326 Project

# 10 STOCKS: AAPL, TSLA, LLY, META, V, XOM, COST, NFLX, PEP, MCD
# Sector: Technology, Automotive, Pharmaceutical, Technology (Social Media), Payment, Energy, Warehouse Store, Media, Food/Beverage, Fast Food
# LOAD FUNCTIONS FIRST

##### Obtain data using the quantmod package #####
library(quantmod)
library(lubridate)
library(MASS)
library(ggplot2)

read.bossa.data <- function(vec.names) {
  p <- length(vec.names)
  n1 <- 20000
  dates <- matrix(99999999, p, n1)
  closes <- matrix(0, p, n1)
  max.n2 <- 0
  
  for (i in 1:p) {
    filename <- paste0(vec.names[i], ".txt")
    tmp <- scan(filename, list(date=numeric(), NULL, NULL, NULL, NULL, NULL, close=numeric()), skip=1, sep="")
    
    n2 <- length(tmp$date)
    max.n2 <- max(n2, max.n2)
    
    dates[i,1:n2] <- tmp$date
    closes[i,1:n2] <- tmp$close
  }
  
  dates <- dates[,1:max.n2]
  closes <- closes[,1:max.n2]
  
  days <- rep(0, n1)
  arranged.closes <- matrix(0, p, n1)
  date.indices <- starting.indices <- rep(1, p)
  already.started <- rep(0, p)
  day <- 1
  
  
  while(max(date.indices) <= max.n2) {
    current.dates <- current.closes <- rep(0, p)
    for (i in 1:p) {
      current.dates[i] <- dates[i,date.indices[i]]
      current.closes[i] <- closes[i,date.indices[i]]
    }
    min.indices <- which(current.dates == min(current.dates))
    days[day] <- current.dates[min.indices[1]]
    arranged.closes[min.indices,day] <- log(current.closes[min.indices])
    arranged.closes[-min.indices,day] <- arranged.closes[-min.indices, max(day-1, 1)]
    already.started[min.indices] <- 1
    starting.indices[-which(already.started == 1)] <- starting.indices[-which(already.started == 1)] + 1
    day <- day + 1
    date.indices[min.indices] <- date.indices[min.indices] + 1
  }
  
  
  days <- days[1:(day-1)]
  arranged.closes <- arranged.closes[,1:(day-1)]
  max.st.ind <- max(starting.indices)
  r <- matrix(0, p, (day-max.st.ind-1))
  
  for (i in 1:p) {
    r[i,] <- diff(arranged.closes[i,max.st.ind:(day-1)])
    r[i,] <- r[i,] / sqrt(var(r[i,]))
    r[i,r[i,]==0] <- rnorm(sum(r[i,]==0))
  }
  
  return(list(dates=dates, closes=closes, days=days, arranged.closes=arranged.closes, starting.indices=starting.indices, r=r))
}

tickers <- c('^GSPC','AAPL','TSLA','LLY','META','V','XOM','COST','NFLX','PEP','MCD')

for (ticker in tickers) {
  getSymbols(ticker)
  
  # remove ^ from '^GSPC'
  ticker <- gsub("^\\^", "", ticker)
  
  ticker.dat <- as.data.frame(get(ticker))
  
  ticker.dat <- cbind(Date = as.numeric(as.Date(rownames(ticker.dat))), ticker.dat)

  # write files
  write.table(ticker.dat, paste0(ticker,".txt"), row.names = FALSE, col.names = FALSE, sep = "\t")
}


#### Read files with bossa.data function ####


clean_tickers <- c('GSPC','AAPL','TSLA','LLY','META','V','XOM','COST','NFLX','PEP','MCD')

ind = read.bossa.data(clean_tickers)


##### Plot log prices #####


colours <- rainbow(length(clean_tickers))
five_ago <- Sys.Date() %m-% years(5)
  
plot(1, type = 'n', 
     xlim = c(five_ago, Sys.Date()), ylim = range(ind$arranged.closes),
     xlab = 'Date', ylab = 'log(Price)', main = 'Log Prices over the Last 5 Years',
     xaxt = 'n')

# plot for last 5 years
for (i in 1:length(clean_tickers)) {
  dates <- as.Date(ind$days)
  log_prices <- ind$arranged.closes[i, ]
  l5_dates <- (dates >= five_ago) & (dates <= Sys.Date())
  
  lines(dates[l5_dates], log_prices[l5_dates], col = colours[i])
}

# format x-axis as points spread 3 months apart
month_seq <- seq(from = floor_date(five_ago, "month"), 
                 to = Sys.Date(), 
                 by = "3 months")

axis(side = 1, at = as.numeric(month_seq),
     labels = format(month_seq, "%b %Y"))
     
legend("bottomright", legend = clean_tickers, col = colours, lty = 1, cex = 0.5)


#### Find optimal lambda from training data ####

mle_lambda <- function(returns) {
  log_likelihood <- function(lambda, returns) {
    n <- length(returns)
    sigma_sq <- numeric(n)
    sigma_sq[1] <- var(returns)  
    
    for (t in 2:n) {
      sigma_sq[t] <- lambda * sigma_sq[t-1] + (1 - lambda) * returns[t-1]^2
    }
    
    ll <- sum(dnorm(returns, mean = 0, sd = sqrt(sigma_sq), log = TRUE))
    return(-ll)  
  }
  
  optim_result <- optim(par = 0.94, fn = log_likelihood, returns = returns, 
                        method = "L-BFGS-B", lower = 0.01, upper = 0.99)
  return(optim_result$par)
}

# consider last 5 years
valid_days <- which(as.Date(ind$days) >= five_ago)

ind$days <- ind$days[valid_days]
ind$r <- ind$r[, (ncol(ind$r) - length(valid_days) + 1):ncol(ind$r), drop = FALSE]

# length to later assign training data
train_length <- floor(0.5 * ncol(ind$r))  
val_length <- floor(0.25 * ncol(ind$r))  
test_length <- ncol(ind$r) - train_length - val_length  

# create separate datasets
train.dat <- list(r = ind$r[, 1:train_length], normal = NA)
val.dat <- list(r = ind$r[, (train_length + 1):(train_length + val_length)], normal = NA)
test.dat <- list(r = ind$r[, (train_length + val_length + 1):ncol(ind$r)], normal = NA)

optimal_lambdas <- numeric(11)
names(optimal_lambdas) <- clean_tickers

for (i in 1:nrow(ind$r)) {
  optimal_lambdas[i] <- mle_lambda(train.dat$r[i, ])
}

optimal_lambdas


#### Estimate volatilities using calculated lambdas ####

vol.exp.sm <- function(x, lambda) {
  
  # Exponential smoothing of x^2 with parameter lambda
  
  sigma2 <- x^2
  n <- length(x)
  
  for (i in 2:n)
    sigma2[i] <- sigma2[i-1] * lambda + x[i-1]^2 * (1-lambda)
  
  sigma <- sqrt(sigma2)
  
  resid <- x/sigma
  resid[is.na(resid)] <- 0
  sq.resid <- resid^2
  
  list(sigma2=sigma2, sigma=sigma, resid = resid, sq.resid = sq.resid)
  
}

volatilities_train <- list()
volatilities_val <- list()
volatilities_test <- list()

for (i in 1:11) {
  lambda <- optimal_lambdas[i]
  vol_train <- vol.exp.sm(train.dat$r[i,], lambda)  # Training set volatility
  vol_val <- vol.exp.sm(val.dat$r[i,], lambda)      # Validation set volatility
  vol_test <- vol.exp.sm(test.dat$r[i,], lambda)    # Test set volatility
  
  volatilities_train[[clean_tickers[i]]] <- vol_train$sigma
  volatilities_val[[clean_tickers[i]]] <- vol_val$sigma
  volatilities_test[[clean_tickers[i]]] <- vol_test$sigma
}

for (i in 1:length(clean_tickers)) {
  plot(volatilities_train[[clean_tickers[i]]], type = 'l',
       main = paste("Volatility for", clean_tickers[i]),
       xlab = "Training Set Time", ylab = "Estimated Volatility")
}


#### Normalise returns with computed lambdas ####


normalise_r <- function(data, volatilities) {
  normal_r <- matrix(NA, nrow = nrow(data$r), ncol = ncol(data$r))
  
  # by row as each row in ind$r is a different stock
  for (i in 1:nrow(data$r)) {
    normal_r[i, ] <- data$r[i, ] / volatilities[[i]]
  }
  
  # new double in ind
  data$normal <- normal_r  
  
  return(data)
}

train.dat <- normalise_r(train.dat, volatilities_train)  # Apply normalization to the training data
val.dat <- normalise_r(val.dat, volatilities_val)        # Apply normalization to the validation data
test.dat <- normalise_r(test.dat, volatilities_test) 


#### Add S&P return elements to datasets ####


datasets <- list(train.dat = train.dat, val.dat = val.dat, test.dat = test.dat)

for (name in names(datasets)) {
  datasets[[name]]$spr <- datasets[[name]]$r[1, , drop = FALSE]  
  datasets[[name]]$spnormal <- datasets[[name]]$normal[1, , drop = FALSE]  
  
  # remove s&p from inital elements
  datasets[[name]]$r <- datasets[[name]]$r[-1, , drop = FALSE]
  datasets[[name]]$normal <- datasets[[name]]$normal[-1, , drop = FALSE]
  
}

# Assign back to original variables
train.dat <- datasets$train.dat
val.dat <- datasets$val.dat
test.dat <- datasets$test.dat


#### An prediction algorithm: invest 1 for +ve returns and -1 for -ve ####


predictive_algorithm <- function(t0, window_length, dataset, name, q_lags = 0) {
  # t0: starting index in the time series for prediction
  # window_length: rolling window length (D)
  # dataset: train.dat, val.dat, or test.dat
  # q_lags: number of lagged S&P returns to include in the model
  
  sp_r <- dataset$spnormal  
  
  if (q_lags > 0) {
    # create lagged S&P returns
    lagged_r <- sapply(1:q_lags, function(lag) {
      sp_r[(1 + lag):(length(sp_r) - (q_lags - lag))]
    })
    
    # adjust dimensions and include lagged returns in $normal
    lagged_matrix <- t(lagged_r)
    dataset$normal <- dataset$normal[, (1 + q_lags):ncol(dataset$normal), drop = FALSE]
    dataset$spr <- dataset$spr[(1 + q_lags):length(sp_r)]
    
    x_n_lag <- rbind(dataset$normal, lagged_matrix)
  } else {
    x_n_lag <- dataset$normal
  }
  
  predicted_r <- numeric(length(dataset$spr))
  
  # loop through the time series using a rolling window
  for (t in seq(max(t0, window_length), ncol(dataset$normal) - 1)) {
    train_indices <- (t - window_length + 1):t
  
    x_window <- t(x_n_lag[, train_indices, drop = FALSE])
    y_window <- dataset$spr[train_indices]
    
    # estimate alpha using least squares
    alpha_hat <- MASS::ginv(t(x_window) %*% x_window) %*% t(x_window) %*% y_window
    
    # predict return at time t+1
    x_next <- x_n_lag[, t + 1, drop = FALSE]
    predicted_r[t] <- sum(alpha_hat * x_next)
  }
  
  # investment strategy: 1 if predicted return > 0, -1 otherwise
  strategy <- ifelse(predicted_r > 0, 1, -1)
  
  # index to times where predictions exist
  strategy_valid <- predicted_r != 0
 
  strategy_r <- strategy[strategy_valid] * dataset$spr[strategy_valid] / 100

  sharpe_ratio <- (mean(strategy_r, na.rm = TRUE) / sd(strategy_r, na.rm = TRUE)) * sqrt(250)
  
  print(paste("Sharpe Ratio:", sharpe_ratio, "for window length D =", window_length))

  plot(strategy_r, type = "l", main = paste(name, "Strategy Returns Over Time ( Window Length:", window_length,")"))
}

# iterate through window lengths for each dataset
window_lengths <- seq(10, 150, 20)

results <- lapply(names(datasets), function(name) {
  cat("Results for", name, "dataset:\n")
  sapply(window_lengths, function(window) {
    predictive_algorithm(t0 = window, window_length = window, dataset = datasets[[name]], name = name, q_lags = 1)
  })
})


#### Algorithm with factors as a tuning parameter ####


pcr_predictive_algorithm <- function(t0, window_length, dataset, num_factors, q_lags = 0) {
  # num_factors: number of principal components to use

  sp_r <- dataset$spnormal  

  if (q_lags > 0) {
    # create lagged S&P returns
    lagged_r <- sapply(1:q_lags, function(lag) {
      sp_r[(1 + lag):(length(sp_r) - (q_lags - lag))]
    })
    
    # adjust dimensions and include lagged returns in $normal
    lagged_matrix <- t(lagged_r)
    dataset$normal <- dataset$normal[, (1 + q_lags):ncol(dataset$normal), drop = FALSE]
    dataset$spr <- dataset$spr[(1 + q_lags):length(sp_r)]
    
    x_n_lag <- rbind(dataset$normal, lagged_matrix)
  } else {
    x_n_lag <- dataset$normal
  }
  
  predicted_r <- numeric(length(dataset$spr))
  
  # loop through the time series using a rolling window
  for (t in seq(max(t0, window_length), ncol(dataset$normal) - 1)) {
    train_indices <- (t - window_length + 1):t
    
    x_window <- t(x_n_lag[, train_indices, drop = FALSE])
    y_window <- dataset$spr[train_indices]
    
    # principal component regression on x_window
    pca_result <- prcomp(x_window, center = TRUE, scale. = TRUE)
    
    # select number of factors
    pc_scores <- pca_result$x[, 1:num_factors, drop = FALSE]
    
    pc_scores_df <- as.data.frame(pc_scores)
    colnames(pc_scores_df) <- paste0("PC", 1:num_factors)
    
    # fit linear model using principal components
    model <- lm(y_window ~ ., data = pc_scores_df)
    
    # predict return at time t+1
    x_next <- x_n_lag[, t + 1, drop = FALSE]
    pc_next <- predict(pca_result, newdata = t(x_next))[1:num_factors]
    pc_next_df <- as.data.frame(t(pc_next))
    colnames(pc_next_df) <- paste0("PC", 1:num_factors)
    
    predicted_r[t] <- predict(model, newdata = pc_next_df)
  }
  
  # investment strategy: 1 if predicted return > 0, -1 otherwise
  strategy <- ifelse(predicted_r > 0, 1, -1)
  
  # index to times where predictions exist
  strategy_valid <- predicted_r != 0
  strategy_r <- strategy[strategy_valid] * dataset$spr[strategy_valid] / 100
  
  sharpe_ratio <- (mean(strategy_r, na.rm = TRUE) / sd(strategy_r, na.rm = TRUE)) * sqrt(250)
  
  print(paste0("Sharpe Ratio: ", sharpe_ratio, ", D = ", window_length, ", Factors: ", num_factors, ", Lags: ", q_lags))
  
  return(sharpe_ratio)
}

num_factors_list <- 1:2
q_lags_list <- 0:1  

# for q = 0,1 & factors = 1,2
results <- lapply(names(datasets), function(name) {
  cat("Results for", name, "dataset:\n")
  result_df <- data.frame()
  
  for (window in window_lengths) {
    for (num_factors in num_factors_list) {
      for (q_lags in q_lags_list) {
        sharpe <- pcr_predictive_algorithm(t0 = window, window_length = window, dataset = datasets[[name]], num_factors = num_factors, q_lags = q_lags)
        result_df <- rbind(result_df, data.frame(Window = window, Factors = num_factors, Lags = q_lags, Sharpe = sharpe))
      }
    }
  }
  
  return(result_df)
})

names(results) <- names(datasets)

# plots across D, lags & factors
for (name in names(datasets)) {
  df <- results[[name]]
  
  p <- ggplot(df, aes(x = Window, y = Sharpe, color = factor(Factors), linetype = factor(Lags), group = interaction(Factors, Lags))) +
    geom_line() +
    geom_point() +
    labs(title = paste("Sharpe Ratios for", name, "dataset"),
         x = "Window (D)",
         y = "Sharpe Ratio",
         color = "No. Factors",
         linetype = "No. Lags") +
    theme_minimal()
  
  print(p)
}
