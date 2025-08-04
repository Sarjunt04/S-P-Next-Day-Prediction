# USmkt-Next-Day-Prediction

An OLS algorithm for predicting next day S&amp;P returns, using recent returns of constituent stocks.

This was a final year project for LSE module ST326 (Financial Statistics), with a goal to predict next day S&P returns using return series of US stocks. My approach used 10 stocks within a variety of sectors, including AAPL, LLY, V and MCD. The overall grade was an 85.

First, I read in data from the last 5 years and inspected the 11 series in a single plot. I used an EWMA model for each return series (case of GARCH(1,1) where alpha + beta = 1), and thus was able to calculate parameters lambda for each stock using MLE. This lambda was used in computing volatilites for the return series, and consequently normalise them.

Finally, I directly applied OLS to forecast next day returns and applied the algorithm to the last 5 years. On the test set, maximum Sharpe was 1.75, attained by using the 10 stocks (no lagged S&P returns) and a 10 period lookback window (last 10 day's returns used in LS calculation).
