# -*- coding: utf-8 -*-
"""
Created on Thu Mar  6 00:15:16 2025

@author: ASUS
"""
# !pip install yfinance

import yfinance as yf
import pandas as pd

tickers = ["TSLA", "AMZN", "NVDA", "SPY"] 

end_date = pd.to_datetime("today").strftime('%Y-%m-%d')
start_date = (pd.to_datetime("today") - pd.DateOffset(years=2)).strftime('%Y-%m-%d')

# download data
data = yf.download(tickers, start=start_date, end=end_date, interval="1d", auto_adjust=True)

#print(data.head())

data.to_csv("C:/Users/ASUS/Desktop/stock_data.csv")


import pandas as pd
import numpy as np
import numpy.linalg as la

# Load data from CSV
file_path = "C:/Users/ASUS/Desktop/stock_data.csv"
data = pd.read_csv(file_path, header=[0, 1], index_col=0, parse_dates=True)

# Extract Close prices only
close_prices = data['Close']

close_prices = close_prices.drop(columns='SPY', errors='ignore')

# Calculate daily log returns
returns = np.log(close_prices / close_prices.shift(1))

# Drop missing values
returns = returns.dropna()

# Calculate mean returns (mu)
mu = returns.mean()

# Calculate covariance matrix (Sigma)
cov_matrix = returns.cov()

# Compute the inverse of the covariance matrix
inv_cov_matrix = la.inv(cov_matrix)

# Create a vector of ones
ones = np.ones(len(mu))

# Calculate optimal portfolio weights using the analytical solution
weights = inv_cov_matrix @ mu / (ones.T @ inv_cov_matrix @ mu)

# Print results
print("Mean daily returns (μ):")
print(mu)
print("\nCovariance matrix (Σ):")
print(cov_matrix)
print("\nOptimal portfolio weights:")
for stock, weight in zip(close_prices.columns, weights):
    print(f"{stock}: {weight:.4f}")

        
    
    
#No Short-Selling   
from scipy.optimize import minimize
import numpy as np


def negative_sharpe(weights, mu, cov_matrix):
    portfolio_return = np.dot(weights, mu)
    portfolio_volatility = np.sqrt(np.dot(weights.T, np.dot(cov_matrix, weights)))
    return -portfolio_return / portfolio_volatility

constraints = ({'type': 'eq', 'fun': lambda weights: np.sum(weights) - 1})

bounds = [(0, 1) for _ in range(len(mu))]

initial_weights = np.array([1 / len(mu)] * len(mu))

result = minimize(negative_sharpe, initial_weights, args=(mu, cov_matrix),
                  method='SLSQP', bounds=bounds, constraints=constraints)

optimal_weights = result.x

print("\nOptimal portfolio weights (no short-selling):")
for stock, weight in zip(close_prices.columns, optimal_weights):
    print(f"{stock}: {weight:.4f}")

portfolio_return = np.dot(optimal_weights, mu)
portfolio_volatility = np.sqrt(np.dot(optimal_weights.T, np.dot(cov_matrix, optimal_weights)))
sharpe_ratio = portfolio_return / portfolio_volatility

print("\nOptimal portfolio expected daily return: {:.4f}".format(portfolio_return))
print("Optimal portfolio daily volatility: {:.4f}".format(portfolio_volatility))
print("Optimal portfolio Sharpe Ratio: {:.4f}".format(sharpe_ratio))



# wtsMarket
data = pd.read_csv(file_path, header=[0, 1], index_col=0, parse_dates=True)

close_prices = data['Close'].drop(columns='SPY', errors='ignore')

market_cap = close_prices.iloc[-1]


total_market_cap = market_cap.sum()
wtsMarket = market_cap / total_market_cap

wtsMarket 

