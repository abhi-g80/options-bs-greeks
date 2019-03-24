#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python script to plot the option value
"""
import scipy
import matplotlib.pyplot as plt
from scipy.stats import norm


def black_scholes_call(S, X, r_f, sigma, T):
    """
    Black scholes formula for implementing Call option pricing.
    
    Arguments
    ---------
    S = Underlying price
    X = Strike price
    r_f = Risk free return (annualized)
    T = Time (in years)
    sigma = Volatility (annualized)
    
    Returns
    -------
    Type Double (4 decimal precision), price of call option.
    """
    d_1 = (scipy.log(S/X) + (r_f + (sigma**2) * 0.5) * T)/(sigma
                                                           * scipy.sqrt(T))
    d_2 = d_1 - sigma * scipy.sqrt(T)
    
    C = S * norm.cdf(d_1) - X * norm.cdf(d_2) * scipy.exp(-r_f * T)
    
    return round(C, 4)


def put_from_call(C, S, X, r_f, T):
    """
    Put option value using put call parity theorem.
    
    Arguments
    ---------
    C = Call option price
    S = Underlying price
    X = Strike price
    r_f = Risk free return (annualized)
    T = Time (in years)

    Returns
    -------
    Type double (4 decimal precision), price of put option.
    """
    P =  C - S + X * scipy.exp(-r_f * T)
    
    return round(P, 4)


S = 43
X = 45
r_f = 0.015
T = 1.5
sigma = 0.38


# Assignment a
call_value = black_scholes_call(S, X, r_f, sigma, T)
print(f"Call option value = {call_value}")
put_value = put_from_call(call_value, S, X, r_f, T)
print(f"Put option value = {put_value}")

# Assignment b
shs_px = list(range(20, 75, 5))

# Calculate intrinsic values
intrinsic_value = [max(S - X, 0) for S in shs_px]

# Calculate option prices for T = 1
call_px = []
T = 1

for S in shs_px:
    call_px.append(black_scholes_call(S, X, r_f, sigma, T))

# print(call_px)

plt.plot(shs_px, call_px, color='green', label='Expiration of 1 year')

# Calculate time value for call prices of exp 1 year
first_time_value = []

for i in map(lambda x,y: x - y, call_px, intrinsic_value):
    first_time_value.append(i)

# Calculate option prices for T = 0.15
call_px = []
T = 0.15

for S in shs_px:
    call_px.append(black_scholes_call(S, X, r_f, sigma, T))
    
# print(call_px)
    
plt.plot(shs_px, call_px, color='orange', label='Expiration of 0.15 year')

# Calculate time value for call prices of exp 0.15 year
second_time_value = []

for i in map(lambda x,y: x - y, call_px, intrinsic_value):
    second_time_value.append(i)

# Plot call prices
plt.xlabel('Share price')
plt.ylabel('Option price')
plt.title('Call option prices for different share price')
plt.legend()
plt.show()
# plt.savefig('call_prices.png')
plt.close('all')

# Plot time value
plt.plot(shs_px, first_time_value, color='blue', label='Time value for 1 year')
plt.plot(shs_px, second_time_value, color='red', 
         label='Time value for 0.15 year')

plt.xlabel('Share price')
plt.ylabel('Time value')
plt.title('Time value of call option for different share price')
plt.legend()
plt.show()
# plt.savefig('time_values.png')
plt.close('all')
