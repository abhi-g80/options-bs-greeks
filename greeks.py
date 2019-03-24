#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculating greeks
"""

import scipy
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
    Type Double (8 decimal precision), price of call option.
    """
    d_1 = (scipy.log(S/X) + (r_f + (sigma**2) * 0.5) * T)/(sigma
                                                           * scipy.sqrt(T))
    d_2 = d_1 - sigma * scipy.sqrt(T)
    
    C = S * norm.cdf(d_1) - X * norm.cdf(d_2) * scipy.exp(-r_f * T)
    
    return round(C, 8)


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
    Type double (8 decimal precision), price of put option.
    """
    P =  C - S + X * scipy.exp(-r_f * T)
    
    return round(P, 8)


def option_greeks(option_type, S, X, r_f, sigma, T):
    """
    Calculation of option greeks based on the formula discussed in
    class.
    
    First we calculate d1 and d2. Then based on the type of option
    we get our greeks.
    """
    d_1 = (scipy.log(S/X) + (r_f + (sigma**2) * 0.5) * T)/(sigma
                                                           * scipy.sqrt(T))
    d_2 = d_1 - sigma * scipy.sqrt(T)

    if option_type.lower() == "call":
        c_delta = norm.cdf(d_1)
        
        c_gamma = scipy.exp(-(d_1**2) * 0.5) / (S * sigma * scipy.sqrt(2 * scipy.pi * T))
        
        c_theta = -S * sigma * scipy.exp(-(d_1**2) * 0.5) / (2 * scipy.sqrt(2 * scipy.pi * T)) - r_f * X * scipy.exp(-r_f * T) * norm.cdf(d_2)

        c_vega = S * scipy.sqrt(T) * scipy.exp(-(d_1**2) * 0.5)/(scipy.sqrt(2 * scipy.pi))
        
        c_rho = X * T * scipy.exp(-r_f * T) * norm.cdf(d_2)

    elif option_type.lower() == "put":
        c_delta = norm.cdf(d_1) - 1
        
        c_gamma = 1 / (S * sigma * scipy.sqrt(2 * scipy.pi * T)) * scipy.exp(-(d_1**2) * 0.5)
        
        c_theta = -S * sigma * 1/(2 * scipy.sqrt(2 * scipy.pi * T)) * scipy.exp(-(d_1**2) * 0.5) + r_f * X * scipy.exp(-r_f * T) * norm.cdf(-d_2)

        c_vega = S * scipy.sqrt(T) * scipy.exp(-(d_1**2) * 0.5)/(scipy.sqrt(2 * scipy.pi))
        
        c_rho = -X * T * scipy.exp(-r_f * T) * norm.cdf(-d_2)
        

    return round(c_delta, 6), round(c_gamma, 6), round(c_theta, 6), round(c_vega, 6), round(c_rho, 6)


def sensitivity_calc(option_px, pct_chng, delta, theta, vega, rho,
                     S, r_f, sigma, T):
    """
    We are calculating sensitivity using the following formula (pg 6, lec 4) :-
    
    px_change = delta * (stock_px_change) + theta * (time_change) 
                + vega * (volatility change) + rho * (rf_change)
                + 1/2 * gamma * (stock_px_change)ˆ2
    
    Since we need to calculate individual sensitivity, when we want
    to calculate, suppose sensitivity towards S (stock price), other 
    changes will be taken as 0. So
    
    stock_px_change = px_change/delta.
    
    Similarly for others.
    """
    new_option_px =  option_px * (1 + pct_chng/100)
    chng_option_px = new_option_px - option_px
    
    # Change in delta to get the required pct change
    S_chng = chng_option_px/delta
    
    # Change in theta to get the required pct change
    T_chng = chng_option_px/theta

    # Change in vega to get the required pct change
    sigma_chng = chng_option_px/vega
    
    # Change in rho to get the required pct change
    rf_chng = chng_option_px/rho

    new_S = round(S + S_chng, 4)
    new_T = round(T + T_chng, 4)
    new_sigma = round(sigma + sigma_chng, 4)
    new_rf = round(r_f + rf_chng, 4)
    
    print(f"Option px = {round(new_option_px, 5)}, S = {new_S}, T = {new_T}, sigma = {new_sigma * 100}%, risk free = {new_rf * 100}%")
    

# Dictionaries to store greeks so that we can rank them later

delta_dic = {}
theta_dic = {}
vega_dic = {}
rho_dic = {}


print("Call 10, 12, 0.02, 0.5, 0.5")
print("----------------------------")
delta, gamma, theta, vega, rho = option_greeks("Call", 10, 12, 0.02, 0.5, 0.5)
print(f"Delta = {delta}\nGamma = {gamma}\nTheta = {theta}\nVega = {vega}\nRho = {rho}")
option_price = black_scholes_call(10, 12, 0.02, 0.5, 0.5)
print(f"Option price = {option_price}")

delta_dic.update({1 : abs(delta)})
theta_dic.update({1 : abs(theta)})
vega_dic.update({1 : abs(vega)})
rho_dic.update({1 : abs(rho)})

print("\nIncrease by 10%")
sensitivity_calc(option_price, 10, delta, theta, vega, rho, 10, 0.02, 0.5, 0.5)

print("\nDecrease by 20%")
sensitivity_calc(option_price, -20, delta, theta, vega, rho, 10, 0.02, 0.5, 0.5)

print("\n")
print("Call 60, 48, 0.02, 0.29, 0.4")
print("----------------------------")
delta, gamma, theta, vega, rho = option_greeks("Call", 60, 48, 0.02, 0.29, 0.4)
print(f"Delta = {delta}\nGamma = {gamma}\nTheta = {theta}\nVega = {vega}\nRho = {rho}")
option_price = black_scholes_call(60, 48, 0.02, 0.29, 0.4)
print(f"Option price = {option_price}")

delta_dic.update({2: abs(delta)})
theta_dic.update({2: abs(theta)})
vega_dic.update({2: abs(vega)})
rho_dic.update({2: abs(rho)})

print("\nIncrease by 10%")
sensitivity_calc(option_price, 10, delta, theta, vega, rho, 60, 0.02, 0.29, 0.4)

print("\nDecrease by 20%")
sensitivity_calc(option_price, -20, delta, theta, vega, rho, 60, 0.02, 0.29, 0.4)

print("\n")
print("Put 50, 60, 0.02, 0.4, 1.2")
print("----------------------------")
delta, gamma, theta, vega, rho = option_greeks("Put", 50, 60, 0.02, 0.4, 1.2)
print(f"Delta = {delta}\nGamma = {gamma}\nTheta = {theta}\nVega = {vega}\nRho = {rho}")
option_price = black_scholes_call(50, 60, 0.02, 0.4, 1.2)
option_price = put_from_call(option_price, 50, 60, 0.02, 1.2)
print(f"Option price = {option_price}")

delta_dic.update({3: abs(delta)})
theta_dic.update({3: abs(theta)})
vega_dic.update({3: abs(vega)})
rho_dic.update({3: abs(rho)})

print("\nIncrease by 10%")
sensitivity_calc(option_price, 10, delta, theta, vega, rho, 50, 0.02, 0.4, 1.2)

print("\nDecrease by 20%")
sensitivity_calc(option_price, -20, delta, theta, vega, rho, 50, 0.02, 0.4, 1.2)

print("\n")
print("Put 40, 40, 0.02, 0.33, 0.75")
print("----------------------------")
delta, gamma, theta, vega, rho = option_greeks("Put", 40, 40, 0.02, 0.33, 0.75)
print(f"Delta = {delta}\nGamma = {gamma}\nTheta = {theta}\nVega = {vega}\nRho = {rho}")
option_price = black_scholes_call(40, 40, 0.02, 0.33, 0.75)
option_price = put_from_call(option_price, 40, 40, 0.02, 0.75)
print(f"Option price = {option_price}")

delta_dic.update({4: abs(delta)})
theta_dic.update({4: abs(theta)})
vega_dic.update({4: abs(vega)})
rho_dic.update({4: abs(rho)})

print("\nIncrease by 10%")
sensitivity_calc(option_price, 10, delta, theta, vega, rho, 40, 0.02, 0.33, 0.75)

print("\nDecrease by 20%")
sensitivity_calc(option_price, -20, delta, theta, vega, rho, 40, 0.02, 0.33, 0.75)


print("\n\nRankings")
print("--------")

print("\nDelta sensitivity")
[print(f"{key} => {value}") for (key, value) in 
                             sorted(delta_dic.items(), 
                                    key=lambda kv: kv[1], reverse=True)]
                             
print("\nTheta sensitivity")
[print(f"{key} => {value}") for (key, value) in 
                             sorted(theta_dic.items(), 
                                    key=lambda kv: kv[1], reverse=True)]

print("\nVega sensitivity")
[print(f"{key} => {value}") for (key, value) in 
                             sorted(vega_dic.items(), 
                                    key=lambda kv: kv[1], reverse=True)]

print("\nRho sensitivity")
[print(f"{key} => {value}") for (key, value) in 
                             sorted(rho_dic.items(), 
                                    key=lambda kv: kv[1], reverse=True)]
