# Author: Mitch Mitchell (jbm8efn@virginia.edu)
# Functions for calculating the price of European options

import math


# S: Current stock price
# K: Strike Price
# n: Number of periods in binomial model
# T: Time to expiration (in years)
# sigma: Standard deviation of stock returns
# r_f: Risk-free rate
def binomial_price(S, K, n, T, sigma, r_f, type={'call', 'put'}):
    u = pow(math.e, sigma * math.sqrt(T / n))
    d = 1 / u
    R = pow(math.e, r_f * T / n)
    q = (R - d) / (u - d)
    q_c = 1 - q
    df = [[-1.0 for i in range(n+1)] for j in range(n+1)]
    temp = S
    height = 0
    for i in range(n+1):
        df[0][i] = temp
        for j in range(height):
            temp *= d * d
            df[j+1][i] = temp
        temp = df[0][i] * u
        height += 1
    if type == 'call':
        for i in range(n+1):
            if df[i][n] > K:
                df[i][n] -= K
            else:
                df[i][n] = 0
        height = n
        for i in reversed(range(n)):
            for j in reversed(range(height)):
                expected = (1 / R) * (q * df[j][i+1] + q_c * df[j+1][i+1])
                df[j][i] = expected
            height -= 1
        return df[0][0]
    elif type == 'put':
        for i in range(n+1):
            if df[i][n] < K:
                df[i][n] = K - df[i][n]
            else:
                df[i][n] = 0
        height = n
        for i in reversed(range(n)):
            for j in reversed(range(height)):
                expected = (1 / R) * (q * df[j][i+1] + q_c * df[j+1][i+1])
                df[j][i] = expected
            height -= 1
        return df[0][0]


# z: z-score
def cdf(z):
    p = 1 + math.erf(z / math.sqrt(2.0))
    return (p / 2)


# S: Current stock price
# K: Strike Price
# T: Time to expiration (in years)
# sigma: Standard deviation of stock returns
# r_f: Risk-free rate
def black_scholes_price(S, K, T, sigma, r_f, type={'call', 'put'}):
    d1 = math.log(S / K) + T * (r_f + (pow(sigma, 2) / 2))
    d1 /= (sigma * math.sqrt(T))
    d2 = d1 - (sigma * math.sqrt(T))
    N1 = cdf(d1)
    N2 = cdf(d2)
    C = (S * N1 - K * pow(math.e, -r_f * T) * N2)
    if type == 'call':
        return C
    elif type == 'put':
        return (C - S + K * pow(math.e, -r_f * T))


# Main Method
if __name__ == '__main__':
    S = 100
    K = 100
    n = 5
    T = 1
    sigma = 0.1
    r_f = 0.01
    type = 'call'
    bin_price = binomial_price(S, K, n, T, sigma, r_f, type)
    bs_price = black_scholes_price(S, K, T, sigma, r_f, type)
    print('European', type, 'option prices:')
    print(n, 'period binomial model: $', bin_price)
    print('Black-Scholes model: $', bs_price)
