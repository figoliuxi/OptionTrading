# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 16:46:37 2018
This is the package that deal with BSM model
@author:
"""

from __future__ import division 
import numpy as np
from scipy.stats import norm
from scipy.optimize import newton
from AssetClass import *

def calcPriceFromVol(EqOption, r,requireGreeks=False):
    """a function to calculate option price based on black-scholes model
    Type is either C or P"""
    callput = EqOption.callPut
    S0 = EqOption.underlying.price
    K = EqOption.K
    T = EqOption.T
    q =  EqOption.underlying.dividend
    sigma = EqOption.vol
    if callput == 'C':
        type = 1.0
    elif callput =='P':
        type = -1.0
    d1 = (np.log(S0/K) + (r-q+0.5*sigma**2)*T) / (sigma*T**0.5)
    d2 = d1 - sigma*T**0.5
    optionValue = type*S0*np.exp(-q*T)*norm.cdf(type*d1) - type*K*np.exp(-r*T)*norm.cdf(type*d2)
    if requireGreeks == True:
        delta = type*np.exp(-q*T)*norm.cdf(type*d1)
        gamma = np.exp(-q*T)*norm.pdf(d1)/S0/sigma/T**0.5
        vega = S0*np.exp(-q*T)*T**0.5*norm.pdf(d1)
        theta = -np.exp(-q*T)*S0*norm.pdf(d1)*sigma/2/T**0.5 - type*r*K*np.exp(-r*T)*norm.cdf(type*d2) + type*q*S0*np.exp(-q*T)*norm.cdf(type*d1)
        rho = type*K*T*np.exp(-r*T)*norm.cdf(type*d2)
        return optionValue,delta,gamma,vega,theta,rho
    return optionValue



def calcVolFromPrice(EqOption, r, price):
    """a function to calculate implied volatility based on Black-Scholes formula"""
    #if intrinsic value is bigger than option value or the price input is NaN
    callput = EqOption.callPut
    S0 = EqOption.underlying.price
    K = EqOption.K
    T = EqOption.T
    q =  EqOption.underlying.dividend
    if (callput == 'C' and S0*np.exp(-q*T) - K*np.exp(-r*T) > price) or price == 'NaN':
        return 'NaN'
    if (callput == 'P' and K*np.exp(-r*T) - S0*np.exp(-q*T) > price) or price == 'NaN':
        return 'NaN'
    else:
        #make a sub-function only has sigma as parameter; easier to use secant method to find the root
        def ObjFun( sigma ):
            return calcPriceFromVol(EqOption, r, sigma, q)-price
        try:
            result = newton(ObjFun,0.2, tol=1e-4,maxiter=100)
            if result == 'NaN':
                print ("cannot find the volatility")
            return result
        except:
            return 'cannot find the volatility'




def main():
    stock1 = Stock('AAPL',100,0.0)
    vol = 0.2
    opt1 = EqOption('opt1','C',stock1,100,1,None,vol)
    price = calcPriceFromVol(opt1,0.02)
    opt1.setPrice(price)
    print (opt1.price)
#    fileObject = open('FOZOptions.csv','rb')
#    optsTableR = csv.reader( fileObject )
#    #skip the header
#    optsTableR.next()
#    for row in optsTableR:
#        if row[0][16] == "C":
#            callput = 'call'
#        elif row[0][16] =="P":
#            callput = 'put'
#        S0 = float(row[2])
#        K = float(row[8])
#        r = float(row[9])/100.0
#        maturity = date(int(row[7].split('/')[2]),int(row[7].split('/')[0]),int(row[7].split('/')[1]))
#        current =  date(int(row[3].split('/')[2].split(' ')[0]),int(row[3].split('/')[0]),int(row[3].split('/')[1]))
#        T = float((maturity - current).days)/365.0
#        #I calculate the implied vol using mid price
#        if row[4] == '#N/A N/A' or row[4] == '#N/A N/A':
#            price = 'NaN'
#        else:
#            price = (float(row[4]) + float(row[5]))/2
#        q = float(row[10])/100.0
#        print bsimpvolsec( callput, S0, K, r, T, price, q, priceTolerance = 0.01, reportCalls = True )      
    
if __name__=='__main__':
   main()
