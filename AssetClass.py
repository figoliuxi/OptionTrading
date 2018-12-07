# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 16:46:37 2018

@author: filiu
Option class
"""

class Stock:
    def __init__(self,name,price,dividend):
        self.name = name
        self.price = price
        self.dividend = dividend

class EqOption:
    def __init__(self,name,callPut,underlying,K,T,price,vol):
        self.name = name
        self.callPut = callPut
        self.underlying = underlying
        self.K=K
        self.T=T
        self.price = price
        self.vol= vol
        
    def setVol(self,vol):
        if self.vol == None:
            self.vol = vol
    
    def setPrice(self,price):
        if self.price == None:
            self.price = price
        
        
        
        
        
