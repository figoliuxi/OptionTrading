# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 16:46:37 2018

@author: filiu
Option class
"""

class Stock:
    def __init__(self,name,price):
        self.name = name
        self.price = price

class EqOption:
    def __init__(self,name,callPut,underlying,K,T):
        self.name = name
        self.callPut = callPut
        self.underlying = underlying
        self.K=K
        self.T=T
        
        
        
        
