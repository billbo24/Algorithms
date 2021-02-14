# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 18:10:14 2021

@author: dingl
"""

def totients(n):
    ans = [i for i in range(n+1)]
    for i in range(2,(n//2)+1):
        if ans[i] == i: #this means we've hit a prime number
            ans[i] -= 1 #totient for p is p-1
            for j in range(2*i,n+1,i):
                ans[j] = (ans[j]*(i-1))//i #if 3 divides a number we know every third number is not relatively prime, so we have to reduce by 2/3
    for i in range((n//2)+1,n+1):
        if ans[i] == i:
            ans[i] -= 1
    return ans


totients(100)


print("egg")
            