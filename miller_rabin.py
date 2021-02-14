# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 23:33:59 2021

@author: dingl

implement miller rabin
"""


def is_prime(n):
    #implementation of miller rabin
    #based on the idea that if n is prime, pow(a,n-1,n)=1 and sqrt(1) = -/+ 1 mod (n)
    #we test these two facts with several bases looking for one that it
    #doesn't hold true.  If we never find one, that means n
    #is likely prime, but if we do find one it is
    #certainly composite
    if n%2 == 0 and n !=2:
        return False
    bases = [2,3,5,7,11,13,17] #if n < 341,550,071,728,321, it is enough to test a = 2, 3, 5, 7, 11, 13, and 17 
    if n in bases:
        return True
    q = n-1
    
    s=0
    while q%2 == 0:
        q = q//2
        s += 1
    
    #n-1 is now factored as 2^s*q
    for b in bases:
        start_num = pow(b,q,n) #compute b^q mod n
        if start_num == 1: #if this is one then all sqrts are one so we don't need to check them
            continue
        
        breaker = 0 #this is how we will tell the loop to skip to the next base
        sqrts = [start_num]
        for i in range(1,s+1):
            print(start_num)
            start_num = (start_num*start_num)%n
            if start_num == n-1: #I believe this means we have -1 as a square root and a^n-1=1, so we don't bother with anything else and move on to next iteration
                breaker = 1
                break
            sqrts.append(start_num)
        
        if breaker == 1:
            continue #go on to the next base
        
        #now we need to iterate backwards over values of s
        flip = sqrts[::-1] #I know we don't need to flip around, but even if n-1 has like 2^100 as a factor this won't really waste time
        for i in flip:
            if (i != 1 and i != n-1): #this means we've found a nontrivial square root of one so our n is definitely composite
                print(flip,b,i)
                return False
    return True #this means all of our bases passed the test, so n islikely prime