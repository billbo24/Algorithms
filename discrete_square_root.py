# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 18:37:50 2021

@author: dingl

my attempt at tonellis hanks discrete square root

"""
        
        

def euler_criterion(n,p):
    return pow(n,(p-1)//2,p)


def get_non_square_mod_p(p):
    for i in range(2,p):
        numm = pow(i,(p-1)//2,p)
        if numm == p - 1:
            return i
    #there will always be a nonsquare residue

def factor_Q_2S(p):
    #we are factyoring p-1 into Q*2^S where Q is odd
    S = 0
    Q = p-1
    while Q % 2 == 0:
        S += 1
        Q = Q // 2
    return Q,S

def discrete_sqrt(n,p):
    a = euler_criterion(n,p)
    if a == p-1:
        return 0 #no square root by eulers criterion
    z = get_non_square_mod_p(p)
    Q,S = factor_Q_2S(p) #factors p-1 into Q*2^S
    
    R = pow(n,(Q+1)//2,p) #this is basically our first guess for the square root of n
    #Note that R^2 = n^(Q+1) = (n^Q)(n) so if n^Q == 1 then we're done.  Otherwise...
    t = pow(n,Q,p) 
    #note that if t != 1 (mod p) then we necessarily have t^(2^(S-1)) = 1 (mod p) because
    # = n^(Q*2^(S-1)) = n^((p-1)//2) = 1 by assumption that n is a quadratic residue
    #also note that z^Q is a (2^(S-1))-th root of -1 by a similar argument
    b = pow(z,Q)
    i = 1 #this basically counts how many times we've interated throuigh our loop.  
    #Note that as we start t is a 2^(S-1)-th root of 1, and we're gonna do stuff
    #to make t' a 2^(S-2)-th root of 1, until eventually t is a 2^0-th root of one ie 1
    
    #now here's the general idea.  We have R^2 = n*t.  Where t is a 2^(S-1) root of 1.  
    #if it's a 2^(S-2)-th root, great we can move on.  If it's not, check this out.  
    #If t raised to the 2^(S-2) isn't 1 then it has to be -1 (since squaring it yields 1)
    #Note this means we have a 2^(S-2)-th root of -1.  Recall that our number z^Q is a 2^S-1
    #of -1.  Well we can square it (z^2Q) and that is now a 2^(S-2) root of -1.  Well then
    #t*(z^(2Q)) is a 2^(S-2)-th root of +1, because both individually are -1.  
    
    #next iteration we can repeat.  Square z^2Q to get z^4Q.  This is a 2^(S-3)-th root of negative one...
    
    while t != 1:
        checker = pow(t,pow(2,S-i-1),p) #this is in effect checking if our 2^(S-1)-th root of 1 is also a 2^(S-2)th root of one
        #this can happen in nontrivial cases.  For example if you repeatedly square and eventurally hit 1,
        #every subsequent square is 1 too, so we can have multiple n-th powers be one
        b2 = (b**2)%p #we will be multiplying t by b**2 and R
        #b will always be a 2^(S-i-1)-th root of -1
        
        #if checker == 1 we don't actually need to do anything extra; however, if...
        print(i,b,b2,R,t,checker)
        if checker != 1:
            t = (t*b2)%p
            R = (R*b)%p
        i += 1
        b = b2 #this is how we repeatedly square our num.  Note that 8th roots become 4th roots become square roots become 1-th roots lol
    print(R)
    return R




