# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 22:47:46 2020

@author: dingl
"""

import math

def get_primes(n):
    """ Returns  a list of primes < n """
    sieve = [True] * (n//2)
    for i in range(3,int(n**0.5)+1,2):
        if sieve[i//2]:
            sieve[i*i//2::i] = [False] * ((n-i*i-1)//(2*i)+1)
    return [2] + [2*i+1 for i in range(1,n//2) if sieve[i]]




def factor(n):
    #factor a give number
    #we will need all primes up to sqrt(n)
    ans = []
    temp = n
    limit = int(n**.5)
    prims = get_primes(limit+1) #in case our number is a prime squared
    for p in prims:
        while temp%p==0:
            temp /= p
            ans.append(p)
        if temp == 1:
            return ans
    return ans + [int(temp)] #one big prime factor greater than sqrt n is very likely.  There can be only one though.  
        


    

def factor_powers(n): #this will use our factor function above
    my_primes = factor(n)
    if len(my_primes) == 1: #this is if our number is prime
        return [[my_primes[0],1]]
    ans = []
    
    num_primes = int(len(my_primes))
    
    pointer = 0
    stored_prime = 0
    while pointer < num_primes:
        #check to see if we have a "new" prime
        if my_primes[pointer] > stored_prime: #factor(n) returns a non dewcreasing tuple
            r = 1 #this will represent our exponent
            stored_prime = my_primes[pointer]
            while 1:
                pointer += 1 #this could push us over the edge and cause an error
                try:
                    if my_primes[pointer] == stored_prime:
                        r += 1 #indicates we've seen another prime of the same size
                    else:
                        ans.append([stored_prime,r])
                        break #indicates we've seen all primes of the current size
                except:
                    ans.append([stored_prime,r])
                    pointer+=1
                    break
                
    return ans



def egcd(a, b):
    #the first item in the returned tuple is the gcd
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)



def CRT(a,b): #a is the tuple of moduli and b are the remainders
    if len(a)==1:#stupid but we'll do it
        return b
    if len(a) == 2:
        ans = egcd(a[0],a[1])
        t = ans[1]*a[0]*b[1]+ans[2]*a[1]*b[0]
        return t%(a[0]*a[1])
    else:
        #simply compute the solution to the first two,m rinse and repeat
        small_a = a[2:]
        small_b = b[2:]
        
        small_ans = CRT(a[:2],b[:2])
        
        new_a = [a[0]*a[1]] + small_a
        new_b = [small_ans] + small_b
        
        egg = CRT(new_a,new_b)
        
        return egg



def prod(n):
    ans = 1
    for i in n:
        ans *= i
    return ans


def nCr(n,r):
    nu_r = r
    if r > n//2:
        nu_r = n - r
    ans = 1
    if r == 0: #n choose 0 is 1
        return ans
    for i in range(1,nu_r+1):
        ans *= (n-i+1) #this effectively decrements the ((n-r)!) factorial by one term
        ans = ans // i #increments the other term
    return ans


def fibo(n):
    ans = 1
    prev_ans = 0
    temp = 0
    for i in range(2,n+1):
       temp = ans
       ans += prev_ans
       prev_ans = temp
       
    return ans



def is_prime(n):
    #implementation of miller rabin
    #based on the idea that if n is prime, pow(a,n-1,n)=1 and sqrt(1) = -/+ 1 mod (n)
    #we test these two facts with several bases looking for one that it
    #doesn't hold true.  If we never find one, that means n
    #is likely prime, but if we do find one it is
    #certainly composite
    if n==1:
        return False
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
                return False
    return True #this means all of our bases passed the test, so n islikely prime






def farey(n): #returns all reduced fractions with denominator at most n between 0 and 1
    #interesting function.  I'll skip the derivation, but it's not too complicated. Relies
    #on the interesting property that if you have to fractions a/b, c/d, then
    #a+c/(c+d) is between those two.  (You can kinda do this process in reverse
    #to start with a/b, and (a+c/b+d) to get the next term in the sequence))
    #I'll admit this mediant function is odd and warrants further explanation
    a = 0
    b = 1
    c = 1
    d = n #1/n is always the next smallest after 0
    
    ans = [[a,b],[c,d]]
    
    while 1:
        first_term = (n + b)//d
        p = first_term*c - a
        q = first_term*d - b
        
        ans.append([p,q])
        
        a,b,c,d = c,d,p,q
        
        if p // q == 1:
            break
    
    return ans


def sign(x):
    if x == 0:
        return 0
    return int(abs(x)//x)

def bisection(f,a,b):
    #returns a root so long as there is at least one between a and b
    
    if f(a) == 0:
        return [a,1]
    if f(b) == 0:
        return [b,1]
    
    if a == 0:
        a = pow(10,-6)
    if b == 0:
        b = pow(10,-6)
    
    
    if sign(f(a))*sign(f(b)) > 0: #indicates they have same sign and algorithm won't work
        return [0,0]
    
    cur_mid = (a+b) / 2
    
    while abs(f(cur_mid))>pow(10,-6): #assuming that this will be a good enough delta.  Shouldn't have loss of significance issues
        
        t = f(cur_mid)
        #print(t)
        if sign(t) == sign(f(a)):
            a = cur_mid
        else:
            b = cur_mid
        cur_mid = (a+b)/2

    return [cur_mid,1]


def fast_factor(n,primes): #fast factoring algorithm that requires precomputation of primes
    if is_prime(n):
        return [[n,1]] #we will be giving the prime followed by the power 
    
    temp = n
    counter = 0
    for i in range(int(len(primes))):
        if n%primes[i] == 0:
            
            while temp%primes[i] == 0:
                temp = temp // primes[i]
                counter += 1
            p = primes[i]    
            break
    
    
    ans = [[p,counter]]
    if temp == 1:
        return ans
    
    leftovers = fast_factor(temp,primes[i+1:])
    return ans + leftovers
    


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
        if checker != 1:
            t = (t*b2)%p
            R = (R*b)%p
        i += 1
        b = b2 #this is how we repeatedly square our num.  Note that 8th roots become 4th roots become square roots become 1-th roots lol
    return R