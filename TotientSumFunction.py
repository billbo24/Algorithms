# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 21:33:32 2021

@author: dingl

One of my all time favorite algorithms.  Going to sum the first n totients.  
Naively, you might say "well surely summing the first n totients is 
an O(n) algorithm at least."  But I don't believe it is! Someone very smart
figured out a sublinear way to do it

Here's the gist of it:
The totient function for n (phi(n)) counts the number of numbers x <= n
where gcd(x,n) = 1.  phi(6) = 2, phi(p) = p-1
for any prime p.  Well a second way to think of this is counting (x,n) lattice points
where gcd(x,n) = 1.  Also note that we can count the total number of 
pairs where gcd(x,n) > 1, and then subtract that from (n)

(NB: x <= y to avoid double counting)

"But Billy, how is that helpful?" Observe: if we want to sum the first n
totients, denoted PHI(n), we can count the number of lattice points (x,y)
where gcd(x,y) > 1, and then subtract that from the total number of lattice points 
(which is n*(n+1)//2).  
"But I still fail to see how that's helpful.  Calculating points where the gcd > 1
seems hard.""  

Well check it out: if gcd(x,y) = 1, then gcd(mx,my) = m.  If we want the number of points
where gcd(x,y) = m (note m can be as large as n, the pair (n,n) has gcd n) and x<=y<=n, 
we can find all points (x',y') where  gcd(x',y')=1
and x'<=y'<=floor(n/m) (because then we could multiply the pair (x',y') by m
to get a new pair (x,y) with gcd(x,y) = m).  And now the clever reader will see 
where the speed up comes in.  
floor(n/m) will take the same value for many m values! Consider an example
where we want to calculate PHI(10) (which equals 32).  When we want to find
all pairs (x,y) that have a gcd = 7, we're looking for (x,y) pairs with
x,y <= 1.  There's only 1 pair, (1,1).  Now repeat for 8.  We see that 
floor(10/8) = 1 also!, same with 6, 9, and 10.  So this means for 
all of these numbers once we find all the (x,y) pairs where x,y <= 1 and gcd(x,y)-1
, we can "step up" this pair to pairs with gcd = 6,7,8,9,10

Now to fully illustrate with n = 10.  Let F(n,m) denote "number of lattice points (x,y)
where x<=y<=n with gcd(x,y) = m"

F(10,1) = PHI(10) (what we're trying to solve for)
F(10,2) = PHI(10//2) = PHI(5)
F(10,3) = PHI(10//3) = PHI(3)
F(10,4) = F(10,5) = PHI(2)
F(10,[6,7,8,9,10]) = PHI(1)

F(5,1) = PHI(5)
F(5,2) = PHI(2)
F(5,[3,4,5]) = PHI(1)

PHI(5) = 5*(5+1)//2 - (PHI(2) + 3*PHI(1)) = 15 - (2+3) = 10

PHI(10) = 10*(10+1)//2 - (PHI(5) + PHI(3) + 2*PHI(2) + 5*PHI(1))
= 55 - (10+4+2*2+5) = 55 - (23) = 32!!!

To calculate PHI(10) we only needed to calculate PHI(5),PHI(3),PHI(2),PHI(1)

If you precompute all the possible floor quotients (which I'll omit from here)
and how many times we observe them (for example, we will always have n//2
numbers k for which n//k = 1) then it really becomes a simple exercise in memoization
and recursion

Floor division notes:
n//k = x. n-k < x*k <= n.  Proof: if n-k >= x*k, then n >= (x+1)k
which violates the definition of floor division

if j<k <= sqrt(n)
then n//j > n//k
Proof:  let n//j = nj, n//k = nk.  Suppose nj = nk.  Then j*nj = j*nk < k*nk. 
Because j < k, j <= k-1. So j*nk <= (k-1)*nk < k*nk.  In fact it is 
nk smaller than k*nk.  But wait! Because k <= sqrt(n), then nk >= sqrt(n),
but that means we have
j*nk <= (k-1)*nk < (k-1)*nk+j < k*nk.  Last step because j < k <= nk.  




"""

import TotientSieve as TS

def compute_floor_quotients(n):
    #Given a number n, give the ranges of numbers k
    #Which gives n//k = m for all possible m values
    
    ans = {}
    ans[n] = [1,1] #The only number which gives n is one
    #now if memory serves we can see every number up to sqrt(n)
    #similarly every number up to the sqrt(n) gives a different floor
    #quotient.  (If they didn't, )
    sqrt = int(n**0.5) #taking a chance there's no floating point shit
    #where like 4**.5 = 1.999999999996
    
    #the indices in the dictionary are floor quotients, and the
    #elements is a tuple with the lower and upper bound of k 
    #such that n//k = floor
    for i in range(2,sqrt+1): #proof above for why these are distinct
        ans[n//i] = [i,i]
    
    
    #Note: 12 illustrates a weird problem: sqrt(12) = 3.  
    #12//1 = 1, 12//2=6,12//3=4, 12//[7-12] = 1, 12//[5-6] = 2
    #...but where is 4???? I believe this is because 12 isn't 
    #a "perfect square" in floor world, whereas 10 is.  10//sqrt = sqrt
    #12//sqrt != sqrt
    
    upper_lim = sqrt
    if n // sqrt > sqrt:
        upper_lim += 1 #This ensures we don't miss the case where
        #our floor is equal to the sqrt
    
    for floor in range(1,upper_lim): 
        #Alright for each floor value there is a maximum k and a 
        #minimum k that will give us the appropriate values
        #the maximum is easy enough.  It's simply n // floor.  
        #nf*f <= n, (nf+1)*f = nf*f + f.  This must be larger than n
        #by the definition of nf.  
        #minimum value is just one larger than the next floors maximum!
        
        ans[floor] = [(n//(floor+1))+1,n//floor]
    
    return ans
        

def PHI(n,memo):
    
    #a couple of base cases
    if n == 1:
        return 1,memo
    elif n == 2:
        return 2,memo
    elif n == 3:
        return 4,memo
    
    #Check to see if we've already seen the value
    try:
        ans = memo[n]
        return ans,memo
    except:
        pass
    
    
    #Pre compute necessary floor values
    floors = compute_floor_quotients(n)
    
    #Recall, we want all floors except n itself
    egg = sorted(floors.keys())
    egg.pop()
    
    ans = 0
    
    for i in egg:
        ans += (floors[i][1] - floors[i][0] + 1) * PHI(i,memo)[0]
        #print(i,floors[i],ans,n)
        
    ans = (n*(n+1))//2 - ans
    memo[n] = ans #memoize the value for later reference
    
    return ans,memo






#Now a quick check to make sure that everything is working as intended

for i in range(10,1000):
    x = PHI(i,{})
    y = TS.totients(i)
    print(f"For n = {i} fast way got {x[0]} and slow way got {sum(y)}")

    
    