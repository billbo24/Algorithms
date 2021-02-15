# -*- coding: utf-8 -*-
"""

FUCK YES
after many hours of coding this function is up and running, and it agrees with my prime sieve function
for 10**8






Created on Thu Dec 24 12:25:22 2020

@author: dingl

take two on implementing a prime counting function.  Found a text book from the two guys who 
made one of the modern algorithms.  They describe what I think is an O(x^2/3) algortithm that's kinda complicated

Recall:
phi(x,a) = number of numbers <= x not divisible by first a primes.  (THIS INCLUDES 1)
P_i(x,a) = number of numbers <= x with i prime factors ALL > p_a

we see phi(x,a) = P_0(x,a) + P_1(x,a) + P_2(x,a) + P_3(x,a) + ...

First clever insight: if we take a = pi(x^1/3) (or more generally pi(x^1/3) <= a <= pi(x^1/2)) this means we're using all primes <= to the cube root of x
Conclusion: P_i(x,a) = 0 for i >= 3.  It's not possible for a number to have 3 or more prime factors
that all excpeed p_a since p_a is the largest prime smaller than the cube root.  

EVERYTHING BELOW HERE IS TAKING A = PI(X^1/3)


P_0(x,a) = 1 kinda by convention (1 is the only number with zero prime factors)
P_1(x,a) = the number of prime numbers greater than p_a
    Note this gives pi(x) - a = P_1(x,a) (all primes less than x equals those greater than a and those less than a)

phi(x,a) = 1 + pi(x)-a + P_2(x,a) ==> pi(x) = phi(x,a) + a - 1 - P_2(x,a)
If we can compute phi(x,a) and P_2(x,a) we're in business. 

First, P_2(x,a).  we can compute this by sieving out all numbers that are a product of two primes > p_a.  

Note added 2/11 upon re reading: it looks like we iterate over all primes p_i where pi(x^1/3) < i <= pi(x^1/2)
Recall we only want to count prime factors strictly greater than p_a, so i would be the first prime larger than a.  
In essence we're sieving with all the primes between the cube root and square root of x.  My initial hunch when approaching
this algorithm was that this sieve would be prohibitively long, but it turns out not to be the case.  WE precompute 
pi(x) for all values up to x^2/3 which means this step is effectively instant (like O(pi(x^1/2))~O(ln(x))).  

Note that sieving up to x^2/3 seems to be what really limits our performance.  if I remember I could reasonably run up to about 
x^12 (sieve up to 100 million) but no higher.  

    pi(x/p) gives all primes which when multiplied by p will have a product <= x
    - pi(p) subtracts off all primes less than p since these will either
        a) already have been counted by a previous iteration of the summation
        b) have a prime <= p_a which we don't want to count
    +1 because we want to count p^2


even trickier now is to compute phi(x,a)
now we have a handy recurrence relation:
    phi(x,a) = phi(x,a-1) - phi(x//p_a,a-1) 
    
We see that we can create a binary tree using this rule.  Naively this tree would have
2^(a-1) leaves at the bottom which is too many.  But some smart SOB said let's institute a truncation rule.  

Truncate the tree at phi(x//n,a) if
i) a == 1 and n <= x^1/3 (i.e. first argument >= x^2/3) (phi(x,1) = floor((x+1)/2)) (ordinary leaf)
ii) n > x^1/3 (special leaf)
    
If I'm being 100% honest I'm not completely grasping the bounds on n here, but that's not terribl;y important.  To actually compute this is 
not a simple recursive function though.  We will first construct the tree and store all of the "leaves", separating them into 
two groups, namely those of type (i) and (ii).  THEN we will use a sieving method to compute the phi(x,a) in reverse order, starting
with a=1, then going to a=2 etc.  

A note about this sieving step: 
    Let's say we have the special leaves [90,1],[85,1],[70,2] and x = 1000.  We initialize a vector with 1's running from 
    1-90 (since we want to effectively count things less than or equal to 90).  First sieve this vector by 2.  
    We can sum this vector from 1-85 to get phi(85,1) and then add 86-90 to that value to get phi(90,1).  Now sieve
    that same vector by 3, and sum it up from 1-70, and that computes phi(70,2).  The one nuance that gets lost here is 


PSEUDOCODE
1) compute all primes up to x^2/3.  Afterward we can pass through the table once and compute pi(n) for 1<= n <= x^2/3 DONE
2) using table from (1),compute a = pi(x^1/3) DONE
3) Evaluate P_2(x,a) by using [1] (O(sqrt(x) terms)) DONE
4) Construct Binary tree for phi(x,a) calculation using recurrence above and truncation rules DONE
5) evaulate the 'ordinary leaves' from (4) using formula  DONE
6) Sort the special leaves by order of increasing b.  Evaluate them using a to-be-discovered sieving method
7) return pi(x) = phi(x,a) - P_2(x,a) + a - 1




"""

import euler_functions as ef
import math 

#putting this here to keep the main function clean
def get_pi(all_primes,upper_limit):
    #given a list of primes and an upper limit on n's, compute
    #pi(x) for 1 <= n <= upper limit
    pi_xs = [0 for i in range(upper_limit+1)] #this will store the value of pi(x) for all x < x^2/3
    
    cur_prime = 0
    cur_func = 0
    
    #this computes pi(n) for all n < x^2/3
    for i in range(2,upper_limit+1):
        #if our n is equal to a prime, then we increment the function value by 1.  
        #if it's not, then we know we simply keep the current function value
        
        
        #things get annoying toward the end.  Basically if we're populating pi(n) and the limit on n isn't prime, 
        #the when it checks for the next prime it's out of range.  Easy fix
        try:
            our_prime = all_primes[cur_prime]
        except:
            our_prime = all_primes[cur_prime-1]
            
        if i == our_prime:
            cur_func += 1
            cur_prime += 1
        pi_xs[i] = cur_func
        
    return pi_xs        


def get_bin(x,a,limit,all_primes,all_leaves): #all primes will be a list such that the k-th element is the k-th prime (no off by one)
    #recall the recurrence relation phi(x,a) = phi(x,a-1) - phi(x//p_a,a-1)
    #truncation rule: if a == 1 and x >= limit^2/3 then we stop and this is an ordinary leave
    #else x < limit ^ 2/3 this is a special leaf
    
    #We also must keep track of whether or not we will be adding or subtracting the leaf
    #I think we can be slick and do this by making the first number positive when it's positive and 
    #negative otherwise
    
    temp_dict = {}
    temp_dict['ordinary leaves'] = [i for i in all_leaves['ordinary leaves']]
    temp_dict['special leaves'] = [i for i in all_leaves['special leaves']]
    
    if a == 1 and int(abs(x)) >= limit:
        egg = all_leaves['ordinary leaves'] + [[x,a]]

        temp_dict['ordinary leaves'] = egg
        return temp_dict
    elif int(abs(x)) < limit:
        egg = all_leaves['special leaves'] + [[x,a]]
        
        temp_dict['special leaves'] = egg
        return temp_dict
    else: #this indicates we must make the recursive call
        left_half = get_bin(x,a-1,limit,all_primes,all_leaves)
        
        #this one is a little weird.  Illustration: when doing 500, one of our dependencies is -phi(71,3)
        #this in turn goes down to -phi(71,2) and phi(14,2).  Because of the floor function and negative numbers
        #if we don't take the absolute value of 71, -71//5 = -15, but that's wrong.  
        #we must take the absolute value for the division, then multiply by the sign again, and then flip that sign
        
        
        right_half = get_bin((abs(x)//all_primes[a])*-1*get_sign(x),a-1,limit,all_primes,all_leaves)

        
        
        temp_dict['ordinary leaves'] = left_half['ordinary leaves'] + right_half['ordinary leaves']
        temp_dict['special leaves'] = left_half['special leaves'] + right_half['special leaves']
        
        return temp_dict



def compute_special_leaves(sorted_special,all_primes):
    #alright now we're getting to the real meat and potatoes of the algorithm
    #This is a guess, but it looks like if we take the maximum value of 
    #of x for the 1's then that should gives us enough room in our sieve
    maxx = 0
    
    #Note: I think the first argument is actually the limit, whereas the 0th argument is the prime index.  
    #I switched them to sort by number of primes
    for i in sorted_special:
        if abs(i[1]) > maxx: #it looks like python might sort these, but just in case I'll actually write the logic
            maxx = abs(i[1])
    
    max_prime_index = sorted_special[-1][0] #I believe the index of the last special leaf tells us the maximum prime we need to sieve with
    Y = [1 for i in range(maxx+1)]
    Y[0] = 0
    
    #I'm going to put the special leaves in a dictionary for ease of access
    special_dict = {}
    
    ans = 0
    for i in range(1,max_prime_index+1):
        #going to store them in a dictionary while I'm here
        special_dict[i] = []
        
        for j in sorted_special:
            if j[0] == i:
                special_dict[i].append([abs(j[1]),get_sign(j[1])]) #my corner cutting ended up costing many hours :(
            elif j[0] > i: #this will save a little time
                break

        p = all_primes[i]
        #Now we're going to sieve our list Y
        Y[p:maxx+1:p] = [0 for i in range(maxx//p)] #we will always have maxx//p elements
        
        #the sieve looks like it's working.  All that remains is to sum Y[1:n+1] For 
        #all of the n in this level of the dict.  This is basically computing
        #phi(K,n) for the n primes we care about and the K values we need
        
        
        run_sum = 0
        sorted_indices = sorted(special_dict[i]) #now they should sort properly 
        #this will keep track of our left end point
        left_pointer = 0
        
        #this should 
        
        #print(sorted_indices)
        for num in sorted_indices:

            #Lets say I give you a series of numbers like 30,40,50,72
            #This loop sums 1-n of our vector (which is a basically random series of 0s and 1s)
            #and takes the sign of the number n into consideration.  
            #It's written so weird because I was trying to 
            #be slick and only pass through the vector Y once.  
            #If I'm not mistaken each iteration of this loop could sum 1-n and multiply by the sign
            #but I was afraid that would be slow
            
            run_sum += sum(Y[left_pointer+1:abs(num[0])+1]) #I think we have to do this to avoid double counting our endpoints
            
            #print(left_pointer,num[0],sum(Y[left_pointer+1:abs(num[0])+1]))
            
            ans += run_sum*num[1] #num[1] is either +1 or -1.  
            #each recursive step subtracts a value, but then when we recurse that, the value that is subtracted 
            #is subtracted, and the net result is a positive
            
            
            left_pointer = abs(num[0]) #this god damn absolute value sign cost me two hours
        
        
    
    return ans
    


def get_sign(x):
    return x//abs(x)


def get_relevant_cube_root(n):
    a = int(n**(1/3))
    if pow(a+1,3) == n:
        return a+1
    return a
    


def pi(x): 
    
    t = get_relevant_cube_root(x)
    
    if t**3 == x: #this means x is a perfect cube and weird shit can happen
        x_two_thirds = t**2
    else:    
        x_two_thirds = int(x**(2/3)) #343 2/3 should be 49, but annoyingly enough it returns. 48.9999999999.  Gonna try adding 1
        
    all_primes = ef.get_primes(x_two_thirds+1)
    
    
    
    
    #get_pi is defined above
    pi_xs = get_pi(all_primes,x_two_thirds)
    
    #alright, now we have pi_xs, where the nth term is equal to pi(n)
    #This initializes the value a we will be using for the rest of the problem
    
    #Floating point errors are causing extremely rare but damaging off by one errors.  343 = 7^3
    #but python thinks 343^(1/3) = 6.99999999999999.  This causes issues because it's the one case
    #where the floor function is wrong.  I will have to build a function to quickly verify that our number 
    #is not a perfect cube.  
    
    a = pi_xs[t] 
    
    #this is a convenience.  I want the nth prime to be on the nth spot
    all_primes = [0] + all_primes
    
    #now we want to compute P_2(x,a)
    #recall the sum above in the green text    
    #for all p with p_a < p <= sqrt(x) sum pi(x/p) - pi(p) + 1
    
    #the upper limit on this loop is a fancy way of sayign "the prime
    #that's closest to the square root of x
    upper_limit = pi_xs[int(x**.5)]
    
    
    P2 = 0
    for i in range(a+1,upper_limit+1): #now this loops will loop through the exact bounds described above
        
        p = all_primes[i]
        
        P2 += pi_xs[x//p] - pi_xs[p] + 1
    
    #alright as far as I can tell I think this compute P2 correctly
    
    #now comes the hard shit.  First we must compute the binary tree for phi(x,a)
    leaves = {}
    leaves['ordinary leaves'] = []
    leaves['special leaves'] = []
    
    lim = (x**.666666)
    tree = get_bin(x,a,lim,all_primes,leaves)
     
    #alright up next is to evaluate the ordinary leaves.  This is the easy step
    phi = 0
    for i in tree['ordinary leaves']:
        phi_component  = get_sign(i[0])*((abs(i[0])+1)//2) #this is a simple formula, won't explain here
        phi +=  phi_component

    
    temp = [[i[1],i[0]] for i in tree['special leaves']]
    
    sorted_special = sorted(temp)
    #print(x,P2,phi,a,tree)
    
    #print("SPACE SPACE SPACE")
    
    
    
    if int(len(sorted_special)) == 0: #this means there are no special leaves and we are done
        return phi - P2 + a - 1
    
    
    else: #this means we must compute the special leaves
        special_summ = compute_special_leaves(sorted_special,all_primes)
        phi += special_summ
        
        #print("sum of special leaves",special_summ)
        
        return phi - P2 + a - 1
    




a=pi(50652) 
b=pi(50653)
    
import time




def compare(L):
    '''
    starter_a = time.time()
    a = ef.get_primes(L+1)
    print("num primes old method",len(a))
    print("time old method",time.time()-starter_a)
    '''
    
    starter_b = time.time()
    b = pi(L)
    print("num primes <= ",L," using new method",b)
    print("time new method",time.time()-starter_b)

   
