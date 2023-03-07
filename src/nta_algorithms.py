"""Module for NTA algorithms

All the algorithms that says if number p is prime return True,
if not, return False.
"""

import math
import random
import os

# print(f"Location: {os.getcwd()}")
# from src import my_timer


def Jacobi_symbol(a, n):
    """Mathmatical funtion to find Jacobi symbol value. It says

    Args:
        a (int): a is natural number.
        n (int): n is natural number that isn't divides by 2.

    Returns:
        int: Jacobi symbol value.
    """

    if n <= 0:
        raise ValueError("ERROR in Jacobi_symbol: n must be natural number.")
    
    if n % 2 == 0:
        raise ValueError("ERROR in Jacobi_symbol: n must be odd.")
    
    rez = 1

    while a != 0:
        while a % 2 == 0:
            if ((n*n - 1) / 8) % 2 == 1: 
                rez *= -1
            a = int(a / 2)

        if a == 1:
            return rez
        
        if a == -1:
            if ((n - 1) / 2) % 2 == 1: 
                rez *= -1
            return rez        

        if ((a - 1) * (n - 1) / 4) % 2 == 1: 
            rez *= -1

        a, n = n % a, a

    return 0

def modular_pow(n: int, pow: int, module: int):
    """Methot that calculate modular power of n
    
    Args:
        n (int): natural number, the base number.
        pow (int): power.
        module (int): module of the power.

    Returns:
        int: modular power
    """
    rez = 1
    for i in range(pow):
        rez *= n

        rez = rez % module

    return rez

def Soloway_Strassen_test(p: int, k: int=10):
    """Solovei-Strassen probability test says if the number
    is really prime.

    Args:
        p (int): p is natural number that we will test if it is
        prime
        k (int): times that algorithm will repeat check (default is 10)

    Returns:
        bool: True if number p is prime, False otherwise
    """
    if p % 2 == 0:
        return False

    for i in range(k):
        x = random.randint(2, p-1)

        x_p_gcd = math.gcd(x, p)
        if x_p_gcd > 1:
            return False
        
        Euler_pseudo_prime = modular_pow(x, int((p-1)/2), p)
        Jacobi_n = Jacobi_symbol(x, p)
        if Jacobi_n == -1:
            Jacobi_n += p
        if Jacobi_n != Euler_pseudo_prime:
            return False
        
    return True

def get_i_th_bit(n: int, i: int):
    """Calculates i-th bit of number n
    
    Args:
        n (int): number
        i (int): index in number, the first (the right one) is on 0 place.
    
    Returns:
        int: i-th bit of number n, start  from 0
    """
    if i < 0:
        raise ValueError("Error: indexes less than 0 aren't resolved at a moment")

    return (n >> i) & 1

def get_sum_ai_prod_ri_mod_m(n: int, m: int):
    """Some formula hard computstions: sum a_i * r_i, from i = 0 to size of n
    
    Args:
        n (int): number from what the sum formula will be computed
        m (int): module for computations
    
    Returns:
        int: rezult of formula
    """
    rez = 0

    r_i = 1
    for i in range(n.bit_length()): 
        rez += (get_i_th_bit(n=n, i=i) * r_i) % m
        r_i = (r_i * 2) % m

    return rez % m

def method_of_trial_divisions(n: int, upper_border: int=47):
    """Try do divide number n with all numbers under upper_border,
    
    Args:
        n (int): number
        upper_border (int): [now is unresolved] upper border for calculations, In recomened 47
    
    Returns:
        bool: False if n is exactly not pirme. True in other way.
    """
    if upper_border == -1:
        upper_border = int(math.sqrt(n))
    elif upper_border < 2:
        raise ValueError(f"Error: upper_border must be greather than 1 but {upper_border} is given")

    if upper_border > int(math.sqrt(n)):
        upper_border = int(math.sqrt(n))

    for m in range(2, upper_border + 1):  # (+1) here because I want uper_border to be in cycle
        n_modul = get_sum_ai_prod_ri_mod_m(n, m)
        if n_modul == 0:
            return False
        
    return True



def main():
    pass


if __name__ == "__main__":
    main()

