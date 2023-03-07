"""Module for NTA algorithms
"""

import math
import random

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
        # print(f"a: {a}, n: {n}, rez: {rez}.")  # for manual tests
        while a % 2 == 0:
            if ((n*n - 1) / 8) % 2 == 1: 
                rez *= -1
            a = int(a / 2)
            # print(f"a: {a}, n: {n}, rez: {rez}.")  # for manual tests

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
        module (int): module of the power
    """
    rez = 1
    for i in range(pow):
        rez *= n

        rez = rez % module

    return rez

def Soloway_Strassen_test(p: int, k: int=10):
    """Solovei-Strassen probability test says if the number
    is prime

    Args:
        p (int): p is natural number that we will test if it is
        prime
        k (int): times that algorithm will repeat check (default is 10)

    Returns:
        bool: True if number p isn't prime, False otherwise
    """
    if p % 2 == 0:
        return False

    for i in range(k):
        # print(f"{i}-th try")
        x = random.randint(2, p-1)

        x_p_gcd = math.gcd(x, p)
        # print(f"x: {x}, x_p_gcd: {x_p_gcd}.")
        if x_p_gcd > 1:
            return False
        
        Euler_pseudo_prime = modular_pow(x, int((p-1)/2), p)  #pow(x, int((p-1)/2)) % p
        Jacobi_n = Jacobi_symbol(x, p)
        if Jacobi_n == -1:
            Jacobi_n += p
        # print(f"x: {x}, Euler_pseudo_prime: {Euler_pseudo_prime}, Jacobi_n: {Jacobi_n}.")
        if Jacobi_n != Euler_pseudo_prime:
            return False
        
    return True

def main():
    print(Soloway_Strassen_test(p=15, k=2))


if __name__ == "__main__":
    main()

