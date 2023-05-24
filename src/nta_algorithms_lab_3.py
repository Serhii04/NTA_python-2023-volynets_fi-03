import math
import random

import nta_algorithms_lab_1 as lab_1
import nta_algorithms_lab_2 as lab_2
import my_timer



# *********************************************
#              Project functions
# *********************************************

# First step of index_calculus
def get_factor_base(n: int):
    __C__ = 3.38;

    base = list()
    base_r = dict()

    max_val = __C__ * math.exp(0.5 * math.sqrt(math.log(n) * math.log(math.log(n))));
    print(f"max_val = {max_val}");

    for p in lab_1.__PRIME_NUMBERS__:
        if(p >= max_val):
            return base, base_r
        
        base.append(p);
        base_r[p] = len(base) - 1

    print("Warning: algorithm might want biger prime numbers")
    
    return base, base_r

def index_calculus(alpha: int, beta: int, n: int) -> int:
    base, base_r = get_factor_base(n=n)
    print(base)
    print(base_r)

    return 0

# **********************************************
#                   Example
# **********************************************

def main():
    alpha = 2
    beta = 3
    n = 12

    x = index_calculus(alpha=alpha, beta=beta, n=n)

    print(f"x = {x}")

    return 0

if __name__ == "__main__":
    main()