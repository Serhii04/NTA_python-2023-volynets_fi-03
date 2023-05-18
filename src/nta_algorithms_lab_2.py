import math
import signal
from contextlib import contextmanager

import nta_algorithms_lab_1 as lab_1


######################################
# help functions to limit time
######################################

@contextmanager
def time_limit(seconds):
    def signal_handler(signum, frame):
        raise UserWarning("Timed out!")
    
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    
    try:
        yield
    finally:
        signal.alarm(0)


######################################
# Brute force function implementation
######################################

def discrete_logarithm_brute_force_timed(alpha: int, beta: int, p: int, timeout: int=10) -> int:
    try:
        with time_limit(timeout):
            return discrete_logarithm_brute_force(alpha=alpha, beta=beta, p=p)
    except UserWarning as e:
        print("Timed out!")
    
    return None

def discrete_logarithm_brute_force(alpha: int, beta: int, p: int) -> int:
    for i in range(p):
        if pow(alpha, i, p) == beta:
            return i
    
    return None

######################################
# S-P-G algorithm
######################################


######################################
# Example
######################################

def main():
    alpha=48
    beta=3
    p=237
    x = discrete_logarithm_brute_force_timed(alpha=alpha, beta=beta, p=p)
    
    if x:
        # print(f"{alpha}^{x} = {beta} (mod {p})")
        print(f"x = {x}")
    else:
        print("No solution")

if __name__ == "__main__":
    main()