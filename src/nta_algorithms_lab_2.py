import math
import signal
from contextlib import contextmanager
from collections import defaultdict

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

# def SPG_timed(alpha: int, beta: int, n: int, timeout: int=10) -> int:
#     try:
#         with time_limit(timeout):
#             return SPG(alpha=alpha, beta=beta, n=n)
#     except UserWarning as e:
#         print("Timed out!")
    
#     return None

def get_canon_degrees(n: int) -> defaultdict:
    canon_n = lab_1.get_canon_number_composition_silent(n)

    if canon_n[0] == None:
        return None

    canon_dict = defaultdict(int)
    for p in canon_n:
        canon_dict[p] += 1
    
    return canon_dict

def reverse(a: int, M: int) -> int:
    if not isinstance(a, int) or not isinstance(M, int):
        raise ValueError("Error: only integer values are allowed")
    
    while a < 0: a += M

    if a >= M: a = a % M
    if a == 0: raise ValueError("Error: zero has no reverse element")
    if a == 1: return 1

    q_vals = []
    r_1 = M
    r_2 = a
    r_3 = 1
    
    while r_3 != 0:
        q_vals.append(int(r_1 / r_2))
        r_3 = r_1 % r_2
        r_1 = r_2
        r_2 = r_3
    
    if r_1 != 1:
        raise ValueError(f"Error: number must have no comon divisors with module while have: {r_1}")
    
    q_vals.pop()
    
    u_vals = [1, 0]
    v_vals = [0, 1]
    for q in q_vals:
        u_vals.append(u_vals[-2] - u_vals[-1] * q)
        v_vals.append(v_vals[-2] - v_vals[-1] * q)
    
    if v_vals[-1] >= 0:
        return v_vals[-1]
    
    return v_vals[-1] + M

def Chinese_remainder_theorem(a_list: list, mod_list: list) -> int:
    M = 1
    for m in mod_list:
        M *= m

    X = 0
    for m_i, a_i in zip(mod_list, a_list):
        M_i = int(M / m_i)
        N_i = reverse(a=M_i, M=m_i)

        X = (X + a_i * M_i * N_i) % M

    return X


def SPG(alpha: int, beta: int, n: int) -> int:
    # 1st step
    canon_n = get_canon_degrees(n)

    # 2nd step
    pre_table = dict()
    for p_i, i in zip(canon_n, range(len(canon_n))):
        pre_table_p_i = dict()
        for j in range(p_i):
            r = pow(alpha, int(n/p_i)*j, n+1)
            pre_table_p_i[r] = j
        
        pre_table[p_i] = pre_table_p_i

    # 3rd -5th steps
    a_list = list()
    mod_list = list()
    for p_i, i in zip(canon_n, range(len(canon_n))):
        cur_alpha = 1
        alpha_neg = reverse(alpha, n+1)
        x_coefs = list()
        for j in range(canon_n[p_i]):
            temp = pow(beta * cur_alpha, int(n / pow(p_i, j+1)), n+1)
            x_coefs.append(pre_table[p_i][temp])
            cur_alpha *= pow(alpha_neg, x_coefs[-1] * pow(p_i, j), n+1)

        m_i = pow(p_i, canon_n[p_i])

        a_i = 0
        for i, x_c in enumerate(x_coefs):
            a_i += x_c * pow(p_i, i)

        a_list.append(a_i)
        mod_list.append(m_i)

    # 6th step
    X = Chinese_remainder_theorem(a_list=a_list, mod_list=mod_list)


    return X


######################################
# Example
######################################

def brute_force_example():
    alpha=48
    beta=3
    p=237
    x = discrete_logarithm_brute_force_timed(alpha=alpha, beta=beta, p=p)
    
    if x:
        # print(f"{alpha}^{x} = {beta} (mod {p})")
        print(f"x = {x}")
    else:
        print("No solution")

def SPG_example():
    alpha = 999559
    beta = 1264341
    p = 4700383

    alpha = 999559
    beta = 1264341
    p = 4700383
    x = SPG(alpha=alpha, beta=beta, n=p-1)
    
    if x:
        print(f"x = {x}")
        # print(f"{alpha}^{x} = {beta} (mod {p})")
    else:
        print("No solution")

def main():
    # brute_force_example()
    SPG_example()

if __name__ == "__main__":
    main()