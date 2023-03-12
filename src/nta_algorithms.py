"""Module for NTA algorithms

All the algorithms that says if number p is prime return True,
if not, return False.
In future updtes that algorithms will return divisor, if number isn't
prime, or 0 if they are prime
"""

import math
import random
import sympy
import itertools
import numpy as np
import os
import decimal

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
    C = 1
    A = n
    for i in range(pow.bit_length()):
        if get_i_th_bit(n=pow, i=i) == 1:
            C = (C * A) % module
        A = (A * A) % module
    
    return C

def modular_pow_l(n: int, pow: int, module: int):
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
    if p == 2:
        return True
    
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
        upper_border (int): [now is unresolved] upper border for calculations
    
    Returns:
        bool: divider if n is exactly not pirme. False in other way.
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
            # print(f">>>m: {m}, n_modul: {n_modul}, n: {n}")
            return m
        
    return False

def _f(x: int, n: int):
    return (x*x + 1) % n

def rho_method_of_Pollard(n: int):
    """Las-Vegas probability algorithm that says if number is fully prime

    Args:
        n (int): number
        # f (function): function from real number
    
    Returns:
        bool: divisor if n is prime, False otherwise.
    """
    while True:
        x = random.randint(a=1, b=n)
        x = _f(x, n)
        y = _f(x, n)
        while x != y:
            x = _f(x, n)
            y = _f(_f(y, n), n)
            # print(f"x: {x}, y: {y}")
            d = math.gcd(x - y, n)
            if d != 1:
                if d == n:
                    return 0
                return d
    
    return False

__PRIME_NUMBERS__ = [
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83,
    89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179,
    181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277,
    281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389,
    397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499,
    503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617,
    619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739,
    743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859,
    863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991,
    997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091,
    1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201,
    1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301,
    1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433,
    1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531,
    1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621,
    1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747,
    1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873,
    1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997,
    1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099,
    2111, 2113, 2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237,
    2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 2333, 2339, 2341,
    2347, 2351, 2357, 2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437, 2441,
    2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539, 2543, 2549, 2551, 2557, 2579, 2591,
    2593, 2609, 2617, 2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 2689, 2693,
    2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741, 2749, 2753, 2767, 2777, 2789, 2791, 2797,
    2801, 2803, 2819, 2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 2909, 2917,
    2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 3001, 3011, 3019, 3023, 3037, 3041, 3049,
    3061, 3067, 3079, 3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, 3187, 3191,
    3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 3259, 3271, 3299, 3301, 3307, 3313, 3319,
    3323, 3329, 3331, 3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, 3433, 3449,
    3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533, 3539, 3541, 3547,
    3557, 3559, 3571, 3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659, 3671,
    3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, 3733, 3739, 3761, 3767, 3769, 3779, 3793,
    3797, 3803, 3821, 3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 3911, 3917,
    3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, 4001, 4003, 4007, 4013, 4019, 4021, 4027,
    4049, 4051, 4057, 4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, 4153, 4157,
    4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, 4241, 4243, 4253, 4259, 4261, 4271, 4273,
    4283, 4289, 4297, 4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409, 4421, 4423,
    4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, 4507, 4513, 4517, 4519, 4523, 4547, 4549,
    4561, 4567, 4583, 4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, 4663, 4673,
    4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, 4759, 4783, 4787, 4789, 4793, 4799, 4801,
    4813, 4817, 4831, 4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, 4943, 4951,
    4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, 5009, 5011, 5021, 5023, 5039, 5051, 5059,
    5077, 5081, 5087, 5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, 5189, 5197,
    5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279, 5281, 5297, 5303, 5309, 5323, 5333, 5347,
    5351, 5381, 5387, 5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, 5449, 5471,
    5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, 5527, 5531, 5557, 5563, 5569, 5573, 5581,
    5591, 5623, 5639, 5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, 5701, 5711,
    5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, 5801, 5807, 5813, 5821, 5827, 5839, 5843,
    5849, 5851, 5857, 5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, 5953, 5981,
    5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, 6067, 6073, 6079, 6089, 6091, 6101, 6113,
    6121, 6131, 6133, 6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, 6229, 6247,
    6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, 6311, 6317, 6323, 6329, 6337, 6343, 6353,
    6359, 6361, 6367, 6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, 6481, 6491,
    6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, 6577, 6581, 6599, 6607, 6619, 6637, 6653,
    6659, 6661, 6673, 6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, 6763, 6779,
    6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, 6841, 6857, 6863, 6869, 6871, 6883, 6899,
]

def get_Brillhart_Morrison_factor_base(n: int, a: float=(1/math.sqrt(2))):
    """Creates Brillhart Morrison factor base

    Args:
        n (int): number that we will use to create base. Shouldnt be divisible by 2
    
    Returns:
        list: list of elements in factor base
    """
    La = math.pow(math.exp(math.pow((math.log2(n) * math.log2(math.log2(n))), 0.5)), a)
    # print(f"La: {La}")

    base = [-1]

    if n % 2 == 1:
        base.append(2)

    for p in __PRIME_NUMBERS__:
        if p >= La:
            return base

        if sympy.legendre_symbol(n, p) == 1:
            base.append(p)

    return base

def get_chain_fraction(n: int, k: int):
    """Creates chain fraction from square root of n

    Args:
        n (int): integer number
        k: size of rezult fraction
    
    Returns:
        list: chain fraction
    """
    v = 1
    alpha = math.sqrt(n)
    a = int(alpha // 1)
    u = a
    # print(f"v: {v}, alpha: {alpha}, a: {a}, u: {u}")

    chain_fraction = [a]
    for i in range(k-1):  # (-1) becuse one had been append
        v = (n - u*u) / (v)
        alpha = (math.sqrt(n) + u)/(v)
        a = int(alpha // 1)
        chain_fraction.append(a)
        u = v*a - u
        # print(f"v: {v}, alpha: {alpha}, a: {a}, u: {u}")
    
    return chain_fraction

def get_B_smooth_list(chain_fraction: list, n: int):
    """Returns b smooth list
    """
    b_prew = 0
    b_cur = 1
    B_smooth_list = list()
    for a in chain_fraction:
        b = (b_cur * a + b_prew) % n
        B_smooth_list.append(b)
        b_prew = b_cur
        b_cur = b
    
    return B_smooth_list

def get_B_smooth_list_square(B_smooth_list, n):
    """Returns b^2 smooth list
    """
    B_list= list(B_smooth_list)
    smooth_list = list(B_list)
    for i in range(len(smooth_list)):
        smooth_list[i] = (smooth_list[i] * smooth_list[i]) % n
    
    half = int(n / 2)
    for i in range(len(smooth_list)):
        if smooth_list[i] > half:
            smooth_list[i] = smooth_list[i] - n
    
    return smooth_list

def convert_B_smooth_list_square_to_vector_list(base: list, B_smooth_list_square: list):
    """convert every number of B_smooth_list_square into vector of coeficients in base

    Args:
        base (list): vector of coeficients
        B_smooth_list_square (list): list of numbers
    
    Returns:
        list of lists: matrix of coeficients
    """
    B_list = list(B_smooth_list_square)
    coef_list = list()
    for i in range(len(B_list)):
        coef_list.append([0 for i in range(len(base))])
        for j in range(len(base)):
            if base[j] != -1:
                while B_list[i] % base[j] == 0:
                    B_list[i] = B_list[i] / base[j]
                    coef_list[i][j] = (coef_list[i][j] + 1) % 2
            else:
                if B_list[i] < 0:
                    B_list[i] = -B_list[i]
                    coef_list[i][j] = 1

        if B_list[i] != 1:  # in case if b_i^2 (mod n) can't be decomposed
            coef_list[i] = None

    return coef_list

def vector_sum_is_null(vectors: list):
    for c in range(len(vectors[0])):
        column_value = 0
        for l in range(len(vectors)):
            column_value = (column_value + vectors[l][c])

        if column_value % 2 == 1:
            return False

    return True

def smolarise_matrix(matrix):
    rez_matrix = [list() for i in range(len(matrix))]
    for c in range(len(matrix[0])):
        column_is_nul = True
        for l in range(len(matrix)):
            if matrix[l][c] != 0:
                column_is_nul = False
        
        if not column_is_nul:
            for l in range(len(matrix)):
                rez_matrix[l].append(matrix[l][c])
    
    return rez_matrix

def transpose_matrix(matrix):
    t_matrix = [list() for old_c in range(len(matrix[0]))]

    for old_l in range(len(matrix)):
        for old_c in range(len(matrix[0])):
            t_matrix[old_c].append(matrix[old_l][old_c])
    
    return t_matrix

def Gaus_matrix(matrix):
    """Returns:
        matrix: matrix with zeros below diagonal line 
    """
    temp_matrix = list(matrix)
    places = [i for i in range(len(matrix))]

    for i in range(min(len(temp_matrix), len(temp_matrix[0]))):  # i - is column id, but will be used as diagonal index
        if temp_matrix[i][i] == 0:
            l = i + 1
            while temp_matrix[i][i] == 0 and l < len(temp_matrix):
                if temp_matrix[l][i] == 1:
                    places[i], places[l] = places[l], places[i]
                    temp_matrix[i], temp_matrix[l] = temp_matrix[l], temp_matrix[i]
                    
                    # for c in range(i, len(temp_matrix[0])):
                    #     temp_matrix[i][c] = (temp_matrix[l][c] + temp_matrix[i][c]) % 2
                
                l += 1
  
        for l in range(i+1, len(temp_matrix)):
            if temp_matrix[l][i] == 1:
                for c in range(i, len(temp_matrix[0])):
                    temp_matrix[l][c] = (temp_matrix[l][c] + temp_matrix[i][c]) % 2

    return temp_matrix, places

def solve_matrix(matrix):
    """Solve linear equation in f(2)
    
    Retruns:
        list: list of lines, with sum{vectors at id in list} = 0
    """
    if len(matrix) == 0:
        return None
    
    # for l in matrix:
    #     print(l)

    s_matrix = smolarise_matrix(matrix=matrix)

    for i in range(1, len(s_matrix)+1):
        for perm in itertools.permutations(range(len(s_matrix)), i):
            for_sum = list()
            for_sum_ids = list()
            sumarise = True
            for i in perm:
                if s_matrix[i] == None:
                    sumarise = False
                for_sum.append(s_matrix[i])
                for_sum_ids.append(i)

            if sumarise:
                if vector_sum_is_null(for_sum):
                    return for_sum_ids

def remove_None_from_matrix(matrix):
    rez_matrix = list()
    
    for line_id in range(len(matrix)):
        if matrix[line_id] is not None:
            rez_matrix.append(matrix[line_id])

    return rez_matrix

def remove_Nones(B_smooth_list: list, B_smooth_list_square: list, S: list):
    short_B_smooth_list = list()
    short_B_smooth_list_square = list()
    short_S = list()
    
    for line_id in range(len(S)):
        if S[line_id] is not None:
            short_S.append(S[line_id])
            short_B_smooth_list.append(B_smooth_list[line_id])
            short_B_smooth_list_square.append(B_smooth_list_square[line_id])

    return short_B_smooth_list, short_B_smooth_list_square, short_S

def Brillhart_Morrison_method(n: int):
    a = 0.001  # 1 / math.sqrt(2)
    while a < 1:
        # print(f"a: {a}")
        base = get_Brillhart_Morrison_factor_base(n=n, a=a)
        # print(f"base: {base}")
        chain_fraction = get_chain_fraction(n=n, k=(len(base) + 1))
        # print(f"chain_fraction: {chain_fraction}")
        B_smooth_list = get_B_smooth_list(n=n, chain_fraction=chain_fraction)
        # print(f"B_smooth_list: {B_smooth_list}")
        B_smooth_list_square = get_B_smooth_list_square(n=n, B_smooth_list=B_smooth_list)
        # print(f"B_smooth_list_square: {B_smooth_list_square}")

        S = convert_B_smooth_list_square_to_vector_list(base=base, B_smooth_list_square=B_smooth_list_square)

        B_smooth_list, B_smooth_list_square, S = remove_Nones(B_smooth_list=B_smooth_list,
                                                              B_smooth_list_square=B_smooth_list_square,
                                                              S=S)
        
        # for l in S:
        #     print(f"\t{l}")

        S_solvation = solve_matrix(matrix=S)

        # print(f"S_solvation: {S_solvation}")
        
        # S_solvation = remove_unnesesary_lines(S_solvation)
        # B_smooth_list, B_smooth_list_square, S_solvation = remove_Nones(B_smooth_list=B_smooth_list,
        #                                                                 B_smooth_list_square=B_smooth_list_square,
        #                                                                 S=S_solvation)

        if S_solvation == None:
            if a < 37:
                a += 0.1
            elif a < 38:
                a = 38
            else:
                a += 1

            continue

        X = 1
        for i in S_solvation:
            X = (X * B_smooth_list[i]) % n

        Y = 1
        for i in S_solvation:
            Y = (Y * B_smooth_list_square[i]) % n
        
        Y = int(math.sqrt(Y))

        first_divisor = math.gcd(X + Y, n)
        second_divisor = math.gcd(X - Y, n)

        if first_divisor != 1 and first_divisor != n:
            return first_divisor, second_divisor

        a += 0.1

    return None, None

def print_prime_numbers(n: int):
    num = 0
    for i in range(3, n):
        i_is_prime = True
        for j in range(2, i):
            if (i % j) == 0:
                # print(f"i: {i}, j: {j}")
                i_is_prime = False
        
        if i_is_prime:
            if num > 13:
                print(" ")
                num = 0

            print(f"{i}, ", end="")
            num += 1
    
    print("")

def brut_forsse(n):
    cur_n = n
    while cur_n != 1:
        for i in range(2, cur_n):
            if cur_n % i == 0:
                print(f"{i}: {cur_n}")
                cur_n = int(cur_n / i)

def get_canon_number_composition(n: int):
    import my_timer

    rez = list()

    timer = my_timer.My_Timer()
    print(f"Start factorisation of number: {n}")

    while n != 1:
        if Soloway_Strassen_test(p=n, k=50):
            rez.append(n)
            print(f"Divisor: {n}, curent time: {timer.now():0.2f}, (spend: {timer.point():0.2f} seconds), Soloway Strassen test")
            n = 1
            return rez
        
        m = method_of_trial_divisions(n=n)
        if m:
            if n % m == 0:
                rez.append(m)
                print(f"Divisor: {m}, curent time: {timer.now():0.2f}, (spend: {timer.point():0.2f} seconds), Trial divisions")
                n = int(decimal.Decimal(n) / m)
                continue

        a = rho_method_of_Pollard(n=n)
        if a:
            rez.append(a)
            print(f"Divisor: {a}, curent time: {timer.now():0.2f}, (spend: {timer.point():0.2f} seconds), Rho method of Pollard")
            n = int(n / a)
            
            if Soloway_Strassen_test(p=n, k=10):
                rez.append(n)
                print(f"Divisor: {n}, curent time: {timer.now():0.2f}, (spend: {timer.point():0.2f} seconds), Soloway Strassen test")
                n = 1
                return rez
            continue
            
        a, b = Brillhart_Morrison_method(n=n)
        if a:
            rez.append(a)
            print(f"Divisor: {a}, curent time: {timer.now():0.2f}, (spend: {timer.point():0.2f} seconds), Brillhart Morrison method")
            n = int(n / a)
            continue

        print("я не можу знайти канонiчний розклад числа :(")
        rez.append(n)
        return rez

__NUMBERS__ = [
    572376191,
    194991569,
    962116741,
    6347083,
    89704823,
    19833303,
    6483919,
    6284849,
    7792637,
    28404193,
    # 
    # 3009182572376191,
    # 1021514194991569,
    # 4000852962116741,
    15196946347083,
    # 499664789704823,
    # 269322119833303,
    # 679321846483919,
    96267366284849,
    # 61333127792637,
    2485021628404193,
]

def main():
    print()
    print(get_canon_number_composition(n=691534156424661573), end="\n\n")
    # print(get_canon_number_composition(n=79120395871928357109483571238610239857103241), end="\n\n")
    # print(get_canon_number_composition(n=823423016272041082880), end="\n\n")
    # print(get_canon_number_composition(n=36954229748537), end="\n\n")
    # print(get_canon_number_composition(n=1449863225586482579), end="\n\n")
    # print(get_canon_number_composition(n=3009182572376191), end="\n\n")
    # print(get_canon_number_composition(n=3*3*3*3*3*3*3*3*3), end="\n\n")
    # print(get_canon_number_composition(n=2*2*2*2), end="\n\n")
    # print(Brillhart_Morrison_method(n=9073))
    # print(Brillhart_Morrison_method(n=1829))

    # import my_timer
    # for n in __NUMBERS__:
    #     timer = my_timer.My_Timer()
    #     print(f"For number {n}: {Brillhart_Morrison_method(n)}: spend: {timer.now()}")


if __name__ == "__main__":
    main()

