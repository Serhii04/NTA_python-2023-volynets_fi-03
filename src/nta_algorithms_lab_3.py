import math
import random
import numpy as np
import itertools
import collections

import nta_algorithms_lab_1 as lab_1
import nta_algorithms_lab_2 as lab_2
import my_timer

# *********************************************
#              Project const values
# *********************************************

__PRIMES__ = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83,
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

# *********************************************
#              Helpful algorithms
# *********************************************

def reverse(a: int, M: int) -> int:
    if not isinstance(a, int) or not isinstance(M, int):
        raise ValueError("Error: only integer values are allowed")
    
    while a < 0:
        a += M

    if a >= M:
        a = a % M

    if a == 0:
        raise ValueError("Error: zero has no reverse element")

    if a == 1:
        return 1

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
    
    # print(f"q_vals: {q_vals}")
    # print(f"u_vals: {u_vals}")
    # print(f"v_vals: {v_vals}")
    
    if v_vals[-1] >= 0:
        return v_vals[-1]
    
    return v_vals[-1] + M

def _gaus_forward(A: np.ndarray, b: np.ndarray, ord: int) -> np.ndarray:
    A_cur = np.array(A, copy=True)
    b_cur = np.array(b, copy=True)
    
    n = len(A_cur)
    m = len(A_cur[0])
    if m > n:
        raise ValueError("Bad matrix size: {n}x{m}, where should be m <= n")

    for j in range(m):
        # if diagonal element is zero
        if A_cur[j][j] == 0:
            big = 0
            k_row = j
        
            for k in range(j + 1, n):
                if abs(A_cur[k][j]) > big:
                    big = abs(A_cur[k][j])
                    k_row = k
            
            for l in range(j, m):
                A_cur[j][l], A_cur[k_row][l] = A_cur[k_row][l], A_cur[j][l]

            b_cur[j], b_cur[k_row] = b_cur[k_row], b_cur[j]

        pivot = A_cur[j][j]

        # error case
        if pivot == 0:
            print_matrix(A=A_cur, text="Singular matrix:")
            raise ValueError("Given matrix is singular")

        # main part
        for i in range(j + 1, n):
            d = math.gcd(A_cur[i][j], A_cur[j][j])
            mult_upper = int(A_cur[i][j] / d)
            mult_below = int(A_cur[j][j] / d)

            for l in range(j, m):
                A_cur[i][l] = (mult_below * A_cur[i][l] - mult_upper * A_cur[j][l]) % ord

            b_cur[i] = (mult_below * b_cur[i] - mult_upper * b_cur[j]) % ord
        
    return A_cur, b_cur

def _gaus_backward(A: np.ndarray, b: np.ndarray, ord: int) -> np.ndarray:
    n = len(A)
    m = len(A[0])
    X = np.zeros((m, 1))

    for i in range(m-1, -1, -1):
        sum = 0

        for j in range(i+1, m):
            sum = sum + X[j] * A[i][j]
        
        sum = (b[i] - sum) % ord
        d = math.gcd(int(A[i][i]), ord)

        if sum % d != 0:
            raise ValueError("AHTUNG!!!")

        try:
            X[i] = (pow(int(A[i][i] / d), -1, int(ord / d)) * (sum/d)) % ord
        except ValueError as e:
            print(f"a = {int(A[i][i])}, ord = {ord}")
            print(f"pow({int(A[i][i] / d)}, -1, {ord / d})")
            print(e)

    # rez_X = list()
    # for i in range(d):
    #     rez_X.append(X + d*i)

    return X

def gaus(A: np.ndarray, b: np.ndarray, ord: int) -> np.ndarray:
    print_matrix(A=A, rez=b, text="gaus:")

    A_c, b_c = _gaus_forward(A=A, b=b, ord=ord)
    print_matrix(A=A_c, rez=b_c, text="gaus 0.5:")

    X = _gaus_backward(A=A_c, b=b_c, ord=ord)
    # print_matrix(A=X, text="gaus rez:")

    return X

def norm_mod(a: int, m: int):
    pass

# *********************************************
#              Print functions
# *********************************************

def print_matrix(A: np.ndarray, rez: np.ndarray=None, text: str=None) -> None:
    if text is not None:
        print(text)
    
    if rez is not None:
        for A_i, b_i in zip(A, rez):
            print(f"{A_i} = {b_i}")
    else:
        for A_i in A:
            print(f"{A_i}")



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

    for p in __PRIMES__:
        if(p >= max_val):
            return base, base_r
        
        base.append(p);
        base_r[p] = len(base) - 1

    print("Warning: algorithm might want biger prime numbers")
    
    return base, base_r

# Second step of index_calculus
def create_equations(alpha: int, beta: int, n: int, base: list, base_r: dict):
    equations = list()
    b_values = list()
    
    expected_len = len(base) + 15
    # k = 1
    while len(equations) < expected_len:
        k = random.randint(0, n-1)

        a = pow(alpha, k, n + 1)
        # print(f">>> {a} = pow({alpha}, {k}, {n + 1})")
        canon_a = None
        try:
            canon_a = lab_1.get_canon_number_composition_silent(a)
        except ZeroDivisionError as e:
            print(e)

        # print(f"{k}) {a} = {canon_a}")
        if canon_a is None:
            continue

        equation = [0 for i in range(len(base))]
        is_smooth = True
        for a_i in canon_a:
            if a_i not in base:
                is_smooth = False
                break

            equation[base_r[a_i]] += 1
        
        if is_smooth:
            equations.append(equation)
            b_values.append(k)
            # print(f"canon_a = {canon_a}")
            # print(f"eq = {equation}")

        # k += 1
    
    return equations, b_values

# Third step of index_calculus
def solve_equations(n: int, base: list, base_r: dict, equations: list, b_values: list) -> list:
    # rez = gaus(A=equations, b=b_values, ord=n)
    # return rez

    canon_n = lab_1.get_canon_number_composition_silent(n=n)

    n_primes = collections.defaultdict(lambda: 1)
    for n_i in canon_n:
        n_primes[n_i] *= n_i
    
    print(f"keys = {n_primes.keys()}")
    print(f"values = {n_primes.values()}")

    partial_moduls = n_primes.values()
    partial_rez = list()
    for cur_mod in partial_moduls:
        partial_rez.append(gaus(A=equations, b=b_values, ord=cur_mod))

    rez = lab_2.Chinese_remainder_theorem(a_list=partial_rez, mod_list=partial_moduls)

    # TODO: write adecwatte algorithmo

    return rez

# Forth step of index_calculus
def find_log(alpha: int, beta: int, n: int, base: list, base_r: dict, logs_values: list) -> int:
    l = 0
    while True:
        is_smooth = True
        a = beta * pow(alpha, l, n + 1) % (n + 1)
        print(f">>> {a} = pow({alpha}, {l}, {n + 1})")
        
        canon_a = lab_1.get_canon_number_composition_silent(a)
        print(f"{l}) {a} = {canon_a}")
        
        if canon_a is None:
            is_smooth = False
        else:
            for a_i in canon_a:
                if a_i not in base:
                    is_smooth = False
            
        if is_smooth:
            rez_sum = -l
            for a_i in canon_a:
                rez_sum = (rez_sum + logs_values[base_r[a_i]]) % n
            
            return rez_sum
    
        l += 1

def index_calculus(alpha: int, beta: int, n: int) -> int:
    print("First step")
    base, base_r = get_factor_base(n=n)
    # print(base)
    # print(base_r)

    print("Second step")
    equations, b_values = create_equations(alpha=alpha, beta=beta, n=n, base=base, base_r=base_r)
    # for eq in equations:
    #     print(eq)
    # print("")

    print("Third step")
    logs_values = solve_equations(n=n, base=base, base_r=base_r, equations=equations, b_values=b_values)
    print(f"logs_values = \n{logs_values}")

    print("Forth step")
    X = find_log(alpha=alpha, beta=beta, n=n, base=base, base_r=base_r, logs_values=logs_values)
    print(f"X = {X}")

    return X

# **********************************************
#                   Example
# **********************************************

def main():
    alpha = 10
    beta = 17
    p = 47

    # alpha = 13 * 5519 * 23
    # p = 3 * 11 * 5521
    # beta = pow(alpha, 189, p)

    alpha = 304
    beta = 615
    p = 977

    x = index_calculus(alpha=alpha, beta=beta, n=p-1)

    print(f"x = {x}")

    return 0

if __name__ == "__main__":
    main()