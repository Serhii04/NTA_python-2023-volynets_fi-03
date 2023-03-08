import pytest

from src import nta_algorithms as nalg

# ---------------------------------------------------------------------
#                       Solowey Strassen tests
# ---------------------------------------------------------------------


def test_Jacobi_symbol():
    assert nalg.Jacobi_symbol(59, 97) == -1
    assert nalg.Jacobi_symbol(3, 7) == -1
    assert nalg.Jacobi_symbol(78695, 23523) == -1

    assert nalg.Jacobi_symbol(1236, 20003) == 1
    assert nalg.Jacobi_symbol(3545, 23523) == 1

    assert nalg.Jacobi_symbol(35, 235) == 0

    with pytest.raises(ValueError):
        nalg.Jacobi_symbol(35, 2352) == 0

def test_modular_pow():
    assert nalg.modular_pow(4, 24, 127) == 64
    assert nalg.modular_pow(47, 98637, 24731) == 12403

def test_Soloway_Strassen_test():
    assert nalg.Soloway_Strassen_test(p=11) == True
    assert nalg.Soloway_Strassen_test(p=31) == True
    assert nalg.Soloway_Strassen_test(p=97) == True
    # assert nalg.Soloway_Strassen_test(p=69026161) == True  # takes a lot of time


# ---------------------------------------------------------------------
#                       test method_of_trial_divisions
# ---------------------------------------------------------------------

def test_get_i_th_bit():
    assert nalg.get_i_th_bit(n=5, i=0) == 1
    assert nalg.get_i_th_bit(n=5, i=1) == 0
    assert nalg.get_i_th_bit(n=5, i=2) == 1
    assert nalg.get_i_th_bit(n=5, i=3) == 0
    assert nalg.get_i_th_bit(n=5, i=4) == 0
    assert nalg.get_i_th_bit(n=6, i=0) == 0
    assert nalg.get_i_th_bit(n=32, i=0) == 0
    assert nalg.get_i_th_bit(n=32, i=1) == 0
    assert nalg.get_i_th_bit(n=32, i=2) == 0
    assert nalg.get_i_th_bit(n=32, i=3) == 0
    assert nalg.get_i_th_bit(n=32, i=4) == 0
    assert nalg.get_i_th_bit(n=32, i=5) == 1
    assert nalg.get_i_th_bit(n=32, i=6) == 0

def test_get_sum_ai_prod_ri_mod_m():
    assert nalg.get_sum_ai_prod_ri_mod_m(n=5, m=3) == 2
    assert nalg.get_sum_ai_prod_ri_mod_m(n=12, m=3) == 0
    assert nalg.get_sum_ai_prod_ri_mod_m(n=42, m=21) == 0
    assert nalg.get_sum_ai_prod_ri_mod_m(n=32, m=8) == 0
    assert nalg.get_sum_ai_prod_ri_mod_m(n=107, m=53) == 1
    assert nalg.get_sum_ai_prod_ri_mod_m(n=6123, m=2) == 1

def test_method_of_trial_divisions():
    assert nalg.method_of_trial_divisions(n=107) ==  False
    assert bool(nalg.method_of_trial_divisions(n=6)) ==  True
    assert bool(nalg.method_of_trial_divisions(n=804)) ==  True
    assert nalg.method_of_trial_divisions(n=(107*107)) ==  False
    assert nalg.method_of_trial_divisions(n=5) ==  False
    assert bool(nalg.method_of_trial_divisions(n=(31*31))) ==  True

# ---------------------------------------------------------------------
#                       test rho_method_of_Pollard
# ---------------------------------------------------------------------

def test_rho_method_of_Pollard():
    assert nalg.rho_method_of_Pollard(n=(2*3*5*7*11*13)) != 0
    assert nalg.rho_method_of_Pollard(n=(31*31)) != 0
    assert nalg.rho_method_of_Pollard(n=(107)) == 0
    assert nalg.rho_method_of_Pollard(n=(107036579817239)) != 0
    assert nalg.rho_method_of_Pollard(n=(9973)) == 0
    assert nalg.rho_method_of_Pollard(n=(3131)) != 0
    assert nalg.rho_method_of_Pollard(n=(3123)) != 0

# ---------------------------------------------------------------------
#               test get_Brillhart_Morrison_factor_base
# ---------------------------------------------------------------------

def test_get_Brillhart_Morrison_factor_base():
    assert nalg.get_Brillhart_Morrison_factor_base(n=7) == [-1, 2, 3]
    assert nalg.get_Brillhart_Morrison_factor_base(n=203) == [-1, 2, 11, 17]
    assert nalg.get_Brillhart_Morrison_factor_base(n=203, a=0.8) == [-1, 2, 11, 17, 41, 43]
    assert nalg.get_Brillhart_Morrison_factor_base(n=203, a=0.2) == [-1, 2]
    assert nalg.get_Brillhart_Morrison_factor_base(n=203, a=0.3) == [-1, 2]
    assert nalg.get_Brillhart_Morrison_factor_base(n=203, a=0.6) == [-1, 2, 11, 17]

def test_get_chain_fraction():
    assert nalg.get_chain_fraction(n=203, k=3) == [14, 4, 28]
    assert nalg.get_chain_fraction(n=203, k=6) == [14, 4, 28, 4, 28, 4]  # TODO: I want it to be ok
    assert nalg.get_chain_fraction(n=203, k=8) == [14, 4, 28, 4, 28, 4, 28, 4]  # TODO: I want it to be ok
    assert nalg.get_chain_fraction(n=2033, k=8) == [45, 11, 3, 1, 4, 1, 7, 2]  # TODO: I want it to be ok

def test_get_B_smooth_list():
    assert nalg.get_B_smooth_list([14, 4, 28], 203) == [14, 57, 189]
    # assert nalg.get_B_smooth_list([14, 4, 28], 203) == [14, 57, -14]
    # assert nalg.get_B_smooth_list([14, 4, 28], 203) == [14, 57, -14]
    # assert nalg.get_B_smooth_list([14, 4, 28], 203) == [14, 57, -14]

def test_get_B_smooth_list_square():
    assert nalg.get_B_smooth_list_square([14, 57, 189], 203) == [-7, 1, -7]
    # assert nalg.get_B_smooth_list_square([14, 57, -14], 203) == [-7, 1, -7]
    # assert nalg.get_B_smooth_list_square([14, 57, -14], 203) == [-7, 1, -7]
    # assert nalg.get_B_smooth_list_square([14, 57, -14], 203) == [-7, 1, -7]

def test_convert_B_smooth_list_to_vector_list():
    assert nalg.convert_B_smooth_list_square_to_vector_list(base=[-1, 11, 17],
                                                     B_smooth_list_square=[196, 1, 196]) == [
            None,
            [0, 0, 0],
            None,
        ]
    assert nalg.convert_B_smooth_list_square_to_vector_list(base=[-1, 2, 3, 7],
                                                     B_smooth_list_square=[-48, 139, -7, 87, -27]) == [
            [1, 0, 1, 0],
            None,
            [1, 0, 0, 1],
            None,
            [1, 0, 1, 0],
        ]

def test_vector_sum_is_null():
    assert nalg.vector_sum_is_null([[1, 1, 1, 0],
                                    [0, 1, 1, 0],
                                    [1, 1, 0, 0],
                                    [0, 1, 0, 0],]) == True
    assert nalg.vector_sum_is_null([[1, 1, 1, 0],
                                    [0, 1, 1, 0],
                                    [1, 1, 0, 1],
                                    [0, 1, 0, 0],]) == False
    assert nalg.vector_sum_is_null([[1, 1, 1, 1],
                                    [1, 1, 1, 1],
                                    [1, 1, 1, 1],
                                    [1, 1, 1, 1],]) == True
    assert nalg.vector_sum_is_null([[0, 0, 0, 0],
                                    [0, 0, 0, 0],
                                    [0, 0, 0, 0],
                                    [0, 0, 0, 0],]) == True
    assert nalg.vector_sum_is_null([[1, 0, 0, 0],
                                    [0, 1, 0, 0],
                                    [0, 0, 1, 0],
                                    [0, 0, 0, 1],]) == False
    assert nalg.vector_sum_is_null([[1, 0, 0, 0],
                                    [0, 0, 0, 0],
                                    [0, 0, 0, 0],
                                    [0, 0, 0, 0],]) == False

def hard_test_solve_matrix():
    assert nalg.solve_matrix(matrix=[  # solve_matrix 1
        [1, 0, 0, 0, 1, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 1, 0, 0],
        [0, 1, 1, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 1, 0],
        [1, 0, 1, 0, 0, 1, 0, 0, 0],
    ]) == None
    assert nalg.solve_matrix(matrix=[  # solve_matrix 2
        [1, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 1, 0],
        [1, 0, 1, 0, 0, 1, 0, 0, 0],
    ]) == [1]
    assert nalg.solve_matrix(matrix=[  # solve_matrix 3
        [0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 1, 0],
        [1, 0, 1, 0, 0, 1, 0, 0, 0],
    ]) == [0]
    assert nalg.solve_matrix(matrix=[  # solve_matrix 4
        [1, 0, 0, 0, 1, 0, 1, 0, 0],
        [1, 0, 0, 0, 1, 0, 1, 0, 0],
        [0, 1, 1, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 1, 0],
        [1, 0, 1, 0, 0, 1, 0, 0, 0],
    ]) == [0, 1]
    assert nalg.solve_matrix(matrix=[  # solve_matrix 5
        [1, 0, 0, 0, 1, 0, 1, 0, 0],
        [1, 0, 0, 1, 1, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 1, 1, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 1, 0],
        [1, 0, 1, 0, 0, 1, 0, 0, 0],
    ]) == [1, 2, 3]
    assert nalg.solve_matrix(matrix=[  # solve_matrix 6
        [1, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 1],
    ]) == None
    assert nalg.solve_matrix(matrix=[  # solve_matrix 7
        [1, 0, 0, 0, 0],
        [1, 1, 0, 0, 0],
        [0, 1, 1, 0, 0],
        [0, 0, 1, 1, 0],
        [0, 0, 0, 1, 0],
    ]) == [0, 1, 2, 3, 4]

def test_remove_Nones():  # TODO: believe in its corectness
    pass


