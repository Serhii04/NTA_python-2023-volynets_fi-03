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

def test_method_of_trial_divisions():
    assert nalg.method_of_trial_divisions(n=107) ==  True
    assert nalg.method_of_trial_divisions(n=6) ==  False
    assert nalg.method_of_trial_divisions(n=804) ==  False
    assert nalg.method_of_trial_divisions(n=(107*107)) ==  True
    assert nalg.method_of_trial_divisions(n=5) ==  True
    assert nalg.method_of_trial_divisions(n=(31*31)) ==  False
