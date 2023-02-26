import pytest

from src import nta_algorithms as nalg

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
