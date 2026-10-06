import numpy as np
from dron_simon.plomes import ploma_gaussiana


def test_zero_contra_vent():
    assert ploma_gaussiana(-10.0, 0.0, 10.0) == 0.0


def test_maxim_a_l_eix():
    assert ploma_gaussiana(100.0, 0.0, 10.0) > ploma_gaussiana(100.0, 20.0, 10.0)


def test_positiva():
    assert np.all(ploma_gaussiana(np.linspace(1, 300, 50), 0.0, 10.0) > 0)
