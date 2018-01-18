import crysfipy.const as C
from pytest import approx

def test_R0():
    assert C.R0 == approx(-5.390841372421595e-15)

def test_ion():
    assert C.ion("Ho").Alpha == approx(-0.0022222222222222222)
    assert C.ion("Ho").Beta == approx(-3.330003330003329e-05)
    assert C.ion("Ho").Gamma == approx(-9.951596826260624e-08)

def test_uB():
    assert C.uB == approx(0.6717138840811654)



