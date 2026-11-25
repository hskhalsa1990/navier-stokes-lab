import numpy as np
from navier_stokes import poiseuille, blasius, taylor_green_spectral

def test_poiseuille():
    _, exact, fd = poiseuille()
    assert np.max(np.abs(exact-fd)) < 1e-10

def test_blasius():
    assert abs(blasius()["fpp0"]-0.3320573) < 1e-5

def test_taylor_green_decay():
    t,e,d=taylor_green_spectral(T=2,nsteps=20)
    assert len(t)==len(e)==len(d)
    assert np.all(np.diff(e)<0) and np.all(e>0)


def test_fourth_order_derivative_formula():
    x = np.linspace(0, 2*np.pi, 101)
    h = x[1]-x[0]
    f = np.sin(x)
    d = (f[:-4]-8*f[1:-3]+8*f[3:-1]-f[4:])/(12*h)
    assert np.max(np.abs(d-np.cos(x[2:-2]))) < 1e-5


def test_cavity_cfd_returns_finite_divergence_free_fields():
    from cavity_cfd import simulate_cavity
    result = simulate_cavity(n=15, reynolds=100, steps=20, dt=0.0015)
    assert result["u"].shape == (15, 15)
    assert np.isfinite(result["speed"]).all()
    assert np.max(np.abs(result["divergence"][1:-1, 1:-1])) < 1e-8
