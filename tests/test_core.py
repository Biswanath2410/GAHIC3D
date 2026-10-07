#sanity tests: python -m pytest tests/  (or python tests/test_core.py)
import os
import sys
import numpy as np
from scipy.optimize import check_grad
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from GAHIC3D import getTransformation, Objective, SinglePointCrossoverer, findEuclideanDist
from RMSD import compare, rmsd


def random_rotation(rng):
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    q = q * np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return q


def test_kabsch_recovers_rotation_and_scale():
    rng = np.random.default_rng(0)
    Y = np.cumsum(rng.normal(size=(50, 3)), axis=0)
    X = 2.5 * Y.dot(random_rotation(rng).T) + 7.0
    Xa, Ya = getTransformation(X, Y)
    assert rmsd(Xa, Ya) < 1e-8


def test_mirror_needs_reflection():
    rng = np.random.default_rng(1)
    Y = np.cumsum(rng.normal(size=(50, 3)), axis=0)
    X = Y * np.array([1, 1, -1])                     # mirror image
    assert rmsd(*getTransformation(X, Y)) > 0.1
    assert rmsd(*getTransformation(X, Y, reflection=True)) < 1e-8


def test_compare_matches_bins_by_index():
    rng = np.random.default_rng(2)
    Y = np.cumsum(rng.normal(size=(40, 3)), axis=0)
    gi = np.arange(3, 43)
    ei = np.arange(0, 45)
    E = np.zeros((45, 3))
    E[3:43] = Y
    E[[0, 1, 2, 43, 44]] = 100.
    r = compare(ei, E, gi, Y)
    assert r['NBINS'] == 40 and r['RMSD'] < 1e-8 and r['PCoeff'] > 0.999999


def test_gradients():
    rng = np.random.default_rng(3)
    C = rng.poisson(3, size=(25, 25)).astype(float)
    C = np.triu(C, 1)
    C = C + C.T
    for kind in ['stress', 'wstress', 'poisson', 'poissonnz']:
        o = Objective(C, 2.0, kind)
        x = rng.normal(size=75)
        err = check_grad(lambda v: o.f_and_grad(v, 25)[0], lambda v: o.f_and_grad(v, 25)[1], x)
        f = o.f_and_grad(x, 25)[0]
        assert err < 1e-4 * max(1., abs(f)), (kind, err)


def test_error_is_scale_invariant():
    rng = np.random.default_rng(4)
    C = rng.poisson(5, size=(20, 20)).astype(float)
    C = np.triu(C, 1)
    C = C + C.T
    X = rng.normal(size=(1, 20, 3))
    for kind in ['stress', 'wstress', 'poisson']:
        o = Objective(C, 2.0, kind)
        e1 = o.error(findEuclideanDist(X, o.I, o.J))
        e2 = o.error(findEuclideanDist(3.7 * X, o.I, o.J))
        assert np.allclose(e1, e2), kind


def test_true_structure_has_low_poisson_error():
    rng = np.random.default_rng(5)
    Y = np.cumsum(rng.normal(size=(60, 3)), axis=0)
    D = np.sqrt(((Y[:, None] - Y[None]) ** 2).sum(-1))
    lam = np.where(D > 0, D, 1) ** -2.0 * 50
    np.fill_diagonal(lam, 0)
    C = rng.poisson(np.triu(lam, 1)).astype(float)
    C = C + C.T
    o = Objective(C, 2.0, 'poisson')
    e_true = o.error(findEuclideanDist(Y[None], o.I, o.J))[0]
    e_rand = o.error(findEuclideanDist(rng.normal(size=(1, 60, 3)), o.I, o.J))[0]
    assert e_true < e_rand


def test_crossover_keeps_chain_connected():
    rng = np.random.default_rng(6)
    a = np.cumsum(rng.normal(size=(30, 3)), axis=0)
    b = np.cumsum(rng.normal(size=(30, 3)), axis=0) + 50
    c = SinglePointCrossoverer(a, b, rng)
    steps = np.sqrt((np.diff(c, axis=0) ** 2).sum(1))
    assert steps.max() < 10


if __name__ == '__main__':
    for name, f in list(globals().items()):
        if name.startswith('test_'):
            f()
            print('ok', name)
