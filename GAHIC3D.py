#GAHIC3D - 3D chromosome reconstruction from Hi-C with a memetic genetic algorithm
#
#usage: python GAHIC3D.py <matrix.txt> <outprefix> [options]
#   matrix.txt : dense symmetric contact matrix (nbins x nbins). Bins without any contact are removed;
#                zero-count pairs between the remaining bins are kept and used in the likelihood
#   defaults   : Poisson likelihood over all pairs (--fitness poisson), contact-distance exponent 2 (--alpha),
#                population 10 (--memetic-pop), 40 local optimizations (--n-relax), 200 L-BFGS-B iterations
#                each (--local-iter), pivot-mutation probability 0.3 (--pivot-rate), --seed 1
#   other modes: --variant multistart  (independent local optimizations, no genetic search)
#                --variant ga          (generational GA without local optimization; --pop, --gen, --pmut)
#                --fitness poissonnz | wstress | stress   (alternative objectives)
#
#outputs : <outprefix>_struct.txt  (bin  x  y  z)
#          <outprefix>_param.txt   (run settings, best error, time)
#          <outprefix>_trace.txt   (gen  best  worst  average)
import sys
import os
import time
import math
import argparse
import numpy as np
from datetime import timedelta
from scipy.optimize import minimize
from iced import normalization
import warnings
warnings.filterwarnings('ignore')


class Population_Gen:

    def __init__(self, gen, population, error):
        self.size = len(population)
        if self.size == 0:
            raise RuntimeError('Empty population')
        self.pop = population
        self.error = error
        self.generation = gen

    def get_best_fitness(self):
        return self.error.min()

    def get_worst_fitness(self):
        return self.error.max()

    def get_average_fitness(self):
        return self.error.mean()

    def get_best_organism(self):
        return self.pop[np.argmin(self.error)]

    def best_so_far(self, elit, elit_gen, elit_indiv):
        best_i = np.argmin(self.error)
        if self.error[best_i] < elit:
            elit = self.error[best_i]
            elit_gen = self.generation
            elit_indiv = self.pop[best_i].copy()
        return [elit_gen, elit, elit_indiv]


def get_ga_command_line_params():
    parser = argparse.ArgumentParser(description='GAHIC3D: genetic algorithm for 3D chromosome structure')
    parser.add_argument('matrix', type=str, help='dense Hi-C contact matrix')
    parser.add_argument('outprefix', type=str, help='output prefix')
    parser.add_argument('--variant', default='memetic', choices=['memetic', 'ga', 'multistart'])
    parser.add_argument('--alpha', type=float, default=2.0, help='contact-distance exponent (counts ~ d^-alpha)')
    parser.add_argument('--pop', dest='population_size', type=int, default=100, help='population size (--variant ga)')
    parser.add_argument('--gen', dest='generations', type=int, default=1000, help='generations (--variant ga)')
    parser.add_argument('--pmut', dest='mutation_rate', type=float, default=0.01, help='mutation rate (--variant ga)')
    parser.add_argument('--pxover', dest='xover_rate', type=float, default=0.9, help='crossover rate (--variant ga)')
    parser.add_argument('--elite', type=int, default=2, help='elite size (--variant ga)')
    parser.add_argument('--local-iter', dest='local_iter', type=int, default=200, help='L-BFGS iterations per relaxation')
    parser.add_argument('--n-relax', dest='n_relax', type=int, default=40, help='local relaxations (memetic, multistart)')
    parser.add_argument('--memetic-pop', dest='memetic_pop', type=int, default=10, help='population size (memetic)')
    parser.add_argument('--pivot-rate', dest='pivot_rate', type=float, default=0.3, help='pivot-mutation probability (memetic)')
    parser.add_argument('--patience', type=int, default=200, help='stop after this many generations without improvement (--variant ga)')
    parser.add_argument('--fitness', default='poisson', choices=['stress', 'wstress', 'poisson', 'poissonnz'])
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--no-ice', dest='ice', action='store_false')
    parser.add_argument('--raw-bias', dest='raw_bias', action='store_true',
                        help='Poisson fitness on raw counts with ICE biases (instead of ICE-normalised counts)')
    return parser.parse_args()


########################## Hi-C -> wish distances ###########################

def makeSymmetric(mat):
    """Make matrix symmetric from upper triangle"""
    up = np.triu(mat, k=1)
    return up + up.T


def remove_nan_col(hic):
    """Drop bins without any contact (all-zero rows); zero-count pairs between the other bins are kept"""
    hic = np.nan_to_num(hic)
    col_sum = np.sum(hic, axis=1)
    idx = np.where(col_sum > 0)[0]
    return hic[np.ix_(idx, idx)], idx


def contactToDist(contactMat, alpha):
    """Convert contact matrix to wish distances (upper triangle, contacted pairs only).
    Returns pair indices I, J and wish distance W, scaled so that median(W) = 1."""
    I, J = np.triu_indices(len(contactMat), k=1)
    cont = contactMat[I, J]
    keep = cont > 0                              # zero counts have no finite wish distance
    I, J, cont = I[keep], J[keep], cont[keep]
    W = cont ** (-1. / alpha)
    W = W / np.median(W)
    return I, J, W


############################ distances / error ##############################

def findEuclideanDist(Population, I, J, chunk=20):
    """Pairwise distances of the listed pairs for every individual -> (pop, npairs)"""
    psize = Population.shape[0]
    EDist = np.empty((psize, len(I)))
    for s in range(0, psize, chunk):
        P = Population[s:s + chunk]
        diff = P[:, I, :] - P[:, J, :]
        EDist[s:s + chunk] = np.sqrt((diff * diff).sum(axis=2))
    return EDist


################################ population #################################

def create_pop3(population_size, nbins, rng):
    """Random-walk initialisation: each bin is a +/- step from the previous"""
    steps = rng.uniform(0.01, 0.025, size=(population_size, nbins, 3))
    sign = rng.choice([-1., 1.], size=(population_size, nbins, 3))
    pop = np.cumsum(steps * sign, axis=1)
    pop[:, 0, :] = rng.uniform(-0.5, 0.5, size=(population_size, 3))
    pop[:, 1:, :] += pop[:, :1, :]
    return pop


def tournamentSelection(population, rng, tournamentSize=2):
    selection_i = rng.integers(0, population.size)
    for i in range(tournamentSize):
        ii = rng.integers(0, population.size)
        if population.error[ii] < population.error[selection_i]:
            selection_i = ii
    return selection_i


def getTransformation(X, Y, centering=False, scaling=True, reflection=False):
    """Kabsch superposition of X onto Y (adapted from http://nghiaho.com/?page_id=671).
    Returns aligned X, Y (each N x 3)."""
    X = X.copy().T
    Y = Y.copy().T
    centroid_X = X.mean(axis=1, keepdims=True)
    centroid_Y = Y.mean(axis=1, keepdims=True)
    X = X - centroid_X
    Y = Y - centroid_Y
    if scaling:
        scale_X = np.sqrt(np.square(X).sum() / X.shape[1])
        scale_Y = np.sqrt(np.square(Y).sum() / Y.shape[1])
        to_scale = scale_Y if type(scaling) is bool else float(scaling)
        X = X / (scale_X / to_scale)
        Y = Y / (scale_Y / to_scale)
    C = np.dot(X, Y.transpose())
    V, S, Wt = np.linalg.svd(C)
    I = np.identity(3)
    if not reflection and np.linalg.det(C) <= 0:
        I[-1, -1] = -1
    U = Wt.transpose().dot(I).dot(V.transpose())
    X = U.dot(X)
    if not centering:
        X = X + centroid_Y
        Y = Y + centroid_Y
    return X.T, Y.T


def SinglePointCrossoverer(ind_a, ind_b, rng, connect=True):
    """Single point crossover along the chromosome (bins < point from a, rest from b).
    connect=True shifts the b-part so the chain stays connected at the junction."""
    size = ind_a.shape[0]
    crossover_point = rng.integers(1, size)
    new_org = np.empty_like(ind_a)
    new_org[:crossover_point] = ind_a[:crossover_point]
    shift = ind_a[crossover_point - 1] - ind_b[crossover_point - 1] if connect else 0.
    new_org[crossover_point:] = ind_b[crossover_point:] + shift
    return new_org


def RandomPointMutator(frequency_of_mutations, indiv, sigma, rng):
    """Gaussian perturbation of single bins"""
    new_org = indiv.copy()
    hit = rng.random(indiv.shape[0]) < frequency_of_mutations
    if hit.any():
        new_org[hit] += rng.normal(0, sigma, size=(hit.sum(), 3))
    return new_org


def random_rotation(angle_max, rng):
    axis = rng.normal(size=3)
    axis /= np.linalg.norm(axis)
    theta = rng.uniform(-angle_max, angle_max)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.identity(3) + math.sin(theta) * K + (1 - math.cos(theta)) * K.dot(K)


def PivotMutator(indiv, angle_max, rng):
    """Polymer pivot move: rotate one arm of the chain around a random bin"""
    new_org = indiv.copy()
    n = indiv.shape[0]
    k = rng.integers(1, n - 1)
    R = random_rotation(angle_max, rng)
    if rng.random() < 0.5:
        new_org[k + 1:] = (indiv[k + 1:] - indiv[k]).dot(R.T) + indiv[k]
    else:
        new_org[:k] = (indiv[:k] - indiv[k]).dot(R.T) + indiv[k]
    return new_org


def sigma_schedule(gen, generations, start_sigma=0.3, end_sigma=0.02):
    """Linearly annealed mutation step (in units of the mean wish distance of neighbours)"""
    return start_sigma + (end_sigma - start_sigma) * min(1.0, gen / float(generations))


############################ local refinement ###############################

class Objective:
    """Fitness of a structure given the Hi-C data.
    poisson   : Poisson negative log-likelihood, c_ij ~ Poisson(beta * d_ij^-alpha), over all pairs (zero counts included)
    poissonnz : the same likelihood over pairs with non-zero counts only
    wstress   : relative stress, sum((s*d - w)/w)^2, on wish distances w = c^(-1/alpha) of contacted pairs
    stress    : unweighted stress, sum(s*d - w)^2
    error() is scale-free (optimal scale / beta in closed form), so structures of any size are comparable."""

    def __init__(self, contactMat, alpha, kind='wstress', bias=None):
        self.kind = kind
        self.alpha = alpha
        self.logw = 0.
        if kind in ('poisson', 'poissonnz'):
            self.I, self.J = np.triu_indices(len(contactMat), k=1)
            self.c = contactMat[self.I, self.J]
            if kind == 'poissonnz':                    # ablation: zero counts ignored (as PASTIS)
                nz = self.c > 0
                self.I, self.J, self.c = self.I[nz], self.J[nz], self.c[nz]
            kind = self.kind = 'poisson'
            if bias is not None:
                # raw counts with per-bin biases: lambda_ij = beta * b_i * b_j * d_ij^-alpha (as PASTIS)
                b = np.asarray(bias, dtype=float).flatten()
                self.logw = np.log(b[self.I]) + np.log(b[self.J])
            self.w_p = np.exp(self.logw) if bias is not None else 1.
            self.c = self.c / self.c[self.c > 0].mean()
            self.clogc = np.where(self.c > 0, self.c * np.log(np.where(self.c > 0, self.c, 1)), 0).sum()
            self.W = None
        else:
            self.I, self.J, self.W = contactToDist(contactMat, alpha)
            self.w = np.ones_like(self.W) if kind == 'stress' else 1. / np.maximum(self.W, 1e-6) ** 2

    def scale_to(self, EDist):
        """Optimal multiplicative scale of each structure"""
        if self.kind == 'poisson':
            # beta fixed to 1: pick s so that sum(lambda) = sum(c)
            return (self.c.sum() / (self.w_p * EDist ** -self.alpha).sum(axis=1)) ** (-1. / self.alpha)
        return (self.w * EDist * self.W).sum(axis=1) / ((self.w * EDist * EDist).sum(axis=1) + 1e-12)

    def error(self, EDist):
        EDist = np.maximum(EDist, 1e-9)
        if self.kind == 'poisson':
            lam = self.w_p * (EDist ** -self.alpha)
            beta = self.c.sum() / lam.sum(axis=1)
            loglam = np.log(beta)[:, None] + self.logw - self.alpha * np.log(EDist)
            # deviance per unit count
            dev = 2 * (self.clogc - (self.c * loglam).sum(axis=1))        # sum(lam)=sum(c) at optimal beta
            return dev / self.c.sum()
        s = self.scale_to(EDist)
        r = s[:, None] * EDist - self.W
        return np.sqrt((self.w * r * r).sum(axis=1) / (self.w * self.W * self.W).sum())

    def f_and_grad(self, x, nbins):
        X = x.reshape(nbins, 3)
        diff = X[self.I] - X[self.J]
        d = np.sqrt((diff * diff).sum(axis=1)) + 1e-9
        if self.kind == 'poisson':
            lam = self.w_p * d ** -self.alpha
            f = (lam - self.c * (self.logw - self.alpha * np.log(d))).sum()
            dfdd = -self.alpha * (lam - self.c) / d
        else:
            r = d - self.W
            f = (self.w * r * r).sum()
            dfdd = 2 * self.w * r
        g = (dfdd / d)[:, None] * diff
        G = np.empty_like(X)
        for k in range(3):
            G[:, k] = np.bincount(self.I, g[:, k], minlength=nbins) - np.bincount(self.J, g[:, k], minlength=nbins)
        return f, G.ravel()

    def rescale(self, Population):
        EDist = findEuclideanDist(Population, self.I, self.J)
        s = self.scale_to(EDist)
        return Population * s[:, None, None]


def local_refinement(indiv, obj, maxiter):
    """L-BFGS on the objective, started from a GA individual (memetic step)"""
    nbins = indiv.shape[0]
    x0 = obj.rescale(indiv[None])[0].ravel()
    res = minimize(obj.f_and_grad, x0, args=(nbins,), jac=True, method='L-BFGS-B',
                   options={'maxiter': maxiter})
    return res.x.reshape(nbins, 3)


def evaluate(Population, obj):
    return obj.error(findEuclideanDist(Population, obj.I, obj.J))


########################### statistics of generation #########################

def stat(current_generation, gen):
    return [gen, current_generation.get_best_fitness(), current_generation.get_worst_fitness(),
            current_generation.get_average_fitness()]


def selector(population, params, gen, obj, sigma, rng):
    """One generation: elitism + tournament selection + aligned crossover + mutation"""
    psize, nbins, _ = population.pop.shape
    NewPopulation = np.empty_like(population.pop)
    order = np.argsort(population.error)
    n_elite = params.elite
    for e in range(n_elite):
        NewPopulation[e] = population.pop[order[e]]
    for p in range(n_elite, psize):
        parent1 = population.pop[tournamentSelection(population, rng)]
        parent2 = population.pop[tournamentSelection(population, rng)]
        if rng.random() < params.xover_rate:
            parent2_align, parent1_ref = getTransformation(parent2, parent1, scaling=False, reflection=True)
            NewOrg = SinglePointCrossoverer(parent1, parent2_align, rng)
        else:
            NewOrg = parent1.copy()
        NewOrgMut = RandomPointMutator(params.mutation_rate, NewOrg, sigma, rng)
        if rng.random() < 0.2:
            NewOrgMut = PivotMutator(NewOrgMut, 0.3, rng)
        NewPopulation[p] = NewOrgMut
    NewPopulation[n_elite:] = obj.rescale(NewPopulation[n_elite:])
    error = evaluate(NewPopulation, obj)
    current_generation = Population_Gen(gen, NewPopulation, error)
    return current_generation, stat(current_generation, gen)


def run_GA(contactMat, params, rng, log=None):
    """Returns best structure (nbins x 3), best error, trace"""
    obj = Objective(contactMat, params.alpha, params.fitness, bias=getattr(params, 'bias', None))
    nbins = len(contactMat)
    psize = params.population_size

    pop_arr = create_pop3(psize, nbins, rng)
    pop_arr = obj.rescale(pop_arr)                  # random walks to the scale of the data

    n_local = params.n_relax                     # number of local relaxations (equal budget for both)

    if params.variant == 'multistart':
        # control: random-walk starts, each relaxed with L-BFGS, keep the best (no evolution)
        if n_local > psize:                          # more starts than random walks: draw extra walks
            pop_arr = np.concatenate([pop_arr, obj.rescale(create_pop3(n_local - psize, nbins, rng))])
        best, best_err = None, math.inf
        trace = []                                   # (relaxation, best so far, this start, best so far)
        for k in range(n_local):
            X = local_refinement(pop_arr[k], obj, params.local_iter)
            e = evaluate(X[None], obj)[0]
            if e < best_err:
                best, best_err = X, e
            trace.append([k + 1, best_err, e, best_err])
        best = local_refinement(best, obj, 1000)
        best_err = evaluate(best[None], obj)[0]
        return best, best_err, trace, obj

    if params.variant == 'memetic':
        # Lamarckian (cut-and-splice) GA, cf. Deaven & Ho 1995: every structure in the population
        # is locally relaxed; children = aligned single-point crossover (+ pivot mutation), relaxed,
        # and replace the worst member if better and not a duplicate.
        P = params.memetic_pop
        pop_m = np.array([local_refinement(pop_arr[k], obj, params.local_iter) for k in range(P)])
        err_m = evaluate(pop_m, obj)
        trace = [[0, err_m.min(), err_m.max(), err_m.mean()]]
        for g in range(1, n_local - P + 1):
            current = Population_Gen(g, pop_m, err_m)
            parent1 = pop_m[tournamentSelection(current, rng)]
            parent2 = pop_m[tournamentSelection(current, rng)]
            parent2_align, _ = getTransformation(parent2, parent1, scaling=True, reflection=True)
            child = SinglePointCrossoverer(parent1, parent2_align, rng)
            if rng.random() < params.pivot_rate:
                child = PivotMutator(child, 0.5, rng)
            child = local_refinement(child, obj, params.local_iter)
            e = evaluate(child[None], obj)[0]
            worst = np.argmax(err_m)
            if e < err_m[worst] and np.min(np.abs(err_m - e)) > 1e-7 * abs(e):
                pop_m[worst] = child
                err_m[worst] = e
            trace.append([g, err_m.min(), err_m.max(), err_m.mean()])
            if log is not None and g % 20 == 0:
                log.write(' {0} child. B: {1:.5f}. W: {2:.5f}. A: {3:.5f}\n'.format(*trace[-1]))
                log.flush()
        best = local_refinement(pop_m[np.argmin(err_m)], obj, 1000)
        best_err = evaluate(best[None], obj)[0]
        return best, best_err, trace, obj

    error = evaluate(pop_arr, obj)
    current_generation = Population_Gen(0, pop_arr, error)
    elit_gen, elit, elit_indiv = current_generation.best_so_far(math.inf, 0, None)
    trace = [stat(current_generation, 0)]

    for gen in range(1, params.generations + 1):
        # mutation step relative to the current size of the best structure
        bond = np.median(np.sqrt(((elit_indiv[1:] - elit_indiv[:-1]) ** 2).sum(axis=1)))
        sigma = sigma_schedule(gen, params.generations) * bond
        current_generation, info = selector(current_generation, params, gen, obj, sigma, rng)

        elit_gen, elit, elit_indiv = current_generation.best_so_far(elit, elit_gen, elit_indiv)
        trace.append(info)
        if log is not None and gen % 100 == 0:
            log.write(' {0} gen. B: {1:.5f}. W: {2:.5f}. A: {3:.5f}\n'.format(*info))
            log.flush()
        if gen - elit_gen > params.patience:
            break

    return elit_indiv, elit, trace, obj


#############################################################
def main():

    params = get_ga_command_line_params()
    rng = np.random.default_rng(params.seed)

    rawmat = np.loadtxt(params.matrix)
    rawmat = makeSymmetric(rawmat)
    rawmat, idx = remove_nan_col(rawmat)
    bias = None
    if params.raw_bias:
        _, bias = normalization.ICE_normalization(rawmat.copy(), output_bias=True)
        contactMat = rawmat
    elif params.ice:
        contactMat = normalization.ICE_normalization(rawmat)
    else:
        contactMat = rawmat
    params.bias = bias

    os.makedirs(os.path.dirname(os.path.abspath(params.outprefix)), exist_ok=True)
    pop = {'memetic': params.memetic_pop, 'multistart': params.n_relax}.get(params.variant, params.population_size)
    start = time.time()
    with open(params.outprefix + '_log.txt', 'w') as log:
        log.write('variant: {0}; fitness: {7}; pop: {1}; gen: {2}; pmut: {3}; alpha: {4}; seed: {5}; bins: {6}\n'.format(
            params.variant, pop, params.generations, params.mutation_rate,
            params.alpha, params.seed, len(idx), params.fitness))
        est, elit, trace, _ = run_GA(contactMat, params, rng, log)
    elapsed = time.time() - start

    np.savetxt(params.outprefix + '_struct.txt', np.column_stack((idx, est)), fmt=['%d', '%.6f', '%.6f', '%.6f'])
    np.savetxt(params.outprefix + '_trace.txt', np.array(trace), fmt='%.6g')
    with open(params.outprefix + '_param.txt', 'w') as rparam:
        rparam.write('variant\tfitness\tpop\tgen_run\tpmut\talpha\tseed\tbins\terror\tseconds\n')
        rparam.write('{0}\t{9}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7:.6f}\t{8:.1f}\n'.format(
            params.variant, pop, len(trace) - 1, params.mutation_rate, params.alpha,
            params.seed, len(idx), elit, elapsed, params.fitness))
    print('error: {0:.5f}  time: {1}'.format(elit, timedelta(seconds=elapsed)))


if __name__ == '__main__':
    main()
