from projection_code.projection_code import ProjectionCode
import numpy as np


def calc_illum_loss(W: np.ndarray, S: np.ndarray, lam: float = 1.0):
    N = len(S)
    loss = (np.linalg.norm(np.matmul(W, S) -1)/N)**2 + lam * 4 * np.matmul(S.T, (1 - S))/N
    return loss


def grad_loss(W: np.ndarray, S: np.ndarray, lam: float = 1.0):
    N=len(S)
    grad = (2*np.matmul(W.T, (np.matmul(W, S) - 1)) +lam*(-8*S + 4))/N
    return grad


def project_onto_constrains(W: np.ndarray, C: float):
    """
    @param W: a vector of size(N,)
    @type W: ndarray
    """
    ones = np.ones_like(W)
    N = len(W)
    W_proj = W + (C - np.sum(W))*ones/N
    return W_proj


def bound_to_constrains(W: np.ndarray, min_val: float = 0.0, max_val: float = 1.0):
    return np.clip(W, min_val, max_val)


def optimize_illum(S_init: np.ndarray, W: np.ndarray, C: float, eta: float = 0.1, lam: float = 0.5, eps: float = 0.001):
    """
    @summary An optimization scheme for mutiplexing code using a qudratic regularization and projection onto constrains.
    @param S_init: an initial multiplexing code, shape of (N_srcs,) where N is number of sources
    @type S_init: numpy array
    @param W: an illumination matrix, describes the weight of each pixel under illumination of a single source, shape of
    (N_cols*N_rows, N_srcs). Each row corresponds to a single detector, each column to a single source.
    @type W: numpy array
    @param C: the maximal sum of W
    @type C: float
    @param eta: step size for optimization scheme
    @type eta: float
    @param lam: controls the quadratic regularzation term
    @type lam: float
    """
    S = []
    loss = []

    S_curr = S_init
    loss.append(calc_illum_loss(W, S_curr))
    S.append(S_curr)
    while True:
        # calc grad
        grad = grad_loss(W, S_curr, lam=lam)
        # project
        S_nxt_unconst = S_curr - grad
        S_nxt_proj = project_onto_constrains(S_nxt_unconst, C)
        d_curr = S_nxt_proj - S_curr
        S_nxt_unbound = S_curr + eta * d_curr
        S_nxt = bound_to_constrains(S_nxt_unbound)
        S_nxt = C * (S_nxt/np.sum(S_nxt))
        S_curr = S_nxt
        # calc loss
        loss.append(calc_illum_loss(W, S_curr))
        S.append(S_curr)
        # check if to stop optimization
        if np.abs(loss[-1] - loss[-2]) < eps:
            break
    return S, loss