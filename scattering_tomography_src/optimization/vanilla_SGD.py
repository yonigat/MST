import numpy as np
from numpy import ndarray


def VSGD(sg, x0: ndarray, nIter: int, stepSize: float = 0.004):
    """
    A method to perform vanilla SGD - uses only the gradient of current iteration
    without any information from previous iterations.
    :param sg: a handle of the gradient calculation method
    :param x0: initial guess for the decision variables
    :param nIter: number of iterations to perform
    :param stepSize: scalar step size
    :return:
    """
    # number of decision variables
    nDecVar = x0.shape[0]

    # allocate output
    xMat = np.zeros((nDecVar, nIter + 1))

    # initial guess
    xMat[:, 0] = x0

    m = np.zeros((nDecVar,))
    v = np.zeros((nDecVar,))

    for ind in np.arange(1, nIter+1):
        # calculate gradient w.r.t objective at the current iteration
        sgCurr = sg(xMat[:, ind-1], ind-1)

        # update
        xMat[:, ind] = np.maximum(xMat[:, ind-1] - stepSize * sgCurr, 0)

        # print("iteration: " + str(ind) + " current dens: " + str(xMat[:, ind]))

    return xMat
