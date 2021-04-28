import numpy as np
from numpy import ndarray
from utils.tf_logger import Logger


def ADAM(sg, x0: ndarray, nIter: int, stepSize = 0.004, beta1: float = 0.9, beta2: float = 0.99, epsilon: float =1.5e-8, constrain: bool=False, constrain_func=None, logger: Logger=None, cfg: dict=None):
    """
    :param sg: a handle of the gradient calculation method
    :param x0: initial guess for the decision variables
    :param nIter: number of iterations to perform
    :param stepSize: scalar step size
    :param beta1: exponential decay rate for the 1st moment estimate
    :param beta2: exponential decay rate for the 2nd moment estimate
    :param epsilon:  back-to-numerical-reality addend, default: `sqrt(eps)` of float type
    :return:
    """
    if type(stepSize) is list:
        stepSize = np.array(stepSize)
        stepSize = np.repeat(stepSize, x0.size/len(stepSize))
    # number of decision variables
    nDecVar = x0.shape[0]

    # allocate output
    xMat = np.zeros((nDecVar, nIter + 1))

    # initial guess
    xMat[:, 0] = x0

    m = np.zeros((nDecVar,))
    v = np.zeros((nDecVar,))

    constrain_schedule = np.unique(np.array(np.geomspace(1,nIter, num=int(nIter*0.3)), dtype=int))

    for ind in np.arange(1, nIter+1):
        # calculate gradient w.r.t objective at the current iteration
        sgCurr = sg(xMat[:, ind-1], ind-1)

        # Update biased 1st moment estimate
        m = beta1 * m + (1 - beta1) * sgCurr
        # Update biased 2nd moment estimate
        v = beta2 * v + (1 - beta2) * (sgCurr**2)

        # Compute bias-corrected 1st moment estimate
        mHat = m / (1 - beta1**ind)
        # Compute bias-corrected 2nd moment estimate
        vHat = v / (1 - beta2**ind)


        # update
        xMat[:, ind] = np.maximum(xMat[:, ind-1] - stepSize * mHat/ (np.sqrt(vHat) + epsilon), 0)
        # xMat[:, ind] = np.maximum(xMat[:, ind - 1] - stepSize * sgCurr,    0)

        # if constrain and ind in constrain_schedule:
        if constrain and ind < 5:
            xMat[:, ind] = constrain_func(xMat[:, ind])

        # stepSize*=0.9

        if logger is not None and cfg is not None:
            diff = xMat[:, ind] - xMat[:, ind-1]
            sgCurrmat = sgCurr
            diff = np.reshape(diff,
                              (cfg["NUM_OF_VOXELS_X"], cfg["NUM_OF_VOXELS_Y"], cfg["NUM_OF_Z_SLICES"], len(cfg["ACTIVE_ELEMENTS"])), order="F")
            sgCurrmat = np.reshape(sgCurrmat,
                              (cfg["NUM_OF_VOXELS_X"], cfg["NUM_OF_VOXELS_Y"], cfg["NUM_OF_Z_SLICES"],
                               len(cfg["ACTIVE_ELEMENTS"])), order="F")
            for i in np.arange(diff.shape[-1]):
                curr_element = cfg["ACTIVE_ELEMENTS"][i]
                logger.log_multiple_figures(f'difference between iterations {curr_element}', ind-1, diff[:, :, cfg['SLICES_TO_PLOT'], i],
                                            colorbar=True, dim_permute=(2, 0, 1))
                logger.log_multiple_figures(f'grad {curr_element}', ind - 1, sgCurrmat[:, :, cfg['SLICES_TO_PLOT'], i],
                                            colorbar=True, dim_permute=(2, 0, 1))
    return xMat
