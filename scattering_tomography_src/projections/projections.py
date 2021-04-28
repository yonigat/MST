import numpy as np
from numpy import ndarray
from scipy import ndimage


class Projections:
    def __init__(self, active_scorer=2):
        self._I_GT = None
        self._I0_GT = None
        self._projections = None
        # self._error = np.ones((1, shots, rows, cols))
        self._error = None
        self._loss = []
        self._active_scorer = active_scorer
        self._flat_const = None
        self.itr = 0
        self._mask = None

    def set_I_GT(self, I: ndarray, built_phantom=True):
        if built_phantom:
            self._I_GT = I
        else:
            self._I0_GT = I
            self._flat_const = np.sum(np.squeeze(I[:, self._active_scorer, :, :])) / np.count_nonzero(np.squeeze(I[:, self._active_scorer, :, :]))

    def avg_n_cols(self, orig_array, N=3):
        sum_arr = np.zeros((orig_array.shape[0]//N +N - orig_array.shape[0] % N, orig_array.shape[1]))
        rep_arr = np.vstack([orig_array, orig_array[:N-orig_array.shape[0]%N, :]])
        for i in range(N):
            sum_arr += rep_arr[i::N, :]
        sum_arr /= N
        return np.repeat(sum_arr, N, axis=0)[:orig_array.shape[0], :]


    def calc_error(self, angles: list = None, I: np.ndarray = None, n_photons_GT: int = 100, n_photons_fward: int = 100):
        curr_I = np.squeeze(I[:, self._active_scorer, :, :])

        GT_I = np.squeeze(self._I_GT[angles, self._active_scorer, :, :])
        if len(curr_I.shape) < 3:
            curr_I = np.expand_dims(curr_I, axis=-1)
            GT_I = np.expand_dims(GT_I, axis=-1)

        curr_J = np.sum(curr_I[self._mask[angles] == 0]) / np.count_nonzero(curr_I[self._mask[angles] == 0])

        # for i in range(curr_I.shape[1]):
        #     curr_I[:,i,:] = np.mean(curr_I[:,i,:], axis=1, keepdims=True)
        #     GT_I[:, i, :] = np.mean(GT_I[:, i, :], axis=1, keepdims=True)

        # for i in range(curr_I.shape[0]):
        #    curr_I[i, :, :] = self.avg_n_cols(curr_I[i, :, :])
        #    GT_I[i, :, :] = self.avg_n_cols(GT_I[i, :, :])

        if self._mask is not None:
            mask = np.squeeze(self._mask[angles, :, :])
            if len(mask.shape) < 3:
                mask = np.expand_dims(mask, axis=-1)
        else:
            mask = np.ones_like(GT_I)
        # weighted_diff = (n_photons_fward * mask * (curr_I - curr_J*(GT_I/self._flat_const)))
        weighted_diff = n_photons_fward * mask * (curr_I - GT_I)
        weighted_diff[curr_I==0] = 0
        # weighted_diff = mask * (curr_I - GT_I) / self._flat_const
        return weighted_diff

    def calc_loss(self, angles: list = None, I: np.ndarray = None, n_photons_fward: int = 100):
        error = self.calc_error(angles=angles, I=I, n_photons_fward=n_photons_fward)
        loss = 0.5 * (np.sum(self._mask[angles] * error ** 2))
        self.itr = self.itr + 1
        return loss

    def get_active_scorer(self):
        return self._active_scorer

    def calc_cone_only_mask(self):
        """
        :return: creates a binary mask which contains only pixels which are inside the projection cones (i.e. where I0>0)
        in case there is random pixels that are not inside the cone and have nonzero values, we apply a median filter over
        the images.
        """
        self._mask = ndimage.median_filter(self._I0_GT[:, self._active_scorer, :, :] > 0, size=[1, 3, 3])

    def calc_direct_ray_mask(self, eps: float):
        I_GT = self._I_GT[:, self._active_scorer, :, :]
        I0_GT = self._I0_GT[:, self._active_scorer, :, :]
        ratio = np.log(I_GT / I0_GT)
        ratio[np.isnan(ratio)] = 0
        self._mask = np.zeros_like(I_GT)
        self._mask[ratio>eps] = 1
        self._mask[ratio<-eps] = 1
        # self._mask = np.bitwise_or()
        I_GT = self._I_GT[:, self._active_scorer, :, :]
        # self._mask[I_GT>0] *= 1/ np.sqrt(I_GT[I_GT>0])
        I0 = self._I0_GT[:,self._active_scorer, :, :]
        self._flat_const = np.mean(I0[self._mask==0])

