from projection_code.projection_code import ProjectionCode
import numpy as np


class NonTrivialProjectionCode(ProjectionCode):

    def __init__(self, num_sources, num_shots, num_projections_in_shot):
        super().__init__(num_sources, num_shots, num_projections_in_shot)
        self._num_sources = num_sources
        self._num_shots = num_shots
        self._num_projection_in_shot = num_sources
        self._projection_code = None
        self._intens_code = None

    def load_full_multiplex_code(self, code: np.ndarray, cfg: dict):
        """This method receives a code of shape (N_shots, N_srcs). Each element in the matrix has a value in the
        range [0,1.0], indicating how many photons this  source will use in a certain shot"""
        if code.shape[0] != self._num_shots:
            raise ValueError('code loaded does not match in number of shots')

        for idx, shot in enumerate(code):
            non_zero_srcs = np.flatnonzero(shot)
            if not idx:  # enters only on 1st iteration
                self._projection_code = np.zeros((self._num_shots, len(non_zero_srcs)), dtype=int)
                self._intens_code = np.zeros((self._num_shots, len(non_zero_srcs)))
            self._projection_code[idx, :] = non_zero_srcs
            self._intens_code[idx, :] = code[idx, non_zero_srcs]

        self._num_projection_in_shot = len(non_zero_srcs)
        cfg['NUM_PROJECTIONS_IN_SHOT'] = int(self._num_projection_in_shot)
        return
