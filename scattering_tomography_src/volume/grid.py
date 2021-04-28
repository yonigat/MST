
class Grid:
    def __init__(self, shape: tuple, voxel_size: tuple):
        """
        Initializes the Grid object
        :param shape: tuple with the grid shape (x_shape: int, y_shape: int, z_shape: int)
        :param voxel_size: tuple with the size of each voxel in cm (x_size: float, y_size: float, z_size: float)
        """
        self._shape = shape
        self._voxel_size = voxel_size

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self._shape == other._shape) and (self._voxel_size == other._voxel_size)
        else:
            return False

    def get_shape(self):
        return self._shape

    def get_voxel_size(self):
        return self._voxel_size

    def get_num_voxels(self):
        return  self._shape[0]*self._shape[1]*self._shape[2]
