import numpy as np
import scipy.io


class TensorTomographyReconstruction:
    """ """

    def __init__(self, path_to_mat_file):
        self.path_to_mat_file = path_to_mat_file
        self.mat = scipy.io.loadmat(self.path_to_mat_file)
        self.mask = self.mat["s"]["mask3D"][0, 0]
        self.theta = self.mat["s"]["theta"][0, 0][0, 0][0]
        self.phi = self.mat["s"]["phi"][0, 0][0, 0][0]
        self.a0 = self.mat["s"]["a"][0, 0][:, 0][0][0]
        self.a2 = self.mat["s"]["a"][0, 0][:, 1][0][0]
        self.a4 = self.mat["s"]["a"][0, 0][:, 2][0][0]
        self.a6 = self.mat["s"]["a"][0, 0][:, 3][0][0]

        self.shape = self.mask.shape  # Shape of the reconstructed volume
        return

    def __len__(self):
        """
        Returns number of voxels
        """
        return len(np.ndarray.flatten(self.mask))

    def get_1D_array(self, ndarray):  # Necessary?
        """
        Returns the 1D array of the reconstructed volume
        """
        return np.ndarray.flatten(ndarray)
