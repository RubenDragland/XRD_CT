import numpy as np
import scipy.io


class TensorTomographyReconstruction:
    """ """

    def __init__(self, path_to_mat_file, dataset=False):

        self.path_to_mat_file = path_to_mat_file
        self.mat = scipy.io.loadmat(self.path_to_mat_file)
        if dataset:
            self.key = "fasit"
        else:
            self.key = "s"
            self.error_data = self.mat["e"]
            self.timing_data = self.mat["timing"]

        self.mask = self.mat[self.key]["mask3D"][0, 0]
        self.theta = self.mat[self.key]["theta"][0, 0][0, 0][0] % (np.pi)
        self.phi = self.mat[self.key]["phi"][0, 0][0, 0][0] % (np.pi)
        self.a0 = self.mat[self.key]["a"][0, 0][:, 0][0][0]
        self.a2 = self.mat[self.key]["a"][0, 0][:, 1][0][0]
        self.a4 = self.mat[self.key]["a"][0, 0][:, 2][0][0]
        self.a6 = self.mat[self.key]["a"][0, 0][:, 3][0][0]

        self.params = np.array(
            [self.a0, self.a2, self.a4, self.a6, self.theta, self.phi]
        )

        self.indexing_dict = {
            "a0": self.a0,
            "a2": self.a2,
            "a4": self.a4,
            "a6": self.a6,
            "theta": self.theta,
            "phi": self.phi,
        }

        self.shape = self.mask.shape  # Shape of the reconstructed volume
        self.format_0_pi()

        return

    def __len__(self):
        """
        Returns number of voxels
        """
        return len(np.ndarray.flatten(self.mask))

    def __getitem__(self, key):
        """
        Returns the value of the parameter given by key
        """
        return self.indexing_dict[key]

    def get_1D_array(self, ndarray):
        """
        Returns the 1D array of the reconstructed volume
        """
        return np.ndarray.flatten(ndarray)

    def format_0_pi(self):
        """
        Formats the theta and phi angles to be between 0 and pi
        """
        self.theta = self.theta % np.pi
        self.phi = self.phi % np.pi

        self.theta[self.theta < 0] += np.pi
        self.phi[self.phi < 0] += np.pi
        return
