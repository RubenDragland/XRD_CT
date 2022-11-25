import numpy as np
import scipy.io


class TensorTomographyReconstruction:
    """ """

    def __init__(self, path_to_mat_file, dataset=False):

        self.path_to_mat_file = path_to_mat_file
        self.mat = scipy.io.loadmat(self.path_to_mat_file)
        if dataset:
            self.key = "fasit"
            self.slice = self.mat["slice"][0]
        else:
            self.key = "s"
            self.error_data = self.mat["e"]
            self.timing_data = np.squeeze(self.mat["timing"][0])
            self.tot_time = np.squeeze(self.mat["tot_time"][0])
            self.loss_curve = np.squeeze(self.mat["Err_hist"])

        self.mask = self.mat[self.key]["mask3D"][0, 0]
        self.theta = self.mat[self.key]["theta"][0, 0][0, 0][0] % (np.pi)
        self.phi = self.mat[self.key]["phi"][0, 0][0, 0][0] % (np.pi)
        self.shape = self.mask.shape
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

    def flatten(self, key):
        """
        Returns the 1D array of the reconstructed volume
        """
        return np.ndarray.flatten(self.indexing_dict[key])

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

    def slice_params(self, fasit_slices):
        """
        Returns the accuracte slice of the reconstructed volume
        """
        self.slice = fasit_slices
        self.sliced_dict = {}
        for key in self.indexing_dict.keys():
            self.sliced_dict[key] = self.indexing_dict[key][
                fasit_slices[0] : fasit_slices[1],
                fasit_slices[0] : fasit_slices[1],
                fasit_slices[0] : fasit_slices[1],
            ]
        return self.sliced_dict


class SH_Reconstruction(TensorTomographyReconstruction):
    """ """

    def __init__(self, path_to_mat_file, dataset=False):
        super().__init__(path_to_mat_file, dataset=dataset)

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


class EXPSIN_Reconstruction(TensorTomographyReconstruction):
    """ """

    def __init__(self, path_to_mat_file, dataset=False):
        super().__init__(path_to_mat_file, dataset=dataset)

        self.A = self.mat[self.key]["A"][0, 0][:, 0][0][0]
        self.B = self.mat[self.key]["B"][0, 0][:, 0][0][0]

        self.params = np.array([self.A, self.B, self.theta, self.phi])

        self.indexing_dict = {
            "A": self.A,
            "B": self.B,
            "theta": self.theta,
            "phi": self.phi,
        }
