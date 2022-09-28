import scipy.io
import numpy as np



def create_smaller_dataset(filename:str, ny:int, nx:int):
    """
    Take raw .mat-file in same folder and save a copy with a new name, and a reduced number of voxels.
    Think it is clever to use adjacent areas.
    """

    assert not(ny%2) and not(nx%2) # Preferably even numbers    

    key = "projection"
    mat_dict = scipy.io.loadmat(filename) #["projection"]#[0,0][2]
    new_dict = mat_dict.copy()

    # Finding slices
    ny_org, nx_org, segments = mat_dict[key][0, 0][2].shape
    y_start, y_end = int( ny_org//2 - ny//2), int( ny_org//2 + ny//2) 
    x_start, x_end = int( nx_org//2 - nx//2), int( nx_org//2 + nx//2) 
    for i in range( mat_dict[key].shape[-1]):
        new_dict[key][0, i][2] = mat_dict[key][0,i][2][y_start:y_end, x_start: x_end, :] # Slicing


    scipy.io.savemat(filename[:-4] + f"ny{ny}nx{nx}" + ".mat", new_dict)
    return

create_smaller_dataset("Data sets/SASTT_carbon_knot_aligned_ASTRA_corrected.mat", 4, 4)