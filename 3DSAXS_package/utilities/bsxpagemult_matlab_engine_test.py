import matlab.engine
import numpy as np
import time
import sys
import torch


def page_multiply(A: torch.tensor, B: torch.tensor):
    # Check if numba compatible. More variables need to be converted to tensor.

    # assert A.size[1] == B.size[0], "A and B are not compatible for multiplication"
    assert len(A.size()) == 3, "A is not a 3D tensor"
    assert len(B.size()) <= len(A.size()), "B is not a 3D or 2D tensor"

    if A.size()[-1] == B.size()[-1] and len(B.size()) == 3:

        C = torch.einsum("ijm,jkm->ikm", A, B)  # Check this.

        # for i in range(A.size()[-1]):
    #
    #    A[:, :, i] = torch.matmul(
    #        A[:, :, i], B[:, :, i]
    #    )  # Check if works compared to creating new Tensor.

    elif len(A.size()) > len(B.size()):

        A = torch.einsum("ijk,jm->imk", A, B)
        # C = torch.zeros(A.size()[1], B.size()[1], A.size()[-1])
        # for i in range(A.size()[-1]):
        # C[:, :, i] = torch.tensordot(
        #    A[:, :, i], B[:, :], dims=([0], [0])
        # )  # torch.matmul(A[:, :, i], B[:, :])
        #    C[:, :, i] = torch.matmul(A[:, :, i], B[:, :])
    return C


# A = torch.rand((3, 3, 64))
# B = torch.rand((3, 8))
A = torch.rand((3, 3, 3))
B = torch.rand((3, 3, 3))

times = 100

tic = time.time()
for i in range(times):
    # C = torch.tensordot(A, B, dims=([1], [0])).reshape((3, 8, 64))
    C = page_multiply(A, B)  # torch.zeros((3, 8, 64))
    pass
tac = time.time()

print("torch.matmul(A,B) time: ", tac - tic)

eng = matlab.engine.start_matlab()
A = matlab.double(A.numpy().tolist())
B = matlab.double(B.numpy().tolist())

tic = time.time()
for i in range(times):
    D = eng.bsxpagemult(A, B)
tac = time.time()

print("bsxpagemult(A,B) time: ", tac - tic)

tic = time.time()
for i in range(times):
    E = eng.pagemtimes(A, B)
tac = time.time()

print("mtimes(A,B) time: ", tac - tic)

if np.allclose(C.numpy(), np.array(D), atol=1e-8):
    print("torch(A,B) is correct")
if np.allclose(C.numpy(), np.array(E), atol=1e-8):
    print("torch(A,B) fits pagemtimes")
if np.allclose(np.array(D), np.array(E), atol=1e-8):
    print("Validation of bsxpagemult and pagemtimes")

print(C.numpy(), np.array(D), np.array(E))
eng.quit()
sys.exit()

# Consider tensordot.

matrix_shape = (6, 4, 4, 4)  # (46,34,34,4,1) # y,x,z,l,m
matrix_shape = (4, 8, 1000)
matrix_shape = (3, 3, 3)
rng = np.random.default_rng()  # (seed=69)


a_temp = np.ones(
    matrix_shape
)  # np.ones(matrix_shape) #rng.normal(size = matrix_shape )
Ylm = np.ones(
    matrix_shape
)  # + 1j*np.ones(matrix_shape) #rng.normal(size = matrix_shape) + 1j * rng.normal(size = matrix_shape)


sys.exit()
times = 100000

# Test speed in python

tic = time.time()
for i in range(times):
    a_Ylm_sum = np.tensordot(a_temp, Ylm, axes=([0, 1], [0, 1]))  # y,x,z, (lxl)
tac = time.time()

print(a_Ylm_sum.shape)

print(f"Python tensordot: {tac-tic} s for {times} runs")

sys.exit()
eng1 = matlab.engine.start_matlab()


a_temp = matlab.double(np.ndarray.tolist(a_temp))
Ylm = matlab.double(np.ndarray.tolist(Ylm), is_complex=True)
print(np.array(Ylm).shape)
# Ylm = eng1.permute(Ylm, matlab.int8([2,1]) )
print(np.array(Ylm).shape)

tic = time.time()
for i in range(times):
    matrix_product = eng1.bsxpagemult(a_temp, Ylm)  # C-written matrix product
tac = time.time()
eng1.quit()

print(f"C bsxpagemult: {tac-tic} s for {times} runs")


print(f"\nPython: \n{a_Ylm_sum}\nC: \n{np.array(matrix_product)}")


"""
Tests:

Matrix multiplication does work. Tensordot significantly faster. 

Matrix times Tensor works. No need to be aware of row*col dimensions. 

Final challenge: Make dimensions align

Speed issue: Check if the speed is due to matlab engine by calling c-function from .m-file. s

Use of matlab engine results in significant increase in computation time for mex-function bsxpagemult. 

















"""
