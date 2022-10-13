import matlab.engine
import numpy as np
import time
import sys


# Consider tensordot.

matrix_shape = (6,4,4,4)#(46,34,34,4,1) # y,x,z,l,m
matrix_shape = (4,8,1000)
matrix_shape = (3,3,3)
rng = np.random.default_rng()#(seed=69)


a_temp = np.ones(matrix_shape ) #np.ones(matrix_shape) #rng.normal(size = matrix_shape )
Ylm = np.ones(matrix_shape) #+ 1j*np.ones(matrix_shape) #rng.normal(size = matrix_shape) + 1j * rng.normal(size = matrix_shape)


sys.exit()
times = 100000

# Test speed in python

tic = time.time()
for i in range(times):
    a_Ylm_sum = np.tensordot( a_temp, Ylm, axes= ( [0,1],[0,1] )) # y,x,z, (lxl) 
tac = time.time()

print(a_Ylm_sum.shape)

print(f"Python tensordot: {tac-tic} s for {times} runs")

sys.exit()
eng1 = matlab.engine.start_matlab()



a_temp = matlab.double(  np.ndarray.tolist(a_temp)  )
Ylm = matlab.double( np.ndarray.tolist( Ylm ), is_complex = True )
print(np.array(Ylm).shape)
#Ylm = eng1.permute(Ylm, matlab.int8([2,1]) )
print(np.array(Ylm).shape)

tic = time.time()
for i in range(times):
    matrix_product = eng1.bsxpagemult(a_temp, Ylm ) # C-written matrix product
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