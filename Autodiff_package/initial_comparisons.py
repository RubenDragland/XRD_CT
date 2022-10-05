from importlib.metadata import requires
import matplotlib.pyplot as plt
import numpy as np
#from numba import njit Numpy version issue
import torch
from autograd import elementwise_grad as egrad
import autograd.numpy as anp

import time, timeit
import sys


# Functions/gradients expression f = ln(x1) + x1*x2 + sin(x2)


def symbolic_f(x1,x2)-> np.array:
    return np.array( [1/x1 + x2 , x1 - np.cos(x2) ] )

def torch_f(x1,x2)-> torch.tensor:
    return torch.log(x1) + x1*x2 - torch.sin(x2)

def numpy_f(x1,x2) -> np.array:
    return np.log(x1) + x1*x2 - np.sin(x2)

def autograd_f(x1,x2) -> anp.array:
    return anp.log(x1) + x1*x2 - anp.sin(x2)



# Functions/gradient expression g = ln(x1*x2)*sin(x2)

def symbolic_g(x1,x2)-> np.array:
    return np.array( [np.sin(x2)/x1, np.sin(x2)/x2 + np.cos(x2)*np.log(x1*x2) ] )

def torch_g(x1,x2)-> torch.tensor:
    return torch.log(x1*x2)*torch.sin(x2)

def numpy_g(x1,x2) -> np.array:
    return np.log(x1*x2)*np.sin(x2)

def autograd_g(x1,x2) -> anp.array:
    return anp.log(x1*x2)*anp.sin(x2)

# Functions/gradients expression h = 64x*(1-x)(1-2x)^2(1-8x+8x^2)^2 (Attempt at slowing down symbolic)

def symbolic_h(x1,x2)-> np.array:
    return 128*x1*(1-x1)*(-8+16*x1)*((1-2*x1)**2)*(1-8*x1+8*x1**2) + 64*(1-x1)*((1-2*x1)**2)*((1-8*x1+8*x1**2)**2) - (64*x1*(1-2*x1)**2)*(1-8*x1+8*x1**2)**2 - 256*x1*(1-x1)*(1-2*x1)*(1-8*x1+8*x1**2)**2

def torch_h(x1,x2)-> torch.tensor:
    return 64*x1*(1-x1)*( (1-2*x1)**2 ) * (1-8*x1+8*x1**2)**2

def numpy_h(x1,x2) -> np.array:
    return 64*x1*(1-x1)*( (1-2*x1)**2 ) * (1-8*x1+8*x1**2)**2

def autograd_h(x1,x2) -> anp.array:
    return 64*x1*(1-x1)*( (1-2*x1)**2 ) * (1-8*x1+8*x1**2)**2


### Symbolic differentiation


def CGD_gradient(): # Has to be defined for all the types of variables.
    return 


def symbolic_diff(f, x1, x2):
    return f(x1,x2)

### Automatic differentiation

def AD_forward_mode():
    return

def AD_reverse_mode(): # If I want to implement it by myself. 
    return


def AD_torch(func, x1: np.array, x2: np.array )-> torch.tensor:
    '''Applies automatic differentiation provided by pytorch'''

    x = torch.tensor(np.array([x1, x2]) , requires_grad=True) 
    y = func(x[0], x[1]) # Assumes expression implemented using torch.
    print(y)
    #y.backward( torch.tensor( np.ones(len(x1)) ) )
    y.backward( )
    return x.grad

def AD_torch_general(func, *args, **kwargs) -> torch.tensor:
    '''
    Applies automatic differentiation provided by pytorch

    Parameters: 
    -----------
    func: function
        The cost function to be called
    *args: list
        The arguments of the cost function
    **kwargs: dict
        The variables for whom one shall differentiate.

    returns: torch.tensor
        The gradients    
    '''
    first_var = list(kwargs)[0]
    x = torch.ones(size = (len(kwargs), *kwargs[first_var].size() ) )# , requires_grad=True)
    for it, (key, elem) in enumerate( kwargs.items() ):
        x[it] = elem
    x.requires_grad_(True)
    y = func( *args, x )

    y.backward( ) 

    return x.grad 



def AD_autograd(func, x1: anp.array, x2: anp.array) -> anp.array:
    return egrad(func, (0,1))(x1,x2) #egrad( func(x1,x2) )


# Numerical

# def central_diff(f, x, h=1e-6):
#     grad = np.zeros(x.shape)
#     grad[0] = ( f( x[0] + h, x[1] ) - f( x[0]-h, x[1] )   )/ (2*h)
#     grad[1] = ( f( x[0], x[1] + h ) - f( x[0], x[1] - h ) )/ (2*h)
#     return grad

def numerical_diff(f, x1: np.array , x2: np.array, h=1e-6):
    grads = np.zeros( (2,len(x1) ) )
    grads[0] = ( f( x1 + h, x2 ) - f( x1-h, x2 )   )/ (2*h)
    grads[1] = ( f( x1, x2 + h ) - f( x1, x2 - h ) )/ (2*h)
    return grads

def numerical_diff_general(f, *args, **kwargs):
    """Generalisation of numerical_diff()"""
    h = 1e-4 # Set in stone
    grads = {}
    for key, elem in kwargs.items():
        grads[key] = ( f( *args, elem + h, ) - f( *args, elem-h) )/ (2*h)
    return grads


def interpolation_test(x1, x2):

    x = torch.zeros(1, requires_grad=True)
    y = torch.zeros(1, requires_grad=True)

    for i in range( len(x1) ):
        if (x1[i]>9) and (x2[i]>9):
            x = x+ x1[i]
            y = y+ x2[i]

    return x**2 + y**2#x1**2 + x2**2#x**2 + y**2


if __name__ == "__main__":

    print(np.arange(1,10,1, dtype=np.float32))

    grad = AD_torch(interpolation_test, np.arange(1,10,1, dtype=np.float32) , np.arange(1,10,1, dtype=np.float32) ) 

    print(grad)


    sys.exit()

    num_elem = int(1e7)

    rng = np.random.default_rng(seed=69)
    x = np.abs( rng.normal(size = num_elem) )
    y = np.abs( rng.normal(size = num_elem) )

    times_dict = {
        "Symbolic": [],
        "AD Torch": [],
        "AD Autograd": [],
        "Numerical" : []
    }

    grad_dict = {
        "Symbolic": [],
        "AD Torch": [],
        "AD Autograd": [],
        "Numerical" : []
    }

    func_dict = {
        "Symbolic": [symbolic_f, symbolic_g, symbolic_h],
        "AD Torch": [torch_f, torch_g, torch_h],
        "AD Autograd": [autograd_f, autograd_g, autograd_h],
        "Numerical" : [numpy_f, numpy_g, numpy_h]
    }

    diff_mode_dict = {
        "Symbolic": symbolic_diff,
        "AD Torch": AD_torch,
        "AD Autograd": AD_autograd,
        "Numerical" : numerical_diff,
    }

    for i, (key, v) in enumerate(times_dict.items() ):

        for func in func_dict[key]:

            tic = time.time()

            grad_dict[key].append( diff_mode_dict[key](func, x, y) )

            tac = time.time()

            times_dict[key].append(tac-tic)

        if num_elem == 1:

            print(f"\n{key}:\ngrad: {grad_dict[key]}\ntime: {times_dict[key]}\ntime/it: {np.array(times_dict[key])/num_elem}\n")
        
        else:
            print(f"\n{key}:\ntime: {times_dict[key]}\ntime/it: {np.array(times_dict[key])/num_elem}\n")


















