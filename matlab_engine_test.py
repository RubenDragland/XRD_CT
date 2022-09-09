import matlab.engine
import numpy as np
import sys
import os

# No importing of matlab packages necessary. Located in matlab environment

eng1 = matlab.engine.start_matlab()

max1 = eng1.SASTT_calculate_measurement_time # Some function

eng1.quit()

print(max1) # Some output
