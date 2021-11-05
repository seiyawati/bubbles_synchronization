import math

nb = 2 # number of bubble
h = 0.5 * math.exp(-9) # time_step_for_computation
Tmax = 1.0 * math.exp(-4) # maximum_time_of_computation
eps = 1.0 * math.exp(-6) # tolerance_of_RKF45_scheme 
Tprt = 5.0 * math.exp(-7) # output_interval_of_pressure_distribution 
Xmin = -1.0 * math.exp(-3) # x_minimum_of_output_interval 
Xmax = 5.0 * math.exp(-3) # x_maximum_of_output_interval
Xint = 1000 # number_of_division_of_output_interval
