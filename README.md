# REIM
This repo is on the MATLAB codes implemented in numerical experiments of the paper "Rational empirical interpolation methods with applications" by Aidi Li and Yuwen Li posted as arXiv:2406.19339

EIM.m is on the poles, interpolation points and the Lebesgue Constants.

REIM_power1_demo.m is on the REIM rational approximation of x^(-s) on [1e-6,1];

REIM_power2_demo.m is on the REIM rational approximation of x^(-s) on [1e-8,1];

REIM_exp_demo.m is on the REIM rational approximation of (exp(-tau x) - 1)/(-tau x) and exp(-tau x) on [1,1e6];

REIM_time_demo.m is on the REIM rational approximation of 1/(x^s+d) on [1e-6,1];

REIM_precon_demo.m is on the REIM rational approximation (x^(-0.5)+Kx^0.5)^-1 on [1e-6,1];
