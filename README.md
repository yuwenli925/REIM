# REIM
This repo is on the MATLAB codes implemented in numerical experiments of the paper "Rational empirical interpolation methods with applications" by Aidi Li and Yuwen Li posted as arXiv:2406.19339

In the folder 'FEM', the functions squaremesh.m, uniformrefine.m, gradbasis.m are borrowed from the iFEM package developed by Prof. Long Chen at UC Irvine.

aaa.m is borrowed from the chebfun package. The AAA algorithm was developed by Yuji Nakatsukasa, Olivier Sète, Lloyd N. Trefethen, see the paper "The AAA algorithm for rational approximation".

EIM.m is on the poles, interpolation points and the Lebesgue Constants, which generates Fig.1 in the paper.

REIM.m is the function used for all other scripts. When the approximation interval and the objective functions change, the dictionary of the negative poles and the interpolation points (bset and xset) should change by users to achieve better approximation results.

REIM_power1_demo.m is on the REIM rational approximation of x^(-s) on [1e-6,1], in which the OGA and the AAA are also used. It generates Fig.2 in the paper.

REIM_power2_demo.m is on the REIM rational approximation of x^(-s) on [1e-8,1], which generates left of Fig.4.

FEM2D_fractional_demo.m solves the fractional Laplace equation on uniform and graded meshes. It outputs Fig. 4 (right). Setting dynamic = 1 in this function produce the graded mesh in Fig. 3 (right).

REIM_time_demo.m is on the REIM rational approximation of 1/(x^s+d) on [1e-6,1], which generates Fig.5.

BDF2_FEM.m is on the L^2 errors of numerical solutions and the accepted/rejected step sizes, which generates Fig.6.

REIM_exp_demo.m is on the REIM rational approximation of (exp(-tau x) - 1)/(-tau x) and exp(-tau x) on [1,1e6], which generates the right of Fig.7.

REIM_precon_demo.m is on the REIM rational approximation (x^(-0.5)+Kx^0.5)^-1 on [1e-6,1], which generates the left of Fig.7.


