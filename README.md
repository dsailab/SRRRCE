# SRRRCE

This repository contains the implementation of the SRRRCE (Sparse Reduced Rank Regression with Covariance Estimation).

We first implemented the algorithm proposed in [1]. In [2], we noted certain issues with the algorithm proposed in [1], where they used the DP-GLASSO algorithm to estimate the inverse covariance matrix, treating it as an unconstrained problem. The algorithm proposed in [2] based on block majorization minimization addresses this issue and improves the efficiency of the algorithm.

- `covsrrr_solver.m`: Implements the efficient algorithm in [1].
- `SRRRCE_solver.m`: Implements the method in [2].
- `dpglasso.m`: Implements the algorithm in [3].

[1] L. Chen and J. Z. Huang, "Sparse Reduced-Rank Regression With Covariance Estimation," Statistical Computing, vol. 26, pp. 461-470, 2016.

[2] F. Li and Z. Zhao, "Efficient Sparse Reduced-Rank Regression With Covariance Estimation," in 2023 IEEE Statistical Signal Processing Workshop (SSP), Hanoi, Vietnam, 2023, pp. 46-50.

[3] R. Mazumder and T. Hastie, "The Graphical Lasso: New Insights and Alternatives," Electronic Journal of Statistics, vol. 6, pp. 2125-2149, 2012.

