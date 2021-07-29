# ConstrainedLeastSquares

This package implements a general least squares method described in the paper *Evaluation of Measurements by the Method of Least Squares* by 
Lars Nielsen (2001).

The input of the least squares method is a set of $m$ measured parameters $z_1 \dots z_m$ and their respective uncertainties and covariances (covariance matrix $\Sigma$), a guess for a set of $k$ unknown parameters $\beta_1 \dots \beta_k$ and a set of $n$ constraints which define the measurement model:

$$
f_1(\beta, \zeta) = 0 \\
f_2(\beta, \zeta) = 0 \\
\dots \\
f_n(\beta, \zeta) = 0
$$

The chi-square function to minimize is:

$$
\chi^2(\vec{\zeta}, \vec{z}) = (\vec{z} - \vec{\zeta})^T\Sigma^{-1}(\vec{z} - \vec{\zeta})
$$

$\vec{\zeta}$ is a set of $m$ parameters that are tuned such that the constraint functions are satisfied while the deviation of $\zeta$ to $z$, weighted by the uncertainty and covariances of $z$, is minimal.
The output of the method is an estimate of optimal parameters $\hat{\beta}$ and $\hat{\zeta}$ and their uncertainties and covariance matrix.

## Usage

...