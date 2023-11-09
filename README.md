# Kalman filter in econometrics
Application of the state-space approach via the discrete-time linear Kalman filter for estimating the class of stochastic volatility models.

# References
Each code (implementation method) includes the exact reference where the particular algorithm was published. If you use these codes in your research, please, cite the corresponding articles mentioned in the codes.

# Remark
The codes have been presented here for their instructional value only. They have been tested with care but are not guaranteed to be free of error and, hence, they should not be relied on as the sole basis to solve problems.

# Steps to reproduce
- `Illustrate*` are the scripts that illustrate the obtained estimates (over time). You can find their call at the end of each script.
- `Fitting*` are the scripts that perform the model calibration by the KF state-space approach and the method of maximum likelihhod. The optimization method starts with various initial values and performs `MC=` number of trials.

# List of the models and estimation methods used
### Stochastic volatility models:
- `Fitting_SV1` is the script for the QML method (Harvey A., Ruiz E., Shephard N., 1994, <a href="https://doi.org/10.2307/2297980">DOI</a>) for estimating the SV models (type 1 parametrization) in a general setting, that is, the unrestricted model.
  - `Fitting_GaussianSV1` is for estimating the SV1 model with Gaussian distribution assumption in the original model 

