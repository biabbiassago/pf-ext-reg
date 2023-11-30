File: `sampling-joints.R` contains functions to sample each parameter.  
File: `bayesian_max_model.R` contains the function to fit the model in question + a couple of auxiliary functions for plots.     
File:: `fit_model1.R` contains the file to run the model with a simulated dataset (in data folder) and save the results.    



### One step model:

Let $s_i$ , $i = 1,2,..$ be the sampled locations, and $y_{ij}$ be the max at location $i$ over time

**Level 1**

$$
\vec{Y}(s_i)  \overset{i.i.d}\sim GEV(\mu(s_i), \sigma, \xi) \quad \xi=0
$$

**Level 2**

$$
\vec{\mu} \sim GP(m,\Sigma)
$$  

$$
m | \alpha_0, \alpha_1, x = \alpha_0 + \alpha_1 \times x
$$  

And the covariance $\Sigma$ is induced by the exponential covariance function:

$$
k(s,s') = \beta_0 \times exp(-\beta_1 \lVert s-s' \rVert)
$$


**Level 3**


$\sigma \sim \text{lnNorm}(0,\frac{1}{2})$  
$\vec{\alpha} \sim MVN_2((0,0),\bar{V})$  
$\beta_0 \sim \text{InvGamma}(3,1)$  
$\beta_1 \sim \text{lnNorm}(-\frac{1}{2},\frac{1}{2})$



