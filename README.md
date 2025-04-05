## Table of Contents
- [Introduction](#introduction)
- [Usage](#usage)
- [Example](#example)

  ## Introduction
Yau-Yau filtering algorithm is designed to handle nonlinear filtering systems in continuous form:

$$
\left\\{
\\begin{array}{ccc}
dX_t &= f(X_t)dt + \sigma_ X dW_t, &X_0 = x_0;\\
\\
\\
dY_t &= h(X_t)dt +dV_t, &Y_0 = 0,
\\end{array}
\right.
$$

where 
- $X_t$, states, $n$-vector; 
- $Y_t$, observation path, $m$-vector;
- $f$, drift term, $n$-vector;
- $\sigma_ X$, diffusion term, $n$-vector;
- $h$, observation term, $m$-vector
- $W_t$, $n$-vector standard Brownian motion;
- $V_t$, $m$-vector standard Browian motion; $\\{W_t, t\ge 0\\}$ and $\\{V_t, t\ge0\\}$ are independent.

Yau-Yau algorithm computes the density function of $E(X_t |\\{Y_s, 0 \le s \le t\\})$ and is **memoryless** and **real-time** in its computation. 
Due to these excellent properties of yauyau filtering, we developed a set of stochastic complex network construction methods based on the Yau-Yau filtering method and named this method as Yau-Yau network. 

  ## Usage
How to apply Yau-Yau network to data analysis?
