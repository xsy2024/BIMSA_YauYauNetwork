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
Characterization based on realistic data without loss of generality, we use a system of stochastic differential equations of the following form to construct the Yau-Yau network:

$$
\left\\{
\\begin{array}{ccc}
dX_t &= L\Phi(X_t)dt + \sigma_ X dW_t, &X_0 = x_0;\\
\\
\\
dY_t &= \sum_{i=1}^{n} h_i(X_{i;t})dt +dV_t, &Y_0 = 0,
\\end{array}
\right.
$$

where $\Phi(X_ t)$ is a functions' libarary of $X_t$, $h_i(X_{i;t})$ are polynomial functions of the $i _ {th}$ term of $X_t$, and the first stochastic differential equation characterizes the population in terms of population mean and variance, the second stochastic differential equation calibrates the specifics of an individual in a population using specific individual observations $Y_t$, each $h_i(X_{i;t})$ captures the strength of the influence of the individual's $i _ {th}$ characteristic on other characteristics.

The yauyauSDE class is mainly implemented to determine the form of equations from data.
After we create an instance of yauyauSDE,
```
yysde = yauyauSDE(name, data)
```
then we can determine the form of the first differential equation with the following line of code:

```
yysde.h()
```
The second diffrential equation can be obtained from the following code:

```
yysde.H_self(no)
```
We can invert the determined equations with this line of code and can visualize the computation:
```
yysde.Euler_gene(time, step, no)
```
At this point we have transformed the problem into a filtering problem, using the file `Yau-Yau filter algorithm` can real state what we really want.

After obtaining the true state, we use the weight coefficients previously calculated to obtain the network graph based on the true state:

```
yysde.w_individual(state=estimate, name , no)
```
For more, you can choose to use `draw_graph` and `draw_static_graph` to draw dynamic networks Video.

##Example
`main.py` gives an example on #64 microbe.

