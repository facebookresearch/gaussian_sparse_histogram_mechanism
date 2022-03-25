# Copyright (c) Meta Platforms, Inc. and affiliates.
# All rights reserved.

# This source code is licensed under the license found in the
# LICENSE file in the root directory of this source tree.

##R Code to Generate Figures for "Exact Privacy Analysis of Gaussian Sparse Histogram Mechanism" paper

delta_bw = function(epsilon, mu){
    delta = pnorm(mu/2 - epsilon/mu) - exp(epsilon)*pnorm(-mu/2 - epsilon/mu)
    return(delta)
}

tau_diff_func = function(delta, c_u, sigma){
    tau_diff = qnorm((1 - delta)^(1/c_u))*sigma
    return(tau_diff)
}

##Necessary optimization functions:

#Root finder for mu:
findrootmu = function(candidate_mu, epsilon, delta){
    m = numeric(100)
    mu = candidate_mu
    for(i in 1:length(m)){
        f = delta_bw(epsilon, mu) - delta
        fp = dnorm(mu/2 - epsilon/mu)*(1/2 + epsilon/mu^2) - delta - exp(epsilon)*dnorm(-mu/2 - epsilon/mu)*(-1/2 + epsilon/mu^2)
        mu = mu - f/fp
        m[i] = mu
    }
    return(mu)
}

#Root finder for epsilon:
findroot = function(candidate_epsilon, mu, delta){
    eps = numeric(100)
    epsilon = candidate_epsilon
    for(i in 1:length(eps)){
        f = delta_bw(epsilon, mu) - delta
        fp = dnorm(mu/2 - epsilon/mu)*(-1/mu) - exp(epsilon)*pnorm(-mu/2 - epsilon/mu) - exp(epsilon)*dnorm(-mu/2 - epsilon/mu)*(-1/mu)
        epsilon = epsilon - f/fp
        eps[i] = epsilon
    }
    return(epsilon)
}

#GSHM:
gshm_delta <- function(c_u, sigma, tau_diff, epsilon){
    mu = sqrt(c_u)/sigma

    delta_arr = delta_arr2 = numeric(c_u)
    c_u_iter = 1:c_u

    for(i in c_u_iter){
        a_eq = i-1
        mu = sqrt(c_u-a_eq)/sigma
        epsilon2 = epsilon - a_eq*log(pnorm(tau_diff/sigma))
        f = delta_bw(epsilon2, mu)
        delta_arr[i] = 1 - pnorm(tau_diff/sigma)^a_eq + pnorm(tau_diff/sigma)^(a_eq)*f
        epsilon3 = epsilon + a_eq*log(pnorm(tau_diff/sigma))
        delta_arr2[i] = delta_bw(epsilon3, mu)

    }

    final_delta = max(c(
        1 - pnorm(tau_diff/sigma)^c_u,

        #Max second term
        max(delta_arr),

        #Max third term
        max(delta_arr2)
    ))

    #Return final epsilon, delta
    return(final_delta)

}

## Set parameters
c_u = 51914
eps_views = 0.349 # Epsilon for the thresholding column
delta = 10^(-5)
k = sqrt(c_u)
mu_arr = 0 # use c() to list mu values for datasets with additional columns
extra_y = 0 # Use to get extra space on the y-axis in figure if needed

##Create sigma, tau curve:
mu = findrootmu(sqrt(eps_views), eps_views, delta)
min_sigma = k/mu

mu = sqrt(mu^2 + sum(mu_arr^2))
epsilon = findroot(.01, mu, delta)

#Cover a range of values from min_sigma to min_sigma + 300
sigma = min_sigma + c(seq(10^(-2), 10^(-1), 10^(-2)), seq(0.1, 300, .01))

#Initialize variables
tau_diff_gshm = tau_diff = numeric(length(sigma))

for(i in 1:length(sigma)){
    mu = k/sigma[i]
    mu = sqrt(mu^2 + sum(mu_arr^2))
    delta_prime = delta_bw(epsilon, mu)
    delta_infinite = delta - delta_prime
    tau_diff[i] = tau_diff_func(delta_infinite, c_u, sigma[i])
    tau_diff_gshm[i] = tau_diff_func(delta, c_u, sigma[i])
}

#Check that value of tau_diff meets target delta for GSHM (within delta^2 of the target delta):
tau_diff_gshm_check = c()
for(i in seq(1, length(sigma), 500)){
    tau_diff_gshm_check = c(
        tau_diff_gshm_check,
        abs(gshm_delta(c_u, sigma[i], tau_diff_gshm[i], epsilon) - delta) < delta^2
    )
}
sum(tau_diff_gshm_check)/length(tau_diff_gshm_check)

min_y = min(c(tau_diff, tau_diff_gshm))
max_y = max(c(tau_diff, tau_diff_gshm))

## Figure A:
plot(sigma, tau_diff, type="l", ylim=c(min_y, max_y+extra_y), ylab=expression(paste(tau, "*", -tau)), xlab=expression(sigma), lwd=2, cex.axis=1.3, cex.lab=2.1, mgp=c(2.5,1,0))

#Thin out the line using seq
lines(sigma[seq(1, length(sigma), 10)], tau_diff_gshm[seq(1, length(sigma), 10)], col="blue", lwd=2, lty=3)
abline(v=round(min_sigma), lty=2, col="gray")

min_sigma_label = substitute(paste(sigma, " = ", min_sigma), list(min_sigma = round(min_sigma)))
legend("topright", legend=c("Add deltas", "GSHM", as.expression(min_sigma_label)), lty=c(1, 3, 2), col=c("black", "blue", "gray"), lwd=c(2, 2, 1), cex=1.8)


###########
#### Epsilon, Delta curves
###########

c_u = 51914
epsilon = 0.349
delta = 10^(-5)
delta_infinite = 10^(-8)
mu = findrootmu(sqrt(epsilon), epsilon, delta)
sigma = sqrt(c_u)/mu

#Find epsilon under delta = delta_infinite
eps_max = findroot(.01, mu, delta_infinite)

# Points for particular delta(epsilon)
delta_plot_vals = 10^(seq(-8, -5, 1))
eps_plot_vals = numeric(length(delta_plot_vals))
for(i in 1:length(delta_plot_vals)){
    eps_plot_vals[i] = findroot(.01, mu, delta_plot_vals[i])
}

# Calculate deltas for points on curve
delta_arr_add_deltas_plot_vals = delta_arr_gshm_plot_vals = numeric(length(eps_plot_vals))
annotations = c()

for(i in 1:length(eps_plot_vals)){
    delta_arr_add_deltas_plot_vals[i] = delta_bw(eps_plot_vals[i], mu) + delta_infinite
    tau_diff = tau_diff_func(delta_infinite, c_u, sigma)
    delta_arr_gshm_plot_vals[i] = gshm_delta(c_u, sigma, tau_diff, eps_plot_vals[i])
    d = delta_arr_gshm_plot_vals[i]
    annotations = c(
        annotations,
        substitute(paste(delta[GSHM], "=", delta_arr_gshm_plot_vals), list(delta_arr_gshm_plot_vals=round( d, ceiling(-log(d, 10)) ) ))
    )
}

# Calculate deltas for the full curve
# delta(epsilon) curves
eps_arr = c(seq(0.1, eps_max, .005), eps_max)
delta_arr_add_deltas = delta_arr_gshm = numeric(length(eps_arr))


for(i in 1:length(eps_arr)){
    delta_arr_add_deltas[i] = delta_bw(eps_arr[i], mu) + delta_infinite
    tau_diff = tau_diff_func(delta_infinite, c_u, sigma)
    delta_arr_gshm[i] = gshm_delta(c_u, sigma, tau_diff, eps_arr[i])
}

###Figure B:
min_val = min(delta_arr_add_deltas/delta_arr_gshm)
max_val = max(delta_arr_add_deltas/delta_arr_gshm)


plot(delta_arr_add_deltas/delta_arr_gshm, eps_arr, type="l", ylab=expression(epsilon), xlab=expression(delta[add]/delta[GSHM]), lwd=2, xlim=c(min_val, max_val+0.06), cex.lab=1.9, las=1, cex.axis=1.3)
points(delta_arr_add_deltas_plot_vals/delta_arr_gshm_plot_vals, eps_plot_vals, pch=20)
text(delta_arr_add_deltas_plot_vals[-1]/delta_arr_gshm_plot_vals[-1], eps_plot_vals[-1]-0.008, as.expression(annotations)[-1], pos=4)
text(delta_arr_add_deltas_plot_vals[1]/delta_arr_gshm_plot_vals[1], eps_plot_vals[1]-0.008, as.expression(annotations)[1], pos=1)
