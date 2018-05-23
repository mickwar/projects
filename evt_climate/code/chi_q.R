# Compute chi(u) for several u
# Note that lim u -> 1 chi(u) = chi, the measure of asymptotic dependence
angle.2.cone = function(theta){
    out = matrix(1, length(theta), 2)
    out[,1] = tan(pi/4 - abs(theta - pi/4))
    ind = which(theta <= pi/4)
    out[ind,] = out[ind, c(2,1)]
    return (out)
    }

# Sample a bunch of angles from the simple Pareto proces
n = 1e5
# V = angle.2.cone(rbeta(n, 1/8, 1/8) * pi/2)
V = angle.2.cone(rbeta(n, 1/8, 1/4) * pi/2)

plot(V[sample(10000, replace = TRUE),], pch = 16, col = rgb(0.0,0.0,0.0,0.1))

plot(density(V[,1]), col = 'red', lwd = 2)
lines(density(V[,2]), col = 'blue', lwd = 2)

mean(V[,1])
mean(V[,2])

# Construct grid for quantiles u
u = seq(0.000, 0.999, by = 0.001)

# Sample from Pareto
Y = 1/runif(n)

# Compute denominators (marginals, from conditional probability)
denom1 = sapply(u, function(x) mean(Y * V[,1] > mean(V[,1]) / (1-x)) )
denom2 = sapply(u, function(x) mean(Y * V[,2] > mean(V[,2]) / (1-x)) )

# Plot marginal
plot(u, denom1, type = 'l', xlim = c(0, 1), ylim = c(0, 1), lwd = 2, col = 'blue')
lines(u, denom2, col = 'red', lty = 2, lwd = 2)
abline(1, -1)
abline(v = 1-mean(V[,1]), lty = 2, col = 'gray50')
abline(h = mean(V[,1]), lty = 2, col = 'gray50')
abline(v = 1-mean(V[,2]), lty = 2, col = 'gray50')
abline(h = mean(V[,2]), lty = 2, col = 'gray50')

# Compute joint probability
Vm1 = cbind(mean(V[,1]) / V[,1], mean(V[,2]) / V[,2])
Vm2 = cbind(V[,1] / mean(V[,1]), V[,2] / mean(V[,2]))
Vmax = apply(Vm1, 1, max)
Vmin = apply(Vm2, 1, min)
numer = sapply(u, function(x) mean(Y > Vmax / (1-x)) )

# The quantity currently given as chi
mean(Vmin)

# The calculation in appendix 3 requires that u be larger than
max(1 - mean(V[,1]), 1 - mean(V[,2]))
# in order for some of the equalities to hold true, which in
# turn allows us to cancel out u (notice the flat portion
# for u greater than the above; it not longer becomes a function
# of u and so the limit goes away)

range(numer)
range(denom1)
range(denom2)

# Plot of each conditional probability
plot(u, numer / denom1, type= 'l', col = 'blue', lwd = 3, ylim = c(0, 1))
lines(u, numer / denom2, type= 'l', col = 'red', lwd = 3, lty = 2)
abline(v = 1-mean(V[,1]), lty = 2, col = 'blue')
abline(v = 1-mean(V[,2]), lty = 2, col = 'red')
abline(h = mean(Vmin), col = 'darkgreen', lwd = 2)
abline(h = mean(Vmin), col = 'darkgreen', lwd = 2)

# When u is large enough, it doesnt matter which variable is being conditioned on.
# In either case, the limit is the same. 

# It would seem that chi is the minimum probability that can be used to describe
# the conditional probability of extreme events (at least as a simple Pareto
# process is concerned)

# Can I get a maximum? Maybe, just consider the function chi(u) when u = 0.
# There are two functions chi(u), depending on which of the two variables is being
# conditioned on. The larger of chi_12(0) and chi_21(0) could serve as a max.

