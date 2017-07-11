library(mwBASE)
library(mwEVT)

### Simulation study:
#   Calculate coverage (can I take 10000 samples, divide into 
#   ten 1000 segments and compute 10 credible bounds? No, need
#   to run the model with new data)
#
#   For several thetas   (3)
#   Along a range of thresholds (12)
#   For several replicates  (3)
#   For several obs / replicate (3)
#   For Ferro/Suveges Bayes
#   3 * 12 * 3 * 3 * 2 = 648 is too many combinations
#
#   Solution:
#   Extreme thetas, 0.15 0.85 (2)
#   One threshold
#   One replicate, 10
#   One obs, 1000
#   Both Ferro/Suveges
#   2 * 1 * 1 * 1 * 2 = 4 combos

library(foreach)
library(doMC)
registerDoMC(4)




doit = function(k, R, n, uq, theta){

    cat("\r", k, "     ")
    ### Simulation: Frechet
    ww = matrix(-1/log(runif(n*R)), n, R)
    y = matrix(0, n, R)
    y[1,] = ww[1,] / theta
    for (i in 2:n)
        y[i,] = apply(rbind((1-theta)*y[i-1,], ww[i,]), 2, max)

    ### For split up sequences
    uu = quantile(y, uq)
    N = sapply(uu, function(x) sum(y > x))

    ferro_ba = theta_hier(y, uu, likelihood = "ferro",
        nburn = 10000, nmcmc = 10000, display = 0)
    tmp1 = c(ferro_ba$keep, R+1)
    fb.v = apply(ferro_ba$mcmc, 2, hpd_mult, force_uni = TRUE)
    fb.r = apply(ferro_ba$mcmc, 2, function(x) sqrt(mean((x-theta)^2)))

    suveges_ba = theta_hier(y, uu, likelihood = "suveges",
        nburn = 10000, nmcmc = 10000, display = 0)
    tmp2 = c(suveges_ba$keep, R+1)
    sb.v = apply(suveges_ba$mcmc, 2, hpd_mult, force_uni = TRUE)
    sb.r = apply(suveges_ba$mcmc, 2, function(x) sqrt(mean((x-theta)^2)))

    out = NA*double(4*(R+1))
    out[tmp1 + 0*(R+1)] = (fb.v[1,] < theta & fb.v[2,] > theta)
    out[tmp1 + 1*(R+1)] = fb.r
    out[tmp2 + 2*(R+1)] = (sb.v[1,] < theta & sb.v[2,] > theta)
    out[tmp2 + 3*(R+1)] = sb.r

    names(out)[(0*(R+1)+1):(1*(R+1))] = paste0("Fe_Int", 1:(R+1))
    names(out)[(1*(R+1)+1):(2*(R+1))] = paste0("Fe_RMSE", 1:(R+1))
    names(out)[(2*(R+1)+1):(3*(R+1))] = paste0("Su_Int", 1:(R+1))
    names(out)[(3*(R+1)+1):(4*(R+1))] = paste0("Su_RMSE", 1:(R+1))

    return (out)
    }

R = 10
n = 1000    # nobs per time-series
B = 100
uq = c(0.95, 0.96, 0.97, 0.98, 0.99)
out = rep(list(NULL), length(uq))
for (j in 1:length(uq)){
    out[[j]] = foreach(k=1:B, .combine=rbind) %dopar% doit(k, n = n, R = R, uq = uq[1], theta = theta)
    cat("\n")
    }

save("out", file = "./comp.RData")






colMeans(out, na.rm = TRUE)[1:11]
colMeans(out, na.rm = TRUE)[23:33]

plot(colMeans(out, na.rm = TRUE)[1:11], colMeans(out, na.rm = TRUE)[23:33])
abline(0, 1)

plot(colMeans(out, na.rm = TRUE)[12:22], colMeans(out, na.rm = TRUE)[34:44])
abline(0, 1)


# Coverage
mean(apply(fb.v, 1, function(x) x[1] < theta && x[2] > theta))
mean(apply(sb.v, 1, function(x) x[1] < theta && x[2] > theta))

plot(0, type = 'n', xlim = c(1, B), ylim = c(0, 1))
segments(x0 = 1:B, x1 = 1:B, y0 = fb.v[,1], y1 = fb.v[,2])
abline(h = theta)

plot(0, type = 'n', xlim = c(1, B), ylim = c(0, 1))
segments(x0 = 1:B, x1 = 1:B, y0 = sb.v[,1], y1 = sb.v[,2])
abline(h = theta)

plot(0, type = 'n', xlim = c(1, B/1), ylim = c(0, 1))
segments(x0 = 1:B, x1 = 1:B, y0 = fb.v[,1], y1 = fb.v[,2], col = rgb(0,0,1,0.5))
points(apply(fb.v, 1, mean), col = rgb(0,0,1,0.5))
segments(x0 = 1:B+0.0, x1 = 1:B+0.0, y0 = sb.v[,1], y1 = sb.v[,2], col = rgb(1,0,0,0.5))
points(apply(sb.v, 1, mean), col = rgb(1,0,0,0.5))
abline(h = theta)

# RMSE
mean(fb.r)
mean(sb.r)

median(fb.r)
median(sb.r)

plot(density(fb.r))
lines(density(sb.r), col = 'red')

plot(density(apply(fb.v, 1, diff)))
lines(density(apply(sb.v, 1, diff)), col = 'red')

plot(apply(fb.v, 1, diff), apply(sb.v, 1, diff))
cor(apply(fb.v, 1, diff), apply(sb.v, 1, diff))
