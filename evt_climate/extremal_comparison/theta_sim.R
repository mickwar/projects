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

R = 10
n = 1000
B = 200


### Simulation: Frechet
set.seed(1)

theta = 0.10

fb.v = matrix(0, B, 2)
fb.r = double(B)
sb.v = matrix(0, B, 2)
sb.r = double(B)

for (k in 1:B){
    ww = matrix(-1/log(runif(n*R)), n, R)
    y = matrix(0, n, R)
    y[1,] = ww[1,] / theta
    for (i in 2:n)
        y[i,] = apply(rbind((1-theta)*y[i-1,], ww[i,]), 2, max)

    ### For split up sequences
    uu = quantile(y, 0.90)
    # K = 12
    # uu = seq(quantile(y, 0.80), quantile(y, 0.97), length = K)
    N = sapply(uu, function(x) sum(y > x))

    # N = round(seq(1000, 10, length = K))
    # uu = sapply(N, function(x) quantile(y, 1 - x/n))

    up = mean(y <= uu)
        
    ferro_ba = theta_hier(y, uu, likelihood = "ferro",
        nburn = 10000, nmcmc = 10000)
    tmp = c(ferro_ba$keep, R+1)
#   fb.v[k,] = quantile(ferro_ba$mcmc[,which(tmp == R+1)], c(0.025, 0.975))
    fb.v[k,] = hpd_mult(ferro_ba$mcmc[,which(tmp == R+1)], force_uni = TRUE)
    

    # RMSE
    fb.r[k] = sqrt(mean((ferro_ba$mcmc[,which(tmp == R+1)] - theta)^2))

    suveges_ba = theta_hier(y, uu, likelihood = "suveges",
        nburn = 10000, nmcmc = 10000)
    tmp = c(suveges_ba$keep, R+1)
#   sb.v[k,] = quantile(suveges_ba$mcmc[,which(tmp == R+1)], c(0.025, 0.975))
    sb.v[k,] = hpd_mult(suveges_ba$mcmc[,which(tmp == R+1)], force_uni = TRUE)

    # RMSE
    sb.r[k] = sqrt(mean((suveges_ba$mcmc[,which(tmp == R+1)] - theta)^2))

    }


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
