library(mwBASE)
library(mwEVT)

### Simulation: Frechet
set.seed(1)
n = 5000
theta = 0.25
ww = -1/log(runif(n))
y = double(n)
y[1] = ww[1] / theta
for (i in 2:n)
    y[i] = max((1-theta)*y[i-1], ww[i])

plot(y, type = 'l')

### For split up sequences
K = 12
uu = seq(quantile(y, 0.80), quantile(y, 0.97), length = K)
N = sapply(uu, function(x) sum(y > x))

# N = round(seq(1000, 10, length = K))
# uu = sapply(N, function(x) quantile(y, 1 - x/n))

fc.m = double(K)
fc.v = matrix(0, K, 2)

fb.m = double(K)
fb.v = matrix(0, K, 2)

sc.m = double(K)
sc.v = matrix(0, K, 2)

sb.m = double(K)
sb.v = matrix(0, K, 2)

u.vec = double(length(uu))


for (j in 1:K){
    up = mean(y <= uu[j])
    u.vec[j] = up
    
    ferro_cl = theta_uni(y, uu[j], likelihood = "ferro",
        method = "classical")
    fc.m[j] = ferro_cl$theta
    fc.v[j,] = quantile(ferro_cl$bootstrap, c(0.025, 0.975))

    ferro_ba = theta_uni(y, uu[j], likelihood = "ferro",
        method = "bayesian", nburn = 40000, nmcmc = 20000)
    fb.m[j] = mean(ferro_ba$mcmc)
    fb.v[j,] = quantile(ferro_ba$mcmc, c(0.025, 0.975))

    suveges_cl = theta_uni(y, uu[j], likelihood = "suveges",
        method = "classical")
    sc.m[j] = suveges_cl$theta
    sc.v[j,] = quantile(suveges_cl$bootstrap, c(0.025, 0.975))

    suveges_ba = theta_uni(y, uu[j], likelihood = "suveges",
        method = "bayesian", nburn = 40000, nmcmc = 20000)
    sb.m[j] = mean(suveges_ba$mcmc)
    sb.v[j,] = quantile(suveges_ba$mcmc, c(0.025, 0.975))
    }

pdf("./figs/frechet_uni.pdf", width = 8, height = 8)
par(mfrow = c(2, 2))
mains = c("Ferro MLE", "Ferro Bayes", "Suveges MLE", "Suveges Bayes")
cols = c("green", "blue", "red", "orange")
for (j in 1:4){
    mm = switch(j,
        "1" = fc.m,
        "2" = fb.m,
        "3" = sc.m,
        "4" = sb.m)
    vv = switch(j,
        "1" = fc.v,
        "2" = fb.v,
        "3" = sc.v,
        "4" = sb.v)
    plot(uu, mm, type = 'b', pch = 15, col = col_mult(cols[j], 'gray75'),
        xlab = "Threshold", ylab = expression(theta), bty = 'n',
        ylim = c(0, 1))
    lines(uu, vv[,1], col = col_fade(cols[j], 0.75))
    lines(uu, vv[,2], col = col_fade(cols[j], 0.75))
    title(main = mains[j], line = 2.5)
    tmp = axTicks(1)
#   axis(3, at = tmp, labels = sapply(tmp, function(x) sum(y > x))) 
    axis(3, at = tmp, labels = round(sapply(tmp, function(x) mean(y < x)), 3) )
    rm(tmp)
    if (exists("theta"))
        abline(h = theta, lty = 2)
    }
par(mfrow = c(1, 1))
dev.off()
