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


par_func = function(k, R, n, uq, theta){

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

theta = 0.9
R = 10
n = 1000    # nobs per time-series
B = 200
uq = c(0.95, 0.96, 0.97, 0.98, 0.99)
out = rep(list(NULL), length(uq))

for (j in length(uq):1){
    out[[j]] = foreach(k=1:B, .combine=rbind) %dopar% par_func(k, n = n, R = R,
        uq = uq[j], theta = theta)
    cat("\n")
    }

save("out", file = "./comp_90.RData")




pdf("./figs/sim_coverage_90.pdf", width = 6, height = 4)
plot(0, type = 'n', xlim = range(uq)+c(-0.0010,0.0010), ylim = c(0, 1), bty = 'n',
    main = "Coverage", ylab = "Probability", xlab = "Threshold quantile")
abline(h = 0.95, lty = 2)
title(main = bquote(paste(theta, " = 0.9")), line = 0.5, cex.main = 1.3)
for (j in 1:R){
    points(uq-0.0005, sapply(out, function(x) mean(x[,j], na.rm = TRUE)),
        col = 'lightblue', pch = 15)
    points(uq+0.0005, sapply(out, function(x) mean(x[,j+2*(R+1)], na.rm = TRUE)),
        col = 'pink', pch = 16)
    }
points(uq-0.0005, sapply(out, function(x) mean(x[,R+1], na.rm = TRUE)),
    col = 'blue', pch = 15)
points(uq+0.0005, sapply(out, function(x) mean(x[,3*(R+1)], na.rm = TRUE)),
    col = 'red', pch = 16)
legend(uq[1], 0.3, legend = c("Ferro", "Suveges"), col = c("blue", "red"),
    pch = c(15, 16), bty = 'n', cex = 1.3)
dev.off()


pdf("./figs/sim_rmse_90.pdf", width = 6, height = 4)
ylim = range(sapply(out, function(x) colMeans(x[,grep("RMSE", colnames(x))], na.rm = TRUE)))
plot(0, type = 'n', xlim = range(uq)+c(-0.0010,0.0010), ylim = ylim, bty = 'n',
    main = "RMSE", ylab = "RMSE", xlab = "Threshold quantile")
title(main = bquote(paste(theta, " = 0.9")), line = 0.5, cex.main = 1.3)
for (j in 1:R){
    points(uq-0.0005, sapply(out, function(x) mean(x[,j+1*(R+1)], na.rm = TRUE)),
        col = 'lightblue', pch = 15)
    points(uq+0.0005, sapply(out, function(x) mean(x[,j+3*(R+1)], na.rm = TRUE)),
        col = 'pink', pch = 16)
    }
points(uq-0.0005, sapply(out, function(x) mean(x[,2*(R+1)], na.rm = TRUE)),
    col = 'blue', pch = 15)
points(uq+0.0005, sapply(out, function(x) mean(x[,4*(R+1)], na.rm = TRUE)),
    col = 'red', pch = 16)
legend(uq[1], ylim[2], legend = c("Ferro", "Suveges"), col = c("blue", "red"),
    pch = c(15, 16), bty = 'n', cex = 1.3)
dev.off()


