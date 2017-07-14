library(mwBASE)
library(mwEVT)
library(MASS)

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

    # Ferro
    ferro_ba = theta_hier(y, uu, likelihood = "ferro",
        nburn = 10000, nmcmc = 10000, display = 0, flag = 1e6)
    tmp1 = c(ferro_ba$keep, R+1)
    fb.v = apply(ferro_ba$mcmc, 2, hpd_mult, force_uni = TRUE)
    fb.r = apply(ferro_ba$mcmc, 2, function(x) sqrt(mean((x-theta)^2)))

    # Suveges K = 1
    suveges_ba = theta_hier(y, uu, likelihood = "suveges", K = 1,
        nburn = 10000, nmcmc = 10000, display = 0, flag = 1e6)
    tmp2 = c(suveges_ba$keep, R+1)
    sb1.v = apply(suveges_ba$mcmc, 2, hpd_mult, force_uni = TRUE)
    sb1.r = apply(suveges_ba$mcmc, 2, function(x) sqrt(mean((x-theta)^2)))

    # Suveges K = 5
    suveges_ba = theta_hier(y, uu, likelihood = "suveges", K = 5,
        nburn = 10000, nmcmc = 10000, display = 0, flag = 1e6)
    tmp3 = c(suveges_ba$keep, R+1)
    sb5.v = apply(suveges_ba$mcmc, 2, hpd_mult, force_uni = TRUE)
    sb5.r = apply(suveges_ba$mcmc, 2, function(x) sqrt(mean((x-theta)^2)))

    out = NA*double(6*(R+1))
    out[tmp1 + 0*(R+1)] = (fb.v[1,] < theta & fb.v[2,] > theta)
    out[tmp1 + 1*(R+1)] = fb.r
    out[tmp2 + 2*(R+1)] = (sb1.v[1,] < theta & sb1.v[2,] > theta)
    out[tmp2 + 3*(R+1)] = sb1.r
    out[tmp3 + 4*(R+1)] = (sb5.v[1,] < theta & sb5.v[2,] > theta)
    out[tmp3 + 5*(R+1)] = sb5.r

    names(out)[(0*(R+1)+1):(1*(R+1))] = paste0("Fe_Int", 1:(R+1))
    names(out)[(1*(R+1)+1):(2*(R+1))] = paste0("Fe_RMSE", 1:(R+1))
    names(out)[(2*(R+1)+1):(3*(R+1))] = paste0("Su1_Int", 1:(R+1))
    names(out)[(3*(R+1)+1):(4*(R+1))] = paste0("Su1_RMSE", 1:(R+1))
    names(out)[(4*(R+1)+1):(5*(R+1))] = paste0("Su5_Int", 1:(R+1))
    names(out)[(5*(R+1)+1):(6*(R+1))] = paste0("Su5_RMSE", 1:(R+1))

    return (out)
    }

for (theta in c(0.5, 0.9)){
    lab = as.character(round(theta*100))
    R = 10
    n = 1000    # nobs per time-series
    B = 500
    uq = c(0.95, 0.96, 0.97, 0.98, 0.99)
    out = rep(list(NULL), length(uq))

    for (j in length(uq):1){
        out[[j]] = foreach(k=1:B, .combine=rbind, .errorhandling = "remove") %dopar%
            par_func(k, n = n, R = R, uq = uq[j], theta = theta)
        cat("\n")
        }

    save("out", file = paste0("./comp_", lab, ".RData"))
    }



pdf(paste0("./figs/sim_coverage_", lab, ".pdf"), width = 6, height = 5)
plot(0, type = 'n', xlim = range(uq)+c(-0.0010,0.0010), ylim = c(0, 1), bty = 'n',
    main = "Coverage", ylab = "Probability", xlab = "Threshold quantile")
abline(h = 0.95, lty = 2)
title(main = bquote(paste(theta, " = ", .(theta))), line = 0.5, cex.main = 1.3)
offset = c(-0.0010, 0, 0.0010)
for (j in 1:R){
    points(uq+offset[1], sapply(out, function(x) mean(x[,j], na.rm = TRUE)),
        col = 'lightblue', pch = 15)
    points(uq+offset[2], sapply(out, function(x) mean(x[,j+2*(R+1)], na.rm = TRUE)),
        col = 'pink', pch = 16)
    points(uq+offset[3], sapply(out, function(x) mean(x[,j+4*(R+1)], na.rm = TRUE)),
        col = 'lightgreen', pch = 16)
    }
points(uq+offset[1], sapply(out, function(x) mean(x[,R+1], na.rm = TRUE)),
    col = 'darkblue', pch = 15)
points(uq+offset[2], sapply(out, function(x) mean(x[,3*(R+1)], na.rm = TRUE)),
    col = 'red', pch = 16)
points(uq+offset[3], sapply(out, function(x) mean(x[,5*(R+1)], na.rm = TRUE)),
    col = 'darkgreen', pch = 16)
legend(uq[1], 0.35, legend = c("Ferro", "Suveges, K=1", "Suveges, K=5"),
    col = c("darkblue", "red", "darkgreen"), pch = c(15, 16, 16), bty = 'n', cex = 1.3)
dev.off()


pdf(paste0("./figs/sim_rmse_", lab, ".pdf"), width = 6, height = 5)
ylim = range(sapply(out, function(x) colMeans(x[,grep("RMSE", colnames(x))], na.rm = TRUE)))
plot(0, type = 'n', xlim = range(uq)+c(-0.0010,0.0010), ylim = ylim, bty = 'n',
    main = "RMSE", ylab = "RMSE", xlab = "Threshold quantile")
title(main = bquote(paste(theta, " = ", .(theta))), line = 0.5, cex.main = 1.3)
offset = c(-0.0010, 0, 0.0010)
for (j in 1:R){
    points(uq+offset[1], sapply(out, function(x) mean(x[,j+1*(R+1)], na.rm = TRUE)),
        col = 'lightblue', pch = 15)
    points(uq+offset[2], sapply(out, function(x) mean(x[,j+3*(R+1)], na.rm = TRUE)),
        col = 'pink', pch = 16)
    points(uq+offset[3], sapply(out, function(x) mean(x[,j+5*(R+1)], na.rm = TRUE)),
        col = 'lightgreen', pch = 16)
    }
points(uq+offset[1], sapply(out, function(x) mean(x[,2*(R+1)], na.rm = TRUE)),
    col = 'darkblue', pch = 15)
points(uq+offset[2], sapply(out, function(x) mean(x[,4*(R+1)], na.rm = TRUE)),
    col = 'red', pch = 16)
points(uq+offset[3], sapply(out, function(x) mean(x[,6*(R+1)], na.rm = TRUE)),
    col = 'darkgreen', pch = 16)
legend(uq[1], ylim[2], legend = c("Ferro", "Suveges, K=1", "Suveges, K=5"),
    col = c("blue", "red", "darkgreen"), pch = c(15, 16, 16), bty = 'n', cex = 1.3)
dev.off()


