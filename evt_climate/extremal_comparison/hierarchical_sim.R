library(mwBASE)
library(mwEVT)

### Simulation study:
#   Calculate coverage (can I take 10000 samples, divide into 
#   ten 1000 segments and compute 10 credible bounds?)
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


### Simulation: Frechet
for (theta in c(0.15, 0.5, 0.85)){
    for (R in c(5, 10, 20)){
        for (n in c(250, 500, 1000)){
            set.seed(1)
            ww = matrix(-1/log(runif(n*R)), n, R)
            y = matrix(0, n, R)
            y[1,] = ww[1,] / theta
            for (i in 2:n)
                y[i,] = apply(rbind((1-theta)*y[i-1,], ww[i,]), 2, max)

            ### For split up sequences
            K = 12
            uu = seq(quantile(y, 0.80), quantile(y, 0.97), length = K)
            N = sapply(uu, function(x) sum(y > x))

            # N = round(seq(1000, 10, length = K))
            # uu = sapply(N, function(x) quantile(y, 1 - x/n))

            fb.m = matrix(0, K, R+1)
            fb.v = matrix(0, K, 2)

            sb.m = matrix(0, K, R+1)
            sb.v = matrix(0, K, 2)

            u.vec = double(length(uu))


            for (j in 1:K){
                up = mean(y <= uu[j])
                u.vec[j] = up
                
                ferro_ba = theta_hier(y, uu[j], likelihood = "ferro",
                    nburn = 40000, nmcmc = 20000)
                tmp = c(ferro_ba$keep, R+1)
                fb.m[j, tmp] = colMeans(ferro_ba$mcmc)
                fb.m[j, -tmp] = NA
                fb.v[j,] = quantile(ferro_ba$mcmc[,which(tmp == R+1)], c(0.025, 0.975))

                suveges_ba = theta_hier(y, uu[j], likelihood = "suveges",
                    nburn = 40000, nmcmc = 20000)
                tmp = c(suveges_ba$keep, R+1)
                sb.m[j, tmp] = colMeans(suveges_ba$mcmc)
                sb.m[j, -tmp] = NA
                sb.v[j,] = quantile(suveges_ba$mcmc[,which(tmp == R+1)], c(0.025, 0.975))
                }

            pdf(paste0("./figs/sim_frechet_hier_", as.character(theta),"_",n*R,"_",R, ".pdf"),
                width = 12, height = 8)
            par(mfrow = c(1, 2))
            mains = c("Ferro Bayes", "Suveges Bayes")
            cols = c("blue", "orange")
            for (j in 1:2){
                mm = switch(j,
                    "1" = fb.m,
                    "2" = sb.m)
                vv = switch(j,
                    "1" = fb.v,
                    "2" = sb.v)
                plot(uu, mm[,R+1], type = 'b', pch = 15, col = col_mult(cols[j], 'gray75'),
                    xlab = "Threshold", ylab = expression(theta), bty = 'n',
                    ylim = c(0, 1))
                matplot(uu, mm[,-(R+1)], add = TRUE, type = 'l',
                    col = col_fade(cols[j], 0.25), lty = 1)
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
            }
        }
    }
