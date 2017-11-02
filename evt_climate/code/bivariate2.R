library(mwBASE)
rgpd = function(n, mu, ksi, sigma){
    if (length(ksi) == 1)
        ksi = rep(ksi, n)
    return (ifelse(ksi == 0, mu - sigma * log(runif(n)),
        mu + sigma * (runif(n)^(-ksi) - 1) / ksi))
    }

# Rescale to [0, 1]
rescale = function(x)
    (x - min(x)) / diff(range(x))

# Empirical CDF transform to (0, 1)
get.ecdf = function(x, n){
    if (missing(n))
        n = length(x)
    ((1:n) / (n + 1))[rank(x)]
    }
angle.2.cone = function(theta){
    out = matrix(1, length(theta), 2)
    out[,1] = tan(pi/4 - abs(theta - pi/4))
    ind = which(theta <= pi/4)
    out[ind,] = out[ind, c(2,1)]
    return (out)
    }
logsum = function(x, log = TRUE){
    w = which.max(x)
    if (log){
        return (x[w] + log(1 + sum(exp(x[-w] - x[w]))))
    } else {
        return (log(x[w]) + log(1 + sum(exp(log(x[-w]) - log(x[w])))))
        }
    }


DIR = "./RData/2017_09_08/"
OUT_DIR = "./figs/2017_09_08/"
files = list.files(DIR)
ind.control = grep("control", files)
ind.decadal = grep("decadal", files)
ind.historical = grep("historical", files)
ind.obs = grep("obs", files)
ind.62 = grep("1962", files)
ind.90 = grep("1990", files)
ind.cal = grep("california", files)
ind.usa = grep("usa", files)
ind.pr = grep("pr", files)
ind.tasmax = grep("tasmax", files)
ind.winter = grep("winter", files)
ind.summer = grep("summer", files)
labs = unlist(strsplit(files, ".RData"))

x = rep(list(NULL), length(labs))
names(x) = labs
for (i in 1:length(labs)){
    load(paste0(DIR, files[i]))
    x[[labs[i]]] = out
    }

lab.year = c("1962", "1990")
lab.reg = c("CA", "USA")
lab.var = c("Pr", "Tasmax")
lab.sea = c("Winter", "Summer")
lab.mod = c("Control", "Decadal", "Historical", "Observations")
lab.shortmod = c("Cont.", "Dec.", "Hist.", "Obs.")
cols = c("blue", "green", "red", "gray50")
# 
# tab.theta.mm = sapply(x, function(y) tail(y$theta.mean, 1))
# tab.theta.qq = t(sapply(x, function(y) c(tail(t(y$theta.ints), 1))))
# tab.threshold = sapply(x, function(y) y$threshold[1])
# tab.T_C = sapply(x, function(y) mean(y$T_C))


wrapper = function(x, obs.ind, sim.ind, zero, one, parts, nburn = 10000, nmcmc = 10000){

    R = ncol(x[[sim.ind]]$varmat)
    n = nrow(x[[obs.ind]]$varmat)

    new.xy = matrix(0, n, R+1)

    ex = x[[obs.ind]]$varmat[,1] - x[[obs.ind]]$threshold
    ksi = mean(x[[obs.ind]]$params[,1])
    sig = mean(x[[obs.ind]]$params[,2])

    new.xy[, 1] = (1 + ksi * ex / sig)^(1/ksi)
    new.xy[is.na(new.xy[,1]),1] = 0

    for (i in 1:R){

        ex = x[[sim.ind]]$varmat[,i] - x[[sim.ind]]$threshold[i]
        ksi = mean(x[[sim.ind]]$params[,i])
        sig = mean(x[[sim.ind]]$params[,i + R])

        new.xy[, i+1] = (1 + ksi * ex / sig)^(1/ksi)
        new.xy[is.na(new.xy[,i+1]),i+1] = 0

        }


    # zero.ind = which(apply(new.xy, 1, function(x) any(x == 0)))
    # ex.ind = which(new.xy[,1] > 1 | new.xy[,2] > 1)

    # ind = setdiff(ex.ind, zero.ind)

    de.ind = rep(list(NULL), R)
    ind = rep(list(NULL), R)
    v = rep(list(NULL), R)
    cone = rep(list(NULL), R)
    phi = rep(list(NULL), R)
    phi.s = rep(list(NULL), R)

    for (i in 1:R){

        # Declustered indices for each margin, joined
        tmp.ind1 = which((x[[obs.ind]]$varmat - x[[obs.ind]]$threshold) %in% x[[obs.ind]]$y)
        tmp.ind2 = which((x[[sim.ind]]$varmat[,i] - x[[sim.ind]]$threshold) %in% x[[sim.ind]]$y[[i]])
        ind[[i]] = sort(unique(c(tmp.ind1, tmp.ind2)))

        # All exceedances
#       ind[[i]] = which(new.xy[,1] > 1 | new.xy[,i+1] > 1)

        v[[i]] = apply(new.xy[ind[[i]], c(1, i+1)], 1, max)
        cone[[i]] = new.xy[ind[[i]], c(1, i+1)] / v[[i]]
        phi[[i]] = atan(new.xy[ind[[i]], i+1] / new.xy[ind[[i]], 1])
        phi.s[[i]] = phi[[i]] * 2 / pi
        }

    phi.s = unlist(phi.s)

    make.post = function(zero = 1, one = 1, parts = 2){
        if (!(zero %in% c(0, 1)))
            stop("zero must be in c(0, 1). (zero = indicator for point mass at zero)")
        if (!(one %in% c(0, 1)))
            stop("one must be in c(0, 1). (one = indicator for point mass at one)")
        if (floor(parts) != parts || parts < 1)
            stop("parts must be a positive integer. (parts = # of beta mixtures)")
        nparam = zero + one + 3*parts - 1
        f = function(dat, params){
            # Beta mixture parameters
            ab.params = exp(head(params, 2*parts))

            # Probability parameters. If false, there is no mixture
            if (nparam > 2*parts){
                prob.params = exp(tail(params, nparam - 2*parts))
                prob.params = c(prob.params, 1 - sum(prob.params))
            } else {
                prob.params = 1
                }

            # Check if (a_i + b_i) < (a_{i+1} + b_{i+1}) to hopefully
            # prevent parameters from swapping labels
            if (parts > 1){
                for (i in 1:(parts-1)){
                    if (sum(ab.params[c(2*(i + 0) - 1, 2*(i + 0))]) <
                        sum(ab.params[c(2*(i + 1) - 1, 2*(i + 1))])){
                        return (-Inf)
                        }
                    }
                }

            # Check the sum-to-one probability constraint
            if (tail(prob.params, 1) < 0)
                return (-Inf)

            # Initialize
            out = 0
            p.ind = 1

            # Remove zeros and ones (for beta parts)
            y = dat[!(dat %in% c(0, 1))]

            # Add zero and one point mass contributions to likelihood
            if (zero == 1){
                zeros = which(dat == 0)
                out = out + length(zeros) * log(prob.params[p.ind])
                p.ind = p.ind + 1
                }
            if (one == 1){
                ones = which(dat == 1)
                out = out + length(ones) * log(prob.params[p.ind])
                p.ind = p.ind + 1
                }

            # Add beta mixture
            tmp = double(length(y))
            for (i in 1:parts){
                tmp = tmp + prob.params[p.ind]*dbeta(y, ab.params[2*i-1], ab.params[2*i])
                p.ind = p.ind + 1
                }
            out = out + sum(log(tmp))


            # Priors on beta parameters
            out = out - sum(log(ab.params))

            # With no ``priors'' for prob.params, this is equivalent
            # to a Dirichlet(1, 1, ..., 1) prior, so long as the support
            # is included

            # Add the Jacobian for the sampler (in this case, everything is log'd)
            out = out + sum(params)

            return (out)
            }
        return (f)
        }

    if (missing(zero))
        zero = any(phi.s == 0)*1
    if (missing(one))
        one = any(phi.s == 1)*1
    if (missing(parts))
        parts = 2

    calc.post = make.post(zero, one, parts)
    nparam = zero + one + 3*parts - 1
    if (nparam == 2*parts){
        mcmc = mcmc_sampler(phi.s, calc.post, nparam = nparam, nburn = nburn, nmcmc = nmcmc,
            chain_init = rep(0, 2*parts))
    } else {
        mcmc = mcmc_sampler(phi.s, calc.post, nparam = nparam, nburn = nburn, nmcmc = nmcmc,
            chain_init = c(rep(0, 2*parts),
                rep(log(1/(zero + one + parts)), zero + one + parts - 1)))
        }

    mcmc$params = exp(mcmc$params)
    if (nparam > 2*parts){
        mcmc$params = cbind(mcmc$params,
            1-apply(as.matrix(mcmc$params[,(2*parts+1):nparam]), 1, sum))
        }

    out = list("xy" = new.xy, "params" = mcmc$params, "phi.s" = phi.s, "v" = v,
        "cone" = cone, "ind" = ind,
        "zero" = zero, "one" = one, "parts" = parts, "nmcmc" = nmcmc)

    return (out)
    }

get.preds = function(x){
    nparam = x$zero + x$one + 3*x$parts
    nprob = nparam - 2*x$parts

    if (nprob == 1)
        return (rbeta(x$nmcmc, x$params[,1], x$params[,2]))

    ind = c(t(apply(x$params[,(2*x$parts+1):NCOL(x$params)], 1,
        function(r) sample(nprob, 1, prob = r))))

    p.ind = 1
    out = double(x$nmcmc)
    if (x$zero){
        out[which(ind == p.ind)] = 0
        p.ind = p.ind + 1
        }
    if (x$one){
        out[which(ind == p.ind)] = 1
        p.ind = p.ind + 1
        }
    for (i in 1:x$parts){
        ind2 = which(ind == p.ind)
        out[ind2] = rbeta(length(ind2), x$params[ind2,2*i-1], x$params[ind2,2*i])
        p.ind = p.ind + 1
        }
    return (out)
    }

get.preds2 = function(x, n){
    nparam = x$zero + x$one + 3*x$parts
    nprob = nparam - 2*x$parts

    ind = apply(x$params[,(2*x$parts+1):NCOL(x$params)], 1,
        function(r) sample(nprob, n, prob = r, replace = TRUE))

    out = matrix(0, x$nmcmc, n)

    for (j in 1:n){
        p.ind = 1
        if (x$zero){
            out[which(ind[j,] == p.ind),j] = 0
            p.ind = p.ind + 1
            }
        if (x$one){
            out[which(ind[j,] == p.ind),j] = 1
            p.ind = p.ind + 1
            }
        for (i in 1:x$parts){
            ind2 = which(ind[j,] == p.ind)
            out[ind2,j] = rbeta(length(ind2), x$params[ind2,2*i-1], x$params[ind2,2*i])
            p.ind = p.ind + 1
            }
        }
        
    return (out)
    }

calc.chi.pred = function(x, n){

    }





nburn = 20000
nmcmc = 80000

R = 10  # Number of replicates
width = 9
height = 9

# rep.num = 1
# comp = c(18, 34)
# comp = seq(1, 47, by = 2)
# comp = c(17, 33)
# comp = c(1, 3)
# n = nrow(x[[comp[1]]]$varmat)

pproc = rep(list(NULL), length(ind.control) + length(ind.decadal) + length(ind.historical))

B5 = list(ind.control, ind.decadal, ind.historical)
for (A1 in 1:2){
    B1 = list(ind.cal, ind.usa)[[A1]]
    for (A2 in 1:2){
        B2 = list(ind.62, ind.90)[[A2]]
        for (A3 in 1:2){
            B3 = list(ind.pr, ind.tasmax)[[A3]]
            for (A4 in 1:2){
                B4 = list(ind.winter, ind.summer)[[A4]]
                for (A5 in 1:3){
                    rr.obs = Reduce(intersect, list(B1, B2, B3, B4, ind.obs))
                    rr.sim = Reduce(intersect, list(B1, B2, B3, B4, B5[[A5]]))

#                   if (rr.sim == 2)
#                       next

                    pproc[[rr.sim]] = wrapper(x, rr.obs, rr.sim, parts = 2,
                        nburn = nburn, nmcmc = nmcmc)

#                   pred = get.preds2(pproc[[rr.sim]], 10)
                    pred = get.preds(pproc[[rr.sim]])
                    pred.cone = angle.2.cone(pred * pi/2)

                    pproc[[rr.sim]]$pred = pred

                    ind = (pred.cone[,1] > 0 & pred.cone[,2] > 0)
                    z1 = pred.cone[ind, 1] / mean(pred.cone[, 1])
                    z2 = pred.cone[ind, 2] / mean(pred.cone[, 2])
                    z = cbind(z1, z2)
                    m1 = apply(z, 1, min)
                    pproc[[rr.sim]]$chi = mean(m1) * mean(ind)

                    hist(pproc[[rr.sim]]$phi.s, breaks = 50, col = 'black', border = 'gray50',
                        freq = FALSE, ylab = "", xlab = expression(phi), main = "Angle",
                        cex.lab = 1.5, cex.main = 2)
                    hist(pred, col = rgb(34/255, 139/255, 34/255, 0.5), breaks = 50, freq = FALSE,
                        border = 'lightgreen', add = TRUE)
                    abline(v = pproc[[rr.sim]]$chi)

                    }
                }
            }
        }
    }

u = seq(0.8, 0.995, length = 20)
chi.u = matrix(0, length(u), 50)
chi.bar.u = matrix(0, length(u), 50)
chi = double(NCOL(chi.u))
for (i in 1:NCOL(chi.u)){

    pred = get.preds(pproc[[rr.sim]])
    pred.cone = angle.2.cone(pred * pi/2)
    pred.v = 1/runif(nrow(pred.cone))
    pred.y = pred.v * pred.cone

    z1 = rank(pred.y[,1]) / nmcmc
    z2 = rank(pred.y[,2]) / nmcmc

    for (j in 1:length(u)){
        ind1 = z1 <= u[j]
        ind2 = z2 <= u[j]
        chi.u[j,i] = 2 - log(mean(ind1 & ind2)) / log(mean(ind2))

        ind1 = z1 > u[j]
        ind2 = z2 > u[j]
        chi.bar.u[j,i] = 2 * log(mean(ind1)) / log(mean(ind1 & ind2)) - 1
        }

    ind = (pred.cone[,1] > 0 & pred.cone[,2] > 0)
    z1 = pred.cone[ind, 1] / mean(pred.cone[, 1])
    z2 = pred.cone[ind, 2] / mean(pred.cone[, 2])
    z = cbind(z1, z2)
    m1 = apply(z, 1, min)
    chi[i] = mean(m1) * mean(ind)
    }

#plot(u, apply(chi.u, 1, mean), ylim = c(0, 1))
plot(u, apply(chi.u, 1, mean), ylim = range(chi.u))
lines(u, apply(chi.u, 1, quantile, 0.025))
lines(u, apply(chi.u, 1, quantile, 0.975))
points(u[length(u)], mean(chi), col = 'red')
points(rep(u[length(u)], 2), quantile(chi, c(0.025, 0.975)), col = 'red')


plot(u, apply(chi.bar.u, 1, mean), ylim = range(chi.bar.u))
lines(u, apply(chi.bar.u, 1, quantile, 0.025))
lines(u, apply(chi.bar.u, 1, quantile, 0.975))







cor(pred.y)

plot(sapply(pproc, function(x) x$chi))
sapply(pproc, function(x) x$chi)

# save(pproc, file = "./biv_pred.RData")



par(mfrow = c(4, 4), mar = c(0,0,0,0), oma = c(6,10,10,6))
pnum = 0
xlim = range(sapply(pproc, function(x) log(x$chi / (1 - x$chi))))
ylim = c(0.5, 3.5)
B5 = list(ind.control, ind.decadal, ind.historical)
for (A1 in 1:2){
    B1 = list(ind.cal, ind.usa)[[A1]]
    for (A2 in 1:2){
        B2 = list(ind.62, ind.90)[[A2]]
        for (A3 in 1:2){
            B3 = list(ind.pr, ind.tasmax)[[A3]]
            for (A4 in 1:2){
                B4 = list(ind.winter, ind.summer)[[A4]]
                plot(0, type='n', xlim = xlim, ylim = ylim, axes = FALSE)
                pnum = pnum + 1
                for (A5 in 1:3){
                    rr = Reduce(intersect, list(B1, B2, B3, B4, B5[[A5]]))
                    points(log(pproc[[rr]]$chi / (1 - pproc[[rr]]$chi)), A5, pch = 15, col = cols[A5])
                    }

                if ((pnum <= 4) && (pnum %% 2 == 1))
                    axis(3, lwd = -1, lwd.ticks = 1)
                if ((pnum >= 13) && (pnum %% 2 == 0))
                    axis(1, lwd = -1, lwd.ticks = 1)
                if (pnum %% 4 == 0)
                    axis(side = 4, labels = lab.shortmod[-4], at = 1:3, las = 1,
                        lwd = -1, lwd.ticks = 1)
                box()

                # Put in thicker lines in the middle
                if ((A1 == 1) && (A2 == 2))
                    axis(1, labels = FALSE, lwd.ticks = 0, lwd = 5,
                        at = xlim + c(-1, 1)*diff(xlim)/5)
                if ((A3 == 1) && (A4 == 2))
                    axis(4, labels = FALSE, lwd.ticks = 0, lwd = 5,
                        at = ylim + c(-1, 1)*diff(ylim)/5)
                }
            }
        }
    }
cex = 1.25
mtext(lab.var, side = 3, outer = TRUE, at = c(0.25, 0.75),
    line = 7, cex = cex)
mtext(lab.sea, side = 3, outer = TRUE, at = c(0.125, 0.375, 0.625, 0.875),
    line = 4, cex = cex)
mtext(lab.reg, side = 2, outer = TRUE, at = rev(c(0.25, 0.75)),
    line = 9, las = 1, cex = cex, adj = 0)
mtext(lab.year, side = 2, outer = TRUE, at = rev(c(0.125, 0.375, 0.625, 0.875)),
    line = 5, las = 1, cex = cex, adj = 0)
mtext(expression(chi), outer = TRUE, side = 3, at = -0.120, line = 4, cex = 2.5, adj = 0)
#par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0,0,0,0))






























## par(mfrow = c(3,4))
## for (i in 1:10){
##     comp = c(17, 33)
## #   comp = c(1, 3)
## #   comp = c(17, 3)
##     out = wrapper(x, comp, rep.num = i, zero = 1, one = 1, parts = 2,
##         nburn = nburn, nmcmc = nmcmc)
## 
##     pred = get.preds(out)
##     pred.cone = angle.2.cone(pred * pi/2)
## 
##     cone = out$cone
##     v = out$v
##     phi.s = out$phi.s
##     new.xy = out$xy
##     ex.ind = which(new.xy[,1] > 1 | new.xy[,2] > 1)
## 
##     hist(phi.s, breaks = 50, col = 'black', border = 'gray50', freq = FALSE, ylab = "",
##         xlab = expression(phi), main = "Angle", cex.lab = 1.5, cex.main = 2)
##     hist(pred, col = rgb(34/255, 139/255, 34/255, 0.5), breaks = 50, freq = FALSE,
##         border = 'lightgreen', add = TRUE)
## 
##     cat(i, out$zero, out$one, "\n")
##     }
## 
## # v = apply(new.xy[ind,], 1, max)
## # cone = new.xy[ind,] / v
## # phi.s = atan(new.xy[ind,2] / new.xy[ind,1]) * 2/pi
## 
## pred = get.preds(out)
## pred.cone = angle.2.cone(pred * pi/2)
## cone = out$cone
## v = out$v
## phi.s = out$phi.s
## new.xy = out$xy
## 
## get.clusters = function(ex.times, runs){
##     cluster = rep(list(NULL), length(ex.times))
##     flag = 1
##     for (i in 1:(length(ex.times)-1)){
##         if (abs(ex.times[i] - ex.times[i+1]) <= runs){
##             cluster[[flag]] = c(cluster[[flag]], ex.times[i])
##         } else {
##             cluster[[flag]] = c(cluster[[flag]], ex.times[i])
##             flag = flag + 1
##             }
##         }
##     if (tail(diff(ex.times), 1) <= runs){
##         cluster[[flag]] = c(cluster[[flag]], tail(ex.times, 1))
##     } else {
##         cluster[[flag+1]] = tail(ex.times, 1)
##         }
##     cluster = cluster[!sapply(cluster, is.null)]
##     return (cluster)
##     }
## 
## ex.times = which(new.xy[,1] > 1 | new.xy[,2] > 1)
## 
## plot(0, type = 'n', xlim = range(f(new.xy[,1])), ylim = range(f(new.xy[,-1])))
## for (i in 2:(R+1))
##     points(f(new.xy[,c(1, i)]), pch = 16, col = rainbow(R)[i])
## 
## cluster = get.clusters(ex.times, runs = 1)
## 
## 
## f = function(x) x^(1/3)
## f = function(x) log(abs(x))
## 
## pdf("~/multi_clusters.pdf", width = width, height = height)
## par(mfrow = c(1,1))
## plot(f(new.xy[ex.times,]), pch = 16, col = 'gray50', xlab = "cube root of obs",
##     ylab = "cube root of sims", main = "some clusters connected with lines")
## flag.ind = sapply(cluster, head, 1)
## for (i in 1:10){
##     if (length(cluster[[i]]) > 1){
##         for (j in 1:(length(cluster[[i]])-1)){
##             segments(x0 = f(new.xy[cluster[[i]][j], 1]), x1 = f(new.xy[cluster[[i]][j+1], 1]),
##                      y0 = f(new.xy[cluster[[i]][j], 2]), y1 = f(new.xy[cluster[[i]][j+1], 2]),
##                 col = i %% 4 + 1, lty = i %% 5 + 1, lwd = 2)
##             }
##         }
##     }
## lines(f(c(0, 1, 1)), f(c(1, 1, 0)), col = 'black', lwd = 2)
## dev.off()
## 
## 
## 
## pdf("~/four_betas.pdf", width = width, height = height)
## hist(unlist(phi.s), breaks = 50, col = 'black', border = 'gray50', freq = FALSE, ylab = "",
##     xlab = expression(phi), main = "Angle", cex.lab = 1.5, cex.main = 2)
## hist(pred, col = rgb(34/255, 139/255, 34/255, 0.5), breaks = 50, freq = FALSE,
##     border = 'lightgreen', add = TRUE)
## for (i in 1:4){
##     g = function(x)
##         dbeta(x, mean(out$params[,2*i-1]), mean(out$params[,2*i]))*mean(out$params[,8+i])
##     curve(g(x), add =TRUE, col = rainbow(4)[i], lwd = 3)
##     }
## legend("topright", col = rainbow(4), legend = 1:4, lwd = 3, lty = 1, box.lty = 0)
## dev.off()
## 
## 0.15 - 0.05*0.15
## 0.05 - 0.05*0.15
## 0.05*0.15
## 
## 
## colMeans(out$params)
## colMeans(out$params[,9:12])
## 
## 
## nu = length(ex.times)
## 
## 1 - (0.85 * 0.95)
## nu / n
## 
## (nu / n)
## 
## h = sample(c(0, 1), n, replace = TRUE, prob = c(0.85*0.95, 1 - (0.85 * 0.95)))
## mean(diff(which(h == 1)) == 1)
## 
## B = 50000
## boot.cluster = matrix(0, B, max(sapply(cluster, length)))
## for (i in 1:B){
##     h = sample(c(0, 1), n, replace = TRUE, prob = c(0.85*0.95, 1 - (0.85 * 0.95)))
##     tmp = get.clusters(which(h == 1), runs = 1)
##     tmp = sapply(tmp, length)
##     for (j in 1:NCOL(boot.cluster)){
##         boot.cluster[i, j] = mean(tmp >= j)
##         }
##     }
## 
## apply(boot.cluster, 2, mean)
## apply(boot.cluster, 2, quantile, c(0.025, 0.975))
## 
## for (j in 1:NCOL(boot.cluster)){
##     plot(density(boot.cluster[,j]), xlim = c(0, 1))
##     abline(v = mean(sapply(cluster, length) >= j), col = 'blue', lwd = 3)
##     readline()
##     }
## 
## mean(sapply(cluster, length) >= 2)
## 
## 1-0.85*0.95
## 
## 
## tmp.ind1 = which((x[[50]]$varmat - x[[50]]$threshold) %in% x[[50]]$y)
## tmp.ind2 = which((x[[2]]$varmat[,1] - x[[2]]$threshold) %in% x[[2]]$y[[1]])
## tmp.ind = sort(unique(c(tmp.ind1, tmp.ind2)))
## par(mfrow = c(1,1))
## plot(f(new.xy[ex.times,]), pch = 16, col = 'gray50', xlab = "cube root of obs",
##     ylab = "cube root of sims", main = "some clusters connected with lines")
## flag.ind = sapply(cluster, head, 1)
## for (i in 1:length(cluster)){
##     if (length(cluster[[i]]) > 1){
##         for (j in 1:(length(cluster[[i]])-1)){
##             segments(x0 = f(new.xy[cluster[[i]][j], 1]), x1 = f(new.xy[cluster[[i]][j+1], 1]),
##                      y0 = f(new.xy[cluster[[i]][j], 2]), y1 = f(new.xy[cluster[[i]][j+1], 2]),
##                 col = i %% 4 + 1, lty = i %% 5 + 1, lwd = 2)
##             }
##         }
##     }
## lines(f(c(0, 1, 1)), f(c(1, 1, 0)), col = 'black', lwd = 2)
## points(f(new.xy[tmp.ind,]), col = col_fade("dodgerblue", 0.5), cex = 1.2, pch = 15)
## 
## sapply(get.clusters(tmp.ind, 1), length)
## 
## intersect(tmp.ind1, tmp.ind2)
## union(tmp.ind1, tmp.ind2)
## 
## x[[50]]$varmat[tmp.ind]
## 
## x[[50]]$varmat[299]
## 
## which.min(abs(x[[50]]$varmat - 59.45))
## 
## x[[50]]$y + x[[50]]$threshold
## x[[50]]$varmat
## 
## # flag.ind = NULL
## # lty.count = 1
## # col.count = 2
## # flag = c(0, 0)
## # for (i in 1:(length(ex.times)-1)){
## #     if (ex.times[i] == ex.times[i+1]-1){
## #         segments(x0 = f(new.xy[ex.times[i], 1]), x1 = f(new.xy[ex.times[i+1], 1]),
## #                  y0 = f(new.xy[ex.times[i], 2]), y1 = f(new.xy[ex.times[i+1], 2]),
## #             col = col.count %% 4 + 1, lty = lty.count %% 5 + 1, lwd = 2)
## #         if (flag[1] == 0){
## #             flag.ind = c(flag.ind, i)
## #             flag[1] = 1
## #             flag[2] = flag[2] + 1
## #             }
## #     } else {
## #         lty.count = lty.count + 1
## #         col.count = col.count + 1
## #         flag[1] = 0
## #         }
## #     }
## points(f(new.xy[flag.ind,]), pch = 15, col = 'dodgerblue', cex = 1.5)
## lines(f(c(0, 1, 1)), f(c(1, 1, 0)), col = 'black', lwd = 2)
## 
## 
## c.max = sapply(cluster, function(x) x[(which.max(new.xy[x,])-1) %% length(x) + 1])
## plot(f(new.xy[ex.times,]), pch = 15, col = 'dodgerblue', cex = 1.5)
## points(f(new.xy[c.max,]), pch = 15, col = 'firebrick', cex = 1.0)
## lines(f(c(0, 1, 1)), f(c(1, 1, 0)), col = 'black', lwd = 2)
## 
## 
## plot(new.xy[c.max,c(2,1)])
## 
## 
## 
## flag[2]                 # Number of exceedances occurring in groups
## sum(diff(ex.times) > 1) # Number of exceedances occurring alone
## flag[2] + sum(diff(ex.times) > 1) # Total Number of cluster exceedance
## 
## 
## 
## par(mfrow = c(1,1))
## pdf("~/anglebeta4.pdf", width = width, height = height)
## hist(unlist(phi.s), breaks = 50, col = 'black', border = 'gray50', freq = FALSE, ylab = "",
##     xlab = expression(phi), main = "Angle", cex.lab = 1.5, cex.main = 2)
## hist(pred, col = rgb(34/255, 139/255, 34/255, 0.5), breaks = 50, freq = FALSE,
##     border = 'lightgreen', add = TRUE)
## title(main = bquote(chi ~ "=" ~ .(round(h1, 3))), line = 0)
## dev.off()
## 
## ind = (pred.cone[,1] > 0 & pred.cone[,2])
## 
## z1 = pred.cone[ind, 1] / mean(pred.cone[, 1])
## z2 = pred.cone[ind, 2] / mean(pred.cone[, 2])
## 
## z = cbind(z1, z2)
## m1 = apply(z, 1, min)
## (h1 = mean(m1) * mean(ind))
## 
## 
## ps = unlist(phi.s)
## 
## 
## ucone = NULL
## for (i in 1:R)
##     ucone = rbind(ucone, cone[[i]])
## 
## ind = (ucone[,1] > 0 & ucone[,2])
## 
## z1 = ucone[ind, 1] / mean(ucone[, 1])
## z2 = ucone[ind, 2] / mean(ucone[, 2])
## 
## z = cbind(z1, z2)
## m1 = apply(z, 1, min)
## (h1 = mean(m1) * mean(ind))
## 
## 
## 
## plot(pred.cone[sample(nmcmc, 10000, replace = TRUE),], bty = 'n',
##     col = rgb(34/255,139/255,34/255,0.01), pch = 15, xlab = "x", ylab = "y",
##     main = "Cone", cex.lab = 1.5, cex.main = 2)
## points(cbind(cone[[1]][cone[[1]][,1]<1,1], cone[[1]][cone[[1]][,1]<1,2]*0.99), pch = 15)
## points(cbind(cone[[1]][cone[[1]][,2]<1,1]*0.99, cone[[1]][cone[[1]][,2]<1,2]), pch = 15)
## 
## par(mfrow = c(3,4))
## for (i in 1:12)
##     plot_hpd(out$params[,i])
## 
## 
## z = cbind(pred.cone[,1] / mean(pred.cone[,1]),
##     pred.cone[,2] / mean(pred.cone[,2]))
## h = mean(apply(z, 1, min))
## 
## 
## pred.v = 1/runif(nmcmc)
## 
## mean(apply(pred.cone, 1, min)) / mean(pred.cone[,2])
## mean(apply(pred.cone, 1, min)) / mean(pred.cone[,1])
## 
## mean(apply(cone, 1, min)) / mean(cone[,2])
## mean(apply(cone, 1, min)) / mean(cone[,1])
## 
## hist(v, col = 'green', border = 'white', breaks = 50, freq = FALSE)
## for (i in 1:1000)
##     hist(1/runif(length(v)), col = rgb(30/255, 144/255, 255/255, 0.01), breaks = 50, freq = FALSE, border = rgb(1,1,1,0), add = TRUE)
## hist(v, col = 'black', border = 'white', breaks = 50, freq = FALSE, add = TRUE)
## curve(1/x^2, add = TRUE, col = 'blue', lwd = 3)
## 
## 
## 
## u = exp(0:8)
## chi.12 = matrix(0, 1000, length(u))
## chi.21 = matrix(0, nrow(chi.12), length(u))
## chi.s1 = matrix(0, nrow(chi.12), length(u))
## chi.s2 = double(nrow(chi.12))
## par(mfrow = c(1, 1))
## plot(log(u), rep(0, length(u)), ylim = c(0, 1), type = 'n')
## for (j in 1:nrow(chi.12)){
##     pred.v = 1/runif(nmcmc)
##     pred = get.preds(out)
## #   pred = rbeta(nmcmc, 1/100, 1/100)
## #   pred = rbeta(nmcmc, 10, 10)
##     pred.cone = angle.2.cone(pred * pi/2)
##     pred.xy = pred.cone * pred.v
## 
##     ind = (pred.cone[,1] > 0 & pred.cone[,2])
## 
##     z1 = pred.cone[ind, 1] / mean(pred.cone[, 1])
##     z2 = pred.cone[ind, 2] / mean(pred.cone[, 2])
## 
##     z = cbind(z1, z2)
##     m1 = apply(z, 1, min)
##     h1 = mean(m1) * mean(ind)
## 
##     
## 
## #   z1 = pred.cone[, 1] / mean(pred.cone[, 1])
## #   z2 = pred.cone[, 2] / mean(pred.cone[, 2])
## #   z = cbind(z1, z2)
## #   m2 = apply(z, 1, min)
## #   h2 = mean(m2)
## 
##     for (i in 1:length(u)){
##         chi.12[j,i] = mean(pred.xy[,1] > u[i] & pred.xy[,2] > u[i]) / mean(pred.xy[,2] > u[i])
##         chi.21[j,i] = mean(pred.xy[,1] > u[i] & pred.xy[,2] > u[i]) / mean(pred.xy[,1] > u[i])
##         chi.s1[j,i] = mean(pred.xy[,1] > u[i] & pred.xy[,2] > u[i])
##         }
## 
##     chi.12[j,] = ifelse(is.na(chi.12[j,]), 0, chi.12[j,])
##     chi.21[j,] = ifelse(is.na(chi.21[j,]), 0, chi.21[j,])
##     chi.s2[j] = mean(apply(pred.cone, 1, min))
## 
## #   mean(pred.xy[,1] > u[i] & pred.xy[,2] > u[i]) / mean(pred.xy[,2] > u[i])
## #   mean(apply(pred.cone, 1, min)) / mean(pred.cone[,2])
## 
## #   mean(pred.xy[,1] > u[i] & pred.xy[,2] > u[i]) / mean(pred.xy[,1] > u[i])
## #   mean(apply(pred.cone, 1, min)) / mean(pred.cone[,1])
## 
##     points(log(u), chi.12[j,], pch = 16, col = rgb(1,0,0,0.1))
##     points(log(u)+0.1, chi.21[j,], pch = 16, col = rgb(0,0,1,0.1))
##     points(log(u)+0.2, chi.s1[j,], pch = 16, col = rgb(1,0,1,0.1))
##     abline(h = mean(chi.s2[1:j]), col = rgb(0,1,0,0.1))
##     abline(h = h1, col = rgb(0,1,1,0.1))
## #   abline(h = h2, col = rgb(1,1,0,0.1))
## #   readline()
##     }
## 
## chi.mm = (chi.12 + chi.21) / 2
## 
## par(mfrow = c(1, 1))
## plot(log(u), rep(0, length(u)), ylim = c(0, 1), type = 'n')
## for (j in 1:nrow(chi.12)){
##     points(log(u), chi.12[j,], pch = 16, col = rgb(1,0,0,0.1))
##     points(log(u)+0.1, chi.21[j,], pch = 16, col = rgb(0,0,1,0.1))
##     }
## lines(log(u), apply(chi.12, 2, mean), col = 'firebrick', lwd = 2)
## lines(log(u)+0.1, apply(chi.21, 2, mean), col = 'dodgerblue', lwd = 2)
## lines(log(u), apply(chi.mm, 2, mean), col = 'forestgreen', lwd = 2)
## 
## apply(chi.mm, 2, mean)
## 
## 
## 
## # Bootstrap the chi?
## B = 1000
## chi.12 = double(B)
## chi.21 = double(B)
## chi = double(B)
## tmp.mins = apply(pred.cone, 1, min)
## tmp.1 = pred.cone[,1]
## tmp.2 = pred.cone[,2]
## for (i in 1:B){
##     b.samp = sample(nmcmc, nmcmc, replace = TRUE)
##     chi.12[i] = mean(tmp.mins[b.samp]) / mean(tmp.2[b.samp])
##     chi.21[i] = mean(tmp.mins[b.samp]) / mean(tmp.1[b.samp])
##     chi[i] = mean(c(chi.12[i], chi.21[i]))
##     }
## 
## hist(chi, col = rgb(1,0,0,0.5), border = rgb(1,1,1,0.5), freq = FALSE,
##     xlim = range(chi.12, chi.21, chi))
## hist(chi.21, col = rgb(0,0,1,0.5), border = rgb(1,1,1,0.5), add = TRUE, freq = FALSE)
## hist(chi.12, col = rgb(0,1,0,0.5), border = rgb(1,1,1,0.5), add = TRUE, freq = FALSE)
## # abline(v = mean(tmp.mins) / mean(tmp.2), col = 'forestgreen', lwd = 3)
## # abline(v = mean(tmp.mins) / mean(tmp.1), col = 'royalblue4', lwd = 3)
## # abline(v = (mean(tmp.mins)/mean(tmp.2) + mean(tmp.mins)/mean(tmp.1))/2, col = 'firebrick', lwd = 3)
## 
## var(chi.12)
## var(chi.21)
## var(chi)
## 
## sd(chi.12)
## sd(chi.21)
## sd(chi)
## 
## 
## 
## 
## preds2 = get.preds2(out, 1000)
## chi1 = apply(preds2, 1,
##     function(v){
##         cone = angle.2.cone(v * pi/2)
##         ind = (cone[,1] > 0 & cone[,2] > 0)
##         z1 = cone[ind, 1] / mean(cone[, 1])
##         z2 = cone[ind, 2] / mean(cone[, 2])
## 
##         if (sum(ind) == 0)
##             return (0)
## 
##         z = cbind(z1, z2)
##         m1 = apply(z, 1, min)
##         return (mean(m1) * mean(ind))
##         })
## chi2 = apply(preds2, 2,
##     function(v){
##         cone = angle.2.cone(v * pi/2)
##         ind = (cone[,1] > 0 & cone[,2])
##         z1 = cone[ind, 1] / mean(cone[, 1])
##         z2 = cone[ind, 2] / mean(cone[, 2])
## 
##         z = cbind(z1, z2)
##         m1 = apply(z, 1, min)
##         (h1 = mean(m1) * mean(ind))
##         })
## 
## mean(chi1 == 0)
## mean(chi2 == 0)
## 
## plot(density(chi1))
## lines(density(chi2))
