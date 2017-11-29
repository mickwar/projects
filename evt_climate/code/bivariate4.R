library(DPpackage)
library(mwBASE)
library(MASS)
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

DIR = "~/files/gdrive/climate_data/2017_09_08/"

#DIR = "./RData/2017_09_08/"
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


get_phi = function(x, obs.ind, sim.ind, cols){

    R = ncol(x[[sim.ind]]$varmat)
    p = length(cols)
    n = nrow(x[[obs.ind]]$varmat)

    new.xy = matrix(0, n, p+1)

    ex = x[[obs.ind]]$varmat[,1] - x[[obs.ind]]$threshold
    ksi = mean(x[[obs.ind]]$params[,1])
    sig = mean(x[[obs.ind]]$params[,2])

    new.xy[, 1] = (1 + ksi * ex / sig)^(1/ksi)
    new.xy[is.na(new.xy[,1]),1] = 0

    for (j in 1:length(cols)){
        i = cols[j]

        ex = x[[sim.ind]]$varmat[,i] - x[[sim.ind]]$threshold[i]
        ksi = mean(x[[sim.ind]]$params[,i])
        sig = mean(x[[sim.ind]]$params[,i + R])

        new.xy[, j+1] = (1 + ksi * ex / sig)^(1/ksi)
        new.xy[is.na(new.xy[,j+1]),j+1] = 0

        }


    # zero.ind = which(apply(new.xy, 1, function(x) any(x == 0)))
    # ex.ind = which(new.xy[,1] > 1 | new.xy[,2] > 1)

    # ind = setdiff(ex.ind, zero.ind)

    de.ind = rep(list(NULL), p)
    ind = rep(list(NULL), p)
    v = rep(list(NULL), p)
    cone = rep(list(NULL), p)
    phi = rep(list(NULL), p)
    phi.s = rep(list(NULL), p)

    for (j in 1:length(cols)){
        i = cols[j]

        # Declustered indices for each margin, joined
        tmp.ind1 = which((x[[obs.ind]]$varmat - x[[obs.ind]]$threshold) %in% x[[obs.ind]]$y)
        tmp.ind2 = which((x[[sim.ind]]$varmat[,i] - x[[sim.ind]]$threshold) %in% x[[sim.ind]]$y[[i]])
        ind[[j]] = sort(unique(c(tmp.ind1, tmp.ind2)))

        # All exceedances
#       ind[[i]] = which(new.xy[,1] > 1 | new.xy[,i+1] > 1)

        v[[j]] = apply(new.xy[ind[[j]], c(1, j+1)], 1, max)
        cone[[j]] = new.xy[ind[[j]], c(1, j+1)] / v[[j]]
        phi[[j]] = atan(new.xy[ind[[j]], j+1] / new.xy[ind[[j]], 1])
        phi.s[[j]] = phi[[j]] * 2 / pi
        }

    return (list("phi.s" = phi.s, "v" = v, "cone" = cone))
    }


wrapper = function(x, obs.ind, sim.ind, cols, zero, one, parts, nburn = 10000, nmcmc = 10000){

    R = ncol(x[[sim.ind]]$varmat)
    p = length(cols)
    n = nrow(x[[obs.ind]]$varmat)

    new.xy = matrix(0, n, p+1)

    ex = x[[obs.ind]]$varmat[,1] - x[[obs.ind]]$threshold
    ksi = mean(x[[obs.ind]]$params[,1])
    sig = mean(x[[obs.ind]]$params[,2])

    new.xy[, 1] = (1 + ksi * ex / sig)^(1/ksi)
    new.xy[is.na(new.xy[,1]),1] = 0

    for (j in 1:length(cols)){
        i = cols[j]

        ex = x[[sim.ind]]$varmat[,i] - x[[sim.ind]]$threshold[i]
        ksi = mean(x[[sim.ind]]$params[,i])
        sig = mean(x[[sim.ind]]$params[,i + R])

        new.xy[, j+1] = (1 + ksi * ex / sig)^(1/ksi)
        new.xy[is.na(new.xy[,j+1]),j+1] = 0

        }


    # zero.ind = which(apply(new.xy, 1, function(x) any(x == 0)))
    # ex.ind = which(new.xy[,1] > 1 | new.xy[,2] > 1)

    # ind = setdiff(ex.ind, zero.ind)

    de.ind = rep(list(NULL), p)
    ind = rep(list(NULL), p)
    v = rep(list(NULL), p)
    cone = rep(list(NULL), p)
    phi = rep(list(NULL), p)
    phi.s = rep(list(NULL), p)

    for (j in 1:length(cols)){
        i = cols[j]

        # Declustered indices for each margin, joined
        tmp.ind1 = which((x[[obs.ind]]$varmat - x[[obs.ind]]$threshold) %in% x[[obs.ind]]$y)
        tmp.ind2 = which((x[[sim.ind]]$varmat[,i] - x[[sim.ind]]$threshold) %in% x[[sim.ind]]$y[[i]])
        ind[[j]] = sort(unique(c(tmp.ind1, tmp.ind2)))

        # All exceedances
#       ind[[i]] = which(new.xy[,1] > 1 | new.xy[,i+1] > 1)

        v[[j]] = apply(new.xy[ind[[j]], c(1, j+1)], 1, max)
        cone[[j]] = new.xy[ind[[j]], c(1, j+1)] / v[[j]]
        phi[[j]] = atan(new.xy[ind[[j]], j+1] / new.xy[ind[[j]], 1])
        phi.s[[j]] = phi[[j]] * 2 / pi
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

            out = ifelse(is.na(out), -Inf, out)

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

    flag = 1
    while (flag){
        calc.post = make.post(zero, one, parts)
        nparam = zero + one + 3*parts - 1
        if (nparam != 2*parts){
            mcmc = try(mcmc_sampler(phi.s, calc.post, nparam = nparam, nburn = nburn, nmcmc = nmcmc,
                chain_init = c(rep(0, 2*parts), rep(log(1/(zero + one + parts)),
                zero + one + parts - 1))), silent = TRUE)
            }
        if (nparam == 2*parts){
            mcmc = mcmc_sampler(phi.s, calc.post, nparam = nparam, nburn = nburn, nmcmc = nmcmc,
                chain_init = rep(0, 2*parts))
            }
        if (class(mcmc) == "try-error"){
            parts = max(parts - 1, 1)
            cat("      MCMC chain broke. Reducing number of beta components to", parts, "\n")
        } else {
            flag = 0
            if (parts > 1){
                for (j in 1:parts){
                    a = mean(mcmc$params[,2*j-1])
                    b = mean(mcmc$params[,2*j])
                    # The 1e-4 is arbitrary, perhaps can have it based on the data?
                    if ( a*b / ( (a+b)^2 * (a+b+1) ) < 1e-6 ){
                        flag = 1
                        }
                    }
                if (flag){
                    parts = max(parts - 1, 1)
                    cat("      A beta component had two low of variance. Reducing number of beta components to", parts, "\n")
                    }
                }
            }
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

calc.chi = function(cone){
    ind = (cone[,1] > 0 & cone[,2] > 0)
    z1 = cone[ind, 1] / mean(cone[, 1])
    z2 = cone[ind, 2] / mean(cone[, 2])
    z = cbind(z1, z2)
    m1 = apply(z, 1, min)
    chi = mean(m1) * mean(ind)
    return (chi)
    }





# MCMC parameters
mcmc = list(nburn = 1000,
    nsave = 1000,
    nskip= 1,
    ndisplay = 100)

# Prior
prior = list(aa0 = 1,
    ab0 = 0.1,
    kmax = 300,
    a0 = 1,
    b0 = 1)
#prior = list(alpha = 0.1,
#    kmax = 100,
#    a0 = 1/2,
#    b0 = 1/2)

B = 100000

R = 10  # Number of replicates
width = 9
height = 9


# rep.num = 1
# comp = c(18, 34)
# comp = seq(1, 47, by = 2)
# comp = c(17, 33)
# comp = c(1, 3)
# n = nrow(x[[comp[1]]]$varmat)

load("~/files/gdrive/climate_data/biv_chi4.RData")
load("~/files/gdrive/climate_data/biv_pred4.RData")

u = seq(0.8, 0.999, length = 30)
# out = rep(list(NULL), length(ind.control) + length(ind.decadal) + length(ind.historical))
# out_samps = rep(list(NULL), length(ind.control) + length(ind.decadal) + length(ind.historical))
#out = matrix(0, length(ind.control) + length(ind.decadal) + length(ind.historical), 11)

par(mfrow = c(4,4))
for (j in 1:3){
    for (i in 1:16){
        hh = 16*(j-1)+i
        hist(out_samps[[hh]]$dat.phi[[11]], freq = FALSE, main = hh)
        lines(out_samps[[hh]]$grid,
            out_samps[[hh]]$fun[,11], col = 'blue')
        }
    readline()
    }

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

                    if (!(rr.sim %in% c(21, 29, 30, 45)))
                        next

                    cat(rr.sim)

#                   out[[rr.sim]]$u = u     # range of u for chi.u and chi.bar.u
#                   out[[rr.sim]]$dat.chi = double(11)  # chi from data
#                   out[[rr.sim]]$dat.chi.u = matrix(0, length(u), 11)
#                   out[[rr.sim]]$dat.chi.bar.u = matrix(0, length(u), 11)
#                   out[[rr.sim]]$pp.chi = double(11)   # chi from pareto process
#                   out[[rr.sim]]$pp.chi.u = matrix(0, length(u), 11)
#                   out[[rr.sim]]$pp.chi.bar.u = matrix(0, length(u), 11)

#                   out_samps[[rr.sim]]$grid = seq(0, 1, length = 1000) # BDPdenity default
#                   out_samps[[rr.sim]]$fun = matrix(0, 1000, 11)
#                   out_samps[[rr.sim]]$probs = matrix(0, 3, 11)
#                   rownames(out_samps[[rr.sim]]$probs) = c("zero", "one", "inside")
#                   out_samps[[rr.sim]]$dat.phi = rep(list(NULL), 11)
                    

                    for (i in 11){
                        cat(paste0("\n", lab.shortmod[A5]), lab.reg[A1], lab.year[A2], lab.var[A3],
                            lab.sea[A4], "Rep =", i, "\n")
                        cols = i
                        if (i == 11) cols = 1:10

                        # Seemed to be getting issues with values too close to the boundary
                        tmp = get_phi(x, rr.obs, rr.sim, cols = cols)
                        phi = unlist(tmp$phi.s)
                        ind.rm.0 = (phi <= 0.005)
                        ind.rm.1 = (phi >= 0.995)
                        ind.rm = (ind.rm.0 | ind.rm.1)
                        cone = do.call(rbind, tmp$cone)
                        v = do.call(c, tmp$v)
                        
                        # Fit Bernstein-Dirichlet prior
                        fit = BDPdensity(y = phi[!ind.rm], prior = prior, mcmc = mcmc,
                            state = NULL, status = TRUE, support = 1)


                        # Get a posterior predictive distribution
                        pred = sample(fit$grid, B, replace = TRUE, prob = fit$fun)

                        # Manually insert 0's and 1's
                        pred.ind = sample(c(0,1,2), B, replace = TRUE,
                            prob = c(mean(ind.rm.0), mean(ind.rm.1), 1-mean(ind.rm)))
                        pred[pred.ind == 0] = 0
                        pred[pred.ind == 1] = 1

                        pred.cone = angle.2.cone(pred*pi/2)


                        out_samps[[rr.sim]]$fun[,i] = fit$fun
#                       out_samps[[rr.sim]]$probs[1,i] = mean(ind.rm.0)
#                       out_samps[[rr.sim]]$probs[2,i] = mean(ind.rm.1)
#                       out_samps[[rr.sim]]$probs[3,i] = 1-mean(ind.rm)
                        out_samps[[rr.sim]]$dat.phi[[i]] = phi


                        out[[rr.sim]]$pp.chi[i] = calc.chi(pred.cone)
                        out[[rr.sim]]$dat.chi[i] = calc.chi(cone)

                        pred.v = 1/runif(nrow(pred.cone))
                        pred.y = pred.v * pred.cone

                        z1 = rank(pred.y[,1]) / nrow(pred.y)
                        z2 = rank(pred.y[,2]) / nrow(pred.y)

                        for (j in 1:length(u)){
                            ind1 = (z1 <= u[j])
                            ind2 = (z2 <= u[j])
                            out[[rr.sim]]$pp.chi.u[j,i]= 2 - log(mean(ind1 & ind2)) / log(mean(ind2))

                            ind1 = (z1 > u[j])
                            ind2 = (z2 > u[j])
                            out[[rr.sim]]$pp.chi.bar.u[j,i] = 2 * log(mean(ind1)) / log(mean(ind1 & ind2)) - 1
                            }

                        pred.y = v * cone

                        z1 = rank(pred.y[,1]) / nrow(pred.y)
                        z2 = rank(pred.y[,2]) / nrow(pred.y)

                        for (j in 1:length(u)){
                            ind1 = (z1 <= u[j])
                            ind2 = (z2 <= u[j])
                            out[[rr.sim]]$dat.chi.u[j,i]= 2 - log(mean(ind1 & ind2)) / log(mean(ind2))

                            ind1 = (z1 > u[j])
                            ind2 = (z2 > u[j])
                            out[[rr.sim]]$dat.chi.bar.u[j,i] = 2 * log(mean(ind1)) / log(mean(ind1 & ind2)) - 1
                            }

                        par(mfrow = c(2, 2))
#                       matplot(u, out[[rr.sim]]$dat.chi.u[,1:i], lty = 1, type = 'l',
#                           ylab = "chi.u", col = c(rainbow(10), "black"), lwd = c(rep(1.5, 10), 2))
#                       points(rep(1, i), out[[rr.sim]]$dat.chi[1:i], pch = 1,
#                           col = c(rainbow(10), "black"), cex = c(rep(1, 10), 1.5))

                        matplot(u, out[[rr.sim]]$pp.chi.u[,1:i], lty = 1, type = 'l',
                            ylab = "chi.u", col = c(rainbow(10), "black"), lwd = c(rep(1.5, 10), 2))
                        points(rep(1.001, i), out[[rr.sim]]$pp.chi[1:i], pch = 15,
                            col = c(rainbow(10), "black"), cex = c(rep(1, 10), 1.5))

#                       matplot(u, out[[rr.sim]]$dat.chi.bar.u[,1:i], lty = 1, type = 'l',
#                           ylab = "chi.bar", col = c(rainbow(10), "black"),
#                           lwd = c(rep(1.5, 10), 2))
                        matplot(u, out[[rr.sim]]$pp.chi.bar.u[,1:i], lty = 1, type = 'l',
                            ylab = "chi.bar", col = c(rainbow(10), "black"),
                            lwd = c(rep(1.5, 10), 2))

                        hist(phi, freq = FALSE, col = 'gray', border = 'white', breaks = 20,
                            main = round(out[[rr.sim]]$pp.chi[i], 3))
                        lines(density(pred), col = c(rainbow(10), "black")[i],
                            lwd = c(rep(1.5, 10), 2.5)[i])

                        }

#                   hist(tmp$phi.s, breaks = 50, col = 'black', border = 'gray50',
#                       freq = FALSE, ylab = "", xlab = expression(phi), main = "Angle",
#                       cex.lab = 1.5, cex.main = 2)
#                   hist(pred, col = rgb(34/255, 139/255, 34/255, 0.5), breaks = 50, freq = FALSE,
#                       border = 'lightgreen', add = TRUE)
#                   abline(v = chi)

                    }
                }
            }
        }
    }

for (i in 1:48){
    out_samps[[i]]$probs = matrix(0, 3, 11)
    rownames(out_samps[[i]]$probs) = c("zero", "one", "inside")
    for (j in 1:11){
        phi = out_samps[[i]]$dat.phi[[j]]
    
        ind.rm.0 = (phi <= 0.005)
        ind.rm.1 = (phi >= 0.995)
        ind.rm = (ind.rm.0 | ind.rm.1)

        out_samps[[i]]$probs[1,j] = mean(ind.rm.0)
        out_samps[[i]]$probs[2,j] = mean(ind.rm.1)
        out_samps[[i]]$probs[3,j] = 1-mean(ind.rm)
        }
    }




save(out, file = "./biv_chi4.RData")
save(out_samps, file = "./biv_pred4.RData")



for (i in 1:length(out_samps)){
    par(mfrow = c(3, 4), mar = c(5.1, 4.1,4.1,2.1), oma = c(0,0,0,0))
    for (j in 1:11){
        hist(out_samps[[i]]$dat[[j]], freq = FALSE, col = 'gray', border = 'white',
            main = j, xlab = "", ylab = "", breaks = 20)
        lines(out_samps[[i]]$grid, out_samps[[i]]$fun[,j], col = 'dodgerblue', lwd = 3)
        legend("topright", legend = round(out[[i]]$pp.chi[j], 3), bty = 'n', cex = 2)
        }
    title(main = i, outer = TRUE, line = -2, cex.main = 2.5)
    readline()
    }


load("~/files/gdrive/climate_data/biv_chi4.RData")


pdf("~/chi4.pdf", height = height, width = width)
par(mfrow = c(4, 4), mar = c(0,0,0,0), oma = c(6,10,10,6))
pnum = 0
#xlim = range(sapply(out, function(x) log(x$pp.chi[11] / (1 - x$pp.chi[11]))))
xlim = range(sapply(out, function(x) c(x$dat.chi, x$pp.chi)))
#xlim = c(0, 1)
ylim = c(0.5, 3.5)
cols = c("blue", "green", "red", "gray50")
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
#                   points(log(out[[rr]]$pp.chi[11] / (1 - out[[rr]]$pp.chi[11])), A5,
#                       pch = 15, col = cols[A5])
#                   points(out[[rr]]$pp.chi[11],
#                       A5, pch = 15, col = cols[A5], cex = 1.5)
#                   points(out[[rr]]$pp.chi, rep(A5, 11), pch = 15,
#                       col = cols[A5], cex = c(rep(1, 10), 1.5))
                    points(out[[rr]]$pp.chi[1:10], rep(A5, 10), pch = 16,
                        col = col_fade(cols[A5], 0.5), cex = 1)
                    points(out[[rr]]$pp.chi[11], A5, pch = 4,
                        col = col_mult(cols[A5], "gray50"), cex = 2.5, lwd = 2.5)

#                   points(out[[rr]]$dat.chi[1:10], rep(A5-0.2, 10), pch = 16,
#                       col = cols[A5], cex = 1)
#                   points(out[[rr]]$dat.chi[11], A5-0.2, pch = 4,
#                       col = col_mult(cols[A5], "gray50"), cex = 2.0, lwd = 2)
                    }

                if ((pnum <= 4) && (pnum %% 2 == 1))
                    axis(3, lwd = -1, lwd.ticks = 1)
                if ((pnum >= 13) && (pnum %% 2 == 0))
                    axis(1, lwd = -1, lwd.ticks = 1)
                if (pnum %% 4 == 0){
                    axis(side = 4, labels = lab.shortmod[-4],
                        at = 1:3, las = 1, lwd = -1, lwd.ticks = 1)
#                   axis(side = 4, labels = paste("PP", lab.shortmod[-4]),
#                       at = 1:3, las = 1, lwd = -1, lwd.ticks = 1)

#                   axis(side = 4, labels = paste("DAT", lab.shortmod[-4]),
#                       at = 1:3 - 0.2, las = 1, lwd = -1, lwd.ticks = 1)
                    }
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
dev.off()












library(DPpackage)

state = NULL

# MCMC parameters
mcmc = list(nburn = 1000,
    nsave = 1000,
    nskip= 1,
    ndisplay = 100)

# Prior
prior = list(aa0 = 1,
    ab0 = 0.1,
    kmax = 100,
    a0 = 1,
    b0 = 1)

fit = BDPdensity(y = y, prior = prior, mcmc = mcmc,
    state = state, status = TRUE, support = 1)


for (j in 1:11){
    hist(out_samps[[i]]$dat[[j]], freq = FALSE, col = 'gray', border = 'white',
        main = j, xlab = "", ylab = "", breaks = 20)
    lines(density(out_samps[[i]]$pred[,j]), col = 'darkgreen', lwd = 3)
    legend("topright", legend = round(out[[i]]$pp.chi[j], 3), bty = 'n', cex = 2)
    }
title(main = i, outer = TRUE, line = -2, cex.main = 2.5)

# Fitting the model
i = 1
j = 3
y = out_samps[[i]]$dat[[j]]
fit = BDPdensity(y = y, prior = prior, mcmc = mcmc,
    state = state, status = TRUE, support = 1)

hist(y, freq = FALSE, col = 'gray', border = 'white',
    main = j, xlab = "", ylab = "", breaks = 20)
lines(fit$grid, fit$fun, col = 'red', lwd = 2)

B = 100000
boot = sample(fit$grid, B, replace = TRUE, prob = fit$fun)
lines(density(boot))

calc.chi(angle.2.cone(boot*pi/2))

object.size(fit)
object.size(cbind(fit$grid, fit$fun))
plot(density(fit$save.state$thetasave[,3]))
