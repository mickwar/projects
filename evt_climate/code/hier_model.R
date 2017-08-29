# data_name = "decadal1989"
# data_dir = "../data/processed/"
# region = "california"
# variable = "tasmax"
# season = "summer"
# date_begin = "1990-01-01"
# date_end = "1999-12-31"
# uq = 0.95
# months = NULL
# min_uq = TRUE
# nburn = 6000
# nmcmc = 800
# window = 500
# anomaly = TRUE

hier_excess = function(data_name, data_dir, region, variable,
    season = c("year", "summer", "winter"), date_begin, date_end,
    r = 0, uq = 0.90, min_uq = TRUE, threshold,
    nburn = 80000, nmcmc = 40000, window = 500,
    months = NULL, anomaly = TRUE){

    require(MASS)
    require(mwBASE)
    require(mwEVT)
#   require(truncnorm)

    output = NULL
    output$model = data_name

    output$region = region
    output$variable = variable

    ### Return the index of a vector of dates (x) which contain months
    subset_data = function(ts, months){
        convert = function(x){
            y = tolower(substr(x, 1, 3))
            out = character(length(y))
            for (i in 1:length(y)){
                out[i] = switch(y[i],
                    "jan" = "01", "feb" = "02", "mar" = "03", "apr" = "04",
                    "may" = "05", "jun" = "06", "jul" = "07", "aug" = "08",
                    "sep" = "09", "oct" = "10", "nov" = "11", "dec" = "12",
                    "1" = "01", "2" = "02", "3" = "03", "4" = "04",
                    "5" = "05", "6" = "06", "7" = "07", "8" = "08",
                    "9" = "09", y[i])
                }
            return (out)
            }
        months = convert(months)
        out = NULL
        for (i in 1:length(months))
            out = c(out, grep(paste0("-", months[i], "-"), ts))
        out = sort(out)
        return (out)
        }

    if (is.null(months)){
        output$season = season[1]
        output$months = switch(output$season,
            "winter" = c("dec", "jan", "feb"),
            "summer" = c("jun", "jul", "aug"),
            c("jan", "feb", "mar", "apr", "may", "jun",
                "jul", "aug", "sep", "oct", "nov", "dec"))
    } else {
#       message(paste0("Note: `months' variable manually set, `season' variable set to `", season[1], "'.\n"))
        output$season = season[1]
        output$months = months
        }
    output$days.in.months = 0
    for (i in 1:length(output$months)){
        output$days.in.months = output$days.in.months + 
            switch(output$months[i],
                "jan"=,"mar"=,"may"=,"jul"=,"aug"=,
                "oct"=,"dec"=31,
                "apr"=,"jun"=,"sep"=,"nov"=30,
                "feb"=28)
        }
    output$anomaly = anomaly

    ### Some settings based on the data used
    output$units = switch(output$variable,
        "pr"        = "m/day",
        "tas"       = "Kelvin",
        "tasmax"    = "Kelvin",
        "tasmin"    = "Kelvin",
        NULL)
    output$units = paste0(output$units, ifelse(output$anomaly, " anomaly", ""))
    output$label = switch(output$variable,
        "pr"        = "Total Precipitation",
        "tas"       = "Average Temperature",
        "tasmax"    = "Average Maximum Temperature",
        "tasmin"    = "Average Minimum Temperature",
        NULL)
    output$main.prefix = switch(substr(output$model, 1, 7),
        "decadal" = sub("d", "D", output$model),
        "histori" = "Historical",
        "picontr" = "piControl",
        "control" = "piControl",
        "Observation")
    output$base.color = switch(substr(output$model, 1, 7),
        "decadal" = "orange",
        "histori" = "firebrick1",
        "picontr" = "dodgerblue",
        "control" = "dodgerblue",
        "gray50")

    ### Load the processed data
    file_name = paste0(data_dir, output$region, "_", output$variable,
        "_", output$model, ".txt")
#   message("Reading data from ", file_name, " ...")
    dat = read.table(file_name, header = TRUE)

    # Remove leap years
    tmp = grep("-02-29", substr(dat[,1], 5, 10))
    if (length(tmp) > 0)
        dat = dat[-tmp,]

    ### DLM anomalies
    if (output$anomaly){
        output$anomaly_type = "dlm"
        source("~/files/repos/r-sandbox/dlm.R")
        require(Matrix)

        # Use a linear trend and 4 seasonal components (first four harmonics)
        Gj = function(j, p)
            matrix(c(cos(2*pi*j/p), -sin(2*pi*j/p), sin(2*pi*j/p), cos(2*pi*j/p)), 2, 2)
#       F.vec = matrix(c(c(1,0), rep(c(1,0), 4)), ncol = 1)
#       G.mat = as.matrix(bdiag(matrix(c(1,0,1,1),2,2),
#           Gj(1, 365), Gj(2, 365), Gj(3, 365), Gj(4, 365)))

        F.vec = matrix(c(c(1,0), rep(c(1,0), 2)), ncol = 1)
        G.mat = as.matrix(bdiag(matrix(c(1,0,1,1),2,2),
            Gj(1, 365), Gj(2, 365)))

        # If available, use earlier and later observations when making the dlm
        dlm_range = c(max(which(date_begin == dat[,1]) - 100, 1),
            min(which(date_end == dat[,1]) + 100, length(dat[,1])))
        dlm_seq = dlm_range[1]:dlm_range[2]

        output$anomaly_shift = matrix(0, NROW(dat), NCOL(dat)-1)


        for (j in 2:NCOL(dat)){
#           eps = 1e-6
#           tmp = list("y" = log(dat[dlm_seq, j] + eps), "F" = F.vec, "G" = G.mat)
            tmp = list("y" = dat[dlm_seq, j], "F" = F.vec, "G" = G.mat)
            n = length(dlm_seq)
            p = length(F.vec)
            m0 = double(p)
            m0[1] = mean(tmp$y)
            C0 = var(tmp$y)*diag(p)
            n0 = 1
            d0 = 100
            params = get.vals(tmp, delta = 0.9999, m0, C0, n0, d0)
            spar = smooth(tmp, params)

            dat[dlm_seq,j] = dat[dlm_seq,j] - spar$ft.s
            output$anomaly_shift[dlm_seq,j-1] = spar$ft.s

#           plot(tmp$y)
#           lines(spar$ft.s, col = 'blue', lwd = 2)
#           plot(tmp$y - spar$ft.s, type = 'l')
            }
        }

    # pdf("~/dlm.pdf", width = 12, height = 6)
    # par(mfrow = c(1, 2))
    # plot(dat[head(dlm_seq, 1000),11] + output$anomaly_shift[head(dlm_seq, 1000),10], 
    #     axes = FALSE, xlab = "Time", ylab = "Average tasmax", type = 'o')
    # axis(2)
    # axis(1, at = seq(1, 1000, length = 7), labels = (dat[head(dlm_seq, 1000), 1])[seq(1, 1000, length = 7)])
    # lines(spar$ft.s, col = 'blue', lwd = 3)
    # abline(v = grep("-06-01", dat[,1]), lty = 2, col = 'green')
    # abline(v = grep("-08-31", dat[,1]), lty = 2, col = 'green')
    # 
    # plot(dat[head(dlm_seq, 1000),11], 
    #     axes = FALSE, xlab = "Time", ylab = "Average tasmax anomaly", type = 'o')
    # axis(2)
    # axis(1, at = seq(1, 1000, length = 7), labels = (dat[head(dlm_seq, 1000), 1])[seq(1, 1000, length = 7)])
    # lines(spar$ft.s, col = 'blue', lwd = 3)
    # abline(v = grep("-06-01", dat[,1]), lty = 2, col = 'green')
    # abline(v = grep("-08-31", dat[,1]), lty = 2, col = 'green')
    # abline(h = 0, col = 'blue', lwd = 3)
    # dev.off()

    ### Deal with those leap years and get things all on the same time scale
#   message(paste0("Indexing data within ", date_begin, " and ", date_end,
#       " for months ", paste(output$months, collapse=" "), " ..."))
    if (substr(output$model, 1, 7) == "control"){
        season.ind = subset_data(dat[,1], output$months)
        dat[,1] = paste0(as.numeric(substr(dat[,1], 1, 4)) +
            as.numeric(substr(date_begin, 1, 4)) - 1, substr(dat[,1], 5, 10))
        dec.ind = c(min(grep(date_begin, dat[,1])), max(grep(date_end, dat[,1])))
        season.ind = season.ind[season.ind >= dec.ind[1] & season.ind <= dec.ind[2]]
        dat = dat[season.ind,]
    } else {
        season.ind = subset_data(dat[,1], output$months)
        ly.ind = grep("-02-29", dat[,1])
        season.ind = season.ind[!(season.ind %in% ly.ind)]
        dec.ind = c(min(grep(date_begin, dat[,1])), max(grep(date_end, dat[,1])))
        season.ind = season.ind[season.ind >= dec.ind[1] & season.ind <= dec.ind[2]]
        dat = dat[season.ind,]
        }
    if (length(output$anomaly_shift) > 0)
        output$anomaly_shift = output$anomaly_shift[season.ind,]


    output$time.dates = as.Date(dat[,1])
    output$varmat = as.matrix(dat[,2:NCOL(dat)])
    nReps = NCOL(output$varmat)
    output$run.index = 1:nReps

    # Number of full years from date_begin to date_end
    output$nyears = length(output$time.dates) / output$days.in.months

    if (missing(threshold)){
        output$uq = uq
        if (min_uq){
            threshold = min(apply(output$varmat, 2, quantile, output$uq))
            output$threshold = rep(threshold, nReps)
        } else {
            if (length(output$uq) == 1)
                output$uq = rep(output$uq, nReps)
            output$threshold = diag(apply(output$varmat, 2, quantile,  output$uq))
            }
    } else {
#       message("Note: value for `threshold' specified, overriding value for `uq'.\n")
        output$uq = NULL
        output$threshold = threshold
        if (length(output$threshold) == 1)
            output$threshold = rep(output$threshold, nReps)
        }


    ### Estimate theta and decluster
#   tmp = hier_theta(output$varmat, output$threshold[1], nburn = 40000, nmcmc = 20000)

    tmp = theta_hier(y = output$varmat, u = output$threshold[1], likelihood = "ferro",
        nburn = 40000, nmcmc = 60000)

    output$y = rep(list(NULL), nReps)
    for (j in 1:nReps)
        output$y[[j]] = tmp$y[[j]] - output$threshold[j]
    output$n_c = sapply(tmp$y, length)
    output$n_u = as.vector(tmp$N)
    output$n = rep(NROW(output$varmat), nReps)

    output$T_C = tmp$T_C
    output$theta = tmp$mcmc
    output$theta.mean = colMeans(tmp$mcmc)
    output$theta.ints = apply(tmp$mcmc, 2, hpd_mult, force_uni = TRUE)


    ### MCMC
    nparam = 3*nReps + 6
    params = matrix(0, nburn + nmcmc, nparam)
    accept = double(nburn + nmcmc)

    cand.sig = diag(0.01, nparam)
# tmp.params = tail(params, 1)

    # priors
       ksi.mean = 0
         ksi.sd = 0.5
     tau2.alpha = 1
      tau2.beta = .1
    sig.alpha.a = 1
    sig.alpha.b = .1
     sig.beta.a = 1
     sig.beta.b = .1
      zeta.mu.a = 1
      zeta.mu.b = 9
     zeta.eta.a = 1
     zeta.eta.b = .1
#    theta.mu.a = 1
#    theta.mu.b = 1
#   theta.eta.a = 1
#   theta.eta.b = .1
    
    params[1,] = 
        c(rep(0, nReps),    # ksi_r,
        rep(1, nReps),      # sig_r
        rep(1-uq, nReps),   # zeta_r
        0, 1, 1, 1, 1-uq, 1)
#   tmp = c(abs(rnorm(1, ksi.mean, ksi.sd)),
#       rgamma(1, tau2.alpha,  tau2.beta),
#       rgamma(1, sig.alpha.a, sig.alpha.b),
#       rgamma(1, sig.beta.a,  sig.beta.b),
#        rbeta(1, zeta.mu.a,   zeta.mu.b),
#       rgamma(1, zeta.eta.a,  zeta.eta.b))
#   params[1,] = 
#       c(abs(rnorm(nReps, tmp[1], sqrt(tmp[2]))),      # ksi_r,
#       rgamma(nReps, tmp[3], tmp[4]),                  # sig_r
#       rbeta(nReps, tmp[5]*tmp[6], (1-tmp[5])*tmp[6]), # zeta_r
#       tmp)
#   params[1,] = tmp.params



    lower = c(rep(-Inf, nReps), rep(0, nReps), rep(0, nReps),
        -Inf, 0, 0, 0, 0, 0)
    upper = c(rep(Inf, nReps), rep(Inf, nReps), rep(1, nReps),
        Inf, Inf, Inf, Inf, 1, Inf)

#   lower = c(rep(-Inf, nReps), rep(0, nReps), rep(0, nReps), rep(0, nReps),
#       -Inf, 0, 0, 0, 0, 0, 0, 0)
#   upper = c(rep(Inf, nReps), rep(Inf, nReps), rep(1, nReps), rep(1, nReps),
#       Inf, Inf, Inf, Inf, 1, Inf, 1, Inf)

#   lower = c(rep(-1, nReps), rep(0, nReps), rep(0, nReps),
#       -1, 0, 0, 0, 0, 0)
#   upper = c(rep(1, nReps), rep(Inf, nReps), rep(1, nReps),
#       1, Inf, Inf, Inf, 1, Inf)

    ldigamma = function(x, alpha, beta)
        alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(x) - beta / x

    calc.post = function(x, param){
        nReps = NCOL(x$varmat)
        ksi.vec = param[1:nReps + nReps*0]
        sig.vec = param[1:nReps + nReps*1]
        zeta.vec = param[1:nReps + nReps*2]
#       theta.vec = param[1:nReps + nReps*3]
        tmp = tail(c(param), 6)
        ksi = tmp[1]
        tau2 = tmp[2]
        sig.alpha = tmp[3]
        sig.beta = tmp[4]
        zeta.mu = tmp[5]
        zeta.eta = tmp[6]
#       theta.mu = tmp[7]
#       theta.eta = tmp[8]

        ind = which(ksi.vec < 0)
        if (any(sapply(x$y[ind], max) >= -sig.vec[ind] / ksi.vec[ind]))
            return (-Inf)

        # Likelihood
        out = 0
        out = out + sum(x$n_u * log(zeta.vec)) + sum((x$n - x$n_u) * log(1 - zeta.vec))
#       out = out + sum(x$n_c * log(theta.vec)) + sum((x$n_u - x$n_c) * log(1 - theta.vec))
        out = out - sum(x$n_c * log(sig.vec))
        for (r in 1:nReps){
            if (ksi.vec[r] != 0){
                out = out - (1 + 1/ksi.vec[r])*sum(log(1+ksi.vec[r]*x$y[[r]]/sig.vec[r]))
            } else {
                out = out - sum(x$y[[r]])/sig.vec[r]
                }
            }

        # Priors
        out = out + sum(dbeta(zeta.vec, zeta.mu*zeta.eta, (1-zeta.mu)*zeta.eta, log = TRUE))
#       out = out + sum(dbeta(theta.vec, theta.mu*theta.eta, (1-theta.mu)*theta.eta, log = TRUE))
        out = out + sum(dnorm(ksi.vec, ksi, sqrt(tau2), log = TRUE))
#       out = out + sum(log(dtruncnorm(ksi.vec, a = -1, b = 1, mean = ksi, sd = sqrt(tau2))))
#       out = out + sum(ldigamma(sig.vec, sig.alpha, sig.beta))
        out = out + sum(dgamma(sig.vec, sig.alpha, sig.beta, log = TRUE))

        out = out +    dnorm(ksi,       ksi.mean,    ksi.sd, log = TRUE)
#       out = out + log(dtruncnorm(ksi, a = -1, b = 1, mean = ksi.mean, sd = ksi.sd))
#       out = out + ldigamma(tau2,      tau2.alpha,  tau2.beta)
        out = out +   dgamma(tau2,      tau2.alpha,  tau2.beta, log = TRUE)
        out = out +   dgamma(sig.alpha, sig.alpha.a, sig.alpha.b, log = TRUE)
        out = out +   dgamma(sig.beta,  sig.beta.a,  sig.beta.b, log = TRUE)
        out = out +    dbeta(zeta.mu,   zeta.mu.a,   zeta.mu.b, log = TRUE)
        out = out +   dgamma(zeta.eta,  zeta.eta.a,  zeta.eta.b, log = TRUE)
#       out = out +  dbeta(theta.mu,  theta.mu.a,  theta.mu.b, log = TRUE)
#       out = out + dgamma(theta.eta, theta.eta.a, theta.eta.b, log = TRUE)
        if (is.na(out))
            out = -Inf
        return (out)
        }

    output$post.vals = double(nburn + nmcmc)
    post = calc.post(output, params[1,])
    if (!is.finite(post)){
        params[1,] = rep(0.5, nparam)
        post = calc.post(output, params[1,])
        }
    output$post.vals[1] = post

    trailing = function(x, digits = 4)
        formatC(x, digits=digits, format="f")
    nice_time = function(seconds){
        # floor() or round() would work as well
        seconds = ceiling(seconds)
        days = seconds %/% (60*60*24)
        seconds = seconds %% (60*60*24)
        hours = seconds %/% (60*60)
        seconds = seconds %% (60*60)
        minutes = seconds %/% (60)
        seconds = seconds %% (60)
        out = ""
        if (days > 0)
            out = paste0(out, days, "d ", hours, "h ", minutes, "m ", seconds, "s")
        if (days == 0 && hours > 0)
            out = paste0(out, hours, "h ", minutes, "m ", seconds, "s")
        if (days == 0 && hours == 0 && minutes > 0)
            out = paste0(out, minutes, "m ", seconds, "s")
        if (days == 0 && hours == 0 && minutes == 0)
            out = paste0(out, seconds, "s")
        return (out)
        }
#   message("Obtaining posterior samples ...")
    display = 1000
    begin_time = as.numeric(Sys.time())
    for (i in 2:(nburn + nmcmc)){
        if (floor(i/display) == i/display && display > 0){
            curr_time = as.numeric(Sys.time()) - begin_time
            cat("\r   ", i, " / ", nburn+nmcmc, " -- ",
                trailing(100*i/(nburn+nmcmc), 2),"% -- Remaining: ",
                nice_time(curr_time*(nmcmc+nburn-i)/(i-1)), "            ", sep = "")
            }
        params[i,] = params[i-1,]
        cand = mvrnorm(1, params[i-1,], cand.sig)
        if (all(cand > lower) && all(cand < upper)){
            cand.post = calc.post(output, cand)
            if (log(runif(1)) <= cand.post - post){
                post = cand.post
                params[i,] = cand
                accept[i] = 1
                }
            }
        output$post.vals[i] = post
        if ((floor(i/window) == i/window) && (i <= nburn))
            cand.sig = autotune(mean(accept[(i-window+1):i]), target = 0.234, k = window/50) *
                (cand.sig + window * var(params[(i-window+1):i,]) / i)
        if (i == (nburn + nmcmc) && display > 0){
            curr_time = as.numeric(Sys.time()) - begin_time
            cat("\r   ", i, " / ", nburn+nmcmc, " -- ",
                trailing(100, 2),"% -- Elapsed: ",
                nice_time(curr_time), "            \n", sep = "")
            }
        }

    output$params = tail(params, nmcmc)
    output$accept = mean(tail(accept, nmcmc))
    output$post.vals = tail(output$post.vals, nmcmc)

#   output$params = params
#   output$accept = accept
#   output$post.vals = output$post.vals

    output$means = c(colMeans(output$params[,1:30]), 0, 0)
    output$hpds = cbind(apply(output$params[,1:30], 2, hpd_mult, force_uni = TRUE), 0, 0)
    output$cand.sig = cand.sig
    output$nmcmc = nmcmc
    output$nburn = nburn

#   output$ksi = cbind(params[,1:nReps], params[,3*nReps + 1])
#   output$sigma = cbind(params[,(nReps+1):(2*nReps)],
#       params[,3*nReps+3] / params[,3*nReps+4])
#   output$zeta = cbind(params[,(2*nReps+1):(3*nReps)], params[,3*nReps+

#   tmp = output$theta
#   output$theta = matrix(0, nmcmc, nReps+1)
#   for (i in 1:(nReps+1))
#       output$theta[,i] = sample(tmp[,i], nmcmc, replace = TRUE)

#   output$rl20 = matrix(0, nmcmc, nReps+1)

#   ksi.star = x[[rr]]$params[,31]
#   sig.star = x[[rr]]$params[,33]/x[[rr]]$params[,34]
#   zeta.star = x[[rr]]$params[,35]
#   param.star = cbind(ksi.star, sig.star, zeta.star)#, theta.star)
#   ndays = x[[rr]]$days.in.months
#   ZZ = apply(param.star, 1,
#       function(y)
#           x[[rr]]$threshold[1] + y[2]/y[1] * (((return.period*ndays)*y[3])^y[1] - 1)
#       )
#   Zq = apply(ZZ, 1, quantile, c(0.025, 0.975))

    return (output)
    }
