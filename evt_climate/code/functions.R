### Obtain posterior parameters for a normal dlm
dlm_get_vals = function(dat, delta, m0, C0, n0, d0){
    y = dat$y
    F.vec = dat$F
    G.mat = dat$G
    n = length(y)
    p = length(m0)

    Rt = array(0, c(n+1, p, p))
    qt = double(n+1)
    At = matrix(0, n+1, p)
    ft = double(n+1)
    et = double(n+1)
    mt = matrix(0, n+1, p)
    nt = double(n+1)
    dt = double(n+1)
    st = double(n+1)
    Ct = array(0, c(n+1, p, p))

    mt[1,] = m0
    Ct[1,,] = C0
    nt[1] = n0
    dt[1] = d0
    st[1] = dt[1] / nt[1]

    for (i in 2:(n+1)){
#       Rt[i,,] = G.mat %*% Ct[i-1,,] %*% t(G.mat) + (1-delta)/delta*Ct[i-1,,]
        Rt[i,,] = 1/delta * G.mat %*% Ct[i-1,,] %*% t(G.mat)
        qt[i] = t(F.vec) %*% Rt[i,,] %*% F.vec + st[i-1]
        At[i,] = (Rt[i,,] %*% F.vec)/qt[i]   
        ft[i] = t(F.vec) %*% G.mat %*% mt[i-1,]
        et[i] = y[i-1] - ft[i]
        mt[i,] = G.mat %*% mt[i-1,] + At[i,] * et[i]
        nt[i] = nt[i-1] + 1
        dt[i] = dt[i-1] + st[i-1]*(et[i]^2) / qt[i]
        st[i] = st[i-1] + st[i-1] / nt[i]*( (et[i]^2) / qt[i] - 1)
        Ct[i,,] = (st[i] / st[i-1])*(Rt[i,,] - At[i,] %*% t(At[i,]) * qt[i])
        }

    
    return (list("Rt"=Rt[-1,,], "qt"=qt[-1], "At"=At[-1,],
        "ft"=ft[-1], "et"=et[-1], "mt"=mt[-1,], "nt"=nt[-1],
        "dt"=dt[-1], "st"=st[-1], "Ct"=Ct[-1,,], "delta"=delta,
        "m0"=m0, "C0"=C0, "n0"=n0, "d0"=d0))
    }

### Smoothing a dlm (conditioning on the last value and going backward)
dlm_smooth = function(dat, params){
    y = dat$y
    F.vec = dat$F
    G.mat = dat$G
    G.inv = solve(G.mat)
    n = length(y)
    p = length(F.vec)
    delta = params$delta
    out.at = matrix(0, n, p)
    out.Rt = array(0, c(n, p, p))
    out.at[n,] = params$mt[n,]
    out.Rt[n,,] = params$Ct[n,,]
    out.ft = double(n)
    out.qt = double(n)
    out.ft[n] = t(F.vec) %*% out.at[n,]
    out.ft[n] = t(F.vec) %*% out.at[n,]
    for (i in (n-1):1){
        out.at[i,] = (1-delta)*params$mt[i,] + delta*(G.inv %*% out.at[i+1,])
        out.Rt[i,,] = (1-delta)*params$Ct[i,,] + delta^2*(G.inv %*% out.Rt[i+1,,] %*% t(G.inv))
        out.ft[i] = t(F.vec) %*% out.at[i,]
        out.qt[i] = t(F.vec) %*% out.Rt[i,,] %*% F.vec + params$st[n]
        }
    return (list("at.s"=out.at, "Rt.s"=out.Rt, "ft.s"=out.ft, "qt.s"=out.qt))
    }


### Return the index of a vector of dates (ts) which contain months
subset_data = function(ts, months){
    # Make things easier to work with in grep
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

### Initialize R object
clim_init = function(data_name, data_dir, region, variable, season, months,
    anomaly, date_begin, date_end){

    output = NULL
    output$model = data_name
    output$data_dir = data_dir

    output$region = region
    output$variable = variable

    # If 'months' not specified, get them based off 'season'
    # season[1] is either 'year' (default) or user input
    if (is.null(months)){
        output$season = season[1]
        output$months = switch(output$season,
            "winter" = c("dec", "jan", "feb"),
            "summer" = c("jun", "jul", "aug"),
            c("jan", "feb", "mar", "apr", "may", "jun",
                "jul", "aug", "sep", "oct", "nov", "dec"))
    } else {
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
    output$date_begin = date_begin
    output$date_end = date_end
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
    return (output)
    }

clim_load_proc = function(data_name, data_dir, region, variable){
    file_name = paste0(data_dir, region, "_", variable, "_", data_name, ".txt")
    dat = read.table(file_name, header = TRUE)

    # Remove leap years
    tmp = grep("-02-29", substr(dat[,1], 5, 10))
    if (length(tmp) > 0)
        dat = dat[-tmp,]

    return (dat)
    }

clim_anomaly_calc = function(dat, date_begin, date_end){
    require(mwBASE)
    require(Matrix)

    # Use a linear trend and 2 seasonal components (first two harmonics)
    Gj = function(j, p)
        matrix(c(cos(2*pi*j/p), -sin(2*pi*j/p), sin(2*pi*j/p), cos(2*pi*j/p)), 2, 2)

    F.vec = matrix(c(c(1,0), rep(c(1,0), 2)), ncol = 1)
    G.mat = as.matrix(bdiag(matrix(c(1,0,1,1),2,2),
        Gj(1, 365), Gj(2, 365)))

    # If available, use earlier and later observations when making the dlm
    dlm_range = c(max(which(date_begin == dat[,1]) - 100, 1),
        min(which(date_end == dat[,1]) + 100, length(dat[,1])))
    dlm_seq = dlm_range[1]:dlm_range[2]

    anomaly_shift = matrix(0, NROW(dat), NCOL(dat)-1)

    for (j in 2:NCOL(dat)){
        tmp = list("y" = dat[dlm_seq, j], "F" = F.vec, "G" = G.mat)
        n = length(dlm_seq)
        p = length(F.vec)
        m0 = double(p)
        m0[1] = mean(tmp$y)
        C0 = var(tmp$y)*diag(p)
        n0 = 1
        d0 = 100
        params = dlm_get_vals(tmp, delta = 0.9999, m0, C0, n0, d0)
        spar = dlm_smooth(tmp, params)

        dat[dlm_seq,j] = dat[dlm_seq,j] - spar$ft.s
        anomaly_shift[dlm_seq,j-1] = spar$ft.s
        }

    return (list("dat" = dat, "anomaly_shift" = anomaly_shift))
    }

clim_same_time = function(model, dat, date_begin, date_end, months){
    ### Deal with those leap years and get things all on the same time scale
    if (substr(model, 1, 7) == "control"){
        season.ind = subset_data(dat[,1], months)
        dat[,1] = paste0(as.numeric(substr(dat[,1], 1, 4)) +
            as.numeric(substr(date_begin, 1, 4)) - 1, substr(dat[,1], 5, 10))
        dec.ind = c(min(grep(date_begin, dat[,1])), max(grep(date_end, dat[,1])))
        season.ind = season.ind[season.ind >= dec.ind[1] & season.ind <= dec.ind[2]]
    } else {
        season.ind = subset_data(dat[,1], months)
        ly.ind = grep("-02-29", dat[,1])
        season.ind = season.ind[!(season.ind %in% ly.ind)]
        dec.ind = c(min(grep(date_begin, dat[,1])), max(grep(date_end, dat[,1])))
        season.ind = season.ind[season.ind >= dec.ind[1] & season.ind <= dec.ind[2]]
        }
    return (season.ind)
    }

wrapper_theta = function(x){
    tmp = theta_hier(y = x$varmat, u = x$threshold[1], likelihood = "ferro",
        nburn = x$nburn, nmcmc = x$nmcmc)

    out = NULL

    out$y = rep(list(NULL), x$R)
    for (j in 1:x$R)
        out$y[[j]] = tmp$y[[j]] - x$threshold[j]

    out$n_c = sapply(tmp$y, length)
    out$n_u = as.vector(tmp$N)
    out$n = rep(NROW(x$varmat), x$R)
    out$ymax = sapply(out$y, max)

    out$T_C = tmp$T_C
    out$theta = tmp$mcmc
    out$theta.mean = colMeans(tmp$mcmc)
    out$theta.ints = apply(tmp$mcmc, 2, hpd_mult, force_uni = TRUE)

    return (out)
    }

wrapper_mcmc = function(x){

    uq = x$uq
    nburn = x$nburn
    nmcmc = x$nmcmc
    nReps = x$R
    out = NULL

    ### Priors
       ksi.mean = 0
         ksi.sd = 0.33
     tau2.alpha = 1
      tau2.beta = .1
    sig.alpha.a = 1
    sig.alpha.b = .1
     sig.beta.a = 1
     sig.beta.b = .1
      zeta.mu.a = (1 - uq)*10
      zeta.mu.b = uq*10
     zeta.eta.a = 1
     zeta.eta.b = .1


    ### Parameters
    par.vec.ksi.sig = matrix(0, nburn + nmcmc, 2*nReps)
    par.vec.zeta = matrix(0, nburn + nmcmc, nReps)
    par.ksi = double(nburn + nmcmc)
    par.tau2 = double(nburn + nmcmc)
    par.alpha.beta = matrix(0, nburn + nmcmc, 2)
    par.zeta.eta = matrix(0, nburn + nmcmc, 2)

    # Initial values
    par.vec.ksi.sig[1,] = c(rep(0, nReps), rep(1, nReps))
    par.vec.zeta[1,] = rep(1-uq, nReps)
    par.ksi[1] = 0
    par.tau2[1] = 1
    par.alpha.beta[1,] = c(1, 1)
    par.zeta.eta[1,] = c(1-uq, 1)

    # Acceptances for walks
    accept.vec.ksi.sig = double(nburn + nmcmc)
    accept.vec.zeta = double(nburn + nmcmc)
    accept.alpha.beta = double(nburn + nmcmc)
    accept.zeta.eta = double(nburn + nmcmc)

    accept.tau2 = double(nburn + nmcmc)

    # Candidate sigmas
    cov.vec.ksi.sig = diag(0.01, 2*nReps)
    cov.vec.zeta = diag(0.01, nReps)
    cov.alpha.beta = diag(0.01, 2)
    cov.zeta.eta = diag(0.01, 2)

    cov.tau2 = 0.01

    chol.vKS = chol(cov.vec.ksi.sig)
    chol.vZ = chol(cov.vec.zeta)
    chol.AB = chol(cov.alpha.beta)
    chol.ZE = chol(cov.zeta.eta)

    chol.T = chol(cov.tau2)

    calc.like = function(y, n_c, ymax, ksi, sig){
        if (missing(n_c))
            n_c = length(y)
        if (missing(ymax))
            ymax = max(y)

        if (ksi < 0 && ymax >= -sig/ksi)
            return (-Inf)
        if (sig < 0)
            return (-Inf)

        if (ksi != 0){
            out = -n_c*log(sig) - (1/ksi + 1) * sum(log(1 + ksi * y / sig))
        } else {
            out = -n_c*log(sig) - sum(y) / sig
            }
        return (out)
        }

    
    curr.like.ksi.sig = double(nReps)
    for (i in 1:nReps){
        curr.like.ksi.sig[i] = calc.like(y = x$y[[i]], n_c = x$n_c[i],
            ymax = x$ymax[i], ksi = par.vec.ksi.sig[1, i],
            sig = par.vec.ksi.sig[1, nReps + i])
        }
    curr.prior.ksi = dnorm(par.vec.ksi.sig[1, 1:nReps], par.ksi[1], sqrt(par.tau2[1]), log = TRUE)
    curr.prior.sig = dgamma(par.vec.ksi.sig[1, nReps + 1:nReps],
        par.alpha.beta[1, 1], par.alpha.beta[1, 2], log = TRUE)

    curr.like.zeta = (x$n-x$n_u)*log(1-par.vec.zeta[1,]) + x$n_u*log(par.vec.zeta[1,])
    curr.prior.zeta = dbeta(par.vec.zeta[1,], par.zeta.eta[1, 1]*par.zeta.eta[1, 2],
        (1-par.zeta.eta[1, 1])*par.zeta.eta[1, 2], log = TRUE)

    cand.like = double(nReps)

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
    display = 1000
    begin_time = as.numeric(Sys.time())
    for (i in 2:(nburn + nmcmc)){
        if (floor(i/display) == i/display && display > 0){
            curr_time = as.numeric(Sys.time()) - begin_time
            cat("\r   ", i, " / ", nburn+nmcmc, " -- ",
                trailing(100*i/(nburn+nmcmc), 2),"% -- Remaining: ",
                nice_time(curr_time*(nmcmc+nburn-i)/(i-1)), "            ", sep = "")
            }

        par.vec.ksi.sig[i,] = par.vec.ksi.sig[i-1,]
        par.vec.zeta[i,] = par.vec.zeta[i-1,]
        par.ksi[i] = par.ksi[i-1]
        par.tau2[i] = par.tau2[i-1]
        par.alpha.beta[i,] = par.alpha.beta[i-1,]
        par.zeta.eta[i,] = par.zeta.eta[i-1,]

        # Update (ksi_i, sig_i)
        cand = rnorm(2*nReps) %*% chol.vKS + par.vec.ksi.sig[i,]
        ind = which(cand[1:nReps] < 0)
        if ((length(ind) == 0) || (x$ymax < -cand[nReps + 1:nReps] / cand[1:nReps])[ind]){
            for (j in 1:nReps)
                cand.like[j] = calc.like(y = x$y[[j]], ksi = cand[j], sig = cand[j + nReps])
            cand.prior.ksi = dnorm(cand[1:nReps], par.ksi[i], sqrt(par.tau2[i]), log = TRUE)
            cand.prior.sig = dgamma(cand[nReps + 1:nReps], par.alpha.beta[i, 1],
                par.alpha.beta[i, 2], log = TRUE)

            if (log(runif(1)) < (sum(cand.like + cand.prior.ksi + cand.prior.sig) - 
                sum(curr.like.ksi.sig + curr.prior.ksi + curr.prior.sig))){
                curr.like.ksi.sig = cand.like
                curr.prior.ksi = cand.prior.ksi
                curr.prior.sig = cand.prior.sig
                accept.vec.ksi.sig[i] = 1
                par.vec.ksi.sig[i,] = cand
                }
            }

        # Update zeta_i
        cand = rnorm(nReps) %*% chol.vZ + par.vec.zeta[i,]
        if (all(cand > 0 & cand < 1)){
            cand.like = (x$n-x$n_u)*log(1-cand) + x$n_u*log(cand)
            cand.prior = dbeta(cand, par.zeta.eta[i, 1]*par.zeta.eta[i, 2],
                (1-par.zeta.eta[i, 1])*par.zeta.eta[i, 2], log = TRUE)
            if (log(runif(1)) < (sum(cand.like + cand.prior) -
                sum(curr.like.zeta + curr.prior.zeta))){
                curr.like.zeta = cand.like
                curr.prior.zeta = cand.prior
                accept.vec.zeta[i] = 1
                par.vec.zeta[i,] = cand
                }
            }

        # Update ksi
        m.star = (par.tau2[i] * ksi.mean + ksi.sd^2*sum(par.vec.ksi.sig[i, 1:nReps])) /
            (par.tau2[i] + ksi.sd^2 * nReps)
        s.star = sqrt((par.tau2[i] * ksi.sd^2) / (par.tau2[i] + ksi.sd^2 * nReps))
        par.ksi[i] = rnorm(1, m.star, s.star)

        curr.prior.ksi = dnorm(par.vec.ksi.sig[i, 1:nReps],
            par.ksi[i], sqrt(par.tau2[i]), log = TRUE)

        # Update tau2 (variance, Ga prior) on log-scale
        cand = rnorm(1) * chol.T + log(par.tau2[i])
        if (is.finite(cand)){
            cand.prior.ksi = dnorm(par.vec.ksi.sig[i, 1:nReps],
                par.ksi[i], sqrt(exp(cand)), log = TRUE)
            if (log(runif(1)) <
                (sum(cand.prior.ksi) + dgamma(exp(cand), tau2.alpha, tau2.beta, log = TRUE) + cand) -
                (sum(curr.prior.ksi) + dgamma(par.tau2[i], tau2.alpha, tau2.beta, log = TRUE) + 
                    log(par.tau2[i]))){
                accept.tau2[i] = 1
                curr.prior.ksi = cand.prior.ksi
                par.tau2[i] = exp(cand)
                }
            }

        # Update (alpha, beta)
        cand = rnorm(2) %*% chol.AB + par.alpha.beta[i,]
        if (all(cand > 0)){
            cand.prior = dgamma(par.vec.ksi.sig[i, nReps + 1:nReps], cand[1], cand[2], log = TRUE)
            if (log(runif(1)) <
                (sum(cand.prior) +
                    dgamma(cand[1], sig.alpha.a, sig.alpha.b, log = TRUE) +
                    dgamma(cand[2], sig.beta.a, sig.beta.b, log = TRUE)) -
                (sum(curr.prior.sig) +
                    dgamma(par.alpha.beta[i, 1], sig.alpha.a, sig.alpha.b, log = TRUE) +
                    dgamma(par.alpha.beta[i, 2], sig.beta.a, sig.beta.b, log = TRUE))){
                curr.prior.sig = cand.prior
                accept.alpha.beta[i] = 1
                par.alpha.beta[i,] = cand
                }
            }

        # Update (zeta, eta)
        cand = rnorm(2) %*% chol.ZE + par.zeta.eta[i,]
        if (all(cand > 0) && cand[1] < 1){
            cand.prior = dbeta(par.vec.zeta[i,], cand[1]*cand[2], (1-cand[1])*cand[2], log = TRUE)
            if (log(runif(1)) <
                (sum(cand.prior) +
                    dbeta(cand[1], zeta.mu.a, zeta.mu.b, log = TRUE) +
                    dgamma(cand[2], zeta.eta.a, zeta.eta.b, log = TRUE)) -
                (sum(curr.prior.zeta) +
                    dbeta(par.zeta.eta[i, 1], zeta.mu.a, zeta.mu.b, log = TRUE) +
                    dgamma(par.zeta.eta[i, 2], zeta.eta.a, zeta.eta.b, log = TRUE))){
                curr.prior.zeta = cand.prior
                accept.zeta.eta[i] = 1
                par.zeta.eta[i,] = cand
                }
            }

        # Change proposal variances
        if ((floor(i/window) == i/window) && (i <= nburn)){
            cov.vec.ksi.sig = autotune(mean(accept.vec.ksi.sig[(i-window+1):i]),
                target = 0.234, k = window/50) *
                (cov.vec.ksi.sig + window * var(par.vec.ksi.sig[(i-window+1):i,]) / i)
            cov.vec.zeta = autotune(mean(accept.vec.zeta[(i-window+1):i]),
                target = 0.234, k = window/50) *
                (cov.vec.zeta + window * var(par.vec.zeta[(i-window+1):i,]) / i)
            cov.alpha.beta = autotune(mean(accept.alpha.beta[(i-window+1):i]),
                target = 0.234, k = window/50) *
                (cov.alpha.beta + window * var(par.alpha.beta[(i-window+1):i,]) / i)
            cov.zeta.eta = autotune(mean(accept.zeta.eta[(i-window+1):i]),
                target = 0.234, k = window/50) *
                (cov.zeta.eta + window * var(par.zeta.eta[(i-window+1):i,]) / i)

            cov.tau2 = autotune(mean(accept.tau2[(i-window+1):i]),
                target = 0.234, k = window/50) *
                (cov.tau2 + window * var(par.tau2[(i-window+1):i]) / i)

            chol.vKS = chol(cov.vec.ksi.sig)
            chol.vZ = chol(cov.vec.zeta)
            chol.AB = chol(cov.alpha.beta)
            chol.ZE = chol(cov.zeta.eta)

            chol.T = chol(cov.tau2)
            }
        if (i == (nburn + nmcmc) && display > 0){
            curr_time = as.numeric(Sys.time()) - begin_time
            cat("\r   ", i, " / ", nburn+nmcmc, " -- ",
                trailing(100, 2),"% -- Elapsed: ",
                nice_time(curr_time), "            \n", sep = "")
            }
        }

    par.vec.ksi.sig = tail(par.vec.ksi.sig, nmcmc)
    par.vec.zeta = tail(par.vec.zeta, nmcmc)
    par.ksi = tail(par.ksi, nmcmc)
    par.tau2 = tail(par.tau2, nmcmc)
    par.alpha.beta = tail(par.alpha.beta, nmcmc)
    par.zeta.eta = tail(par.zeta.eta, nmcmc)

    accept.vec.ksi.sig = tail(accept.vec.ksi.sig, nmcmc)
    accept.vec.zeta = tail(accept.vec.zeta, nmcmc)
    accept.alpha.beta = tail(accept.alpha.beta, nmcmc)
    accept.zeta.eta = tail(accept.zeta.eta, nmcmc)

    accept.tau2 = tail(accept.tau2, nmcmc)

    colnames(par.vec.ksi.sig) = c(paste0("ksi", 1:nReps), paste0("sig", 1:nReps))
    colnames(par.vec.zeta) = paste0("zeta", 1:nReps)
    colnames(par.alpha.beta) = c("alpha", "beta")
    colnames(par.zeta.eta) = c("zeta", "eta")
    
    params = cbind(par.vec.ksi.sig, par.vec.zeta, "ksi_mean" = par.ksi,
        "tau2" = par.tau2, par.alpha.beta, par.zeta.eta)
    accept = c("vec.ksi.sig" = mean(accept.vec.ksi.sig), "vec.zeta" = mean(accept.vec.zeta),
        "tau2" = mean(accept.tau2), "alpha.beta" = mean(accept.alpha.beta),
        "zeta.eta" = mean(accept.zeta.eta))

    return (list("params" = params, "accept" = accept))
    }

hier_excess = function(data_name, data_dir, region, variable,
    season = c("year", "summer", "winter"), date_begin, date_end,
    r = 0, uq = 0.90, min_uq = TRUE, threshold,
    nburn = 80000, nmcmc = 40000, window = 500,
    months = NULL, anomaly = TRUE){

    require(MASS)
    require(mwBASE)
    require(mwEVT)

    ### Initialize output object
    output = clim_init(data_name, region, variable, season, months, anomaly)

    ### Load the processed data
    dat = clim_load_proc(data_name, data_dir, region, variable)

    ### Compute the anomalies using a DLM
    if (output$anomaly){
        output$anomaly_type = "dlm"
        tmp = clim_anomaly_calc(dat, date_begin, date_end)
        dat = tmp$dat
        output$anomaly_shift = tmp$anomaly_shift
        rm(tmp)
        }

    ### Get appropriate time index
    ind = clim_same_time(output$model, dat, output$date_begin,
        output$date_end, output$months)

    output$time.dates = as.Date(dat[ind,1])
    output$varmat = as.matrix(dat[ind,-1])
    output$anomaly_shift = output$anomaly_shift[ind,]

    ### Threshold
    # Number of full years from date_begin to date_end
    output$R = NCOL(output$varmat)
    output$nyears = length(output$time.dates) / output$days.in.months

    if (missing(threshold)){
        output$uq = uq
        if (min_uq){
            threshold = min(apply(output$varmat, 2, quantile, output$uq))
            output$threshold = rep(threshold, output$R)
        } else {
            if (length(output$uq) == 1)
                output$uq = rep(output$uq, output$R)
            output$threshold = diag(apply(output$varmat, 2, quantile, output$uq))
            }
    } else {
        output$uq = NULL
        output$threshold = threshold
        if (length(output$threshold) == 1)
            output$threshold = rep(output$threshold, output$R)
        }

    ### Extremal index and declustering (retain index of exceedances for bivariate declustering)
    output$nburn = nburn
    output$nmcmc = nmcmc

    tmp = wrapper_theta(output)
    output = c(output, tmp)
    rm(tmp)

    ### MCMC (hierarchical model)
    tmp = wrapper_mcmc(output)
    output = c(output, tmp)
    rm(tmp)


    ### Compute 20- and 50- year return levels


    ### Compute within replicate Bhattacharyya distance (bootstrap this maybe?)


    ### Output (make sure things are all good for bivariate/multivariate analysis)

    }
