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
        }

    ### Get appropriate time index
    ind = clim_same_time(output$model, dat, output$date_begin,
        output$date_end, output$months)

    output$time.dates = as.Date(dat[ind,1])
    output$varmat = as.matrix(dat[ind,-1])
    output$anomaly_shift = output$anomaly_shift[ind,]
