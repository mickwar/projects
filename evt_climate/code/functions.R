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

### Deal with those leap years and get things all on the same time scale
clim_same_time = function(data_name, dat, date_begin, date_end, m)
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
        clim_anomaly_calc(dat, date_begin, date_end)    # makes changes to dat, need to make assignments
        }
