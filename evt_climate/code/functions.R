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
