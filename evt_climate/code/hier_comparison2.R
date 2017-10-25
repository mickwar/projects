library(mwBASE)
rgpd = function(n, mu, ksi, sigma){
    if (length(ksi) == 1)
        ksi = rep(ksi, n)
    return (ifelse(ksi == 0, mu - sigma * log(runif(n)),
        mu + sigma * (runif(n)^(-ksi) - 1) / ksi))
    }
bhatta.dist = function(x, y, ...){
    # Bhattacharyya distance for two continuous random variables x and y
    left = max(min(x), min(y))
    right = min(max(x), max(y))

    if (left > right)
        return (Inf)

    dx = density(x, from = left, to = right, ...)
    dy = density(y, from = left, to = right, ...)

    xx = dx$x
    p = dx$y
    q = dy$y

    BC = sum(mean(diff(xx)) * sqrt(p*q))
    BC = ifelse(BC < 0, 0, BC)
    BC = ifelse(BC >= 1, 1, BC)

    DB = -log(BC)
    
    return (DB)
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

tab.theta.mm = sapply(x, function(y) tail(y$theta.mean, 1))
tab.theta.qq = t(sapply(x, function(y) c(tail(t(y$theta.ints), 1))))
tab.threshold = sapply(x, function(y) y$threshold[1])
tab.T_C = sapply(x, function(y) mean(y$T_C))

R = 10  # Number of replicates
width = 9
height = 9


### Ksi
pdf(paste0(OUT_DIR, "shape.pdf"), height = height, width = width)

ylim = c(0, 5)
xlim = range(sapply(x[-ind.obs], function(y) hpd_mult(y$params[,31], force_uni = TRUE)))
xlim = range(xlim, range(sapply(x[ind.obs], function(y) hpd_mult(y$params[,1], force_uni = TRUE))))

if (!exists("bh.list.ksi"))
    bh.list.ksi = rep(list(matrix(0, R+2, 3)), 16)
tmp.mean = double(3)

pnum = 0
lwd = 3
cex = 1.5
par(mfrow = c(4, 4), mar = c(0,0,0,0), oma = c(6,10,10,6))
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
                for (A5 in 1:4){
                    B5 = list(ind.control, ind.decadal, ind.historical, ind.obs)[[A5]]
                    rr = Reduce(intersect, list(B1, B2, B3, B4, B5))
                    if (length(rr) == 0) next
                    if (A5 != 4){
                        p = x[[rr]]$params[,31]
                        if (all(bh.list.ksi[[pnum]][,A5] == 0)){
                            bh.list.ksi[[pnum]][2, A5] = rr
                            bh.list.ksi[[pnum]][3:(R+2),A5] = apply(x[[rr]]$params[,1:10], 2,
                                function(y) bhatta.dist(y, p, n = 2^12))
                            }
                        tmp.mean[A5] = mean(p)
                    } else {
                        p = x[[rr]]$params[,1]
                        if (all(bh.list.ksi[[pnum]][1,] == 0)){
                            bh.list.ksi[[pnum]][1, 1] = bhatta.dist(p,
                                x[[bh.list.ksi[[pnum]][2, 1]]]$params[,31], n = 2^12)
                            bh.list.ksi[[pnum]][1, 2] = bhatta.dist(p,
                                x[[bh.list.ksi[[pnum]][2, 2]]]$params[,31], n = 2^12)
                            bh.list.ksi[[pnum]][1, 3] = bhatta.dist(p,
                                x[[bh.list.ksi[[pnum]][2, 3]]]$params[,31], n = 2^12)
                            }
                        }


                    z = hpd_mult(p, force_uni = TRUE)

                    lines(z, rep(A5, 2), lwd = lwd, col = cols[A5])
                    points(mean(p), A5, col = cols[A5], pch = 16, cex = cex)
#                   text(mean(p), A5+0.35, label = round(tab.threshold[rr], 3))
#                   text(mean(p), A5+0.35, label = round(tab.threshold[rr], 3))
                    }
                tmp = 1*apply(bh.list.ksi[[pnum]], 2, function(y)
                    (y[1] >= min(tail(y, R))) && (y[1] <= max(tail(y, R))))
                text(tmp.mean, (1:3)+0.35, label = tmp)
                if ((pnum <= 4) && (pnum %% 2 == 1))
                    axis(3, lwd = -1, lwd.ticks = 1)
                if ((pnum >= 13) && (pnum %% 2 == 0))
                    axis(1, lwd = -1, lwd.ticks = 1)
                if (pnum %% 4 == 0)
                    axis(side = 4, labels = lab.shortmod, at = 1:4, las = 1,
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
mtext(expression(xi), outer = TRUE, side = 3, at = -0.120, line = 4, cex = 2.5, adj = 0)
#par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0,0,0,0))
dev.off()

# par(mfrow = c(1,3), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0,0,0,0))
# pnum = 4
# bh.cont = apply(x[[0 + pnum]]$params[,1:10], 2,
#     function(y) bhatta.dist(y, x[[0 + pnum]]$params[,31], n = 2^12))
# bh.dec = apply(x[[16 + pnum]]$params[,1:10], 2,
#     function(y) bhatta.dist(y, x[[16 + pnum]]$params[,31], n = 2^12))
# bh.hist = apply(x[[32 + pnum]]$params[,1:10], 2, 
#     function(y) bhatta.dist(y, x[[32 + pnum]]$params[,31], n = 2^12))
# bh.obs.cont = bhatta.dist(x[[48 + pnum]]$params[,1], x[[0 + pnum]]$params[,31], n = 2^12)
# bh.obs.dec = bhatta.dist(x[[48 + pnum]]$params[,1], x[[16 + pnum]]$params[,31], n = 2^12)
# bh.obs.hist = bhatta.dist(x[[48 + pnum]]$params[,1], x[[32 + pnum]]$params[,31], n = 2^12)
# 
# plot(c(bh.cont, bh.obs.cont))
# plot(c(bh.dec, bh.obs.dec))
# plot(c(bh.hist, bh.obs.hist))




### log Sigma
pdf(paste0(OUT_DIR, "log_sigma.pdf"), height = height, width = width)

ylim = c(0, 5)
xlim = range(sapply(x[-ind.obs], function(y)
    hpd_mult(log(y$params[,33]/y$params[,34]), force_uni = TRUE)))
xlim = range(xlim, range(sapply(x[ind.obs], function(y) hpd_mult(log(y$params[,2]), force_uni = TRUE))))

if (!exists("bh.list.sigma"))
    bh.list.sigma = rep(list(matrix(0, R+2, 3)), 16)
tmp.mean = double(3)

pnum = 0
lwd = 3
cex = 1.5
par(mfrow = c(4, 4), mar = c(0,0,0,0), oma = c(6,10,10,6))
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
                for (A5 in 1:4){
                    B5 = list(ind.control, ind.decadal, ind.historical, ind.obs)[[A5]]
                    rr = Reduce(intersect, list(B1, B2, B3, B4, B5))
                    if (length(rr) == 0) next
                    if (A5 != 4){
                        p = log(x[[rr]]$params[,33] / x[[rr]]$params[,34])
                        if (all(bh.list.sigma[[pnum]][,A5] == 0)){
                            bh.list.sigma[[pnum]][2, A5] = rr
                            bh.list.sigma[[pnum]][3:(R+2),A5] = apply(log(x[[rr]]$params[,11:20]), 2,
                                function(y) bhatta.dist(y, p, n = 2^12))
                            }
                        tmp.mean[A5] = mean(p)
                    } else {
                        p = log(x[[rr]]$params[,2])
                        if (all(bh.list.sigma[[pnum]][1,] == 0)){
                            tmp.rr = bh.list.sigma[[pnum]][2, 1]
                            bh.list.sigma[[pnum]][1, 1] = bhatta.dist(p,
                                log(x[[tmp.rr]]$params[,33] / x[[tmp.rr]]$params[,34]), n = 2^12)
                            tmp.rr = bh.list.sigma[[pnum]][2, 2]
                            bh.list.sigma[[pnum]][1, 2] = bhatta.dist(p,
                                log(x[[tmp.rr]]$params[,33] / x[[tmp.rr]]$params[,34]), n = 2^12)
                            tmp.rr = bh.list.sigma[[pnum]][2, 3]
                            bh.list.sigma[[pnum]][1, 3] = bhatta.dist(p,
                                log(x[[tmp.rr]]$params[,33] / x[[tmp.rr]]$params[,34]), n = 2^12)
                            }
                        }

                    z = hpd_mult(p, force_uni = TRUE)

                    lines(z, rep(A5, 2), lwd = lwd, col = cols[A5])
                    points(mean(p), A5, col = cols[A5], pch = 16, cex = cex)
                    }
                tmp = 1*apply(bh.list.sigma[[pnum]], 2, function(y)
                    (y[1] >= min(tail(y, R))) && (y[1] <= max(tail(y, R))))
                text(tmp.mean, (1:3)+0.35, label = tmp)

                if ((pnum <= 4) && (pnum %% 2 == 1))
                    axis(3, lwd = -1, lwd.ticks = 1)
                if ((pnum >= 13) && (pnum %% 2 == 0))
                    axis(1, lwd = -1, lwd.ticks = 1)
                if (pnum %% 4 == 0)
                    axis(side = 4, labels = lab.shortmod, at = 1:4, las = 1,
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
mtext(expression(log ~ sigma), outer = TRUE, side = 3, at = -0.120, line = 4, cex = 2.5, adj = 0)
#par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0,0,0,0))
dev.off()



### zeta
pdf(paste0(OUT_DIR, "zeta.pdf"), height = height, width = width)

ylim = c(0, 5)
xlim = range(sapply(x[-ind.obs], function(y) hpd_mult(y$params[,35], force_uni = TRUE)))
xlim = range(xlim, range(sapply(x[ind.obs], function(y) hpd_mult(y$params[,3], force_uni = TRUE))))

if (!exists("bh.list.zeta"))
    bh.list.zeta = rep(list(matrix(0, R+2, 3)), 16)
tmp.mean = double(3)

pnum = 0
lwd = 3
cex = 1.5
par(mfrow = c(4, 4), mar = c(0,0,0,0), oma = c(6,10,10,6))
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
                for (A5 in 1:4){
                    B5 = list(ind.control, ind.decadal, ind.historical, ind.obs)[[A5]]
                    rr = Reduce(intersect, list(B1, B2, B3, B4, B5))
                    if (length(rr) == 0) next
                    if (A5 != 4){
                        p = x[[rr]]$params[,35]
                        if (all(bh.list.zeta[[pnum]][,A5] == 0)){
                            bh.list.zeta[[pnum]][2, A5] = rr
                            bh.list.zeta[[pnum]][3:(R+2),A5] = apply(x[[rr]]$params[,21:30], 2,
                                function(y) bhatta.dist(y, p, n = 2^12))
                            }
                        tmp.mean[A5] = mean(p)
                    } else {
                        p = x[[rr]]$params[,3]
                        if (all(bh.list.zeta[[pnum]][1,] == 0)){
                            tmp.rr = bh.list.zeta[[pnum]][2, 1]
                            bh.list.zeta[[pnum]][1, 1] = bhatta.dist(p,
                                x[[tmp.rr]]$params[,35], n = 2^12)
                            tmp.rr = bh.list.zeta[[pnum]][2, 2]
                            bh.list.zeta[[pnum]][1, 2] = bhatta.dist(p,
                                x[[tmp.rr]]$params[,35], n = 2^12)
                            tmp.rr = bh.list.zeta[[pnum]][2, 3]
                            bh.list.zeta[[pnum]][1, 3] = bhatta.dist(p,
                                x[[tmp.rr]]$params[,35], n = 2^12)
                            }
                        }

                    z = hpd_mult(p, force_uni = TRUE)

                    lines(z, rep(A5, 2), lwd = lwd, col = cols[A5])
                    points(mean(p), A5, col = cols[A5], pch = 16, cex = cex)
                    }
                tmp = 1*apply(bh.list.zeta[[pnum]], 2, function(y)
                    (y[1] >= min(tail(y, R))) && (y[1] <= max(tail(y, R))))
                text(tmp.mean, (1:3)+0.35, label = tmp)

                if ((pnum <= 4) && (pnum %% 2 == 1))
                    axis(3, lwd = -1, lwd.ticks = 1)
                if ((pnum >= 13) && (pnum %% 2 == 0))
                    axis(1, lwd = -1, lwd.ticks = 1)
                if (pnum %% 4 == 0)
                    axis(side = 4, labels = lab.shortmod, at = 1:4, las = 1,
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
mtext(expression(zeta), outer = TRUE, side = 3, at = -0.120, line = 4, cex = 2.5, adj = 0)
#par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0,0,0,0))
dev.off()


### Return level set up
#return.period = sort(unique(c(20, 30, 50, exp(seq(log(0.1), log(100), length = 100)))))
# rl20.bw = 0.50
# rl50.bw = 0.50
return.period = c(20, 50)
rl.all = array(0, dim = c(length(return.period), R+1, x[[1]]$nmcmc, length(files)))
# rl.20 = rep(list(NULL), length(files))
# rl.50 = rep(list(NULL), length(files))

if (!exists("bh.list.rl20"))
    bh.list.rl20 = rep(list(matrix(0, R+2, 3)), 16)
if (!exists("bh.list.rl50"))
    bh.list.rl50 = rep(list(matrix(0, R+2, 3)), 16)


# ylim.20 = rep(list(c(0,0)), 16)
# xlim.20 = rep(list(NULL), 16)
# ylim.50 = rep(list(c(0,0)), 16)
# xlim.50 = rep(list(NULL), 16)
#xlim = range(return.period)
pnum = 0
for (A1 in 1:2){
    B1 = list(ind.cal, ind.usa)[[A1]]
    for (A2 in 1:2){
        B2 = list(ind.62, ind.90)[[A2]]
        for (A3 in 1:2){
            B3 = list(ind.pr, ind.tasmax)[[A3]]
            for (A4 in 1:2){
                B4 = list(ind.winter, ind.summer)[[A4]]
                pnum = pnum + 1
                cat(pnum, "/", 16, "\r")
                for (A5 in 1:4){
                    B5 = list(ind.control, ind.decadal, ind.historical, ind.obs)[[A5]]
                    rr = Reduce(intersect, list(B1, B2, B3, B4, B5))
                    if (length(rr) == 0) next
                    theta = tail(x[[rr]]$theta.mean, 1)
                    ndays = x[[rr]]$days.in.months
                    if (A5 != 4){
                        bh.list.rl20[[pnum]][2, A5] = rr
                        bh.list.rl50[[pnum]][2, A5] = rr
                        ksi.star = x[[rr]]$params[,31]
                        sig.star = x[[rr]]$params[,33]/x[[rr]]$params[,34]
                        zeta.star = x[[rr]]$params[,35]
                        param.star = cbind(ksi.star, sig.star, zeta.star)#, theta.star)
                        ZZ = apply(param.star, 1,
                            function(y)
                                x[[rr]]$threshold[1] + y[2]/y[1] *
                                    (((theta*return.period*ndays)*y[3])^y[1] - 1)
                            )
                        rl.all[,R+1,,rr] = ZZ
                        for (j in 1:R){
                            theta = x[[rr]]$theta.mean[j]
                            ksi.star = x[[rr]]$params[,j]
                            sig.star = x[[rr]]$params[,R+j]
                            zeta.star = x[[rr]]$params[,2*R+j]
                            param.star = cbind(ksi.star, sig.star, zeta.star)#, theta.star)
                            tmp.ZZ = apply(param.star, 1,
                                function(y)
                                    x[[rr]]$threshold[1] + y[2]/y[1] *
                                        (((theta*return.period*ndays)*y[3])^y[1] - 1)
                                )
                            rl.all[,j,,rr] = tmp.ZZ
                            }
                    } else {
                        bh.list.rl20[[pnum]][1,] = rr
                        bh.list.rl50[[pnum]][1,] = rr
                        ZZ = apply(x[[rr]]$params, 1,
                            function(y)
                                x[[rr]]$threshold[1] + y[2]/y[1] *
                                    (((theta*return.period*ndays)*y[3])^y[1] - 1)
                            )
                        rl.all[,R+1,,rr] = ZZ
                        }

                    }
                }
            }
        }
    }

for (i in 1:16){
    rr = bh.list.rl20[[i]][1, 1]
    p3 = rl.all[1,R+1,,rr]

    rr = bh.list.rl20[[i]][2, 1]
    p1 = rl.all[1,-(R+1),,rr]
    p2 = rl.all[1,R+1,,rr]
    bh.list.rl20[[i]][3:(R+2), 1] = apply(p1, 1, function(x) bhatta.dist(x, p2, n = 2^12))
    bh.list.rl20[[i]][1, 1] = bhatta.dist(p2, p3, n = 2^12)
    
    rr = bh.list.rl20[[i]][2, 2]
    p1 = rl.all[1,-(R+1),,rr]
    p2 = rl.all[1,R+1,,rr]
    bh.list.rl20[[i]][3:(R+2), 2] = apply(p1, 1, function(x) bhatta.dist(x, p2, n = 2^12))
    bh.list.rl20[[i]][1, 2] = bhatta.dist(p2, p3, n = 2^12)

    rr = bh.list.rl20[[i]][2, 3]
    p1 = rl.all[1,-(R+1),,rr]
    p2 = rl.all[1,R+1,,rr]
    bh.list.rl20[[i]][3:(R+2), 3] = apply(p1, 1, function(x) bhatta.dist(x, p2, n = 2^12))
    bh.list.rl20[[i]][1, 3] = bhatta.dist(p2, p3, n = 2^12)
    }
for (i in 1:16){
    rr = bh.list.rl50[[i]][1, 1]
    p3 = rl.all[2,R+1,,rr]

    rr = bh.list.rl50[[i]][2, 1]
    p1 = rl.all[2,-(R+1),,rr]
    p2 = rl.all[2,R+1,,rr]
    bh.list.rl50[[i]][3:(R+2), 1] = apply(p1, 1, function(x) bhatta.dist(x, p2, n = 2^12))
    bh.list.rl50[[i]][1, 1] = bhatta.dist(p2, p3, n = 2^12)
    
    rr = bh.list.rl50[[i]][2, 2]
    p1 = rl.all[2,-(R+1),,rr]
    p2 = rl.all[2,R+1,,rr]
    bh.list.rl50[[i]][3:(R+2), 2] = apply(p1, 1, function(x) bhatta.dist(x, p2, n = 2^12))
    bh.list.rl50[[i]][1, 2] = bhatta.dist(p2, p3, n = 2^12)

    rr = bh.list.rl50[[i]][2, 3]
    p1 = rl.all[2,-(R+1),,rr]
    p2 = rl.all[2,R+1,,rr]
    bh.list.rl50[[i]][3:(R+2), 3] = apply(p1, 1, function(x) bhatta.dist(x, p2, n = 2^12))
    bh.list.rl50[[i]][1, 3] = bhatta.dist(p2, p3, n = 2^12)
    }


### Return level - 20 years
# xlim = range(sapply(rl.20[-ind.obs], function(y) hpd_mult(y, force_uni = TRUE)))
# xlim = range(xlim, range(sapply(rl.20[ind.obs], function(y) hpd_mult(y, force_uni = TRUE))))

tmp.x = matrix(0, 2, 2)
tmp.x[1,] = range(sapply(ind.pr, function(x) hpd_mult(rl.all[1,R+1,,x], force_uni = TRUE)))
tmp.x[1,] = range(tmp.x[1,], range(sapply(ind.pr, function(x) hpd_mult(rl.all[2,R+1,,x], force_uni = TRUE))))
tmp.x[2,] = range(sapply(ind.tasmax, function(x) hpd_mult(rl.all[1,R+1,,x], force_uni = TRUE)))
tmp.x[2,] = range(tmp.x[2,], range(sapply(ind.tasmax, function(x) hpd_mult(rl.all[1,R+1,,x], force_uni = TRUE))))

pdf(paste0(OUT_DIR, "rl20.pdf"), height = height, width = width)
ylim = c(0, 5)
tmp.mean = double(3)

pnum = 0
lwd = 3
cex = 1.5
par(mfrow = c(4, 4), mar = c(0,0,0,0), oma = c(6,10,10,6))
for (A1 in 1:2){
    B1 = list(ind.cal, ind.usa)[[A1]]
    for (A2 in 1:2){
        B2 = list(ind.62, ind.90)[[A2]]
        for (A3 in 1:2){
            B3 = list(ind.pr, ind.tasmax)[[A3]]
            for (A4 in 1:2){
                B4 = list(ind.winter, ind.summer)[[A4]]
                pnum = pnum + 1
                plot(0, type='n', xlim = tmp.x[A3,], ylim = ylim, axes = FALSE)
                for (A5 in 1:4){
                    B5 = list(ind.control, ind.decadal, ind.historical, ind.obs)[[A5]]
                    rr = Reduce(intersect, list(B1, B2, B3, B4, B5))
                    if (length(rr) == 0) next
                    if (A5 != 4){
                        p = rl.all[1,R+1,,rr]
                        tmp.mean[A5] = mean(p)
                    } else {
                        p = rl.all[1,R+1,,rr]
                        }

                    z = hpd_mult(p, force_uni = TRUE)

                    lines(z, rep(A5, 2), lwd = lwd, col = cols[A5])
                    points(mean(p), A5, col = cols[A5], pch = 16, cex = cex)
                    }
                tmp = 1*apply(bh.list.rl20[[pnum]], 2, function(y)
                    (y[1] >= min(tail(y, R))) && (y[1] <= max(tail(y, R))))
                text(tmp.mean, (1:3)+0.35, label = tmp)

                if ((pnum <= 4) && (pnum %% 2 == 1))
                    axis(3, lwd = -1, lwd.ticks = 1)
                if ((pnum >= 13) && (pnum %% 2 == 0))
                    axis(1, lwd = -1, lwd.ticks = 1)
                if (pnum %% 4 == 0)
                    axis(side = 4, labels = lab.shortmod, at = 1:4, las = 1,
                        lwd = -1, lwd.ticks = 1)
                box()
                if ((A1 == 1) && (A2 == 2))
                    axis(1, labels = FALSE, lwd.ticks = 0, lwd = 5,
                        at = tmp.x[A3,] + c(-1, 1)*diff(tmp.x[A3,])/5)
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
mtext("RL-20", outer = TRUE, side = 3, at = -0.120, line = 4, cex = 2.5, adj = 0)
#par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0,0,0,0))
dev.off()




### Return level - 50 years
pdf(paste0(OUT_DIR, "rl50.pdf"), height = height, width = width)

ylim = c(0, 5)
# xlim = range(sapply(rl.50[-ind.obs], function(y) hpd_mult(y, force_uni = TRUE)))
# xlim = range(xlim, range(sapply(rl.50[ind.obs], function(y) hpd_mult(y, force_uni = TRUE))))

# tmp.x = matrix(0, 2, 2)
# tmp.x[1,] = range(sapply(ind.pr, function(x) hpd_mult(rl.50[[x]], force_uni = TRUE)))
# tmp.x[2,] = range(sapply(ind.tasmax, function(x) hpd_mult(rl.50[[x]], force_uni = TRUE)))

pnum = 0
lwd = 3
cex = 1.5
par(mfrow = c(4, 4), mar = c(0,0,0,0), oma = c(6,10,10,6))
for (A1 in 1:2){
    B1 = list(ind.cal, ind.usa)[[A1]]
    for (A2 in 1:2){
        B2 = list(ind.62, ind.90)[[A2]]
        for (A3 in 1:2){
            B3 = list(ind.pr, ind.tasmax)[[A3]]
            for (A4 in 1:2){
                B4 = list(ind.winter, ind.summer)[[A4]]
                pnum = pnum + 1
                plot(0, type='n', xlim = tmp.x[A3,], ylim = ylim, axes = FALSE)
                for (A5 in 1:4){
                    B5 = list(ind.control, ind.decadal, ind.historical, ind.obs)[[A5]]
                    rr = Reduce(intersect, list(B1, B2, B3, B4, B5))
                    if (length(rr) == 0) next
                    if (A5 != 4){
                        p = rl.all[2,R+1,,rr]
                        tmp.mean[A5] = mean(p)
                    } else {
                        p = rl.all[2,R+1,,rr]
                        }

                    z = hpd_mult(p, force_uni = TRUE)

                    lines(z, rep(A5, 2), lwd = lwd, col = cols[A5])
                    points(mean(p), A5, col = cols[A5], pch = 16, cex = cex)
                    }
                tmp = 1*apply(bh.list.rl50[[pnum]], 2, function(y)
                    (y[1] >= min(tail(y, R))) && (y[1] <= max(tail(y, R))))
                text(tmp.mean, (1:3)+0.35, label = tmp)
                if ((pnum <= 4) && (pnum %% 2 == 1))
                    axis(3, lwd = -1, lwd.ticks = 1)
                if ((pnum >= 13) && (pnum %% 2 == 0))
                    axis(1, lwd = -1, lwd.ticks = 1)
                if (pnum %% 4 == 0)
                    axis(side = 4, labels = lab.shortmod, at = 1:4, las = 1,
                        lwd = -1, lwd.ticks = 1)
                box()
                if ((A1 == 1) && (A2 == 2))
                    axis(1, labels = FALSE, lwd.ticks = 0, lwd = 5,
                        at = tmp.x[A3,] + c(-1, 1)*diff(tmp.x[A3,])/5)
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
mtext("RL-50", outer = TRUE, side = 3, at = -0.120, line = 4, cex = 2.5, adj = 0)
dev.off()




# ### RL-20 and RL-50 together
# ylim = c(0, 5)
# # xlim = range(sapply(rl.50[-ind.obs], function(y) hpd_mult(y, force_uni = TRUE)))
# # xlim = range(xlim, range(sapply(rl.50[ind.obs], function(y) hpd_mult(y, force_uni = TRUE))))
# 
# tmp.x = matrix(0, 2, 2)
# tmp.x[1,] = range(sapply(ind.pr, function(x) hpd_mult(rl.50[[x]], force_uni = TRUE)))
# tmp.x[2,] = range(sapply(ind.tasmax, function(x) hpd_mult(rl.50[[x]], force_uni = TRUE)))
# 
# pnum = 0
# lwd = 3
# cex = 1.5
# par(mfrow = c(4, 4), mar = c(0,0,0,0), oma = c(6,10,10,6))
# for (A1 in 1:2){
#     B1 = list(ind.cal, ind.usa)[[A1]]
#     for (A2 in 1:2){
#         B2 = list(ind.62, ind.90)[[A2]]
#         for (A3 in 1:2){
#             B3 = list(ind.pr, ind.tasmax)[[A3]]
#             for (A4 in 1:2){
#                 B4 = list(ind.winter, ind.summer)[[A4]]
#                 pnum = pnum + 1
#                 plot(0, type='n', xlim = tmp.x[A3,], ylim = ylim, axes = FALSE)
#                 for (A5 in 1:4){
#                     B5 = list(ind.control, ind.decadal, ind.historical, ind.obs)[[A5]]
#                     rr = Reduce(intersect, list(B1, B2, B3, B4, B5))
# 
#                     if (length(rr) == 0) next
# 
#                     if (A5 != 4){
#                         p = rl.20[[rr]]
#                     } else {
#                         p = rl.20[[rr]]
#                         }
#                     z = hpd_mult(p, force_uni = TRUE)
#                     lines(z, rep(A5-0.25, 2), lwd = lwd, col = cols[A5])
#                     points(mean(p), A5-0.25, col = cols[A5], pch = 16, cex = cex)
# 
#                     if (A5 != 4){
#                         p = rl.50[[rr]]
#                     } else {
#                         p = rl.50[[rr]]
#                         }
#                     z = hpd_mult(p, force_uni = TRUE)
#                     lines(z, rep(A5, 2), lwd = lwd, col = cols[A5])
#                     points(mean(p), A5, col = cols[A5], pch = 16, cex = cex)
# 
#                     }
#                 if ((pnum <= 4) && (pnum %% 2 == 1))
#                     axis(3, lwd = -1, lwd.ticks = 1)
#                 if ((pnum >= 13) && (pnum %% 2 == 0))
#                     axis(1, lwd = -1, lwd.ticks = 1)
#                 if (pnum %% 4 == 0){
#                     axis(side = 4, labels = paste(lab.shortmod, "RL-50"), at = 1:4, las = 1,
#                         lwd = -1, lwd.ticks = 1)
#                     axis(side = 4, labels = paste(lab.shortmod, "RL-20"), at = 1:4-0.25, las = 1,
#                         lwd = -1, lwd.ticks = 1)
#                     }
#                 box()
#                 }
#             }
#         }
#     }
# cex = 1.25
# mtext(lab.var, side = 3, outer = TRUE, at = c(0.25, 0.75),
#     line = 7, cex = cex)
# mtext(lab.sea, side = 3, outer = TRUE, at = c(0.125, 0.375, 0.625, 0.875),
#     line = 4, cex = cex)
# mtext(lab.reg, side = 2, outer = TRUE, at = rev(c(0.25, 0.75)),
#     line = 9, las = 1, cex = cex, adj = 0)
# mtext(lab.year, side = 2, outer = TRUE, at = rev(c(0.125, 0.375, 0.625, 0.875)),
#     line = 5, las = 1, cex = cex, adj = 0)
# mtext("20-50", outer = TRUE, side = 3, at = -0.120, line = 4, cex = 2.5, adj = 0)



### Extremal index
pdf(paste0(OUT_DIR, "theta.pdf"), height = height, width = width)

ylim = c(0, 5)
xlim = c(0, 1)
pnum = 0
lwd = 3
cex = 1.5
par(mfrow = c(4, 4), mar = c(0,0,0,0), oma = c(6,10,10,6))
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
                for (A5 in 1:4){
                    B5 = list(ind.control, ind.decadal, ind.historical, ind.obs)[[A5]]
                    rr = Reduce(intersect, list(B1, B2, B3, B4, B5))

                    lines(tab.theta.qq[rr,], rep(A5, 2), lwd = lwd, col = cols[A5])
                    points(tab.theta.mm[rr], A5, col = cols[A5], pch = 16, cex = cex)
#                   text(tab.theta.mm[rr], A5-0.5, label = tab.T_C[rr])
                    }

                if ((pnum <= 4) && (pnum %% 2 == 1))
                    axis(3, lwd = -1, lwd.ticks = 1)
                if ((pnum >= 13) && (pnum %% 2 == 0))
                    axis(1, lwd = -1, lwd.ticks = 1)
                if (pnum %% 4 == 0)
                    axis(side = 4, labels = lab.shortmod, at = 1:4, las = 1,
                        lwd = -1, lwd.ticks = 1)
                box()
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
mtext(expression(theta), outer = TRUE, side = 3, at = -0.120, line = 4, cex = 2.5, adj = 0)
dev.off()
                    



if (!exists("bh.list.tail"))
    bh.list.tail = rep(list(matrix(0, R+2, 3)), 16)

tails = rep(list(NULL), length(files))
pnum = 0
for (A1 in 1:2){
    B1 = list(ind.cal, ind.usa)[[A1]]
    for (A2 in 1:2){
        B2 = list(ind.62, ind.90)[[A2]]
        for (A3 in 1:2){
            B3 = list(ind.pr, ind.tasmax)[[A3]]
            for (A4 in 1:2){
                B4 = list(ind.winter, ind.summer)[[A4]]
                pnum = pnum + 1
                for (A5 in 1:4){
                    B5 = list(ind.control, ind.decadal, ind.historical, ind.obs)[[A5]]
                    rr = Reduce(intersect, list(B1, B2, B3, B4, B5))
                    if (length(rr) == 0) next
                    if (A5 != 4){
                        tails[[rr]] = rgpd(nrow(x[[rr]]$params), x[[rr]]$threshold[1],
                            x[[rr]]$params[,31], x[[rr]]$params[,33] / x[[rr]]$params[,34])
                        if (all(bh.list.tail[[pnum]][,A5] == 0)){
                            bh.list.tail[[pnum]][2, A5] = rr
                            for (i in 1:R){
                                tmp = rgpd(nrow(x[[rr]]$params), x[[rr]]$threshold[1],
                                    x[[rr]]$params[,i], x[[rr]]$params[,R+i])
                                bh.list.tail[[pnum]][2+i, A5] =
                                    bhatta.dist(tails[[rr]], tmp, n = 2^12)
                                }
                            }
                    } else {
                        tails[[rr]] = rgpd(nrow(x[[rr]]$params), x[[rr]]$threshold[1],
                            x[[rr]]$params[,1], x[[rr]]$params[,2])
                        if (all(bh.list.tail[[pnum]][1,] == 0)){
                            tmp.rr = bh.list.tail[[pnum]][2, 1]
                            bh.list.tail[[pnum]][1, 1] = 
                                bhatta.dist(tails[[tmp.rr]], tails[[rr]], n = 2^12)

                            tmp.rr = bh.list.tail[[pnum]][2, 2]
                            bh.list.tail[[pnum]][1, 2] =
                                bhatta.dist(tails[[tmp.rr]], tails[[rr]], n = 2^12)

                            tmp.rr = bh.list.tail[[pnum]][2, 3]
                            bh.list.tail[[pnum]][1, 3] =
                                bhatta.dist(tails[[tmp.rr]], tails[[rr]], n = 2^12)
                            }
                        }
                    }
                }
            }
        }
    }

tmp_mean = sapply(tails, mean)
tmp_hpd = sapply(tails, hpd_mult, force_uni = TRUE)





pdf(paste0(OUT_DIR, "tail.pdf"), height = height, width = width)
ylim = c(0, 5)
# xlim = range(sapply(rl.50[-ind.obs], function(y) hpd_mult(y, force_uni = TRUE)))
# xlim = range(xlim, range(sapply(rl.50[ind.obs], function(y) hpd_mult(y, force_uni = TRUE))))


tmp.x = matrix(0, 2, 2)
tmp.x[1,] = range(tmp_hpd[,ind.pr])
tmp.x[2,] = range(tmp_hpd[,ind.tasmax])

pnum = 0
lwd = 3
cex = 1.5
par(mfrow = c(4, 4), mar = c(0,0,0,0), oma = c(6,10,10,6))
for (A1 in 1:2){
    B1 = list(ind.cal, ind.usa)[[A1]]
    for (A2 in 1:2){
        B2 = list(ind.62, ind.90)[[A2]]
        for (A3 in 1:2){
            B3 = list(ind.pr, ind.tasmax)[[A3]]
            for (A4 in 1:2){
                B4 = list(ind.winter, ind.summer)[[A4]]
                pnum = pnum + 1
                plot(0, type='n', xlim = tmp.x[A3,], ylim = ylim, axes = FALSE)
                for (A5 in 1:4){
                    B5 = list(ind.control, ind.decadal, ind.historical, ind.obs)[[A5]]
                    rr = Reduce(intersect, list(B1, B2, B3, B4, B5))

                    if (length(rr) == 0) next

                    lines(tmp_hpd[,rr], rep(A5, 2), lwd = lwd, col = cols[A5])
                    points(tmp_mean[rr], A5, col = cols[A5], pch = 16, cex = cex)

                    }
                tmp = 1*apply(bh.list.tail[[pnum]], 2, function(y)
                    (y[1] >= min(tail(y, R))) && (y[1] <= max(tail(y, R))))
                tmp.mean = tmp_mean[bh.list.tail[[pnum]][2,]]
                text(tmp.mean, (1:3)+0.35, label = tmp)

                if ((pnum <= 4) && (pnum %% 2 == 1))
                    axis(3, lwd = -1, lwd.ticks = 1)
                if ((pnum >= 13) && (pnum %% 2 == 0))
                    axis(1, lwd = -1, lwd.ticks = 1)
                if (pnum %% 4 == 0){
                    axis(side = 4, labels = lab.shortmod, at = 1:4-0.25, las = 1,
                        lwd = -1, lwd.ticks = 1)
                    }
                box()
                if ((A1 == 1) && (A2 == 2))
                    axis(1, labels = FALSE, lwd.ticks = 0, lwd = 5,
                        at = tmp.x[A3,] + c(-1, 1)*diff(tmp.x[A3,])/5)
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
mtext("Tail", outer = TRUE, side = 3, at = -0.120, line = 4, cex = 2.5, adj = 0)
dev.off()










### Table for extremal indexes, thresholds, T_C
# tab.out = cbind(tab.theta.mm, tab.theta.qq, tab.threshold, tab.T_C)
# rownames(tab.out) = NULL
# xtable(tab.out, label = "",
#     align = c("l", "r", "r", "r", "r", "r"))


# ### Checking trace plots
# for (j in 5:length(x)){
#     cat(j, "\r")
#     if (j <= 48){
#         par(mfrow = c(6, 6), mar = c(0,0,0,0), oma = c(0,0,0,0))
#         for (i in 1:36)
#             plot(x[[j]]$params[,i])
#     } else {
#         par(mfrow = c(3, 1), mar = c(0,0,0,0), oma = c(0,0,0,0))
#         for (i in 1:3)
#             plot(x[[j]]$params[,i])
#         }
#     readline()
#     }
# 
# 
# 
# # File size by acceptance rate
# file.info(paste0(DIR, files[1]))$size / 1024 / 1024
# mean(x[[1]]$accept)
# 
# tx = sapply(files[-ind.obs], function(y) file.info(paste0(DIR, y))$size / 1024 / 1024)
# ty = sapply(x[-ind.obs], function(y) y$accept)
# 
# tx = matrix(rep(tx, each = nrow(ty)), nrow = nrow(ty))
# 
# col = double(48)+1
# bad = which(tx[1,] > 11)
# 
# col[which(tx[1,] >10)] = 2
# 
# 
# par(mfrow = c(1, 1))
# plot(0, type = 'n', xlim = range(ty), ylim = range(tx))
# for (i in 1:nrow(ty))
#     points(ty[i,], tx[i,], col = col, pch = i)
# 
# 
# range(sapply(x[ind.obs], function(y) y$accept))
# 
# plot(x[[41]]$params[,31], type = 'l')
# 
# plot_hpd(x[[41]]$params[,5])
# plot(x[[41]]$params[,3], type = 'l')
# 
# plot(log(x[[2]]$params[,32]), type = 'l')
# 
# par(mfrow = c(2, 2), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0,0,0,0))
# for (i in (1:length(x)[-ind.obs])){
#     plot_hpd(x[[i]]$params[,1])
#     plot_hpd(x[[i]]$params[,11])
#     plot(log(x[[i]]$params[,32]), type = 'l')
#     plot_hpd(x[[i]]$params[,31], col1 = 'dodgerblue')
#     title(main = c(i, tx[1,i]), sub = paste0(range(ty[,i]), collapse = " "), outer = TRUE, line = -2)
#     readline()
#     }
# 
