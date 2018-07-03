rm(list=ls())
options(stringsAsFactors=FALSE, digits=3)
setwd("~/dropbox/composite/")

source("bgdp_functions.R")

library("xtable")
library("texreg")
library("rjags")
library("stringi")
library("countrycode")
library("segmented")
library("plm")

### data ###
############
# clio data with missing obs
da = read.csv("allcliodata.csv")
da = da[da$decade < 2010, ]

# clio data with log-linear imputations
di = read.csv('allcliodata_imputed50ymax.csv')
di = di[di$decade < 2010, ]

pop = read.csv("regionpop.csv")
# nb: pop and pop2 in da not identical

# wvs life satisfaction data
sat = read.csv('wvslifesat.csv')
sat$decade = trunc(sat$year / 10) * 10
da = merge(da, sat[, c('sat', 'decade', 'ccode')], by=c('ccode', 'decade'), all.x=T, all.y=F)

# regional gdp figures
gdp = gdp_rgnl = di
gdp = ww_averages(dat=gdp, vrb='gdp', pop='pop', region='region', decade='decade')
gdp$gdp = scale(gdp$wm)

gdp_rgnl = aggregate(gdp ~ region + decade, data=gdp_rgnl, mean)
gdp_rgnl$gdp = scale(gdp_rgnl$gdp)

### consts and vrbs ###
#######################
vrbs = c('lab', 'hgt', 'lif', 'edu', 'ine', 'pol', 'hom')
grps = c( 5,     3,     2,     1,     6,     4,    7)
grpr = c( 1e-5,  100,   100,   1e-5,  1e-5,  1e-5, 1e-5)
ngrp = max(grps)
round(cor(da[, vrbs], use='pairwise.complete'), 2)

lngvrbs = c('Lab. real wage', 'Height', 'Life exp.', 'Av. years edu.', 'Inequality', 'Polity', 'Homicide rate')
names(lngvrbs) = vrbs

vrbs_gdp = c(vrbs, 'gdp')
lngvrbs_gdp = c(lngvrbs, 'GDPpc')
names(lngvrbs_gdp) = vrbs_gdp

### sumstats ###
################
N    = sapply(da[, vrbs_gdp], function(x) sum(!is.na(x)))
mns  = sapply(da[, vrbs_gdp], mean, na.rm=T)
sds  = sapply(da[, vrbs_gdp], sd, na.rm=T)
mins = sapply(da[, vrbs_gdp], min, na.rm=T)
maxs = sapply(da[, vrbs_gdp], max, na.rm=T)
miss = sapply(da[, vrbs_gdp], function(x) sum(is.na(x)))
prds = Ncml = NULL
for (vrb in 1:length(vrbs_gdp)){
    Ncml[vrb] = sum(complete.cases(da[, vrbs_gdp[1:vrb]]))
    prds[vrb] = paste(range(da$decade[!is.na(da[, vrbs_gdp[vrb]])]), collapse='-')
}

out = data.frame(mns, sds, mins, maxs, N, row.names=lngvrbs_gdp[vrbs_gdp])
out[, -ncol(out)] = round(out[, -ncol(out)], 1)

# print(file='sumstats.tex', placement='h',
#     xtable(out, digits=1, label='tab:sumstats', 
#     caption='Summary statistics of wellbeing indicators, 1820–2000'))

### global trends ###
#####################
pdf('worldaverages.pdf', width=8, height=4)
par(mfrow=c(2, 4), mar=mar, font.main=1)
for (vrb in vrbs_gdp){
    dat = ww_averages(da, vrb=vrb)
  if (vrb=='hgt') dat = dat[dat$decade <= 1980, ]
  dat = dat[dat$region=='World', ]
  plot(dat[, 'decade'], dat[, 'wm'], 
    ylab=lngvrbs[vrb], xlab='', type='b', col=red, bty='l')
  title(main=lngvrbs[vrb], line=-0.9)
}
dev.off()

### correlations ###
####################
dcds = unique(da$decade)
cormat = matrix(NA, nrow=length(dcds), ncol=length(vrbs))
dimnames(cormat) = list(dcds, vrbs)
nmat = lowmat = hghmat = cormat

oiloutl = c(634, 414, 784)
dc = da[!da$ccode %in% oiloutl, ]
dc = dc[!is.na(dc$gdp), ]

longcomps = list()
for (vrb in vrbs){
    longcomps[[vrb]] = dc$ccode[dc$decade < 1850 & !is.na(dc$gdp) & !is.na(dc[, vrb])]
}

for (dcd in as.character(dcds)){
    for (vrb in vrbs){
        ds = dc[dc$decade == dcd, c("gdp", vrb)]
        nobs = sum(complete.cases(ds))
        nmat[dcd, vrb] = nobs
        if (nobs > 2){
            co = cor.test(ds[, 1], ds[, 2], use="pairwise.complete")
            cormat[dcd, vrb] = co$estimate
            if (nobs > 3){
                lowmat[dcd, vrb] = co$conf.int[1]
                hghmat[dcd, vrb] = co$conf.int[2]
            }
        }
    }
}

pdf('correlationpanel.pdf', width=8, height=4)
par(mfrow=c(2, 4), mar=mar, font.main=1)
for (vrb in vrbs){
    corm = cormat[!is.na(cormat[, vrb]), ]
    lowm = lowmat[!is.na(lowmat[, vrb]), ]
    hghm = hghmat[!is.na(hghmat[, vrb]), ]
    plot(rownames(corm), corm[, vrb], type='l', col=red,
        lwd=2, ylim=c(-1, 1), bty='l', 
        ylab='corr. coef.', xlab='decade')
    title(main=lngvrbs[vrb], line=-0.8)
    lines(rownames(lowm), lowm[, vrb], col=pin)
    lines(rownames(hghm), hghm[, vrb], col=pin)
}
dev.off()

### lifesat regs ###
####################
dcs = dc
dcs$loggdp = log(dcs$gdp)
dcs[, vrbs_gdp] = lapply(dcs[, vrbs_gdp], scale)

f = paste('sat ~ ' , paste(vrbs_gdp, collapse=' + '), '+ factor(decade) - 1')
msat = lm(f, data=dcs)
summary(msat)

mlist = list()
for (vrb in vrbs){
    f = paste('sat ~ loggdp + factor(decade) - 1 +',  vrb)
    mlist[[vrb]] = lm(f, data=dcs)
}

# texreg(mlist, stars = c(0.01, 0.05, 0.1), 
#   caption='Regression of country-average life satisfaction responses in the WVS on wellbeing indicators',
#   label='tab:satreg',
#   float.pos='h',
#   custom.coef.names=c('Log(GDPpc)', '1981-4', '1990-3', '1999-2001', lngvrbs[vrbs]))
#   file='satreg.tex', fontsize='footnotesize')

### latent variable model ###
#############################
y = da
y = y[order(y$region, y$decade, y$ccode), ]
y[vrbs] = scale(y[vrbs])
lvls = as.factor(paste(y$region, y$decade, sep='_'))

dl = list(y = y[vrbs], N = nrow(y), K = length(vrbs),
    level = lvls, M = length(unique(lvls)),           # levels in multil. model
    g0 = c(0, 0), G0 = diag(1e-7, 2),                 # means and vcov for mnorm
    ngrp = ngrp, grpr = grpr, grps = grps)            # intergroup cor. allowance

inits = make_inits(y, vrbs)

il = list(beta = inits$gamma, tau = 1/ (inits$omega^2), x = inits$xi,
    .RNG.name = 'base::Mersenne-Twister', .RNG.seed = 1438996)

# m = jags.model("fm_hier.bug", data = dl, inits = il)
# update(m, 4e3)
# postr = coda.samples(m, n.iter = 1e4, thin = 10,
#   variable.names = c('country', 'load', 'intercept', 'region', 'y'))
# takes about 10-15m to run

# check convergence by running multiple chains
# 5e4 for final run seems sufficient
il = list(list(beta = inits$gamma, tau = 1/ (inits$omega^2), x = inits$xi,
    .RNG.name = 'base::Mersenne-Twister', .RNG.seed=473826),
          list(beta = inits$gamma, tau = 1/ (inits$omega^2), x = inits$xi,
    .RNG.name = 'base::Mersenne-Twister', .RNG.seed=195749),
          list(beta = inits$gamma, tau = 1/ (inits$omega^2), x = inits$xi,
    .RNG.name = 'base::Mersenne-Twister', .RNG.seed=404312))
m = jags.model("fm_hier.bug", data = dl, inits = il, n.chains=3)
update(m, 4e3)
postr = coda.samples(m, n.iter = 1e4, thin = 10,
  variable.names = c('country', 'load', 'intercept', 'region', 'y'))
# # 
gelman.diag(postr[, grep('load', colnames(postr[[1]]))])
# gelman.diag(postr[, grep('intercept', colnames(postr[[1]]))])
# gelman.diag(postr[, grep('region', colnames(postr[[1]]))])
# # etc.
# gelman.diag(postr[, grep('country\\[1', colnames(postr[[1]]))])

# plot(postr[, grep('load', colnames(postr[[1]]))])
# plot(postr[, grep('region', colnames(postr[[1]]))])
# plot(postr[, grep('country', colnames(postr[[1]]))])
# plot(postr[, grep('intercept', colnames(postr[[1]]))])
# # etc. 

save(postr, file = "~/downloads/data/bgdp_data.rda")

loads = extract(postr, "load")
loadsbar = sumstatsDF(loads)
loadsbar$vrb = vrbs
rownames(loadsbar) = lngvrbs[vrbs]

intercepts = extract(postr, "intercept")
interceptsbar = sumstatsDF(intercepts)
interceptsbar$vrb = vrbs
rownames(interceptsbar) = lngvrbs[vrbs]

countries = extract(postr, "country")
countriesbar = sumstatsDF(countries)
countriesbar$ccode = y$ccode
countriesbar$pop = y$pop
countriesbar$region = y$region
countriesbar$decade = y$decade
countriesbar$gdp = y$gdp
countriesbar$gdp_std = scale(countriesbar$gdp)

write.csv(countriesbar, "composite_country_estimates.csv", row.names = F, na = "")

regions = extract(postr, "region")
regionsbar = sumstatsDF(regions)
regionsbar$region = stringi::stri_split_fixed(unique(lvls), '_', simplify = TRUE)[, 1]
regionsbar$decade = stringi::stri_split_fixed(unique(lvls), '_', simplify = TRUE)[, 2]
regionsbar$decade = as.numeric(regionsbar$decade)
regionsbar = regionsbar[order(regionsbar$region, regionsbar$decade), ]

y_imp = extract(postr, "^y")
y_imp = sumstatsDF(y_imp)
y_imp$vrb = rep(vrbs, each = nrow(y))
y_imp$ccode = rep(y$ccode, times = length(vrbs))
y_imp$decade = rep(y$decade, times = length(vrbs))
y_imp = reshape(y_imp, direction = "wide", timevar = "vrb", idvar = c("ccode", "decade"))
y_imp = merge(y_imp, da[, c("pop", "decade", "ccode", "region")], 
    by = c("ccode", "decade"), all.x = TRUE, all.y = FALSE)
y_imp = y_imp[order(y_imp$region, y_imp$decade, y_imp$ccode), ]
all.equal(countriesbar$ccode, y_imp$ccode)


pdf('regions_nogdp.pdf', height=6, width=10)
par(mfrow=c(2, 4), mar=mar, font.main=1)
for (region in unique(regionsbar$region)){
    toplot = regionsbar[regionsbar$region == region, ]
    plot(c(1820, 2000), c(-2, 3), type='n', bty='l', 
        xlab='decade', ylab='Comp. Ind.; GDP/c')
    title(main = region, line = -0.7, cex.main = 0.9)
    plot_bg(dat = regionsbar, yvrb = 'q50', xvrb = 'decade', byvrb = 'region')
    lines(q50 ~ decade, data = toplot, col = red, lwd = 2)
    # lines(gdp ~ decade, data = toplot, col = pin, lwd = 2)
    segments(x0 = toplot$decade, y0 = toplot$q05, y1 = toplot$q95, col = red)
}
dev.off()

pdf('regions_comparegdp.pdf', height=6, width=10)
par(mfrow=c(2, 4), mar=mar, font.main=1)
for (region in unique(regionsbar$region)){
    toplot = regionsbar[regionsbar$region == region, ]
    plot(c(1820, 2000), c(-2, 3), type='n', bty='l', 
        xlab='decade', ylab='Comp. Ind.; GDP/c')
    title(main = region, line = -0.7, cex.main = 0.9)
    plot_bg(dat = regionsbar, yvrb = 'q50', xvrb = 'decade', byvrb = 'region')
    lines(q50 ~ decade, data = toplot, col = red, lwd = 2)
    lines(gdp ~ decade, data = gdp_rgnl[gdp_rgnl$region == region, ], col = pin, lwd = 2)
    segments(x0 = toplot$decade, y0 = toplot$q05, y1 = toplot$q95, col = red)
}
legend('bottomright', c('Comp. Ind', 'Std. GDP/c'), fill=c(red, pin))
dev.off()

pdf('greatdiv.pdf', height=5)
par(mfrow=c(1, 1))
plot(c(1820, 2000), c(-3, 5), bty='l', type='n', xlab='decade', ylab='Std. pcGDP, Composite indicator')
zaf = countriesbar[countriesbar$ccode==710, ]
lines(gdp_std ~ decade, data=zaf, col=1)
arrows(zaf$decade, zaf$gdp_std, zaf$decade, zaf$q50, length=0.1, col=1)
ind = countriesbar[countriesbar$ccode==356, ]
lines(gdp_std ~ decade, data=ind, col=gre)
arrows(ind$decade, ind$gdp_std, ind$decade, ind$q50, length=0.1, col=gre)
gbr = countriesbar[countriesbar$ccode==826, ]
lines(gdp_std ~ decade, data=gbr, col=red)
arrows(gbr$decade, gbr$gdp_std, gbr$decade, gbr$q50, length=0.1, col=red)
legend('topleft', c('South Africa', 'India', 'Great Britain'), fill=c(1, gre, red))
dev.off()

xl = range(countriesbar[, grep('q', names(countriesbar))])
countriesbar$iso3 = countrycode::countrycode(countriesbar$ccode, 'iso3n', 'iso3c')
countriesbar$big = countriesbar$ccode %in% unique(da$ccode[da$pop > 1000])
countriesbar$clio = countriesbar$ccode %in% cliocodes

pdf("countrysteps.pdf")
par(mfrow=c(2, 2), mar=c(3, 0.5, 1.5, 0.5), font.main=1)
for (decade in c(1850, 1900, 1950, 2000)){
    dat = countriesbar[countriesbar$decade == decade, ]
    dat = dat[order(dat$q50), ]
    N = nrow(dat)
    dat$N = 1:nrow(dat)
    plot(xl, c(1, N), type='n', axes=F, ylab='', xlab='Cmp. Index', main=decade)
    axis(1)
    segments(x0=0, y0=1, y1=N, col='gray')
    points(dat$q50, dat$N, pch=19, cex=0.8, col=red)
    segments(x0=dat$q05, x1=dat$q95, y0=dat$N, col=pin)

    # dat_clio = dat[dat$clio, ]
    # points(dat_clio$q50, dat_clio$N, pch=19, cex=0.8, col=red)
    # segments(x0=dat_clio$q05, x1=dat_clio$q95, y0=dat_clio$N, col=pin)
    dat$iso3[!dat$clio] = ""
    place_text(dat, dat$iso3, cex=0.8)
}
dev.off()

### segmented relation ###
##########################

linmod = lm(q50 ~ gdp, data=countriesbar)
linmod2fe = lm(q50 ~ gdp + factor(decade) + factor(ccode), data=countriesbar)
segm = segmented::segmented(linmod, ~ gdp, psi=7000)
summary(segm)
slope(segm)

pdf("segmented.pdf", height=5)
par(mfrow=c(1, 2))
plot(q50 ~ gdp, data=countriesbar, col="gray", bty='l',
    xlab='GDP per capita', ylab='Composite Indicator', xlim=c(200, 35000))
plot(segm, add=T, col=red, lty=1, conf.level=0.95)
legend('topleft', legend=c('Fitted values', '95% confidence interval'),
                  col=red, lty=c(1:2))
plot(q50 ~ log(gdp), data=countriesbar, col="gray", bty='l',
  xlab='log(GDP per capita)', ylab='Composite Indicator')
abline(lm(q50 ~ log(gdp), data=countriesbar), col=red)
dev.off()

### tradeoffs ###
#################
# verify that country scores are beta0 + beta1*data
yhat = t(interceptsbar$q50 + loadsbar$q50 * t(y_imp[, paste0('q50.', vrbs)]))
yhat = scale(rowSums(yhat))
plot(yhat, countriesbar$q50)
curve(1*x, add=T, col=2)
lm(rowSums(yshat) ~ countriesbar$q50)


derivs = loadsbar$q50 / sapply(da[vrbs], sd, na.rm=T) # not on imputed data because all scaled
tomat = matrix(NA, ncol=length(vrbs), nrow=length(vrbs))
dimnames(tomat) = list(vrbs, vrbs)
for (vrb in vrbs){
    tomat[vrb, ] = derivs[vrb] / derivs
}
print(file='tradeoffs.tex', table.placement='!hbtp', size='small',
    x=xtable::xtable(tomat, digits=2, label='tab:tradeoffs',
    caption='Tradeoffs in the composite indicator'))

### technology ###
##################
countriesbar$loggdp = log(countriesbar$gdp)
countriesbar$loggdps = scale(countriesbar$loggdp)

cbp = plm::pdata.frame(countriesbar, index=c("ccode", "decade"))

mpl = plm::plm(q50 ~ loggdps, data=cbp, model='pooling')
mti = plm::plm(q50 ~ loggdps, data=cbp, effect='time')
mbo = plm::plm(q50 ~ loggdps, data=cbp, effect='twoway')

ctreff = summary(plm::fixef(mbo, effect='individual'))
ctreff = data.frame(cf = ctreff[, 1], se=ctreff[, 2])
ctreff$lw = ctreff$cf - 2*ctreff$se
ctreff$up = ctreff$cf + 2*ctreff$se
ctreff$iso3c = countrycode(rownames(ctreff), 'iso3n', 'iso3c')
ctreff['230', 'iso3c'] = 'ETH'
ctreff['890', 'iso3c'] = 'YUG'
ctreff['810', 'iso3c'] = 'SUN'
ctreff[is.na(ctreff$iso3c),]
# use twoway, want "global technology" without polution from
# country specific effects; want country characteristics (country technology)
# without distortion from global trends

clioctreff = ctreff[rownames(ctreff) %in% c(cliocodes, 810), ]
clioctreff = clioctreff[order(clioctreff$cf), ]

pdf('countryperformance.pdf')
dotchart(clioctreff$cf, labels=clioctreff$iso3, xlim=range(clioctreff[, c('up', 'lw')]))
arrows(x0=clioctreff$lw, x1=clioctreff$up, y0=1:length(clioctreff$cf), 
    angle=90, code=3, length=0.05)
dev.off()

timeff = summary(plm::fixef(mbo, effect="time"))
timeff = data.frame(cf = timeff[, 1], se=timeff[, 2])
timeff$decade = as.numeric(rownames(timeff))
timeff$lw = timeff$cf - 2*timeff$se
timeff$up = timeff$cf + 2*timeff$se

pdf('timeperformance.pdf')
plot(cf ~ decade, data=timeff, type="l", col=red, lwd=2, 
    ylim=range(timeff[, c('up', 'lw')]),
    xlab='decade', ylab='Dec. deviation')
lines(up ~ decade, data=timeff, col=pin)
lines(lw ~ decade, data=timeff, col=pin)
dev.off()

### sensitivity to weights ###
##############################
wgts = rep(1, times = length(vrbs))
names(wgts) = vrbs
wgts[c("ine", "hom")] = -1
wgts[c("hgt", "lif")] = 0.5
wgtmat = matrix(wgts, nrow=nrow(y_imp), ncol=length(vrbs), byrow=TRUE)

comblist = list()
for (i in 1:length(vrbs)){
    comblist[[i]] = combn(1:length(vrbs), i)
}

cimat = matrix(NA, ncol=sum(sapply(comblist, ncol)), nrow=nrow(y_imp))
col = 1
for (i in 1:length(comblist)){
    for (j in 1:ncol(comblist[[i]])){
        wgt_shft = wgtmat
        wgt_shft[, comblist[[i]][, j]] = wgtmat[, comblist[[i]][, j]] * 0.25
        cimat[, col] = scale(rowSums(y_imp[, paste0('q50.', vrbs)] * wgt_shft))
        col = col + 1
    }
}

cir = aggregate.data.frame(cimat, by=list(region = y_imp$region, decade = y_imp$decade), mean)
cir[, 3:ncol(cir)] <- sapply(cir[, 3:ncol(cir)], scale)

trgray <- rgb(90, 90, 90, 20, maxColorValue=255)

pdf('shifweights.pdf', height=6, width=10)
par(mfrow=c(2, 4), mar=mar, font.main=1)
for (region in unique(y_imp$region)){
    plot(x=c(1820, 2010), y=c(-3, 3.3), type='n', bty='l', xlab='decade', ylab='Comp. Ind.')
    title(main = region, line = -0.7, cex.main=0.9)
    lines(q50 ~ decade, data=regionsbar[regionsbar$region == region, ], col=1)
    for (i in 3:ncol(cir)){
        lines(cir$decade[cir$region == region], cir[cir$region == region, i], col=trgray)
    }
}
dev.off()
