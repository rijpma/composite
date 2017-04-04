mar = c(4, 3, 3.5, 1)
red = '#F52A01'# red2 = '#C10033'
pin = '#9A0000'# pin2 = '#A10065'
yel = '#FFCD00'# yel2 = '#F5D300'
ora = '#F08000'# ora2 = '#DC931A'
blu = '#094D8E'# blu2 = '#222D80'
lbl = '#36A2C9'# lbl2 = '#5DB4E5'
pur = '#791D72'# pur2 = '#5F247D'
lgr = '#FEF200'# lgr2 = '#CCCE1E'
gre = '#419702'# gre2 = '#118F40'

clioctrs = c('GBR', 'NLD', 'FRA', 'DEU', 'DDR', 'ITA', 'ESP', 'SWE',
             'POL', 'RUS', 'SUN',               "AUS", "CAN", "USA", 
             "MEX", "BRA", "ARG",                      'EGY', 'TUR', 
             'KEN', 'NGA', 'ZAF',                      'CHN', 'JPN',
             'IND', 'IDN', 'THA')
cliocodes = countrycode::countrycode(clioctrs, 'iso3c', 'iso3n')
cliocodes[clioctrs=='SUN'] = 810
cliocodes[clioctrs=='DDR'] = 278
cliocodes = c(cliocodes, 280)

ww_averages = function(dat, vrb, 
  pop='pop', region='region', decade='decade', threshold=0.4){

    allpop = aggregate(list(cov=dat[, pop]), 
        by=list(region=dat[, region], decade=dat[, decade]), sum, na.rm=T)
    allpop_wrld = aggregate(list(cov=dat[, pop]), 
        by=list(decade=dat[, decade]), sum, na.rm=T)

    dat = dat[!is.na(dat[, vrb]) & !is.na(dat[, pop]), ]
    dat$vrbxpop = dat[, vrb] * dat[, pop]

    reg_vrb = aggregate(list(vrbxpop=dat[, 'vrbxpop'], pop=dat[, pop]), 
        by=list(region=dat[, region], decade=dat[, decade]), sum, na.rm=T)
    
    wld_vrb = aggregate(list(vrbxpop=dat[, 'vrbxpop'], pop=dat[, pop]), 
        by=list(decade=dat[, decade]), sum, na.rm=T)
    allpop_wrld$region = wld_vrb$region = 'World'

    out = rbind(reg_vrb, wld_vrb)
    out$wm = out$vrbxpop / out$pop

    allpop = rbind(allpop, allpop_wrld)
    out = merge(out, allpop, by=c(region, decade), all=TRUE)

    out$cov = out$pop / out$cov
    out$wm[out$cov < threshold] = NA

    return(out)
}

make_inits = function(y, vrs){
    # initialising values for factor loadings (gamma)
    # variance of factor scores (omega)
    # and factors scores (xi)
    # multiple combinations of variables for gamma because 
    # factanal only works on complete cases

    y$init = NA
    for (i in 3:length(vrs)){ # 3=min for fantanal
        combinations = combn(vrs, i)
        for (j in 1:ncol(combinations)){
            variables = combinations[, j]
            dfa = y[complete.cases(y[, variables]), variables]
            try(fa <- factanal(dfa, 1, scores='Bartlet'))
            y$init[match(rownames(fa$scores), rownames(y))] = fa$scores
        }
  }

  initxi = y$init
  initgamma = t(apply(y[, variables], 2, function(x) lm(x ~ y$init)$coef))
  initomega = apply(y[, variables], 2, function(x) summary(lm(x ~ y$init))$sigma)
  
  return(list(xi=initxi, gamma=initgamma, omega=initomega))
}

extract = function(coda, varname){
    keep = grep(varname, colnames(coda[[1]]))
    out = coda[, keep, drop=FALSE]
    return(as.matrix(out))
}

sumstats = function(x){
  out = c(mean(x), quantile(x, c(0.05, 0.5, 0.95)))
  names(out) = c('mean', 'q05', 'q50', 'q95')
  return(out)
}

sumstatsDF = function(coda){
  out = apply(coda, 2, sumstats)
  out = t(out)
  out = as.data.frame(out)
  return(out)
}

vstrsplit = function(strs, split){
  spltlist = strsplit(strs, split)
  spltvec  = do.call(rbind, spltlist)
  return(spltvec)
}

plot_bg = function(dat, yvrb, xvrb, byvrb, col="gray"){
    for (i in unique(dat[, byvrb])){
        dat_sub = dat[dat[, byvrb] == i, c(yvrb, xvrb)]
        dat_sub = dat_sub[!is.na(dat_sub[, yvrb]), ]

        lines(x = dat_sub[, xvrb], y = dat_sub[, yvrb], col=col)
    }
}

place_text <- function(dat, countries, col=1, cex=0.6){
    N <- nrow(dat)

    y_left <- seq(from=1, to=N, by=2)
    x_left <- dat$q05[y_left] 
    text_left <- countries[y_left]
    text(x_left, y_left, text_left, cex=cex, pos=2, col=col)

    y_right <- seq(from=2, to=N, by=2)
    x_right <- dat$q95[y_right]
    text_right <- countries[y_right]
    text(x_right, y_right, text_right, cex=cex, pos=4, col=col)
}