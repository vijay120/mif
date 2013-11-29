
compute.cov <- function(data, startDate, nHistDays = 252, datecol = 'date') {

  ## make sure data is sorted
  data <- data[ match(data[[datecol]], sort(data[[datecol]])),]

  ## find the valid daterange
  eIdx <- max(which( data[[datecol]] < startDate-1  ))
  sIdx <- pmax( eIdx - nHistDays, 1)

  ## subset down to the desired daterange
  data <- data[sIdx:eIdx,]

  ## grab the columns that include returns
  retCols <- names(data)[grep('.ret$', names(data))]

  ## actually compute covariance
  cov( data[,retCols] )
}


preprocess.yf.csv <- function(csv) {

  dat <- read.csv(csv)
  ## clean up the Y! Finance data to be useful R format
  dat$Date <- as.POSIXct(as.character(dat$Date))
  dat <- dat[match(dat$Date, sort(dat$Date)),]

  ## add returns on Adj.Close (includes dividend, stock split).
  ## Note that the as-of date depends on the daterange you download,
  ## so data should be internally consistent but may change if you
  ## redownload on a different day
  dat$ret <- c(NaN,dat$Adj.Close[-1] / dat$Adj.Close[ -1*nrow(dat)])-1

  dat
}


# creates a single data.frame with Date and returns columns for all the assets
combine.yf.csv <- function(csvs) {
  ## csvs - named list of csv files, name of each list item is ticker

  
  ## read all the csvs
  datas <- lapply( csvs, preprocess.yf.csv)

  ## find out what days are in each
  days <- sapply( datas, function(x){ as.character(x[['Date']]) })

  ## get the intersection of all days (i.e. only have an output row when all inputs are present
  uniqueDays <- as.character(days[[1]]);
  for( day in days[-1] ) {
    uniqueDays <- intersect( uniqueDays, as.character(day) )
  }
  uniqueDays <- as.POSIXct(uniqueDays)

  ## sort by date
  res <- data.frame( date = sort(uniqueDays) )

  ## now actually subset each per-stock data.frame down and add its column
  for( ticker in names(csvs) ) {
    dataSubset <- datas[[ticker]];

    dataSubset <- dataSubset[dataSubset[['Date']] %in% uniqueDays,]
    ## ensure still sorted
    dataSubset <- dataSubset[ match(dataSubset[['Date']], uniqueDays),]

    res[[paste(ticker,'.ret',sep = '')]] <- dataSubset[['ret']]
  }
  res
}

## wgt vec is named vector w/names = ticker and contents = portfolio weights
compute.portfolio.pnl <- function(data, sdate, edate, wgtVec, dateCol = 'date', norm = TRUE )  {

  ## force weights to sum to 1 - blows up if you have a 0 nmv portfolio
  if(norm) { wgtVec <- wgtVec/sum(wgtVec) }

  ## find the date range
  sIdx <- match( sdate, data[[dateCol]] )
  eIdx <- match( edate, data[[dateCol]] )

  ## find the returns
  retCols <- paste( names(wgtVec), '.ret', sep = '')
  returns <- data[sIdx:eIdx, retCols]

  ## compute the actual pnl per date
  pnlTs <- as.numeric( as.matrix(returns) %*% wgtVec )

  ## format the output
  names(pnlTs) <- as.character(data[[dateCol]][ sIdx:eIdx ])

  pnlTs
  
}


compute.pnl.rolling.cov <- function( data, alpVec, daysPer = 252 ) {

  ## figure out what periods we're going to use to compute the covariance
  breakDates <- data[seq(1, nrow(data), daysPer ),'date']

  pnlTs <- rep(0, nrow(data))
  ## for each date in the daterange, get the covariance, and figure out the
  ## target weights from alpha and covariance. Grab the pnl for those dates
  ## with the target weights.
  for( idx in 1:(length(breakDates)-2) ) {

    ## get the cov
    covSdate <- breakDates[idx]
    covEdate <- breakDates[idx+1]
    covMat <- compute.cov(data, covEdate, nHist = daysPer  )

    ## figure out weights
    portfolio.wgts <- as.vector(alpVec %*% solve(covMat))
    print(portfolio.wgts)
    names(portfolio.wgts) <- names(alpVec)

    ## populate the output for the appropriate date range
    outIdx <-  seq(match(breakDates[idx+1], data[['date']]),match(breakDates[idx+2], data[['date']]) ,1)
    pnlTs[outIdx] <- compute.portfolio.pnl( data, breakDates[idx+1], breakDates[idx+2], portfolio.wgts)
  }
  idx <- length(breakDates)-1
  outIdx <- seq(match(breakDates[idx+1], data[['date']]), nrow(data),1)
  pnlTs[outIdx] <- compute.portfolio.pnl( data, breakDates[idx+1], max(data[['date']]), portfolio.wgts)

  ## omit the first n zeros that we had from above
  pnlTs[-1:-(daysPer+1)]
  
}

sharpe <- function(vec){
  r.mn(vec)/r.sd(vec)*sqrt(252)
}

gen.summary <- function() {

  ## save the first year of data so we can compute the historical covariance matrix without lookahead
  oos.sdate <- as.POSIXct("2005-11-18")

  ## load up the data
  csvs <- list( spy = '~/Downloads/spy.csv', tlt = '~/Downloads/tlt.csv', gld  = '~/Downloads/gld.csv')
  retDf <- combine.yf.csv( csvs )
  retDf <- retDf[-1,]


  data <- retDf
  data.oos <- data[ data[['date']] > oos.sdate, ]

  ########### Check performance under different weighting schemes

  ## Arbitrary "typical" asset allocation
  wgtVec <- c(0.6, 0.3, 0.1)
  names(wgtVec) <- c('spy','tlt','gld')
  pnl.typical <- compute.portfolio.pnl( data.oos, min(data.oos$date), max(data.oos$date), wgtVec = wgtVec )

  ## what if we actually knew what returns were going to be?
  ## How good can we be if we know returns and can use them to get the true realized covariance
  ## (note this is cheating in two ways)
  alpVec.lookahead <- c(r.mn(data.oos[['spy.ret']]), r.mn(data.oos[['tlt.ret']]), r.mn(data.oos[['gld.ret']]))
  names(alpVec.lookahead) <- c('spy','tlt','gld')
  cov <- compute.cov( data.oos, max(data.oos$date), 5000000 )
  lookahead.cov.wgt <- as.vector(alpVec.lookahead %*% solve(cov))
  names(lookahead.cov.wgt) <- c('spy','tlt','gld')
  pnl.lookaheadCovAlp <- compute.portfolio.pnl( data.oos, min(data.oos$date), max(data.oos$date), wgtVec = lookahead.cov.wgt )

  ## what if we just use weights proportional to what returns are going to be (i.e. ignore risk)
  pnl.lookaheadAlp <- compute.portfolio.pnl( data.oos, min(data.oos$date), max(data.oos$date), wgtVec = alpVec.lookahead )

  ## Russ' made up alpha going forward
  alpVec <- c(0.08, 0.03, 0.02)
  names(alpVec) <- c('spy','tlt','gld')
  lookahead.cov.wgt <- as.vector(alpVec %*% solve(cov))
  names(lookahead.cov.wgt) <- c('spy','tlt','gld')
  pnl.lookaheadCovOnly <- compute.portfolio.pnl( data.oos, min(data.oos$date), max(data.oos$date), wgtVec = lookahead.cov.wgt )

  ## What if we use the true future returns from above but use the realized covariance as our risk model?
  rolling.cov.lookahead <- compute.pnl.rolling.cov( data, alpVec.lookahead )

  ## No look ahead at all, use a priori alphas and rolling covariance
  ## this one is actually not cheating
  rolling.cov <- compute.pnl.rolling.cov( data, alpVec )

  ## single asset performance
  wgt.spy <- c(1)
  names(wgt.spy) <- 'spy'
  pnl.spy <- compute.portfolio.pnl( data.oos, min(data.oos$date), max(data.oos$date), wgtVec = wgt.spy )

  wgt.gld <- c(1)
  names(wgt.gld) <- 'gld'
  pnl.gld <- compute.portfolio.pnl( data.oos, min(data.oos$date), max(data.oos$date), wgtVec = wgt.gld )
  
  wgt.tlt <- c(1)
  names(wgt.tlt) <- 'tlt'
  pnl.tlt <- compute.portfolio.pnl( data.oos, min(data.oos$date), max(data.oos$date), wgtVec = wgt.tlt )

  pnl.tses <- list( spy=pnl.spy, gld=pnl.gld, tlt=pnl.tlt, typical=pnl.typical, risk=rolling.cov.lookahead, alp=pnl.lookaheadAlp)
  

  ## set up plots
  par(mfrow = c(1,2))
  ys <- lapply( pnl.tses, cumsum)
  ylim <- 1.05*c( min(sapply(ys, min)), max(sapply(ys, max)) )

  plot(y=cumsum(pnl.spy), x=data.oos$date, type ='l', main = 'PnL for Various Portfolios', xlab = '', ylab = '', ylim = ylim)
  points(y=cumsum(pnl.gld), x=data.oos$date, type ='l',col=2);
  points(y=cumsum(pnl.tlt), x=data.oos$date, type ='l',col=3)
  points(y=cumsum(pnl.lookaheadAlp), x=data.oos$date, type ='l',col=4);
  points(y=cumsum(rolling.cov.lookahead ), x=data.oos$date, type ='l',col=5);
  points(y=cumsum(pnl.typical ), x=data.oos$date, type ='l',col=6);

  legend('topleft', c('spy','gld','tlt','alp','risk','typical'), text.col = 1:6)

  ## make risk-adjusted plots
  scaling <- 1/sapply( pnl.tses, r.sd)
  ys <- lapply( seq_along(pnl.tses), function(idx){ cumsum(pnl.tses[[idx]]*scaling[idx])})
  ylim <- 1.05*c( min(sapply(ys, min)), max(sapply(ys, max)) )

  plot(y=cumsum(pnl.spy)*scaling[['spy']], x=data.oos$date, type ='l', main = 'Risk-Adjusted PnL for Various Portfolios', xlab = '', ylab = '', ylim = ylim)
  points(y=cumsum(pnl.gld)*scaling[['gld']], x=data.oos$date, type ='l',col=2);
  points(y=cumsum(pnl.tlt)*scaling[['tlt']], x=data.oos$date, type ='l',col=3)
  points(y=cumsum(pnl.lookaheadAlp)*scaling[['alp']], x=data.oos$date, type ='l',col=4);
  points(y=cumsum(rolling.cov.lookahead )*scaling[['risk']], x=data.oos$date, type ='l',col=5);
  points(y=cumsum(pnl.typical )*scaling[['typical']], x=data.oos$date, type ='l',col=6);
  legend('topleft', c('spy','gld','tlt','alp','risk','typical'), text.col = 1:6)
  
}
