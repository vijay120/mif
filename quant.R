library(stockPortfolio)
library(quadprog)
library(corpcor)
library(ggplot2)

stocks <- c('AAPL',
'ARMH',
'BND',
'CHL',
'COST',
'DBC',
'EEMS',
'GLD',
'IEFA',
'IEMG',
'IWC',
'IWR',
'MO',
'SCZ',
'TIP',
'VAL',
'VGIT',
'VNQ',
'VNQI',
'VONE',
'VTWO')

currentWeights <- c(0.080796759,
0.072231243,
0.052856677,
0.040520345,
0.130094855,
0.011853133,
0.009644557,
0.011834489,
0.050605374,
0.031440556,
0.011838684,
0.060550259,
0.057999404,
0.034692438,
0.094000985,
0.025449511,
0.048586814,
0.011876283,
0.01252464,
0.093530837,
0.05707216)


eff.frontier <- function (returns, short="no", max.allocation=NULL, risk.premium.up=.5, risk.increment=.005){
	
	# return argument should be a m x n matrix with one column per security
	# short argument is whether short-selling is allowed; default is no (short selling prohibited)
	# max.allocation is the maximum % allowed for any one security (reduces concentration)
	# risk.premium.up is the upper limit of the risk premium modeled (see for loop below)
	# risk.increment is the increment (by) value used in the for loop
	
	#covariance matrix was not positive definite
	covariance <- cov(returns)	
	covariance <- make.positive.definite(covariance)	

	n <- ncol(covariance)
	
	# Create initial Amat and bvec assuming only equality constraint (short-selling is allowed, no allocation constraints)
	Amat <- matrix (1, nrow=n)
	bvec <- 1
	meq <- 1
	
	# Then modify the Amat and bvec if short-selling is prohibited
	if(short=="no"){
		Amat <- cbind(1, diag(n))
		bvec <- c(bvec, rep(0, n))
	}
	# And modify Amat and bvec if a max allocation (concentration) is specified
	if(!is.null(max.allocation)){
		if(max.allocation > 1 | max.allocation <0){
		stop("max.allocation must be greater than 0 and less than 1")
		}
		if(max.allocation * n < 1){
		stop("Need to set max.allocation higher; not enough assets to add to 1")
		}
		Amat <- cbind(Amat, -diag(n))
		bvec <- c(bvec, rep(-max.allocation, n))
	}
	
	# Calculate the number of loops based on how high to vary the risk premium and by what increment
	loops <- risk.premium.up / risk.increment + 1
	loop <- 1
	
	# Initialize a matrix to contain allocation and statistics
	# This is not necessary, but speeds up processing and uses less memory
	eff <- matrix(nrow=loops, ncol=n+3)
	# Now I need to give the matrix column names
	colnames(eff) <- c(colnames(returns), "Std.Dev", "Exp.Return", "sharpe")
	
	# Loop through the quadratic program solver
	for (i in seq(from=0, to=risk.premium.up, by=risk.increment)){
		dvec <- colMeans(returns) * i # This moves the solution up along the efficient frontier
		sol <- solve.QP(covariance, dvec=dvec, Amat=Amat, bvec=bvec, meq=meq)
		eff[loop,"Std.Dev"] <- sqrt(sum(sol$solution *colSums((covariance * sol$solution))))
		eff[loop,"Exp.Return"] <- as.numeric(sol$solution %*% colMeans(returns))
		eff[loop,"sharpe"] <- eff[loop,"Exp.Return"] / eff[loop,"Std.Dev"]
		eff[loop,1:n] <- sol$solution
		loop <- loop+1
	}
	
	return(as.data.frame(eff))
}


returns <- getReturns(stocks, freq="month")
covariance <- cov(returns$R)
	
#print(covariance)
covariance <- make.positive.definite(covariance)

expReturn <- as.numeric(currentWeights %*% colMeans(returns$R))
print(expReturn)
std <- sqrt(sum(currentWeights * colSums((covariance * currentWeights))))
print(std)
sharpe <- expReturn/std
print(sharpe)

result <- aes(std, expReturn, sharpe)

eff <- eff.frontier(returns=returns$R, short="no", max.allocation=NULL, risk.premium.up=.5, risk.increment=.001)

eff.optimal.point <- eff[eff$sharpe==max(eff$sharpe),]
eff.optimal.point

# Color Scheme
ealred  <- "#7D110C"
ealtan  <- "#CDC4B6"
eallighttan <- "#F7F6F0"
ealdark  <- "#423C30"
ggplot(eff, aes(x=Std.Dev, y=Exp.Return)) + geom_point(alpha=.1, color=ealdark) +
 geom_point(data=eff.optimal.point, aes(x=Std.Dev, y=Exp.Return, label=sharpe), color=ealred, size=5) +
geom_point(result) + annotate(geom="text", x=std, y=expReturn, label=paste("Risk: ", round(std, digits=3),"\nReturn: ",
 round(expReturn*100, digits=4),"%\nSharpe: ",
 round(sharpe*100, digits=2), "%", sep=""), hjust=0, vjust=1.2) +
 annotate(geom="text", x=eff.optimal.point$Std.Dev, y=eff.optimal.point$Exp.Return,
 label=paste("Risk: ", round(eff.optimal.point$Std.Dev, digits=3),"\nReturn: ",
 round(eff.optimal.point$Exp.Return*100, digits=4),"%\nSharpe: ",
 round(eff.optimal.point$sharpe*100, digits=2), "%", sep=""), hjust=0, vjust=1.2) +
 ggtitle("Efficient Frontier\nand Optimal Portfolio") + labs(x="Risk (standard deviation of portfolio variance)", y="Return") +
 theme(panel.background=element_rect(fill=eallighttan), text=element_text(color=ealdark),
 plot.title=element_text(size=24, color=ealred))
