library("rstan")
options(mc.cores = parallel::detectCores())

########################################
#######De novo only
DNextTADA <- "
data {
    int<lower=1> NN; //Number of genes
int<lower=1> K; //Number of classes
int<lower=1> NCdn; //Number of de novo classes
int Ndn[NCdn]; //Number of trios

int dataDN[NN, NCdn]; //denovo data: Kdn classes
real mutRate[NN, NCdn]; //mutation rates: Kdn classes
real rankPercentile[NN]; //gene expression rank percentile
real betaPars[3]; //These parameters are used to adjust beta values as a non-linear function of gamma means

//These below parameters should be default
int<lower=0> NCdnNCdn; //number of pis
real<lower=0> upperPi0;
real<lower=0> lowerPi0;
real<lower=0> lowerHyperGamma; //Low limit for mean of relative risks
real<lower=0> lowerGamma; //Low limit for relative risks
real<lower=0> lowerBeta;
real<lower=0> hyperBetaDN0[NN, NCdn];
int<lower=0> adjustHyperBeta; //Option to adjust betas or not; if this option is used (!=0) then beta is a non-linear function of gamma means
real<lower=0> hyper2GammaMeanDN[NCdn]; //hyperGammaMeanDN ~ gamma(hyper2GammaMeanDN, hyper2BetaDN)
real<lower=0> hyper2BetaDN[NCdn]; //hyperGammaMeanDN ~ gamma(hyper2GammaMeanDN, hyper2BetaDN)

}

parameters {
simplex[NCdnNCdn] pi0; //Proportion of risk genes
real<lower=0> percentileA[NCdn];
real<lower=0> percentileB[NCdn];
real<lower=0, upper=100> percentileC[NCdn];
//real<lower=lowerHyperGamma> hyperGammaMeanDN[NCdn]; //Hyper parameters for de novo relative risks
real<lower=lowerGamma> gammaMeanDN[NN, NCdn]; //parameters (in the sampling process) for de novo relative risks

}

transformed parameters {
real hyperBetaDN[NN, NCdn];
if (adjustHyperBeta != 0) {
for (i in 1:NN) {
for (j in 1:NCdn){
hyperBetaDN[i, j] = exp(betaPars[1]*(percentileA[j]/(1+exp(-percentileB[j]*(rankPercentile[i]-percentileC[j]))) + 1)^(betaPars[2]) + betaPars[3]);
}
}
}
else {
hyperBetaDN = hyperBetaDN0;
}
}

model {
int newIndex;
real ps[NCdnNCdn];
real sumDN[2];
int index[NCdn];
int indexx;
real hyperGammaMeanDN[NN, NCdn]; //Hyper parameters for de novo relative risks
percentileA ~ gamma(0.001, 0.001); // very week prior
percentileB ~ gamma(0.001, 0.001);
percentileC ~ gamma(0.001, 0.001);
pi0 ~ dirichlet(rep_vector(2.0, NCdnNCdn)); //prior for the proportion of risk genes

for (i in 1:NN) {    
for (j in 1:NCdn) {
hyperGammaMeanDN[i, j] = percentileA[j]/(1+exp(-percentileB[j]*(rankPercentile[i]-percentileC[j]))) + 1;
}
}
//De novo data: sample for hyper priors (NPdn populations and Kdn categories)
for (im in 1:NN) {
for (ip in 1:NCdn){
gammaMeanDN[im, ip] ~ gamma(hyperGammaMeanDN[im, ip]*hyperBetaDN[im, ip], hyperBetaDN[im, ip]);
}
}

////Main program
//Loop through data points
//// 
for (ii in 1:NN){
for (jj in 1:NCdnNCdn) {
ps[jj] = log(pi0[jj]);
}
//For de novo data
for (jj in 1:NCdnNCdn){
indexx = jj;
for (kk in 1:NCdn) {
if (indexx % 2 == 0) {
ps[jj] = ps[jj] + poisson_lpmf(dataDN[ii, kk] | Ndn[kk]*2*mutRate[ii, kk]); 
} else {
ps[jj] = ps[jj] + poisson_lpmf(dataDN[ii, kk] | Ndn[kk]*2*mutRate[ii, kk]*gammaMeanDN[ii, kk]);
}
indexx = (indexx - indexx % 2) / 2;
}
}
target += log_sum_exp(ps);
//increment_log_prob(log_sum_exp(ps));
}
}

"

################


library(locfit)
loc2plot <- function(x,y,cprob=0.5, xlim,gxlim,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,scale=c(sc1,sc2), maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	else fit <- locfit(~x+y,scale=c(sc1,sc2), xlim=xlim,maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	lev <- sort(fitted(fit))[floor(cprob*length(x))]
	plot(fit,lev=lev,m=100,label=paste(as.character(100*(1-cprob)),"%",sep=""),
	xlim=gxlim,...)
}

# finds the mode for a bivariate density
loc2mode <- function(x,y,alpha=0.5,xlim,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2), maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),xlim=xlim,maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	tt <- max(fitted(fit))
	wt <- fitted(fit) == tt
	c(x[wt][1],y[wt][1])
}

# this works for univariate data; gives a vector with
# mode (global), hpd1_low, hpd1_high, hpd2_low, hpd2_hi, etc
#The reason for multiple hpd limits is if you have multimodal data.
#prob = prob *outside* the limit; i.e for a normal distribution 0.05 is expected to give
#     the 0.025 and 0.975 quantiles.
# this won't work for weighted data, use loc1statsx instead.
# xlim is optional - use it to define the support of the density.
loc1stats <- function(x,xlim,prob=0.05,...)
{
	if(missing(xlim)){
		fit <- locfit(~x)
	}
	else {
		fit <- locfit(~x,xlim=xlim)
	}
	fx <- fitted(fit)
	x.modef <- max(fx)
	x.mode <- x[fx == x.modef]
	if(!missing(xlim)){
		if(predict(fit,xlim[1]) > x.modef){
			x.modef <- predict(fit,xlim[1])
			x.mode <- xlim[1]
		}
		if(predict(fit,xlim[2]) > x.modef){
			x.modef <- predict(fit,xlim[2])
			x.mode <- xlim[2]
		}
	}

	if(length(x.mode)>1)x.mode <- x.mode[1]
	lev <- sort(fx)[floor(prob*length(x))]
#	print("log difference from max is ")
#	print(log(x.modef)-log(lev))
	l1 <- list()
	l1[[1]] <- x.mode
	indx <- order(x)
	ii <- 2
	flip <- TRUE
	for(j in 2:length(x)){
		if(flip && fx[indx[j]] > lev){
			l1[[ii]] <- x[indx[j-1]]
			if(j==2 && !missing(xlim)){
				if(predict(fit,xlim[1]) >= lev)l1[[ii]] <- xlim[1]
			}
			flip <- FALSE
			ii <- ii+1
		}
		else if(!flip && fx[indx[j]] < lev){
			l1[[ii]] <- x[indx[j]]
			flip <- TRUE
			ii <- ii+1
		}
		if(!flip && !missing(xlim) && j == length(x)){
			l1[[ii]] <- xlim[2]
			flip <- TRUE
		}
	}
	if(!flip)stop("HPD interval not closed")
	as.numeric(l1)
}








loc2plot.old <- function(x,y,cprob=0.5,alpha=0.5,xlim,gxlim,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2), maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2), xlim=xlim,maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	lev <- sort(fitted(fit))[floor(cprob*length(x))]
	plot(fit,lev=lev,m=100,label=paste(as.character(100*(1-cprob)),"%",sep=""),
	xlim=gxlim,...)
}


#gives the HPD value for an observation px,py in a density constructed from x, y.
gethpdprob2 <- function(x,y,px,py,alpha=0.5,xlim,gxlim,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2), maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2), xlim=xlim,maxk=maxk,mint=100,maxit=100)
#,cut=0.8
#	d1 <- (x-px)^2+(y-py)^2
#	best <- d1 == min(d1)
#	lev <- mean(fitted(fit)[best])
	lev <- predict.locfit(fit,list(px,py))
	slev <- sort(fitted(fit))
	indic <- slev <= lev
	sum(indic)/length(x)
}

