\name{GFMR.cv}
\alias{GFMR.cv}

\title{Tuning parameter selection using validation likelihood for GFMR.}
\description{
This routine implements K fold cross validation for group fused multinomial regression.  
}
\usage{
GFMR.cv(Y,X,lamb,sampID,H,n.cores=1,rho=10^-8,...)
}
\arguments{
  \item{Y}{ A matrix of response category counts where the columns represent the categories and rows represent the
  observations.  Currently supported for n=1.}
  \item{X}{A matrix of predictor variables.  The columns represent predictors and rows represent observations.}
  \item{lamb}{tuning parameter for fusion penalty}
  \item{sampID}{An identified or the sampleID for the cross validation routine.  Should take values 1:k and be user supplied.}
  \item{H}{An indicator matrix representing the edge set of the penalty set.  The matrix is square and symmetric with dimension number of response categories, and if two categories are in the penalty set a 1 should be in the row column combination.}
  \item{n.cores}{The number of cores for the mclapply function}
  \item{rho}{Step Size parameter of ADMM}
  \item{...}{Other arguments for Group fused multinomial regression}
}
\details{
This routine implements the validation likelihood approach proposed by Price et. al. to obtain the tuning parameter for Group Fused Multinomial Regression
which automatically combines response categories in multinomial regression. Show to do well with regard to predicting the true category probabilities.  
}
\value{
A list is returned with elements 
\item{vl}{The validation likelihood for all tuning parameters}
\item{vl.sd}{The standard deviation of the validation likelihood}
\item{lambda}{tuning parameter used}
\item{vl.mat}{Validation likelihood for each}
}
\references{
Price, B.S, Geyer, C.J. and Rothman, A.J. "Automatic Response Category Combination in Multinomial Logistic Regression."  https://arxiv.org/abs/1705.03594. 
}

\author{ Brad Price, \email{brad.price@mail.wvu.edu}.}


\examples{

\dontrun{data(nes96)
attach(nes96)
Response=matrix(0,944,7)
for(i in 1:944){
  if(PID[i]=="strRep"){Response[i,1]=1}
  if(PID[i]=="weakRep"){Response[i,2]=1}
  if(PID[i]=="indRep"){Response[i,3]=1}
  if(PID[i]=="indind"){Response[i,4]=1}
  if(PID[i]=="indDem"){Response[i,5]=1}
  if(PID[i]=="weakDem"){Response[i,6]=1}
  if(PID[i]=="strDem"){Response[i,7]=1}
}

Hmat=matrix(1,dim(Response)[2],dim(Response)[2])
diag(Hmat)=0
ModMat<-lm(popul~age,x=TRUE)$x

X=cbind(ModMat[,1],apply(ModMat[,-1],2,scale))

set.seed(1010)
n=dim(Response)[1]
sampID=rep(5,n)
samps=sample(1:n)
mine=floor(n/5)
for(j in 1:4){
  sampID[samps[((j-1)*mine+1):(j*mine)]]=j
}

o1<-GFMR.cv(Response,X,lamb = 2^seq(4.2,4.3,.1),H=Hmat2,sampID = sampID,n.cores =5)

which(o1$vl==max(o1$vl))
}
}
\keyword{regression}% at least one, from doc/KEYWORDS

