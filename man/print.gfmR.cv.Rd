\name{print.gfmR.cv}
\alias{print.gfmR.cv}

\title{ print method for group fused multinomial logistic regression validation likelihood tuning parameter selection }
\description{
This routine fits the group fused multinomial logistic regression model, which uses fusion shrinkage
to automatically combine response categories.  This specifically focuses on tuning parameter selection with validation likelihood.  
}
\usage{
\method{print}{gfmR.cv}(x,...)
}
\arguments{
  \item{x}{ A gfmr.cv object which specifically is the output from the GroupFusedMulti function.}
  \item{...}{Other arguments}
}
\details{
print method for gfmR.cv objects.
}
\value{
A readable printout of cross validation 
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
o1
}

}
\keyword{regression}% at least one, from doc/KEYWORDS

