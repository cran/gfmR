\name{GroupFusedMulti}
\alias{GroupFusedMulti}

\title{ Group Fused Multinomial Logistic Regression  }
\description{
This routine fits the group fused multinomial logistic regression model, which uses fusion shrinkage
to automatically combine response categories.  
}
\usage{
GroupFusedMulti(Y,X,lambda,H,tol1=10^-7,tol2=10^-7,TD=2,rho=10^-8,tau1=10^-9,iter=1e3)
}
\arguments{
  \item{Y}{ A matrix of response category counts where the columns represent the categories and rows represent the
  observations.  Currently supported for n=1.}
  \item{X}{A matrix of predictor variables.  The columns represent predictors and rows represent observations.}
  \item{lambda}{tuning parameter for fusion penalty}
  \item{H}{An indicator matrix representing the edge set of the penalty set.  The matrix is square and symmetric with dimension number of response categories, and if two categories are in the penalty set a 1 should be in the row column combination.}
  \item{tol1}{Convergence tolerance for ADMM}
  \item{tol2}{Convergence tolerance of ADMM}
  \item{TD}{Step size Adjustment for iterative step size}
  \item{rho}{ADMM step-size parameter, iterative implementation}
  \item{tau1}{The threshold parameter to 0 for final estimates}
  \item{iter}{Maximum number of iterations of the algorithm}
}
\details{
Implements the ADMM algorithm for the group fused multinomial regression estimates proposed by Price et. al.  
}
\value{
The function returns a list:
\item{Coeff}{Final coefficient estimates for the gfmr model}
\item{Approx}{Final iterates of beta before threshold}
\item{Z}{Final iterate of Z in the ADMM algorithm}
\item{lambda}{Tuning Parameter}
\item{Converge}{Indicator of algorithm convergence}
\item{NGroups}{Number of groups in final estimates}
\item{Groups}{The response category groups in the final estimates}
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
mod<-GroupFusedMulti(Response,X,lambda=2^4.3,H=Hmat2,rho=10^2,iter=50,tol1=10^-4,tol2=10^-4)
predict(mod,X[1,])
}

}
\keyword{regression}% at least one, from doc/KEYWORDS

