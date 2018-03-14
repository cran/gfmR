#' CalcW
#'@param Beta matrix regression coefficients
#'@param X matrix of predictor valeus that creates probabilities
#'@param Interest which response vector is that of interest
#'@return List of the W matrix in the hessian update and probability for class of interest
#'@author Brad Price <brad.price@mail.wvu.edu>
## Needed for numerically stable probability calculation
## Hidden from user

CalcW<-function(Beta,X,Interest){
  LO=matrix(0,dim(X)[1],dim(Beta)[2])
  for(m in 1:dim(Beta)[2]){
    LO[,m]=X%*%Beta[,m]
  }
  
  A=apply(LO,1,max)
  Bottom=log(exp(-A)+apply(exp(LO-A),1,sum))
  Prob=exp(LO[,Interest]-A-Bottom)
  
  W2=rep(0,dim(X)[1])
  for( i in 1:dim(X)[1]){
    W2[i]=(1-Prob[i])*Prob[i]
  }
  W=diag(W2,dim(X)[1],dim(X)[1])
  return(list(W=W,Probs=Prob))
}


#' CalcW
#'@param Beta matrix regression coefficients
#'@param X matrix of predictor valeus that creates probabilities
#'@return Corresponding probabilities for each vector
#'@author Brad Price <brad.price@mail.wvu.edu>
## Numerically stable probability calculation
## Hidden from User
CalcProbs<-function(Beta,X){
  Probs=matrix(0,dim(Beta))
  LO=matrix(0,dim(X)[1],dim(Beta)[2])
  for(m in 1:dim(Beta)[2]){
    LO[,m]=X%*%Beta[,m]
  }
  A=apply(LO,1,max)
  Bottom=log(exp(-A)+apply(exp(LO-A),1,sum))
  Probs=exp(LO-A-Bottom)
  Probs=cbind(Probs,1-apply(Probs,1,sum))
  return(Probs)
}	

#' CalcW
#'@param Beta matrix regression coefficients
#'@param X matrix of predictor valeus that creates probabilities
#'@return denominators for probabilties and scalars for each observation
#'@author Brad Price <brad.price@mail.wvu.edu>
#'Hidden from User

CalcBottom<-function(X,Beta){
  LO=matrix(0,dim(X)[1],dim(Beta)[2])
  for(m in 1:dim(Beta)[2]){
    LO[,m]=X%*%Beta[,m]
  }
  A=apply(LO,1,max)
  Bottom=log(exp(-A)+apply(exp(LO-A),1,sum))
  VM=list(Bottom=Bottom,A=A)
  return(VM)
}

#' OneD
#'@param C number of response categories
#'@param p number of variables in X
#'@param H matrix defining relationship of the penalty set
#'@param interest1 variable number of first variable
#'@param interest2 variable number 2 that we are interested in
#'@return A Matrix for a single pair in the penalty
#'@author Brad Price <brad.price@mail.wvu.edu>
## Creates an A matrix for a single pair of fusion variables

OneD<-function(C,p,H,interest1, interest2){
  D1=NULL
  for(m in 1:C-1){
    if((m!=interest1 && m!=interest2 || H[interest1,interest2]==0)){
      D1=cbind(D1,diag(0,p,p))
    }
    if(m==interest1 && H[interest1,interest2]==1){
      D1=cbind(D1,diag(1,p,p))
    }
    if(m==interest2 && interest2!=C && H[interest1,interest2]==1){
      D1=cbind(D1,diag(-1,p,p))
    }
  }
  D1=D1[,-(1:p)]
  return(D1)
}


#' AllD
#'@param C number of response categories
#'@param p number of variables in X
#'@param H matrix defining relationship of the penalty set
#'@author Brad Price <brad.price@mail.wvu.edu>
#'@reuturn the A matrix in the ADMM
AllD<-function(C,p,H){
  D2=NULL
  for(k in 1:(C-1)){
    for(m in (k+1):C){	
      D2=rbind(D2,OneD(C,p,H,k,m))
    }
  }	
  return(D2)
}

#' UpdateBeta
#'@param Y2 matrix of response variables with category counts in the columns
#'@param X matrix of predictor values
#'@param Z the Dual variables
#'@param U the Lagrangian variables (Note this is done is a scaled form)
#'@param rho the step-size parameter
#'@param H the matrix defining the penalty set
#'@param tol the convergence tolerance
#'@param Inits Is there an intialization
#'@param iters the maximum number of iterations 
#'@return Beta the regression coefficients from the update in the ADMM
#'@author Brad Price <brad.price@mail.wvu.edu>


UpdateBeta<-function(Y2,X,Z,U,rho,H,tol=10^-7,Inits=FALSE,iters=10^3){
  if(Inits==FALSE){
    p=dim(X)[2]
    C=dim(Y2)[2]
    BetaN=matrix(0,p,C)}else{
      BetaN=Inits
    }
  
  BetaN1=BetaN+1
  C=dim(BetaN)[2]
  it=0
  # while(sum(abs(BetaN-BetaN1)^2)/sum(abs(BetaN1)^2)>tol && it<=iters){
  while(sum(abs(BetaN-BetaN1))>tol && it<=iters){
    BetaN1=BetaN
    for(k in 1:(C-1)){
      Beta2=BetaN[,k]
      W2=CalcW(BetaN[,-C],X,k)
      Front=solve(t(X)%*%W2$W%*%X+diag(2*sum(H[k,])*rho,p,p))
      One=t(X)%*%(Y2[,k]-W2$Probs)
      Two=rep(0,p)
      for(j in 1:C){
        Two=Two+(BetaN[,k]-BetaN[,j]-Z[[k]][,j]+U[[k]][,j])*H[k,j]
      }
      
      BetaN[,k]=BetaN[,k]+Front%*%(One-2*rho*Two)
      Change=sum(abs(BetaN[,k]-Beta2))
      ##iter=iter+1
    }
    ##print(BetaN)
    it=it+1
  }
  Out=list(Beta=BetaN)
  return(Out)
}


#' UpdateZU
#' @param BetaN current iterate of Beta Matrix
#' @param Z a place holder for Z variables
#' @param U the current iterates of the U variables
#' @param D the A matrix of ADMM corresponding
#' @param lambda the tuning parameter
#' @param rho step-size parameter
#' @param H The matrix defining the penalty set
#' @param tol1 Tolerance on residuals
#' @param tol2 Tolerence on residuals
#' @param mu stepzise adjustment
#' @param TD stepsize adjustment
#' @return Updated versions of Z, U, and rho. 
#' @author Brad Price <brad.price@mail.wvu.edu>


UpdateZU<-function(BetaN,Z,U,D,lambda,rho,H,tol=10^-4,tol2=10^-4,mu=10,TD=2){
  C=dim(BetaN)[2]
  p=dim(BetaN)[1]
  U1=U
  Z1=Z
  Z1M=NULL
  for(k in 1:(C-1)){
    for(m in (k+1):C){
      Z1M=c(Z1M,Z[[k]][,m])
    }
  }
  Tol1=0
  Tol2=0
  Tol3=0
  for(k in 1:(C-1)){
    for(m in k:C){
      if(H[k,m]==1){
        Acm=BetaN[,k]-BetaN[,m]+U[[k]][,m]
        Z[[k]][,m]=Acm*pmax(1-(lambda/rho)*(1/sqrt(sum(Acm^2))),0)
        Z[[m]][,k]=-Z[[k]][,m]
        U[[k]][,m]=U[[k]][,m]+(BetaN[,k]-BetaN[,m]-Z[[k]][,m])
        U[[m]][,k]=U[[m]][,k]+(BetaN[,m]-BetaN[,k]-Z[[m]][,k])
        Tol1=Tol1+sum((U1[[k]][,m]-U[[k]][,m])^2)
        Tol2=Tol2+sum((BetaN[,k]-BetaN[,m]-Z[[k]][,m])^2)
      }
    }
  }
  Tol2=sqrt(Tol2)
  Tol1=sqrt(Tol1)
  
  
  ZM=NULL
  U1M=NULL
  for(k in 1:(C-1)){
    for(m in (k+1):C){
      ZM=c(ZM,Z[[k]][,m])
      U1M=c(U1M,U[[k]][,m])
    }
  }
  Tol3=sqrt(sum((rho*t(D)%*%(ZM-Z1M))^2))
  Dual=sqrt(length(U1M))*tol2+tol2*sqrt(sum((t(D)%*%U1M)^2))
  Pri=sqrt(length(as.vector(BetaN[,-C])))*tol2+tol2*max(sqrt(sum((D%*%as.vector(BetaN[,-C]))^2)),sqrt(sum(Z1M^2)),sqrt(sum(U1M^2)))
  ##print(Tol1)
  ##print(Tol2)
  ##print(Tol3)
  Converge=FALSE
  if(Tol1<=tol && Tol2<=Pri && Tol3<=Dual){
    Converge=TRUE
  }
  
  if(Tol2>mu*Tol3){
    rho=TD*rho
    for(k in 1:C){
      U[[k]]=U[[k]]/TD
    }
  }
  if(Tol3>mu*Tol2){
    rho=rho/TD
    for(k in 1:C){
      U[[k]]=U[[k]]*TD
    }
  }
  
  Out=list(Z=Z,U=U,Converge=Converge,rho=rho)
  return(Out)
}

#' GroupFusedMulti
#' @param Y Matrix of category counts for each Y 
#' @param X Matrix of predictor variables of each X
#' @param lambda tuning parameter for GFMR
#' @param H Matrix representing set in L
#' @param tol1 convergence tolerance for ADMM
#' @param tol2 convergence tolderance for ADMM
#' @param  TD update parameter for ADMM
#' @param rho step-size parameter for ADMM
#' @param tau1 ADMM adjustment of rho
#' @param iter Max iterations of ADMM for convergence
#' @result These report the GFMR estimates with tuning parammeter lambda, number of gorups, Z and U, and an indicator of converegence
#' @author Brad Price <brad.price@mail.wvu.edu>
## This is the workhorse function of GFMR and returns much information that a user will want to see. 

GroupFusedMulti<-function(Y,X,lambda,H,tol1=10^-7,tol2=10^-7,TD=2,rho=10^-8,tau1=10^-9,iter=1e3){
  C=dim(Y)[2]
  p=dim(X)[2]	
  diag(H)=rep(0,C)
  mu=10
  U=vector("list",length=C)
  for(m in 1:C){
    U[[m]]=matrix(0,p,C)
  }
  Z=U
  D=AllD(C,p,H)
  it2=0
  CONV=FALSE
  while(CONV==FALSE && it2<iter){
    Out1=UpdateBeta(Y,X,Z,U,rho=rho,H=H,tol=tol1,iters=100)$Beta
    Out2=UpdateZU(Out1,Z,U,D,lambda,rho=rho,H=H,tol=tol1,tol2,TD=TD)
    Z=Out2$Z
    U=Out2$U
    rho=Out2$rho	
    it2=it2+1
    CONV=Out2$Converge
  }
  
  BetaN=Out1
  BetaRes=BetaN
  for(k in 1:(C-1)){
    BetaRep=BetaN[,k]
    lot=0
    for(j in 1:C){
      if(j!=k){
        if(sum(Z[[k]][,j]^2)==0 && H[k,j]==1){
          BetaRep=cbind(BetaRep,BetaN[,j])
          lot=lot+1
        }
      }
    }
    if(lot>0){
      BetaRes[,k]=apply(BetaRep,1,mean)
    }else{
      BetaRes[,k]=BetaRep
    }
  }
  
  for(k in 1:(C-1)){
    if(sum(abs(BetaRes[,k]))<tau1*p || Z[[k]][1,C]==0 && H[k,j]==1){
      BetaRes[,k]=rep(0,p)
    }
  }
  Groups=vector(mode="list",0)
  Set=1:C
  while(length(Set)>0){
    LG=length(Groups)
    M=which((BetaRes[,Set[1]]==BetaRes)[1,]==TRUE)
    for(i in 1:length(M)){
      OutNow=which(Set==M[i])
      Set=Set[-OutNow]
    }
    Groups[[LG+1]]=M
  }		
  NGroups=length(Groups)		
  
  
  
  
  Output=list(Coeff=BetaRes, Approx=BetaN, Z=Z, lambda=lambda,Converge=CONV,NGroups=NGroups,Groups=Groups)
  class(Output)="gfmR"
  return(Output)
}


### Generic print function to return 

print.gfmR<-function(x,...){
  obj=x
  cat("Number of groups")
  print(obj$NGroups)
  cat("The corresponding groups are")
  print(obj$Groups)
}


predict.gfmR<-function(object,newdata,type="probs",...){
  obj=object
  if(type=="probs"){
    return(CalcProbs(Beta = obj$Coeff[,-dim(obj$Coeff)[2]],X = cbind(1,newdata)))
  }
  if(type=="response"){
    return(cbind(1,as.vector(newdata))%*%obj$Coeff)
  }
}



#' GroupFusedMultiL 
#' @param Set is a set containing all elements of interest of GroupFusedMulti function
#' @return a list function of GroupFusedMulti

GroupFusedMultiL<-function(Set){
  
  AM=GroupFusedMulti(Set$Y,Set$X,lambda=Set$lambda,H=Set$H,rho=Set$rho,iter=50,tol1=10^-4,tol2=10^-4)
  Output=list(Groups=AM$Groups,NGroups=AM$NGroups,Coeff=AM$Coeff)
}

#' GFMR.cv
#' @param Y matrix of response category combinations 
#' @param X matrix of predictor category combinations
#' @param lamb the lambda tuning parameters that would be needed for prediction
#' @param sampID the identifiers for the cross validation sets for the cross validation function
#' @param H the matrix identifying the penalty set
#' @param n.cores if parallel needed specify the number of cores
#' @return cross validaiton scores and standard deviations along with optimal results
#'@author Brad Price <brad.price@mail.wvu.edu>

## The workhorse function of the cross valdiation process. 

GFMR.cv<-function(Y,X,lamb,sampID,H,n.cores=1,rho=10^-8,...){
  br=matrix(0,max(sampID),length(lamb))
  for(m in 1:length(lamb)){
    cat(m,"out of", length(lamb),"values")
    res=NULL
    print(m)
    S1<-vector("list",max(sampID))
    for(j in 1:max(sampID)){
      set1=which(sampID==j)
      S1[[j]]=list(Y=Y[-set1,],X=X[-set1,],lambda=lamb[m],H=H,rho=rho)
    }
    
    
    L2<-mclapply(S1,GroupFusedMultiL,mc.cores=n.cores)
    for(j in 1:max(sampID)){
      res=c(res,sum(log(CalcProbs(X=X[sampID==j,],Beta=L2[[j]]$Coeff[,-dim(Y)[2]])^Y[sampID==j,])))
    }
    br[,m]=res
  }
  out=list(vl=apply(br,2,mean),vl.sd=apply(br,2,sd),lambda=lamb,vl.mat=br)
  class(out)<-"gfmR.cv"
  return(out)
}

print.gfmR.cv<-function(x,...){
  obj=x
  cat("Best Tuning Parameter")
  print(obj$lambda[which(x$vl==max(x$vl)[1])])
  cat("Validation Likelihood")
  print(x$vl)
  cat("Standard Deviation of Validation Likelihood")
  print(x$vl.sd)
}

#'FusedMultiAIC
#' @param Y matrix of response category combinations 
#' @param X matrix of predictor category combinations
#' @param H the matrix identifying the penalty set
#' @param Start the lambda tuning parameters that would be needed for prediction
#' @param tau the update for the ADMM
#' @param rho the step size parameter for the ADMM. 
#' @return reuslting statitics and group structures of candidate values
## Employs a gird search until at least one group structure with each possible number of groups is found.
## Contains the AIC computaiton off the multinomial function under the assumption of equal probabilities.  Allows for a comparison under a different
## number of response categories.  
FusedMultiAIC<-function(Y,X,H,Start,tau=10^-9,rho1=1){
  BICVals<-rep(0,dim(Y)[2])
  Lambdas<-rep(0,dim(Y)[2])
  GroupL=BICVals
  C=dim(Y)[2]
  GroupsL=vector("list",C)
  Target=1
  UL=Start
  LL=0
  New=0
  while(Target<=C){
    if(Lambdas[Target]==0){
      FusedM<-GroupFusedMulti(Y,X,Start,H,tau1=tau,rho=rho1)
      NG=FusedM$NGroups
      GR=FusedM$Groups
      if(NG==Target){
        YR=matrix(0,dim(Y),NG)
        for(m in 1:NG){
          if(length(GR[[m]])>1){
            YR[,m]=apply(Y[,GR[[m]]],1,sum)
          }
          else{
            YR[,m]=Y[,GR[[m]]]
          }
        }
        if(NG>1){
          invisible(ModOut<-multinom(YR~X[,-1],trace=FALSE))
          
          BICVals[Target]=ModOut$dev+2*sum(apply(YR,2,sum)*log(unlist(lapply(GR,length))))+2*dim(X)[2]*(NG-1)
          ##ModOut$dev
          Lambdas[Target]=Start
          GroupsL[[Target]]=GR
          
          
        }else{
          BICVals[Target]=dim(Y)[1]*log(C)*2
          Lambdas[Target]=Start
          GroupsL[[Target]]=GR
        }	
      }
      if(NG>Target){
        LL=Start
        
        if(Lambdas[1]==0){
          Start=Start+1
        }else{
          Start=mean(c(Start,UL))
          
        }	
      }
      if(NG<Target){
        UL=Start
        Start=mean(c(LL,Start))
        
      }
      
    }else{
      Target=Target+1	
      New=0
      Start=Start/2
      LL=0
      UL=Lambdas[Target-1]
    }
  }	
  
  Avgs=BICVals
  Best=which(Avgs==min(Avgs))
  
  Output<-list(Dev=Avgs,Lambdas=Lambdas,Optimal=Lambdas[Best[1]],GroupsStruc=GroupsL, OptGroups=Best)
  return(Output)
}
