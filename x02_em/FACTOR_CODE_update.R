require(MCMCpack)
require(partitions)
require(nloptr)

############# FACTOR ANALYSIS WITH AUTOMATIC ROTATIONS TO SPARSITY #####################


# Y.......n x G 		matrix of n observations on G continuous responses
# lambda0    		spike penalty
# lambda1   			slab penalty
# start           		list consisting of 		B 		   	: G x K matrix of loading values
#													sigma 	: K x 1 vector of residual standard deviations
#												    theta	: K x 1 vector of prior inclusion probabilities                        
# K						factor dimension (threshold for the stick-breaking approximation), 
# epsilon      		margin for convergence
# alpha					tuning parameter for the IBP prior, should be of the order 1/G
# PX           			TRUE/ FALSE     if TRUE then parameter expansion used, when FALSE plain EM will be deployed
# approximate 		TRUE/FALSE   if TRUE then fast approximate M-step is used, beneficial for large data
# maxiter_PX         	maximal number of steps for PXEM, if not converged switch to monotone EM
# varimax    			0: no varimax steps used, c>0: varimax rotation performed after every c-th iteration



# ****** E-step *************************************************


E_step_omega<-function(Y,B_k,sigma_k,A_0){
	

	Matrix<-t(B_k/sigma_k)%*%B_k+solve(A_0)
	inv_Matrix<-solve(Matrix)
	list(
		omega=t(inv_Matrix%*%t(B_k/sigma_k)%*%(t(Y))),
	    Matrix=Matrix,
		inv_Matrix=inv_Matrix
	    )
        
}



laplace<-function(x,lambda){
		lambda/2*exp(-abs(x)*lambda)
		 		}


E_step_gamma<-function(B_k,theta_k,lambda0,lambda1){
	
			 		
	K<-ncol(B_k)
	G<-nrow(B_k)
	a<-cbind(apply(B_k,2,laplace,lambda=lambda1))
    b<-cbind(apply(B_k,2,laplace,lambda=lambda0))     
	theta<-matrix(rep(theta_k,each=G),G,K)
    P_star<-(theta*a)/(theta*a+(1-theta)*b)
    P_star
}


# ***************** M-step *******************************


# Exact M-step with glmnet

M_beta_lasso<-function(Y,X,sigma_k,penalties){

     require(glmnet)
     n<-nrow(X)
     Xtilde<-scale(X,center=F,scale=penalties)
	 result<-glmnet(Xtilde, Y, family=c("gaussian"),alpha = 1,lambda=2*sigma_k^2/n,
    	standardize = FALSE, intercept=FALSE, thresh = 1e-07, standardize.response=FALSE)
	as.numeric(result$beta)/penalties    
     
     }




M_beta<-function(Y,Omega,P_star,M,sigma_k,lambda1,lambda0){

     K<-ncol(Omega)
     G<-ncol(Y)
 	   n<-nrow(Y)
     beta<-matrix(0,G,K)	      	     	
     Ystar<-rbind(Y,matrix(0,K,G))
     L<-t(chol(M))
 	 Xstar<-rbind(Omega,sqrt(n)*L) 
	 
     for (i in (1:G)){		
	   		penalties<-P_star[i,]*lambda1+(1-P_star[i,])*lambda0
	   		beta[i,]<-M_beta_lasso(Ystar[,i],Xstar,sigma_k[i],penalties)
				  }
                             
      beta

}



# Approximate M-step with closed form back-substitution

lasso_back<-function(Ystar,XXstar,Astar,sigma_k,penalties){

	n<-nrow(XXstar)
	K<-ncol(XXstar)
	result<-numeric(K)
	Z<-t(XXstar)%*%Ystar/(n-K)
	result[K]<-max(abs(Z[K])-Astar[K,K]*sigma_k^2*penalties[K]/(n-K),0)*sign(Z[K])
	for ( i in (K-1):1)
	{
		index<-(i+1):K
		add<-(result[index]%*%Astar[i,index])/Astar[i,i]
		val<-Z[i]+add
		result[i]<-max(abs(val)-Astar[i,i]*sigma_k^2*penalties[i]/(n-K),0)*sign(val)-add
	}
	result
}


M_beta2<-function(Y,Omega,P_star,M,sigma_k,lambda1,lambda0){

     K<-ncol(Omega)
     G<-ncol(Y)
     n<-nrow(Y)
     beta<-matrix(0,G,K)	      	     	
     Ystar<-rbind(Y,matrix(0,K,G))
     L<-t(chol(M))
   
     Xstar<-rbind(Omega,sqrt(n)*L)
	   ch<-chol(1/n*t(Xstar)%*%Xstar)
	   Astar<-solve(ch)
	   XXstar<-Xstar%*%Astar
	 


     for (i in (1:G)){		
	   		penalties<-P_star[i,]*lambda1+(1-P_star[i,])*lambda0
	   		beta[i,]<-lasso_back(Ystar[,i],XXstar,diag(K),sigma_k[i],penalties)
				  }
				  
     beta				  
}



# Update the residual standard deviations

M_sigma<-function(Y,Omega,B_k,eta,lambda){
  
    
   sqrt((
    apply(Y-cbind(Omega)%*%t(cbind(B_k)),2,
	  function(x){sum(x^2)})+eta*lambda)/
    (nrow(Y)+eta))
}


# Update the stick-breaking fractions

M_theta_NLP<-function(P_star,alpha){
	
	require(nloptr)
	coefs<-apply(P_star,2,sum)
	K<-ncol(P_star)
	N<-nrow(P_star)

    eval_jac_g<-function(x){
	Matrix<--1*diag(K)
	for (i in (1:K-1)){
	Matrix[i,i+1]<-1}
	Matrix[K,K]<-0
	Matrix[K,1]<--1
	return(Matrix)}


	eval_g_ineq<-function(x){
	text<-paste("-x[",1:(K-1),"]+x[",(2:K),"]",sep="",collapse=",")
	text<-paste("c(",text,",-x[1])",sep="")
	eval(parse(text=text))}


	eval_grad_f<-function(x){
	paste1<-paste("-coefs[",1:(K-1),"]/","x[",1:(K-1),"]","+(N-coefs[",1:(K-1),"])/","(1-x[",1:(K-1),"])",sep="",collapse=",")
	paste2<-paste("-(alpha-1+coefs[",K,"])/","x[",K,"]","+(N-coefs[",K,"])/","(1-x[",K,"])",sep="",collapse=",")
	text<-paste(c(paste1,paste2),collapse=",")
	text<-paste("c(",text,")",sep="")
	eval(parse(text=text))}
	
	eval_f<-function(x){
	paste1<-paste("-coefs[",1:K,"]*","log(x[",1:K,"])",sep="",collapse="+")
	paste2<-paste("(N-coefs[",1:K,"])*","log(1-x[",1:K,"])",sep="",collapse="-")
	paste3<-paste("(alpha-1)*log(x[",K,"])",sep="")
	text<-paste(c(paste1,paste2,paste3),collapse="-")
	eval(parse(text=text))}


 	opts<-list("algorithm"="NLOPT_LD_MMA","check_derivatives"=F,"xtol_rel"=10^-10)
    x0<-sort(rbeta(K,1,1),decreasing=TRUE)
    res<-nloptr(x0=x0,eval_f=eval_f,eval_grad_f=eval_grad_f,eval_g_ineq=eval_g_ineq,
    		eval_jac_g_ineq=eval_jac_g,
            opts=opts,lb=rep(0,K),ub=rep(1,K))
    res$solution
	
}




# ************ The Main Program *************************************



FACTOR_ROTATE<-function(Y,lambda0,lambda1,start,K,epsilon,alpha,PX,approximate,stop,varimax,plot=TRUE){
	
	G				<-ncol(Y)
	n				<-nrow(Y)

	# initialize the procedure
 	
 	B_k			<-start$B
	B_new		<-B_k
	sigma_k  	<-start$sigma
	theta_k		<-start$theta
	
	## Jiguang: why overwriting?
	sigma_k		<-rep(1,G)
	theta_k		<-rep(0.5,K)	

	Omega		<-NULL
	Gamma_new<-matrix(rbinom(G*K,1,0.5),G,K) # draw bern 0.5 for each entries 


	wait			<-0
	niter			<-0
	eps			<-epsilon+1
	ch				<-NULL
	test			<-NULL
	A_0			<-diag(K) # identity matrix


	if(plot){myImagePlot(abs(B_k),F)  }


while(eps>epsilon){

	niter<-niter+1
	
	o<-1:K
	
	if(niter==stop){print("Switching to EM")}

	# E-step

	E_step1		    <-E_step_omega(Y,B_k,sigma_k,A_0)
	Omega			<-E_step1$omega
	M			    <-E_step1$inv_Matrix

  P_star			<-E_step_gamma(B_k,theta_k,lambda0,lambda1)     

    
 	Gamma_k			<-B_k!=0
   
   
    # M-step

	sigma_k			<-M_sigma(Y,Omega,B_k,eta=1,lambda=1) 
	theta_k			<-M_theta_NLP(P_star,alpha) # stick breaking

	# parameter expansion

   ## Jiguang: when approximate = T, PX must be false?
	 if (approximate==T){    
     B_k<-M_beta2(Y,Omega,P_star,M,sigma_k,lambda1,lambda0)
       
   			if(PX==F){
   	  				S<-M+1/n*t(Omega)%*%Omega	
      				ch<-chol(S)
  	  				B_k<-B_k%*%solve(t(ch)) # solve is matrix inverse, t is transpose
      						}
    								} else{
     B_k<-M_beta(Y,Omega,P_star,M,sigma_k,lambda1,lambda0)
    	
    			if(PX==T&niter<=stop){
   	  				S<-M+1/n*t(Omega)%*%Omega	
      				ch<-chol(S)
  	  				B_k<-B_k%*%(t(ch))
      						}
   											}
    
     if(varimax>0){
     	
         B_k<-varimax(B_k+0.0000001)$loadings
         wait<-wait+1
         if(wait==varimax){wait=0}
         				
         				  }
  	
  	if(plot){myImagePlot(abs(B_k),F)  }

  	Gamma_k			<-B_k!=0
    
     	
	Gamma_new	<-Gamma_k
 	eps					<-max(abs(B_new[,o]-B_k))
	B_new				<-B_k
	print(eps)


}
 
 	
 	B_k<-M_beta(Y,Omega,P_star,M,sigma_k,lambda1,lambda0) # last iteration using exact EM step
 	
  
 	if(plot){myImagePlot(abs(B_k),F)  }
 
  
	list(Y=Y,B=B_k,sigma=sigma_k,theta=theta_k,P_star=P_star,
	     niter=niter,Omega=Omega,rotate=ch)
}

# ****** Evaluation Run  ************************ 


# log-IBP prior

ibp<-function(Gamma,alpha){

	N<-nrow(Gamma)
	K<-ncol(Gamma)
	binaries<-apply(Gamma,2,todec)
	o<-order(binaries,decreasing=T)
	newGamma<-Gamma[,o]
	ms<-apply(newGamma,2,sum)
	Kplus<-sum(ms>0)
	Hn<-sum(1/(1:N))
	KH<-as.numeric(table(binaries))
	gamma1s<-ms+alpha/K
	gamma2s<-N-ms+1
	lfactorial(K)-sum(lfactorial(KH))+K*(log(alpha/K)-lgamma(N+1+alpha/K))+sum(lgamma(gamma1s)+
	lgamma(gamma2s))
	
}


logprior_B<-function(B,Gamma){
	dim<-sum(Gamma)
	dim*log(lambda1/2)-lambda1*sum(abs(B))
}

logprior_sigma<-function(sigma){
	-sum(0.5*log(sigma^2)+1/(2*sigma^2))
}

FACTOR_EVALUATION<-function(Y,result,lambda1,K,epsilon,alpha,PX,approximate,stop){

    Gamma<-apply(result$B,2,function(x){as.numeric(x!=0)})

    B_k<-result$B
	B_k[Gamma==0]<-0
	B_new<-B_k
	criterion<-0

	Omega<-result$Omega

	sigma_k<-result$sigma
	
	niter<-0
	eps<-epsilon+1
	ch<-NULL
	
while(eps>epsilon&niter<=stop){

	niter<-niter+1

	# E-step

 
	E_step1<-E_step_omega(Y,B_k,sigma_k,diag(K))
	Omega<-E_step1$omega
	M<-E_step1$inv_Matrix


 	
	# parameter expansion

 
	 if (approximate==T){    
     B_k<-M_beta2(Y,Omega,Gamma,M,sigma_k,lambda1,lambda0=10000)
       
   			if(PX==F){
   	  				S<-M+1/n*t(Omega)%*%Omega	
      				ch<-chol(S)
  	  				B_k<-B_k%*%solve(t(ch))
      						}
    								} else{
     B_k<-M_beta(Y,Omega,Gamma,M,sigma_k,lambda1,lambda0=10000)
    	
    			if(PX==T){
   	  				S<-M+1/n*t(Omega)%*%Omega	
      				ch<-chol(S)
  	  				B_k<-B_k%*%(t(ch))
      						}
   											}
    
  
 	sigma_k<-M_sigma(Y,Omega,B_k,eta=1,lambda=1)

	eps<-max(abs(B_new-B_k))
	B_new<-B_k
	#print(eps)

}

	B_new<-M_beta(Y,Omega,Gamma,M,sigma_k,lambda1,lambda0=10000)
 

	covariance<-B_new%*%t(B_new)+diag(sigma_k)

	criterion<-0
	criterion<-sum(dmvnorm(Y,numeric(G),covariance,log=TRUE))
	criterion<-criterion+ibp(Gamma,alpha)+logprior_B(B_new,Gamma)+logprior_sigma(sigma_k)
	list(CRITERION=criterion,B=B_new,sigma=sigma_k)
	
}



myImagePlot<-function(x,lines){
	 x<-x+min(x)
     min <- min(x)
     max <- max(x)
     yLabels <- rownames(x)
     xLabels <- colnames(x)
     title <-c()
 

 # Red and green range from 0 to 1 while Blue ranges from 1 to 0
 ColorRamp <- rgb( seq(1,0,length=100),  # Red
                   seq(1,0,length=100),  # Green
                   seq(1,1,length=100))  # Blue
# ColorRamp <- rainbow(20)

 ColorLevels <- seq(min, max, length=length(ColorRamp))
 ColorLevels <- ColorRamp
 


 # Data Map

 image(1:ncol(x),1:nrow(x),t(x)[,nrow(x):1], col=ColorLevels, xlab="",
 ylab="", axes=F, zlim=c(min,max))

 par(xpd=T)
 
 text(1:ncol(x),-0.2,labels=xLabels,srt = 45, adj = 1)
 text(0.1,1:nrow(x),labels=yLabels)
 
# text(p/2,-6,labels="Regression coefficients",cex=0.9)
# text(d/2,-6,labels="Factor loadings",cex=0.9)


 par(xpd=F)
 
# abline(v=p+0.5,col="black",lwd=2,lty=2)
 
if(lines){ 
 for (i in (1:ncol(x))){
		abline(v=0.5+i,lty=2,col="grey")
	}	}
 box(lwd=2) 
 


}


