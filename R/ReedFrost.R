#' A function to simulate epidemics using the Reed-Frost model
#'
#' This function allows you to simulate epidemics.
#' @param Nsim.
#' @keywords Reed-Frost
#' @export
#' @examples
#' simulateReedFrost()

simulateReedFrost <- function(Nsim=100000, Npop=10000, p=0.00015, I0=1, Nit=50, fig=T)
{
  S0=Npop-I0
  I=t(matrix(rep(c(I0,rep(0,Nit-1)),Nsim),ncol = Nsim))
  S=t(matrix(rep(c(S0,rep(0,Nit-1)),Nsim),ncol = Nsim))
  for(t in 2:Nit)
  {
    I[,t]=rbinom(Nsim,size=S[,t-1],1-(1-p)^I[,t-1])
    S[,t]=S[,t-1]-I[,t]
  }

  if(fig)
    hist(apply(I,1,max),breaks = 100, freq = T)

  return(list(S=S,I=I))
}

RF_with_obs<-function(Npop=10000, p=0.00015, Np = 1, p_obs = 0.05, size = 2)
{
  initial_pop= Npop -rpois(n= Np, lambda = 5)#Npop-rpois(n = 1, lambda = 20)
  initial_I=Npop- initial_pop #rpois(n=1, lambda = 5)
  simu_epi=simulateReedFrost(Nsim=1,Npop=initial_pop, p=p, I0=initial_I, fig=F)
  X<-cbind(simu_epi$S[1,],simu_epi$I[1,]) #the "hidden" Markov chain
  Y<-rnbinom(n=length(X[,1]),mu=p_obs*X[,2], size=size) #the observered process

  return(list(X=X,Y=Y))
}

RF_SIR <- function(Y, Np=1000, Npop=10000, p=0.00015, p_obs = 0.05, size = 2)
{
  N=length(Y)

  #Sequential important sampling
  Xp<-array(rep(0,2*N*Np),dim=c(Np,N,2))      #particles
  gammap<-matrix(rep(0,N*Np),ncol=N)  #partial likelihood
  wp<-matrix(rep(0,N*Np),ncol=N)      #unormalised importance weights
  Wp<-matrix(rep(0,N*Np),ncol=N)      #normalised importance weights
  #Xrp<-matrix(rep(0,N*Np),ncol=N)     #resampled particles
  A<-matrix(rep(0,N*Np),ncol=N) #Ancestry of particle

  #Initialisation of the algorithm
  Xp[,1,1]<-Npop -rpois(n= Np, lambda = 5)#-rpois(n = Np, lambda = 20) #initialisation of the size of susceptible population for each particle
  Xp[,1,2]<-Npop- Xp[,1,1] # rpois(n= Np, lambda = 5) #initialisation of the size of the infectious population for each particle
  gammap[,1]<-dpois(Npop-Xp[,1,1], lambda = 20)*dpois(Xp[,1,2],lambda=5)*dnbinom(Y[1],mu=p_obs*Xp[,1,2],size=size)

  wp[,1]<-gammap[,1]/(dpois(Npop-Xp[,1,1], lambda = 20)*dpois(Xp[,1,2],lambda=5))
  Wp[,1]<-wp[,1]/sum(wp[,1])

  #Sequential calculation of particles and importance weights
  for(i in 2:N)
  {
    #Resampling step -using systematic resampling
    U1<-runif(1,min = 0,max = 1/Np)
    cumWj<-0
    lU<-U1
    index<-1
    for(j in 1:Np)
    {
      cumWjminus<-cumWj
      cumWj<-cumWj+Wp[j,i-1]
      if(lU<cumWj) #test if at least one Uk is between two corresponding cumulative Wp
      {
        Nji<-1+floor((cumWj-lU)*Np)
        A[index:(index+Nji-1),i]<-j
        lU<-lU+Nji/Np
      }
    }

    #IS step
    Xp[,i,2]<-rbinom(Np,size=Xp[,i-1,1],1-(1-p)^Xp[,i-1,2])
    Xp[,i,1]<-Xp[,i-1,1]-Xp[,i,2]
    wp[,i]<-dbinom(Xp[,i,2],size=Xp[,i-1,1],1-(1-p)^Xp[,i-1,2])*dnbinom(Y[i],mu=p_obs*Xp[,i,2],size=size)/dbinom(Xp[,i,2],size=Xp[,i-1,1],1-(1-p)^Xp[,i-1,2])
    Wp[,i]<-wp[,i]/sum(wp[,i])
  }
  return(list(Xp=Xp, wp=wp, A=A))
}

computeMeanPathRF<-function(particleSet)
{
  N=length(particleSet$Xp[1,,1])

  #Compute the mean particles path
  mu_Xp<-matrix(rep(0,3*N),ncol=3)
  ESS<-rep(0,N)
  low_Xp<-rep(0,N)
  up_Xp<-rep(0,N)
  for(i in 1:N)
  {
    #ESS[i]<-1/sum(particleSet$Wp[,i]^2)
    mu_Xp[i,1]<-sum(particleSet$Xp[,i,2]*particleSet$wp[,i])/sum(particleSet$wp[,i])
    #mu_Xp[i,1]<-mean(particleSet$Xp[rmultinom(N,N,prob=particleSet$wp[,i]),i,2])
    mu_Xp[i,2:3]<-quantile(sample(particleSet$Xp[,i,2],N,replace=T,prob=particleSet$wp[,i]),c(0.025,0.975))
  }

  return(list(ESS=ESS,mu_Xp=mu_Xp))
}

RF_PMMH <- function(Y, Np=1000, Nchain=1000, Npop=10000, p0=0.0001, PBar = T)
{
  marginal_likelihood <- rep(0, Nchain)
  p_estimate <- rep(0, Nchain)

  p_estimate[1]<-p0

  filterParticles<-RF_SIR(simuRF$Y, Np=Np, p=p_estimate[1])
  marginal_likelihood[1]<-prod(colSums(filterParticles$wp)/Np)

  if(PBar) pb <- txtProgressBar(min = 0, max = Nchain, style=3)

  for(i in 2:Nchain)
  {
    if(PBar) setTxtProgressBar(pb,i)

    p_prop <- p_estimate[i-1] + rnorm(1,sd=0.01/Np)

    filterParticles<-RF_SIR(simuRF$Y, Np=Np, p= p_prop)
    marginal_likelihood_prop<-prod(colSums(filterParticles$wp)/Np)

    if(runif(n=1) < marginal_likelihood_prop/marginal_likelihood[i-1])
    {
      p_estimate[i]<-p_prop
      marginal_likelihood[i]<-marginal_likelihood_prop
    } else {
      p_estimate[i]<-p_estimate[i-1]
      marginal_likelihood[i]<-marginal_likelihood[i-1]
    }
  }
  return(list(LL=marginal_likelihood,p=p_estimate))
}
