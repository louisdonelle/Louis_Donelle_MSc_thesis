#### This code was not meant to be shared. That being said, anyone can use it, but the code was written to meet computational restrictions, and can therefore be painful to read/follow from time to time. Feel free to request clarifications or a user-friendly code by e-mailing me at louis.donelle@gmail.com ###

rm(list=ls())  
library(geoR)
library(nlme)
library(lme4)
library(glmm)
library(dplyr)
library(plyr)
library(matrixStats)
library(ggplot2)
library(boot)
library(MASS)
library("foreach")
library("doParallel")
registerDoParallel(24)

#home-made diversity functions
simpson<- function(abun){
  return(if(sum(abun)>0){1-sum((abun/sum(abun))^2)}else{NA})
}
shannon<- function(abun){
  return(if(sum(abun)>0){-sum((abun/sum(abun))*log((abun/sum(abun))))}else{NA})
}

#generate an autocorrelated environment
autocor_env_mvnorm <- function(size, alpha, C0) {
  range<-alpha
  spatDat <- expand.grid(1:size,1:size)
  x<-spatDat[,1]
  y<-spatDat[,2]
  spatDat <- data.frame(x=x,y=y)
  geo.distance <- as.matrix(dist(spatDat))
  
  Sigma<-cov.spatial(geo.distance,cov.model="spherical",cov.pars=c((1-C0)^1,range))
  diag(Sigma)<-diag(Sigma)+C0^1
  #Sigma<-Sigma+C0
  sim<-mvrnorm(n = 1, mu = rnorm(size^2,0,0), Sigma)
  return(sim)
}


#parameters
sp<-1000
k<-5000
site<-50^2
local.heter<-.1
disp.rate<-.3
range<-5
strength<-1
r.max<-5


simulation<-function(range,strength,disp.rate,disturbance,sp=250,k=1000,site=40^2,local.heter=.1,r.max=5,generation=500,init=4){
  #initial settings
  #generate spatially structured landscape
  env<- autocor_env_mvnorm(sqrt(site),range,1-strength) 
  env<-scale(env)
  env.long<-cbind(1:site,expand.grid((1:sqrt(site))-1,(1:sqrt(site))-1))
  env.long[,4]<-as.vector(env)
  colnames(env.long)<-c("site","x","y","env")
  
  #generate species pool
  sp.pool<- data.frame(
    1:sp,
    rnorm(sp),
    rlnorm(sp,-2.75,0.75)
  )
  colnames(sp.pool)<-c("ID","opt","tol")
  
  #generate a matrix giving the site ID of the 8 neighbours for each site. This consider the landscape as a torus and speeds up computation
  neighbor.matrix<- matrix(0,site,8)
  neighbor.matrix[,1]<-((0:(site-1)) - sqrt(site))%%site +1
  neighbor.matrix[,2]<- sqrt(site)*floor((((0:(site-1)) - sqrt(site))%%site)/sqrt(site))+ ((0:(site-1)) - sqrt(site)+1)%%sqrt(site) +1
  neighbor.matrix[,3]<- sqrt(site)*floor((((0:(site-1)) )%%site)/sqrt(site))+ ((0:(site-1)) +1)%%sqrt(site) +1
  neighbor.matrix[,4]<- sqrt(site)*floor((((0:(site-1)) + sqrt(site))%%site)/sqrt(site))+ ((0:(site-1)) + sqrt(site)+1)%%sqrt(site) +1
  neighbor.matrix[,5]<-((0:(site-1)) + sqrt(site))%%site +1
  neighbor.matrix[,6]<- sqrt(site)*floor((((0:(site-1)) + sqrt(site))%%site)/sqrt(site))+ ((0:(site-1)) + sqrt(site)-1)%%sqrt(site) +1
  neighbor.matrix[,7]<- sqrt(site)*floor((((0:(site-1)) )%%site)/sqrt(site))+ ((0:(site-1)) -1)%%sqrt(site) +1
  neighbor.matrix[,8]<- sqrt(site)*floor((((0:(site-1)) - sqrt(site))%%site)/sqrt(site))+ ((0:(site-1)) - sqrt(site)-1)%%sqrt(site) +1
  neighbor.matrix<-matrix(as.integer(neighbor.matrix),nrow(neighbor.matrix),ncol(neighbor.matrix))
  # same as neighbor.matrix but gives local environment value of site instead of site ID
  neighbor.env<- matrix(env.long$env[neighbor.matrix],nrow(neighbor.matrix),ncol(neighbor.matrix))
  
  
  distribution<- expand.grid(1:site,1:sp)
  colnames(distribution)<- c("site","sp")
  distribution$abun <- as.integer(rpois(nrow(distribution),k/sp/2*init)*rbinom(nrow(distribution),1,1/init)) #generate random initil abundances
  
  ##
  result<-matrix(NA,sp*site,generation/50)
  for(i in 1:generation){
    
    #### Dispersal process ####
    
    #combined probability of dispersing and surviving first step
    dispersal<-matrix(rbinom(length(distribution$abun)*8,distribution$abun,disp.rate*exp(-((neighbor.env[distribution$site,]-sp.pool$opt[distribution$sp])^2)/(2*sp.pool$tol[distribution$sp]^2))/8),length(distribution$abun),8) # if target patch is random
    
    #because probabilities of dispersing in each of the 8 neighbouring site are modeled as independant eventhough they are not, it is not impossible that we have more dispersers than individuals. In which case, we redraw until it is not the case.
    while(sum(rowSums(dispersal)>distribution$abun)!=0){
      dispersal[(rowSums(dispersal)>distribution$abun)!=0,]<-matrix(rbinom(sum(rowSums(dispersal)>distribution$abun)*8,distribution$abun[(rowSums(dispersal)>distribution$abun)!=0],disp.rate*exp(-((neighbor.env[distribution$site,]-sp.pool$opt[distribution$sp])^2)/(2*sp.pool$tol[distribution$sp]^2))[(rowSums(dispersal)>distribution$abun)!=0,]/8),sum(rowSums(dispersal)>distribution$abun),8) # if target patch is random
    }
    
    #assigned dispersers from first dispersal step to new home
    immigration<-data.frame(rep(distribution$sp,8)[dispersal>0],neighbor.matrix[distribution$site,][dispersal>0],dispersal[dispersal>0])
    colnames(immigration)<-c("sp","site","abun")
    immigration<-dplyr::summarise(group_by((immigration[order(immigration[,1], immigration[,2]),]),sp,site),sum(abun))
    colnames(immigration)<-c("sp","site","abun")
    immigration.site<-(immigration$sp-1)*site+immigration$site
    distribution.imm<-distribution
    distribution.imm$abun<-0
    distribution.imm$abun[immigration.site]<-immigration$abun
    distribution$abun[immigration.site]<-distribution$abun[immigration.site]+immigration$abun
    distribution$abun<- distribution$abun-rowSums(dispersal)
    
    #2nd dispersal step
    dispersal<-matrix(rbinom(length(distribution.imm$abun)*8,distribution.imm$abun,disp.rate*exp(-((neighbor.env[distribution.imm$site,]-sp.pool$opt[distribution.imm$sp])^2)/(2*sp.pool$tol[distribution.imm$sp]^2))/8),length(distribution.imm$abun),8) # if target patch is random
    
    while(sum(rowSums(dispersal)>distribution.imm$abun)!=0){
      dispersal[(rowSums(dispersal)>distribution.imm$abun)!=0,]<-matrix(rbinom(sum(rowSums(dispersal)>distribution.imm$abun)*8,distribution.imm$abun[(rowSums(dispersal)>distribution.imm$abun)!=0],disp.rate*exp(-((neighbor.env[distribution.imm$site,]-sp.pool$opt[distribution.imm$sp])^2)/(2*sp.pool$tol[distribution.imm$sp]^2))[(rowSums(dispersal)>distribution.imm$abun)!=0,]/8),sum(rowSums(dispersal)>distribution.imm$abun),8) # if target patch is random
    }
    
    #so on so forth until all have stop their dispersal
    while(sum(dispersal)>0){
      immigration<-data.frame(rep(distribution$sp,8)[dispersal>0],neighbor.matrix[distribution$site,][dispersal>0],dispersal[dispersal>0])
      colnames(immigration)<-c("sp","site","abun")
      immigration<-dplyr::summarise(group_by((immigration[order(immigration[,1], immigration[,2]),]),sp,site),sum(abun))
      colnames(immigration)<-c("sp","site","abun")
      immigration.site<-(immigration$sp-1)*site+immigration$site
      distribution.imm$abun<-0
      distribution.imm$abun[immigration.site]<-immigration$abun
      distribution$abun[immigration.site]<-distribution$abun[immigration.site]+immigration$abun
      dispersal<-matrix(rbinom(length(distribution.imm$abun)*8,distribution.imm$abun,disp.rate*exp(-((neighbor.env[distribution.imm$site,]-sp.pool$opt[distribution.imm$sp])^2)/(2*sp.pool$tol[distribution.imm$sp]^2))/8),length(distribution.imm$abun),8) # if target patch is random
      
      while(sum(rowSums(dispersal)>distribution.imm$abun)!=0){
        dispersal[(rowSums(dispersal)>distribution.imm$abun)!=0,]<-matrix(rbinom(sum(rowSums(dispersal)>distribution.imm$abun)*8,distribution.imm$abun[(rowSums(dispersal)>distribution.imm$abun)!=0],disp.rate*exp(-((neighbor.env[distribution.imm$site,]-sp.pool$opt[distribution.imm$sp])^2)/(2*sp.pool$tol[distribution.imm$sp]^2))[(rowSums(dispersal)>distribution.imm$abun)!=0,]/8),sum(rowSums(dispersal)>distribution.imm$abun),8) # if target patch is random
      }
      
    }
    
    
    #calculate fitness of each population in each patch based on the matching of local environment and species niche
    env.mat.stoch<-data.frame(distribution$site[distribution$abun>0],distribution$abun[distribution$abun>0]*
                                matrix(dnorm((matrix(env.long$env[distribution$site][distribution$abun>0],sum(distribution$abun>0),10)+matrix(rnorm(10*sum(distribution$abun>0),0,local.heter),sum(distribution$abun>0),10, byrow = TRUE)),
                                             matrix(sp.pool$opt[distribution$sp][distribution$abun>0],sum(distribution$abun>0),10),matrix(sp.pool$tol[distribution$sp][distribution$abun>0],sum(distribution$abun>0),10)),sum(distribution$abun>0),10))
    colnames(env.mat.stoch)[1]<-"site"
    #local growth rate (based on relative fitness)
    r<-k/10*rowSums(env.mat.stoch[,2:11]/(dplyr::mutate_at(group_by(env.mat.stoch,site),.vars = paste(rep("X",10),1:10,sep = ""),sum)[,2:11]))/distribution$abun[distribution$abun>0]
    r[is.na(r)]<-0
    #population at t+1 based on growth rate, population at t and demographic stochasticity (Poisson draw)
    distribution$abun[distribution$abun>0]<-rpois(sum(distribution$abun>0),distribution$abun[distribution$abun>0]*r*(r.max/(r.max+r)))
    distribution$abun[is.element(distribution$site,c(1:site)[rbinom(site,1,disturbance)==1])]<-0
    rm(env.mat.stoch)
    
    #save species distribution every 50 generation
    if(i%%50==0){
      result[,i/50]<-as.integer(distribution$abun)
    }
    #print(i)
    #distribution$abun[immigration.site]<-distribution$abun[immigration.site]+immigration$`sum(abun)`
  }
  return(list(result,distribution,sp.pool,env.long))
}



setwd(dir = "~/data/ibmdata/simulpop/")
foreach(rep=1:30,.combine = "c", .inorder = F, .multicombine = F) %:% 
  foreach(range=c(3,10,20),.combine = "c", .inorder = F, .multicombine = F) %:%
  foreach(disp.rate=c(.1,0.25,.5),.combine = "c", .inorder = F, .multicombine = F) %:%
  foreach(disturbance=c(0,.01),.combine = "c", .inorder = F, .multicombine = F) %:%
  foreach(init=c(1,4),.combine = "c", .inorder = F, .multicombine = F) %:%
  foreach(strength=c(1,.666,.333),.combine = "c", .inorder = F, .multicombine = F)  %dopar%{
    filename<-paste("simulpop5",range,strength,disp.rate,disturbance,init,rep,"Rdata",sep = ".")
    if(is.element(filename,list.files(path = "./", pattern= "simulpop5"))==F){ #this allow to restart simulation (after computer shutdown), without overwritng progress
      result<-list(simulation(range,strength,disp.rate,disturbance,init=init),c(range,strength,disp.rate,disturbance,init,rep))
      save(result,file = filename)
    }else{
      print(filename)
    }
  }
