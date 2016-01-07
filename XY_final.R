############################################################################
#--------------------------------------------------------------------------#
#----------------------------Projet SMC -----------------------------------#
#--------------------------------------------------------------------------#
#-------- Gautier Appert & Lionel Riou Durand & Guillaume Salha------------#
############################################################################

#### DIMENSION DU LATTICE :
d = 8

#------------------------#
#-----Packages requis----#
#------------------------#

library(rotations)
library("circular")
library("sqldf")


#----------------------------#
#-----Matrice d'adjacence----#
#----------------------------#


Adjacent=function(x,y,z,w,d){ifelse(((abs(x-z)%in%c(1,d-1))&(abs(y-w)==0))||((abs(y-w)%in%c(1,d-1))&(abs(x-z)==0)),1,0)}
M_Adj=array(NA,dim=c(d,d,d,d))
for (x in 1:d){
  for (y in 1:d){
    for (z in 1:d){
      for (w in 1:d){
        M_Adj[x,y,z,w]=Adjacent(x,y,z,w,d)
      } 
    }
  }
}


#----------------------------------------#
#-----Fonction qui trouve les voisins----#
#----------------------------------------#

find_Neighbours = function(list_vertices,M_adj){
  mat =  M_Adj[list_vertices[1,1],list_vertices[1,2],,]  
  Neighbours =  which(mat!=0,arr.ind = T)
  for (v in 2:nrow(list_vertices)){
    mat =  M_Adj[list_vertices[v,1],list_vertices[v,2],,]  
    Neighbours = rbind(Neighbours, which(mat!=0,arr.ind = T))
  }
  Neighbours = Neighbours[!duplicated(Neighbours), ]
  
  a1=as.data.frame(Neighbours)
  a2=as.data.frame(list_vertices)
  Neighbours=sqldf('SELECT * FROM a1 EXCEPT SELECT * FROM a2')
  return(as.matrix(Neighbours))
}


###########################################
#-----------------------------------------#
#-----Fonction qui trouve les liaisons----#
#-----------------------------------------#
###########################################


find_Liaisons = function(list_vertices,M_adj){
  n=nrow(list_vertices)
  mat =  M_Adj[list_vertices[n,1],list_vertices[n,2],,]  
  Neighbours =  which(mat!=0,arr.ind = T)
  
  a1=as.data.frame(Neighbours)
  a2=as.data.frame(list_vertices)[-nrow(list_vertices),]
  Liaisons=sqldf('SELECT * FROM a1 INTERSECT SELECT * FROM a2')
  return(as.matrix(Liaisons))
}


##########################################################################
#------------------------------------------------------------------------#
#-----Fonction qui calcule la somme des COSINUS (sum phaser property)----#
#------------------------------------------------------------------------#
##########################################################################



phaser=function(A1,A2,theta1,theta2){
  B=sqrt(A1^2+A2^2+2*A1*A2*cos(theta1-theta2))
  #B=sqrt((A1*cos(theta1)+A2*cos(theta2))^2+(A1*sin(theta1)+A2*sin(theta2))^2)
  delta=atan2(A1*sin(theta1)+A2*sin(theta2),A1*cos(theta1)+A2*cos(theta2))
  return(list(B=B,delta=delta))
}

#################################################################
#---------------------------------------------------------------#
#-----Calcul des parametres KAPPA et MU (a l'aide du phaser)----#
#---------------------------------------------------------------#
#################################################################



parameter=function(beta,valeurs){
  n=length(valeurs)
  if(n==1){kappa=beta;mu=valeurs}
  else{B=1;delta=-valeurs[1];
  for (j in 2:n){
    faz=phaser(B,1,delta,-valeurs[j])
    B=faz$B;delta=faz$delta;
  }
  kappa=beta*B;mu=-delta;
  }
  return(list(kappa=kappa,mu=mu))
}

####################################################
#--------------------------------------------------#
#-----Calcul du proposal von MISES distribution----#
#--------------------------------------------------#
####################################################

proposal=function(Liaisons,X,N,beta){
  valeurs=rep(NA,nrow(Liaisons))
  proposal=rep(NA,N)
  kappa=rep(NA,N);mu=rep(NA,N)
  for (i in 1:N){
    for (j in 1:nrow(Liaisons)){
      valeurs[j]=X[Liaisons[j,1],Liaisons[j,2],i]
    }
    param=parameter(beta,valeurs)
    kappa[i]=param$kappa
    mu[i]=param$mu
    proposal[i]=(rvmises(1,kappa[i])+pi+mu[i])%%(2*pi)-pi
    #proposal[i]=(as.numeric(rvonmises(1,mu[i],kappa[i]))+pi)%%(2*pi)-pi
    #proposal[i]=ifelse(proposal[i]<pi,proposal[i],proposal[i]-2*pi)
  }
  return(list(proposal=proposal,kappa=kappa))
}


Random_Neighbours=function(d,N,M_Adj){
  
  beta=1.1
  X=array(NA,dim=c(d,d,N))
  v1=sample(1:d,2,replace=T)
  lnu = rep(0,N)
  X[v1[1],v1[2],]=runif(N,-pi,pi)
  lZ=log(2*pi)
  Vt=t(matrix(v1))
  W=rep(1/N,N)
  mat =  M_Adj[Vt[1,1],Vt[1,2],,]  
  Neighbours =  which(mat!=0,arr.ind = T)
  Vt=rbind(Vt,Neighbours[sample(nrow(Neighbours),1),])
  Liaisons=find_Liaisons(Vt,M_adj)
  
  PROPE = proposal(Liaisons,X,N,beta)
  X[Vt[2,1],Vt[2,2],]= PROPE$proposal
  lnu = lnu + log(2*pi*besselI(PROPE$kappa,nu=0))
  nu = exp(lnu- max(lnu))
  lZ = lZ + max(lnu) + log(sum(W*nu))
  W=W*nu/sum(W*nu);ESS=1/(sum(W^2));
  lnu=rep(0,N)
  
  for(t in 3:(d*d)){
    
    Neighbours=find_Neighbours(Vt)
    Vt=rbind(Vt,Neighbours[sample(nrow(Neighbours),1),])
    Liaisons=find_Liaisons(Vt,M_adj)
    
    
    PROPE = proposal(Liaisons,X,N,beta)
    X[Vt[t,1],Vt[t,2],]= PROPE$proposal
    lnu = lnu + log(2*pi*besselI(PROPE$kappa,nu=0))
    nu = exp(lnu- max(lnu))
    lZ = lZ + max(lnu) + log(sum(W*nu))
    W=W*nu/sum(W*nu);ESS=1/(sum(W^2))
    #if (ESS<N/2){
    #  X=X[,,sample(N,size=N,replace=T,prob=W)];W=rep(1/N,N)
    #}
    lnu=rep(0,N)
  }
  
  return(lZ)
}


#1st variance measure

N=100
log_Z=rep(NA,100)
for (niter in 1:100){
  log_Z[niter]=Random_Neighbours(8,N,M_Adj)
  print(niter)
}

var(log_Z)
hist(log_Z)



###############################
####Deterministic visiting#####
###############################


V=matrix(NA,nr=d*d,nc=2)
t=0
for (i in 1:d){
  for (j in 1:d){
    t=t+1
    V[t,]=c(i,j)
  }
}

Diag_order=matrix(NA,nr=d*d,nc=2)
C=rep(NA,2)
for(s in 1:(d*d)){
  C=rbind(C,V[V[,1]+V[,2]==s,])
}
Diag_order=C[-1,]

Spiral_order=matrix(NA,nr=d*d,nc=2)
#bon courage...

LR_order=V


####################################
####Random Neighbour visiting#######
####################################

RN_order=matrix(NA,nr=d*d,nc=2)

v1=sample(1:d,2,replace=T)
Vt=t(matrix(v1))
mat =  M_Adj[Vt[1,1],Vt[1,2],,]  
Neighbours =  which(mat!=0,arr.ind = T)
Vt=rbind(Vt,Neighbours[sample(nrow(Neighbours),1),])
for(t in 3:(d*d)){
  Neighbours=find_Neighbours(Vt)
  Vt=rbind(Vt,Neighbours[sample(nrow(Neighbours),1),])
}

RN_order=Vt


############################################################
###  general XY (function of some order of vertices)    !###
############################################################


XY=function(d,N,M_Adj,order,NT){
  ESS.list = rep(NA,d*d)
  beta=1.1
  X=array(NA,dim=c(d,d,N))
  v1=order[1,]
  lnu = rep(0,N)
  X[v1[1],v1[2],]=runif(N,-pi,pi)
  lZ=log(2*pi)
  Vt=t(matrix(v1))
  W=rep(1/N,N)
  ESS.list[1]=100
  mat =  M_Adj[Vt[1,1],Vt[1,2],,]  
  Vt=rbind(Vt,order[2,])
  Liaisons=find_Liaisons(Vt,M_adj)
  
  PROPE = proposal(Liaisons,X,N,beta)
  X[Vt[2,1],Vt[2,2],]= PROPE$proposal
  lnu = lnu + log(2*pi*besselI(PROPE$kappa,nu=0))
  nu = exp(lnu- max(lnu))
  lZ = lZ + max(lnu) + log(sum(W*nu))
  W=W*nu/sum(W*nu);ESS=(100/N)*1/(sum(W^2));
  ESS.list[2]=ESS
  lnu=rep(0,N)
  
  for(t in 3:(d*d)){
    
    Neighbours=find_Neighbours(Vt)
    Vt=rbind(Vt,order[t,])
    Liaisons=find_Liaisons(Vt,M_adj)
    
    
    PROPE = proposal(Liaisons,X,N,beta)
    X[Vt[t,1],Vt[t,2],]= PROPE$proposal
    lnu = lnu + log(2*pi*besselI(PROPE$kappa,nu=0))
    nu = exp(lnu- max(lnu))
    #lZ = lZ + max(lnu) + log(sum(nu)) - log(N)
    lZ = lZ + max(lnu) + log(sum(W*nu))
    W=W*nu/sum(W*nu);ESS=(100/N)*1/(sum(W^2));
    ESS.list[t]=ESS
    if (ESS<=NT*100){
      X=X[,,sample(N,size=N,replace=T,prob=W)]
      W=rep(1/N,N)
    }
    lnu=rep(0,N)
  }
  
  return(list(lZ=lZ,ESS=ESS.list,X=X))
}


d=8

Adjacent=function(x,y,z,w,d){ifelse(((abs(x-z)%in%c(1,d-1))&(abs(y-w)==0))||((abs(y-w)%in%c(1,d-1))&(abs(x-z)==0)),1,0)}
M_Adj=array(NA,dim=c(d,d,d,d))
for (x in 1:d){
  for (y in 1:d){
    for (z in 1:d){
      for (w in 1:d){
        M_Adj[x,y,z,w]=Adjacent(x,y,z,w,d)
      } 
    }
  }
}



#####################
###      TEST     ###
#####################

NT=0.5
N = 1000
XY_DIAG = XY(d,N,M_Adj,Diag_order,NT)
XY_LR = XY(d,N,M_Adj,LR_order,NT)
XY_RN = XY(d,N,M_Adj,RN_order,NT)




#######################################################
### Impact du re-sampling sur le MSE (Left-Right)   ###
#######################################################





MSE_SIS = list()
MSE_RES = list()
MSE_ADA = list()


t=0
n_iter = 100
for(j in 3:10){
  N=2^(j)
  
  t=t+1
  
  NT=0
  XY_SIS = XY(d,N,M_Adj,LR_order,NT)
  MSE_SIS[[t]] = c(XY_SIS$lZ)
  
  NT=1
  XY_RES = XY(d,N,M_Adj,LR_order,NT)
  MSE_RES[[t]] = c(XY_RES$lZ)
  
  NT=0.5
  XY_ADA = XY(d,N,M_Adj,LR_order,NT)
  MSE_ADA[[t]] = c(XY_ADA$lZ)
  
  
  for(k in 2:n_iter){
    message("###")
    message("iteration    ",j,"_",k)
    NT=0
    XY_SIS = XY(d,N,M_Adj,LR_order,NT)
    MSE_SIS[[t]] = c(MSE_SIS[[t]],XY_SIS$lZ)
    
    NT=1
    XY_RES = XY(d,N,M_Adj,LR_order,NT)
    MSE_RES[[t]] = c(MSE_RES[[t]],XY_RES$lZ)
    
    NT=0.5
    XY_ADA = XY(d,N,M_Adj,LR_order,NT)
    MSE_ADA[[t]] = c(MSE_ADA[[t]],XY_ADA$lZ)
    
  }
  
}


# Avoiding numerical overflow in the variance computation
log_mse_z=function(x){log(var(exp(x-160)))+320}


mse_sis = unlist(lapply(MSE_SIS,FUN=log_mse_z))
mse_res = unlist(lapply(MSE_RES,FUN=log_mse_z))
mse_ada = unlist(lapply(MSE_ADA,FUN=log_mse_z))

results_mse_ess=cbind(rbind(mse_sis,mse_res,mse_ada),save)

setwd("C:/Users/lriou.durand/Documents/probabilistic graphical models/PROJET")
write.csv2(results_mse_ess,file="results_mse_ess.csv")

mse_sis=results_mse_ess[1,]
mse_res=results_mse_ess[2,]
mse_ada=results_mse_ess[3,]

#mse_sis = unlist(lapply(MSE_SIS,var))
#mse_res = unlist(lapply(MSE_RES,var))
#mse_ada = unlist(lapply(MSE_ADA,var))

plot(2^(3:10),mse_res,log="x",col='blue',lwd=2,type='l',ylab="log-MSE(Z)",ylim=c(328,334),xlab="N : number of particles (scale log2)")
lines(2^(3:10),mse_ada,log="x",col='red',lwd=2)
legend("topright",c("SMC sampler","Adaptative SMC"),lwd=2,col=c("blue","red"),cex=1.3)


#####################
### PLOT DE l'ESS ###
#####################

plot(XY_DIAG$ESS,type='l',main="", xlab="iteration", ylab="ESS %",col = 'blue', lwd=2)
points(XY_DIAG$ESS,pch=19)


plot(XY_LR$ESS,type='l',main="", xlab="Iterations (dimension=8X8=64)", ylab="ESS %",col = 'blue', lwd=2,ylim=c(0,100))
points(XY_LR$ESS,pch=19)
abline(a=50,b=0,col='red')
legend("bottomright",c("ESS %","Threshold"),lwd=2,col=c("blue","red"),cex=1.1)


plot(XY_RN$ESS,type='l',main="", xlab="Iterations (dimension=8X8=64)", ylab="ESS %",col = 'blue', lwd=2)
points(XY_RN$ESS,pch=19)


###########
### MSE ###
###########

MSE_DIAG = list()
MSE_LR = list()
MSE_RN = list()


t=0
NT=0.5
n_iter = 100
for(j in 5:10){
  N=2^(j)
  
  t=t+1
  
  XY_DIAG = XY(d,N,M_Adj,Diag_order,NT)
  MSE_DIAG[[t]] = c(XY_DIAG$lZ)
  
  XY_LR = XY(d,N,M_Adj,LR_order,NT)
  MSE_LR[[t]] = c(XY_LR$lZ)
  
  
  XY_RN = XY(d,N,M_Adj,RN_order,NT)
  MSE_RN[[t]] = c(XY_RN$lZ)
  
  
  for(k in 2:n_iter){
    message("###")
    message("iteration    ",j,"_",k)
    XY_DIAG = XY(d,N,M_Adj,Diag_order,NT)
    MSE_DIAG[[t]] = c(MSE_DIAG[[t]],XY_DIAG$lZ)
    
    
    XY_LR = XY(d,N,M_Adj,LR_order,NT)
    MSE_LR[[t]] = c(MSE_LR[[t]],XY_LR$lZ)
    
    
    XY_RN = XY(d,N,M_Adj,RN_order,NT)
    MSE_RN[[t]] = c(MSE_RN[[t]],XY_RN$lZ)
    
  }
  
}



mse_diag = unlist(lapply(MSE_DIAG,FUN=log_mse_z))
mse_rn = unlist(lapply(MSE_RN,FUN=log_mse_z))
mse_lr = unlist(lapply(MSE_LR,FUN=log_mse_z))


results_mse_seq=rbind(mse_diag,mse_rn,mse_lr)
write.csv2(results_mse_seq,file="results_mse_seq.csv")

df = data.frame(cbind(mse_diag,mse_rn,mse_lr))





#########################
### Graphical results ###
#########################



plot(mse_diag,type = 'l',col='blue',ylab="MSE",xlab='iteration')
points(mse_diag,pch=19)

x  <- 1:3
g <- ggplot(df, aes(x))
g <- g + geom_line(aes(y=mse_diag), colour="red")
g <- g + geom_line(aes(y=mse_rn), colour="green")
g <- g + geom_line(aes(y=mse_lr), colour="blue")
g <- g + ylab("MSE") + xlab("N particules")
g




plot(2^(5:10),mse_rn,log="x",col='green',type='l',ylim=c(328,337),ylab="log-MSE(Z)",lwd=2,xlab="N : number of particles (scale log2)")
lines(2^(5:10),mse_lr,log="x",col='blue',lwd=2)
lines(2^(5:10),mse_diag,log="x",col='red',lwd=2)
legend("topright",c("Random Neighbours","Left-Right","Diagonal"),lwd=2,col=c('green',"blue","red"),cex=1.15)





