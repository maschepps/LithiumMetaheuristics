#PFIM 4.0 
#April 2014
#Copyright (C) PFIM 4.0 - Universite Paris Diderot and INSERM.
#------------------------------------------------------------------------------------                                                                                       

#Test to accept stdin file of older PFIM versions & default values
#-----------------------------------------------------------------------------
#names.datax and names.datay objects
ngraph<-2
vec<-c("x","y")
err1<-tryCatch(names.data_test<-lapply(1:ngraph,function(ngraph,x,vec) get(x), x=paste("names.data",vec[ngraph],sep=""), vec=vec), error=function(e) 4)
 if(length(err1)==1 && err1==4)
 {
 names.datax<-rep("Time",nr)
 names.datay<-rep("Concentration",nr)
 }
 
#option object
err2<-tryCatch(get("option"), error=function(e) 4)
if(err2==4) option<-1

#covariate.model
err3<-tryCatch(get("covariate.model"), error=function(e) 4)
if(err3==4) covariate.model<-F

#IOV
err4<-tryCatch(get("n_occ"), error=function(e) -4)
if(err4==-4) n_occ<-0

#covariate_occ.model
err5<-tryCatch(get("covariate_occ.model"), error=function(e) 4)
if(err5==4) covariate_occ.model<-F
    
#numerical derivatives for user-defined models(standard function writing)
err6<-tryCatch(get("NUM"), error=function(e) 4)
if(err6==4) NUM<-F 

#choice of population, idividual or bayesian Fisher information matrix
err7<-tryCatch(get("FIM"), error=function(e) 4)
if(err7==4) FIM<-"P"

#saving the FIM, taking into account previous information
err8<-tryCatch(get("outputFIM"), error=function(e) 4)
if(err8==4) outputFIM<-""
err9<-tryCatch(get("previous.FIM"), error=function(e) 4)
if(err9==4) previous.FIM<-""

#To display only  graphs of models and/or sensitivity functions before computing the Fisher Information matrix
err10<-tryCatch(get("graph.only"), error=function(e) 4)
if(err10==4) graph.only<-F 

#Some parameters may be not estimated (not estimated = T, estimated = F) 
err11<-tryCatch(get("beta.fixed")[length(beta)], error=function(e) -1)
if(is.na(err11)==T) stop("The length of the vector 'beta.fixed' must be equal to the length of the vectors 'beta' and 'parameters")
if(err11==-1) beta.fixed<-rep(F,length(beta)) 

#graphical representation of sensitivity functions (Yes=T, No=F)
err12<-tryCatch(get("graphsensi.logical"), error=function(e) 4)
if(err12==4) graphsensi.logical<-F 

  
#-----------------------------------------------------------------------------
#---------------------- PRE - COMPUTATION ------------------------------------
#----------------------------------------------------------------------------- 

theta<-beta
#lf<-length(form)
p<-length(theta)

#ld<-length(protPD)                   #meme nombre de protocoles elementaires
#prots<-list(prot,protPD)

prots<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("prot",LETTERS[1:nr],sep=""))
sigmainter<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("sig.inter",LETTERS[1:nr],sep=""))
sigmaslope<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("sig.slope",LETTERS[1:nr],sep=""))
#bornes<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("borne",LETTERS[1:nr],sep=""))

erro<-tryCatch(graphinf<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("graph.inf",LETTERS[1:nr],sep="")),error=function(e) 1)          
#e=0 si pas d'erreur et 1 sinon#
if (length(which(unlist(erro)==1))>=1)
{
vecg<-replicate(nr,0)
graphinf<- as.vector(vecg,mode="list")
graphsup<-lapply(1:nr,function(nr,x) max(unlist(get(x[nr]))),x=paste("prot",LETTERS[1:nr],sep=""))
} else
{
graphinf<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("graph.inf",LETTERS[1:nr],sep=""))
graphsup<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("graph.sup",LETTERS[1:nr],sep=""))
}

#test si y.range exist y.rangeA<-NULL # default range
erro1<-tryCatch(y.range<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("y.range",LETTERS[1:nr],sep="")),error=function(e) 1)          
if (length(which(unlist(erro1)==1))>=1)
{
vecg1<-replicate(nr,NULL)
y.range<- as.vector(vecg1,mode="list")
} else
{
y.range<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("y.range",LETTERS[1:nr],sep=""))
}

lp<-length(prots[[1]])

if (condinit.identical==T) condinit<-rep(condinit[1],lp) else NULL

if (length(which(omega==0))==length(omega)) 
{
cat("All elements of the vector omega are equal to 0, no random effect is considered","\n")
#ls<-length(subjects); subjects<-rep(1,ls)
}

#-----------------------------------------------------------------------------
#------------------------ COMPUTE --------------------------------------------


#Function used in Population Fisher information matrix  function 
#------------------------------------
pts<-function(nr,i,prots){
protk<-c()
lprotk<-c()
for (ki in 1:nr){
protk <-c(protk,prots[[ki]][i])  
}
return(protk)
}

#Function used in sensibilite function
#------------------------------------

inter<-function(theta,kit,l)
{
	
	for (i in 1:p) {assign(parameters[i],theta[i]) }
	cond<-eval(condinit[kit])
  res<-lsoda(cond,times=c(time.condinit,theta[p+1]),formED,theta,rtol=RtolEQ,atol=AtolEQ,hmax=Hmax)
  if (nr==1)res<-res[,dim(res)[2]][2]
  else res<-res[,((dim(res)[2]-nr+1):(dim(res)[2]))][,l] [-1]
}


#Sensitivity matrix :
#--------------------
sensibilite<-function(nr,protk,kit,lprotk,condinit)
{	

	for (i in 1:p) {assign(parameters[i],theta[i]) }
	mod<-list()
	dd<-list()
	dd2<-list()
	dtot<-list()
  s<-list()
	s1<-matrix(nrow=p,ncol=0)
  s2dp<-list()
	cond<-eval(condinit[kit])	
	for (l in 1:nr){
	dd2[[l]]<-list()
	dd[[l]]<-list()
	if (lprotk[[l]]==0){
	s<-matrix(nrow=p,ncol=0)
  } else
  {
  mod[[l]]<-lsoda(cond,c(time.condinit,sort(protk[[l]])),formED,theta,rtol=RtolEQ,atol=AtolEQ,hmax=Hmax)
  if (nr==1) mod[[l]]<-mod[[l]][,dim(mod[[l]])[2]][-1]      #single response
  else { mod[[l]]<-mod[[l]][,((dim(mod[[l]])[2]-nr+1):(dim(mod[[l]])[2]))][,l][-1]}
  
  lfd<-lprotk[[l]]
  #computation of gradient and hessian  in one list
  dd[[l]]<-lapply(1:lfd,function(lfd,theta,protk,kit,l) fdHess(c(theta,protk[lfd]),inter,kit=kit,l=l)$gradient,theta=theta,protk=protk[[l]],kit=kit,l=l)
  if (lprotk[[l]]!=0) dd2[[l]]<-lapply(1:lfd,function(lfd,theta,protk,kit,l) hessian(inter,c(theta,protk[lfd]),method="Richardson",kit=kit,l=l),theta=theta,protk=protk[[l]],kit=kit,l=l)
  
  s<-t(matrix(unlist(dd[[l]]),ncol=lprotk[[l]]))[,1:p]
  
  if(lprotk[[l]]==1) s<-matrix(s,nrow=p,ncol=1) else s<-t(s)
  
  }
  if (Trand==2) {
  s<-s*theta
  }
  s1<-cbind(s1,s)
  
  if (lprotk[[l]]!=0) { 
  #matrice hessienne 
  s2d<-NULL
  s2dp[[l]]<-list()
       for (j1 in 1:p)
        {
           for (j in 1:lprotk[[l]])
            {
             if (Trand==2)
             {
             tran<-matrix(dd2[[l]][[j]][,j1][1:p]*theta,ncol=1)
             tran[j1]<-tran[j1]+s[j1,j]/theta[j1]
             s2d<-cbind(s2d,tran)
             }
             else
             s2d<-cbind(s2d,matrix(dd2[[l]][[j]][,j1][1:p],ncol=1))
            
            }
         s2dp[[l]][[j1]]<-s2d
         s2d<-NULL
        }
  }
  }
  #joindre les morceaux de chaque reponse par parametres#
     s3dp<-list()
     for (j4 in 1:p)
     {
     mo1<-NULL
     for (j3 in 1:length(s2dp))
     {
     mo<-NULL
     mo<-s2dp[[j3]][[j4]]
     mo1<-cbind(mo1,mo)
     }
     s3dp[[j4]]<-mo1
     }    
	l1<-list(unlist(mod),s1,s3dp)
 
	return(l1) 
}


#population matrix
bloc_PFIM<-function(protk,sensi,var,invvar,dvar,mod,lprotk,bout4)
{
	resB<-matrix(nrow=p+nr*2,ncol=p+nr*2)     #matrice aleatoire #
	resC<-matrix(rep(0,(p+nr*2)*p),nrow=p+nr*2,ncol=p)     #matrice incomplete option 1 covariance==0#
	resA2<-matrix(ncol=p,nrow=p)
	if (Trand==2) sensia<-sensi/theta 	else sensia<-sensi 
	resA1<-(sensia %*% invvar) %*% t(sensia)  #matrice parametres fixe#
	l1<-list()
	li<-list()
	l2<-list()
	#premiere reponse
	c<-cumsum(unlist(lprotk))
	long<-c[nr]-c[1]
  a1<-2*diag(c(bout4[[1]],rep(0,long)),length(unlist(protk)))
	ai<-invvar%*%a1%*%invvar
	resB[p+1,p+1]<-sum(diag(ai%*%a1))
	a2<-a1*mod
	resB[p+2,p+2]<-sum(diag(a2%*%invvar%*%a2%*%invvar))
	resB[p+1,p+2]<-sum(diag(ai%*%a2))
	resB[p+2,p+1]<-resB[p+1,p+2]
	l1[[1]]<-a1
	li[[1]]<-ai
	l2[[1]]<-a2
	
	if (nr!=1){
	  for (i in 1:(nr-1))
   {
	    long<-c[i]
      vec<-rep(0,length(unlist(protk)))
      if(length(protk[[i+1]])!=0){
        vec[long+1:(lprotk[[i+1]])]=bout4[[i+1]]
      }
      else{
        vec<-vec
      } 
  
	   a1<-2*diag(vec,length(unlist(protk)))
	   ai<-invvar%*%a1%*%invvar
	   resB[p+i*2+1,p+i*2+1]<-sum(diag(ai%*%a1))
	   a2<-a1*mod
      resB[p+i*2+2,p+i*2+2]<-sum(diag(a2%*%invvar%*%a2%*%invvar))
      resB[p+(i*2)+1,p+i*2+2]<-sum(diag(ai%*%a2))
      resB[p+i*2+2,p+(i*2)+1]<-resB[p+(i*2)+1,p+i*2+2]
      l1[[i+1]]<-a1
	     li[[i+1]]<-ai
	     l2[[i+1]]<-a2
  }
 
  for (i in 1: (nr)){
    for (j in  2 : (nr)){

	     resB[(p+j*2-1),p+i*2-1]<-sum(diag(li[[j]]%*%l1[[i]])) 
       resB[p+i*2-1,(p+j*2-1)]<-resB[(p+j*2-1),p+i*2-1] #resB[p+1,p+3]<-resB[p+3,p+1]
	     resB[p+j*2,p+i*2-1]<-sum(diag(li[[i]]%*%l2[[j]])) 	#resB[p+1,p+4]<-resB[p+4,p+1]
       resB[p+i*2-1,p+j*2]<-resB[p+j*2,p+i*2-1]
	     resB[p+j*2-1,p+i*2]<-sum(diag(li[[j]]%*%l2[[i]]))  #resB[p+2,p+3]<-resB[p+3,p+2]
	     resB[p+i*2,p+j*2-1]<-resB[p+j*2-1,p+i*2]
	     resB[p+j*2,p+i*2]<-sum(diag(invvar%*%l2[[j]]%*%invvar%*%l2[[i]]))
       resB[p+i*2,p+j*2]<-resB[p+j*2,p+i*2] #resB[p+2,p+4]<-resB[p+4,p+2]
	}
	}
	
	for (i in 1:p)
	{
		resA3<-invvar%*%dvar[[i]]%*%invvar
		resA2[i,]<-sapply(1:p,function(dvar,resA3,p) sum(diag(resA3%*%dvar[[p]])),dvar=dvar,resA3=resA3)
    resB1<-sensi[i,]%*%t(sensi[i,])
		resB2<-invvar%*%resB1%*%invvar
		resB[i,1:p]<-sapply(1:p,function(resB2,sensi,p) sum(diag(resB2%*%sensi[p,]%*%t(sensi[p,]))),resB2=resB2,sensi=sensi)
		resB[i,p+1]<-sum(diag(l1[[1]]%*%resB2))
		resB[i,p+2]<-sum(diag(l2[[1]]%*%resB2))

    for ( j in 2:nr) 
    {
		  resB[i,p+j*2-1]<-sum(diag(l1[[j]]%*%resB2))
		  resB[i,p+j*2]<-sum(diag(l2[[j]]%*%resB2)) 	
    }
      resB[(p+1):(p+(2*nr)),i]<-resB[i,(p+1):(p+(2*nr))]
  
  #BLOCKc
    for( j in 1:p)
    {
      resC1<-invvar%*%dvar[[j]]%*%invvar
		  resC[i,j]<-sum(diag(resC1%*%resB1))
		  resC[p+1,j]<-sum(diag(resC1%*%l1[[1]]))
		  resC[p+2,j]<-sum(diag(resC1%*%l2[[1]]))
		  for ( j2 in 2:nr) {  
		    resC[p+j2*2-1,j]<-sum(diag(resC1%*%l1[[j2]]))
		    resC[p+j2*2,j]<-sum(diag(resC1%*%l2[[j2]]))
      }
    }
  
  }
  } 
  else
  {
        #cas une seule reponse
        for (i in 1:p)
      	{
   	      resA3<-invvar%*%dvar[[i]]%*%invvar
          resA2[i,]<-sapply(1:p,function(dvar,resA3,p) sum(diag(resA3%*%dvar[[p]])),dvar=dvar,resA3=resA3)
          resB1<-sensi[i,]%*%t(sensi[i,])
          resB2<-invvar%*%resB1%*%invvar
		      resB[i,1:p]<-sapply(1:p,function(resB2,sensi,p) sum(diag(resB2%*%sensi[p,]%*%t(sensi[p,]))),resB2=resB2,sensi=sensi)
		      resB[i,p+1]<-sum(diag(l1[[1]]%*%resB2))
		      resB[i,p+2]<-sum(diag(l2[[1]]%*%resB2))
		      resB[(p+1):(p+2),i]<-resB[i,(p+1):(p+2)]	
	       #matrice bloc C
	         for(j in 1:p)
          {
            resC1<-invvar%*%dvar[[j]]%*%invvar
		        resC[i,j]<-sum(diag(resC1%*%resB1))
		        resC[p+1,j]<-sum(diag(resC1%*%l1[[1]]))
		        resC[p+2,j]<-sum(diag(resC1%*%l2[[1]]))
		       
		      }
		    }
  }

 	resA<-resA1+0.5*resA2

	resB<-0.5*resB

	resC<-0.5*resC
	
	h<-cbind(rbind(resA,resC),rbind(t(resC),resB))
  

	return(h) 
}

#individual matrix
bloc_IFIM<-function(protk,sensi,var,invvar,dvar,mod,lprotk,bout4)
{
	
  resB<-matrix(nrow=nr*2,ncol=nr*2)     #matrice aleatoire (no random effect) #
	resC<-matrix(rep(0,(nr*2)*p),nrow=nr*2,ncol=p)     #matrice incomplete option 1 covariance==0#
	resA2<-matrix(ncol=p,nrow=p)
	if (Trand==2) sensia<-sensi/theta 	else sensia<-sensi 
	resA1<-(sensia %*% invvar) %*% t(sensia)  #matrice parametres fixe#
	l1<-list()
	li<-list()
	l2<-list()
	#premiere reponse
	c<-cumsum(unlist(lprotk))
	long<-c[nr]-c[1]
  a1<-2*diag(c(bout4[[1]],rep(0,long)),length(unlist(protk)))
	ai<-invvar%*%a1%*%invvar
	resB[1,1]<-sum(diag(ai%*%a1))
	a2<-a1*mod
	resB[2,2]<-sum(diag(a2%*%invvar%*%a2%*%invvar))
	resB[1,2]<-sum(diag(ai%*%a2))
	resB[2,1]<-resB[1,2]
	l1[[1]]<-a1
	li[[1]]<-ai
	l2[[1]]<-a2
	
	if (nr!=1){
	  for (i in 1:(nr-1))
   {
	    long<-c[i]
      vec<-rep(0,length(unlist(protk)))
      if(length(protk[[i+1]])!=0){
        vec[long+1:(lprotk[[i+1]])]=bout4[[i+1]]
      }
      else{
        vec<-vec
      } 
  
	   a1<-2*diag(vec,length(unlist(protk)))
	   ai<-invvar%*%a1%*%invvar
	   resB[i*2+1,i*2+1]<-sum(diag(ai%*%a1))
	   a2<-a1*mod
      resB[i*2+2,i*2+2]<-sum(diag(a2%*%invvar%*%a2%*%invvar))
      resB[(i*2)+1,i*2+2]<-sum(diag(ai%*%a2))
      resB[i*2+2,(i*2)+1]<-resB[(i*2)+1,i*2+2]
      l1[[i+1]]<-a1
	     li[[i+1]]<-ai
	     l2[[i+1]]<-a2
  }
 
  for (i in 1: (nr)){
    for (j in  2 : (nr)){

	     resB[(j*2-1),i*2-1]<-sum(diag(li[[j]]%*%l1[[i]])) 
       resB[i*2-1,(j*2-1)]<-resB[(j*2-1),i*2-1] 
	     resB[j*2,i*2-1]<-sum(diag(li[[i]]%*%l2[[j]])) 	
       resB[i*2-1,j*2]<-resB[j*2,p+i*2-1]
	     resB[j*2-1,i*2]<-sum(diag(li[[j]]%*%l2[[i]])) 
	     resB[i*2,j*2-1]<-resB[j*2-1,p+i*2]
	     resB[j*2,i*2]<-sum(diag(invvar%*%l2[[j]]%*%invvar%*%l2[[i]]))
       resB[i*2,j*2]<-resB[j*2,i*2] 
	}
	}
	
	for (i in 1:p)
	{
		resA3<-invvar%*%dvar[[i]]%*%invvar
		resA2[i,]<-sapply(1:p,function(dvar,resA3,p) sum(diag(resA3%*%dvar[[p]])),dvar=dvar,resA3=resA3)
    #resB1<-sensi[i,]%*%t(sensi[i,])
		#resB2<-invvar%*%resB1%*%invvar
		#resB[i,1:p]<-sapply(1:p,function(resB2,sensi,p) sum(diag(resB2%*%sensi[p,]%*%t(sensi[p,]))),resB2=resB2,sensi=sensi)
		#resB[i,p+1]<-sum(diag(l1[[1]]%*%resB2))
		#resB[i,p+2]<-sum(diag(l2[[1]]%*%resB2))

    #for ( j in 2:nr) 
    #{
		#  resB[i,p+j*2-1]<-sum(diag(l1[[j]]%*%resB2))
		#  resB[i,p+j*2]<-sum(diag(l2[[j]]%*%resB2)) 	
    #}
    #  resB[(p+1):(p+(2*nr)),i]<-resB[i,(p+1):(p+(2*nr))]
  
  #BLOCKc
    for( j in 1:p)
    {
      resC1<-invvar%*%dvar[[j]]%*%invvar
		  #resC[i,j]<-sum(diag(resC1%*%resB1))
		  resC[1,j]<-sum(diag(resC1%*%l1[[1]]))
		  resC[2,j]<-sum(diag(resC1%*%l2[[1]]))
		  for ( j2 in 2:nr) {  
		    resC[j2*2-1,j]<-sum(diag(resC1%*%l1[[j2]]))
		    resC[j2*2,j]<-sum(diag(resC1%*%l2[[j2]]))
      }
    }
  
  }
  
  } 
  else
  {
        #cas une seule reponse
        for (i in 1:p)
      	{
   	      resA3<-invvar%*%dvar[[i]]%*%invvar
          resA2[i,]<-sapply(1:p,function(dvar,resA3,p) sum(diag(resA3%*%dvar[[p]])),dvar=dvar,resA3=resA3)
          #resB1<-sensi[i,]%*%t(sensi[i,])
          #resB2<-invvar%*%resB1%*%invvar
		      #resB[i,1:p]<-sapply(1:p,function(resB2,sensi,p) sum(diag(resB2%*%sensi[p,]%*%t(sensi[p,]))),resB2=resB2,sensi=sensi)
		      #resB[i,p+1]<-sum(diag(l1[[1]]%*%resB2))
		      #resB[i,p+2]<-sum(diag(l2[[1]]%*%resB2))
		      #resB[(p+1):(p+2),i]<-resB[i,(p+1):(p+2)]	
	       
         #matrice bloc C
	         for(j in 1:p)
          {
            resC1<-invvar%*%dvar[[j]]%*%invvar
		        #resC[i,j]<-sum(diag(resC1%*%resB1))
		        resC[1,j]<-sum(diag(resC1%*%l1[[1]]))
		        resC[2,j]<-sum(diag(resC1%*%l2[[1]]))
		       
		      }
		    }
  }

 	resA<-resA1+0.5*resA2

	resB<-0.5*resB

	resC<-0.5*resC
	
	h<-cbind(rbind(resA,resC),rbind(t(resC),resB))

	return(h) 
}

#bayesian matrix
#-------------------------------------------------------------------------------------------------------
bloc_BFIM<-function(sensi,invvar,dvar,omega)  
{

#Block A only
#-------
  resA1<-matrix(nrow=p,ncol=p)
  #taking into account a priori knowlegde on random effect distribution
 	omega1<-as.matrix(as.matrix(omega[,diag(omega)!=0])[diag(omega)!=0,])
 	if (Trand==2) sensia<-sensi/theta 	else sensia<-sensi
  if (Trand==2) resA1<-((sensia %*% invvar) %*% t(sensia))[diag(omega)!=0,diag(omega)!=0] +solve(diag(beta[diag(omega)!=0])%*%omega1%*%diag(beta[diag(omega)!=0]))
  else resA1<-((sensia %*% invvar) %*% t(sensia))[diag(omega)!=0,diag(omega)!=0] +solve(omega1)
 	
  resA2<-matrix(ncol=p,nrow=p)
  for (i in 1:p){
		resA3<-invvar%*%dvar[[i]]%*%invvar
		resA2[i,]<-sapply(1:p,function(dvar,resA3,p) sum(diag(resA3%*%dvar[[p]])),dvar=dvar,resA3=resA3)
  } 


 	resA<-resA1+0.5*resA2[diag(omega)!=0,diag(omega)!=0]
 	
	return(resA) 
}


#Population Fisher information matrix
#------------------------------------
fisher<-function()
{	 
	pp<-2*p+2*nr #2 parameters for the residual variance by response#
	if (FIM=="I"){pp<-p+2*nr}
	if (FIM=="B"){pp<-p-length(which(diag(omega)==0))}
  somme<-matrix(c(rep(0,pp*pp)),ncol=pp) 
	se<-numeric() 
	cv<-numeric() 
	
	for(i in 1: length(prots[[1]]))  #pour chaque protocole meme nombre de protocole elementaire en pk et pd#
	{
	    protk<-pts(nr,i,prots)
	    lprotk<-lapply(1:nr,function(nr,protk) length(protk[[nr]]),protk=protk )
	    lcumprotk<-c(1,cumsum(unlist(lprotk)))
	    if (sum(lcumprotk[-1])>0){
			sensibili<-sensibilite(nr,protk,i,lprotk,condinit)  #i the elementar protocol
			mod<-sensibili[[1]] 
			sensi<-sensibili[[2]] 
      dsensi2<-sensibili[[3]]
			bout1<-omega%*%sensi
			bout4<-list()
			t1<-0
			t<-0
			eltdvar<-list()
			for (n in 1:(nr)){
			t<-1+t1
			t1<-t-1+lprotk[[n]]
			if(lprotk[[n]]==0){
			bout4[[n]]<-NULL 
      if(n==nr){bout4[[n+1]]=Inf}
      } else {
			bout3<-sigmainter[[n]]+sigmaslope[[n]]*mod[t:t1]
			bout4[[n]]<-bout3
			}
			if(lprotk[[n]]!=0) eltdvar[[n]]<-lapply(1:p, function(p,sensi,t,t1,sigmaslope,n,theta) sigmaslope[[n]]*sensi[p,t:t1]/theta[p],sensi=sensi,t=t,t1=t1,sigmaslope=sigmaslope,n=n,theta=theta)
      }
			t<-NULL
			t1<-NULL
			
			diagsig<-diag(unlist(bout4[1:nr])^2,length(unlist(protk))) 
      invdiagsig<-solve(as.matrix(diagsig))
      
			#element dans dvar  --> diag(c(sigmaslope*(sensi[p,1:lprotk]/theta[p]),sigmaslopePD*(sensi[p,(c+1):(c+c1)]/theta[p])),length(prots))
      #pour reconstituer une liste sur les parametres et non sur le nombre de reponses
      eltdvarf<-list()
      for (j5 in 1:p)
      {
          mo1<-NULL
            for (j6 in 1:length(eltdvar))
            {
              mo<-NULL
              mo<-eltdvar[[j6]][[j5]]
              mo1<-c(mo1,mo)
            }
          eltdvarf[[j5]]<-mo1
      }
      
			
			if (FIM=="I") {
      dvar<-lapply(1:p,function(bout4,sigmaslope,theta,p,protk)
			 			2*diag(unlist(bout4[1:nr]),length(unlist(protk)))*diag(eltdvarf[[p]],length(unlist(protk))),
						bout4=bout4,sigmaslope=sigmaslope,theta=theta,protk=protk)
      fish<-bloc_IFIM(protk,sensi,diagsig,invdiagsig,dvar,mod,lprotk,bout4)}
      
      if (FIM=="B") {
      dvar<-lapply(1:p,function(bout4,sigmaslope,theta,p,protk)
			 			2*diag(unlist(bout4[1:nr]),length(unlist(protk)))*diag(eltdvarf[[p]],length(unlist(protk))),
						bout4=bout4,sigmaslope=sigmaslope,theta=theta,protk=protk)
      fish<-bloc_BFIM(sensi,invdiagsig,dvar,omega) 
      subjects<-c(1) 
      }
      
      if (FIM=="P") {
      dvar<-lapply(1:p,function(dsensi2,bout1,bout4,sigmaslope,sensi,theta,p,protk)
			 			t(dsensi2[[p]])%*%bout1+t(bout1)%*%dsensi2[[p]]+2*diag(unlist(bout4[1:nr]),length(unlist(protk)))*diag(eltdvarf[[p]],length(unlist(protk))),
							dsensi2=dsensi2,bout1=bout1,bout4=bout4,sigmaslope=sigmaslope,sensi=sensi,theta=theta,protk=protk)
			
      var<-t(bout1) %*% sensi + diag(unlist(bout4[1:nr])^2,length(unlist(protk))) 
			if (length(which(dim(var)==0))!=2) {invvar<-solve(as.matrix(var))
      fish<-bloc_PFIM(protk,sensi,var,invvar,dvar,mod,lprotk,bout4) }
      else {fish<-matrix(c(rep(0,pp*pp)),ncol=pp)}
      
      }
      }
      
      else{fish<-matrix(c(rep(0,pp*pp)),ncol=pp)} 
			somme<-somme+subjects[i]*fish 
     
		}
	
	
	#SIG<<-c(sig.inter,sig.slope,sig.interPD,sig.slopePD)
	sigmaslopen<-lapply(1:nr,function(nr,x) x[nr],x=paste("sig.slope",LETTERS[1:nr],sep=""))
  sigmaintern<-lapply(1:nr,function(nr,x) x[nr],x=paste("sig.inter",LETTERS[1:nr],sep=""))
  xn<-rbind(sigmaintern,sigmaslopen)
  xv<-rbind(sigmainter,sigmaslope)
  SIG<<-unlist(c(xv))
  SIGN<<-unlist(c(xn))
	vec2<<-SIG!=0
	SIG<<-SIG[vec2]	
	SIGN<-SIGN[vec2]
	lSIG<<-length(SIG)	
	
	
 #fixed parameters	if beta.fixed	== T
  if (FIM=="B") {
      beta.estim<-!beta.fixed[diag(omega)!=0] 
      theta<-theta[diag(omega)!=0][beta.estim]
      p<-p-length(which((beta.fixed+(diag(omega)==0))>0))
  }
  else{
  beta.estim<-!beta.fixed 
  theta<-theta[beta.estim]
  ncovfixe<-0
  p<-length(beta.estim)-length(which(beta.estim==F))
  }
      
   omega<-as.matrix(omega)
	vec<<-diag(omega)!=0 
	omega1<<-as.matrix(as.matrix(omega[,vec])[vec,]) 
	
  #shrinkage
    sh<-c()
	  
 
  if (FIM=="I"){vec<-c(beta.estim,vec2)}
	if (FIM=="P"){vec<-c(beta.estim,vec,vec2)}
  if (FIM=="B"){
    #shrinkage matrix before removing the fixed parameters
    I_W<-solve(as.matrix(somme))*solve(as.matrix(omega1))
    vec<-c(beta.estim)
    #shrinkage after removing fixed parameters
	  sh<-diag(as.matrix(I_W[vec,vec]))*100
    }  

#remove the row i and the col i of somme if omega[i,i]=0   
	somme<-somme[vec,vec]

	 
     #taking into account previous information
    if (previous.FIM!=  ""){
    previous_information<-read.table(file.path(directory,previous.FIM),header=FALSE,sep="")
    #On teste si la previous information matrix a la meme taille que la matrice actuelle. Si ce n'est pas le cas, on arrete l'execution et on affiche un message d'erreur : 
    if (nrow(previous_information) != nrow(somme) | ncol(previous_information) != ncol(somme) ) {stop ("The previous FIM has not the same size as the current FIM","\n") }
    somme<-as.matrix(previous_information)+somme
    }

    # Ecriture de la MF dans un fichier :	
    if (outputFIM !=  ""){
    #Correction
  if (ncol(somme) >= 4)
  {sommecomplete<-somme
 first.line<-c("# Project: ",project,                        '        Date: ',date(),  rep(" ",(ncol(somme)-4)))  }
 else
{ 
sommecomplete<-matrix("",nrow=nrow(somme),ncol=(4-ncol(somme)))
sommecomplete<-cbind(somme,sommecomplete)
 first.line<-c("# Project: ",project,                        '        Date: ',date() ) }

    sommeoutput<-rbind(first.line,sommecomplete)
    write.table(sommeoutput,file.path(directory,outputFIM),quote=F,append=F,row.names=F,col.names=F)
    }

	
	
pp<-dim(somme)[1] 

det<-det(as.matrix(somme)) 


	inv<-try(solve(as.matrix(somme))) 
	if(!is.null(attributes(inv)$class))
	 {
		se<-rep(NA,pp)
		cv<-se	
		crit<-(det)^(1/pp) 
		} else
     {
		se<-sqrt(diag(inv)) 
		cv1<-se[1:p]/theta*100 
		if (FIM=="P") {
		cv2<-se[(p+1):(pp-lSIG)]/diag(as.matrix(omega1))*100 
		cv3<-se[(pp-(lSIG-1)):pp]/SIG*100 
		cv<-abs(c(cv1,cv2,cv3))}
    if (FIM=="I") {
      cv3<-se[(pp-(lSIG-1)):pp]/SIG*100 
      cv<-abs(c(cv1,cv3))
    }	
    if (FIM=="B") cv<-abs(c(cv1))
		crit<-(det)^(1/pp) 
	}

return(list(formED=formED,somme=somme,subjects=subjects,prots=prots,se=se,cv=cv,det=det,crit=crit,p=p,pp=pp,xn=xn,sensibili=sensibili,sh=sh)) 

}


#GRAPHES


graph<-function(formED,prots,graphinf,graphsup,names.datax,names.datay,y.range)
{
		for (i in 1:p) {assign(parameters[i],theta[i]) }
		tt<-list()
		ff1<-list
    e<-list()
    e1<-list()
    tt<-lapply(1:nr, function(nr,graphinf,graphsup) seq(graphinf[[nr]],graphsup[[nr]],(graphsup[[nr]]-graphinf[[nr]])/1000),graphinf=graphinf,graphsup=graphsup)
	
		
		if (condinit.identical==T)
		{
				par(mfrow=c(1,nr))
				cond<-eval(condinit[1])

        e<-lapply(1:nr,function(nr,time.condinit,tt) c(time.condinit,tt[[nr]]),time.condinit=time.condinit,tt=tt)
       
				ff1<-list()
        for (li in 1 : nr){
				  ff1[[li]]<-lsoda(cond,e[[li]],formED,theta,rtol=RtolEQ,atol=AtolEQ,hmax=Hmax)
          if (nr==1) ff1[[li]]<-ff1[[li]][,dim(ff1[[li]])[2]] [-1]
          else  ff1[[li]]<-ff1[[li]][,((dim(ff1[[li]])[2]-nr+1):(dim(ff1[[li]])[2]))][,li] [-1]
         
        }
      

				for (l in 1:nr){
          if (log.logical!=F) plot(tt[[l]],ff1[[l]],type='l',xlab=paste(names.datax[l]),ylab=paste(names.datay[l]),log=log.logical,ylim=c(y.range[[l]],range(ff1[[l]]))[1:2]) else  plot(tt[[l]],ff1[[l]],type='l',xlab=paste(names.datax[l]),ylab=paste(names.datay[l]),ylim=c(y.range[[l]],range(ff1[[l]]))[1:2])
      	
         if(graph.only==F){   
          ff2<-list()
          
         for (j in 1:length(prots[[l]]))
				{     if (length(prots[[l]][[j]]) !=0){
				      ff2[[l]]<-lsoda(cond,c(time.condinit,sort(prots[[l]][[j]])),formED,theta,rtol=RtolEQ,atol=AtolEQ,hmax=Hmax)
              if (nr==1) ff2[[l]]<-ff2[[l]][,dim(ff2[[l]])[2]][-1] else ff2[[l]]<-ff2[[l]][,((dim(ff2[[l]])[2]-nr+1):(dim(ff2[[l]])[2]))][,l][-1]
               
				       points(prots[[l]][[j]],ff2[[l]],pch=paste(j),cex=2,col=paste(j))}
				       
				       
         }
             
         }
			      title(paste( names.datay[l],"model ","\n Initial Conditions = ", condinit[1],sep=" "))  
        }

		}
		else
		{
			par(mfrow=c(nr,length(prots[[1]])))
		  e<-list()
      ff1<-list()
			ff2<-list()	
			for (l in 1:nr){
			 for (ki in 1:length(prots[[l]]))
			 {
        cond<-eval(condinit[ki])
        e<-lapply(1:nr,function(nr,time.condinit,tt) c(time.condinit,tt[[nr]]),time.condinit=time.condinit,tt=tt)
        
				ff1<-list()
        for (li in 1 : nr){
				  ff1[[li]]<-lsoda(cond,e[[li]],formED,theta,rtol=RtolEQ,atol=AtolEQ,hmax=Hmax)
          if (nr==1) ff1[[li]]<-ff1[[li]][,dim(ff1[[li]])[2]][-1]
          else  ff1[[li]]<-ff1[[li]][,((dim(ff1[[li]])[2]-nr+1):(dim(ff1[[li]])[2]))][,li] [-1]
        }
        if (log.logical!=F) plot(tt[[l]],ff1[[l]],type='l',xlab=paste(names.datax[l]),ylab=paste(names.datay[l]),log=log.logical,ylim=c(y.range[[l]],range(ff1[[l]]))[1:2]) else  plot(tt[[l]],ff1[[l]],type='l',xlab=paste(names.datax[l]),ylab=paste(names.datay[l]),ylim=c(y.range[[l]],range(ff1[[l]]))[1:2])
		  
      if(graph.only==F){
         if (length(prots[[l]][[ki]]) !=0){
         ff2[[l]]<-lsoda(cond,c(time.condinit,sort(prots[[l]][[ki]])),formED,theta,rtol=RtolEQ,atol=AtolEQ,hmax=Hmax)
         if (nr==1) ff2[[l]]<-ff2[[l]][,dim(ff2[[l]])[2]][-1] 
         else  ff2[[l]]<-ff2[[l]][,((dim(ff2[[l]])[2]-nr+1):(dim(ff2[[l]])[2]))][,l][-1]
         points(prots[[l]][[ki]],ff2[[l]],pch=paste(ki),cex=2)}
	       
         }
         title(paste(names.datay[l],"model  ","\n Initial Conditions ",  condinit[ki],sep=" "))
			 }
		  }
      } 
}




#graphs for sensitivity functions


graphsensi<-function(formED,prots,graphinf,graphsup,names.datax,names.datay,y.range)
{
		for (i in 1:p) {assign(parameters[i],theta[i]) }
		tt<-list()
		df1<-list()
		dd<-list()
    e<-list()
    e1<-list()
    tt<-lapply(1:nr, function(nr,graphinf,graphsup) seq(graphinf[[nr]],graphsup[[nr]],(graphsup[[nr]]-graphinf[[nr]])/100),graphinf=graphinf,graphsup=graphsup)
	   		
     par(mfrow=c(nr,p))
				
				cond<-eval(condinit[1])
        
        e<-lapply(1:nr,function(nr,time.condinit,tt) c(time.condinit,tt[[nr]]),time.condinit=time.condinit,tt=tt)
        
        for (li in 1 : nr){
          ti<-tt[[li]]
				  lfd<-length(e[[li]])
          ki<-1
          #inter<-function(theta,kit,l)
          #inter<-function(thetah,kit,t,l,condinit)
          dd[[li]]<-lapply(1:lfd,function(lfd,theta,protk,kit,l) fdHess(c(theta,e[[li]][lfd]),inter,kit=ki,l=li)$gradient,theta=theta)
          df1[[li]]<-t(matrix(unlist(dd[[li]]),ncol=lfd))[-1,1:p] 
          
        }
         
        for (l in 1:nr){
				  for (i in 1:p){
             plot(tt[[l]],df1[[l]][,i],type='l',xlab=paste(names.datax[l]),ylab=paste("d f/d",paste(parameters[i])),ylim=range(unlist(df1[[l]])),main=paste("Response",LETTERS[l]))

          }
			       
        }
	
}

#------------------------------------------------------------------------------
#------------------------------ OUTPUT -----------------------------------

out<-function()
{
  ###GRAPHS##################################################################################################

	if(graph.logical==T) {graph(formED,prots,graphinf,graphsup,names.datax,names.datay,y.range)} else NULL 
	if(graphsensi.logical==T) {
  graphsensi(formED,prots,graphinf,graphsup,names.datax,names.datay,y.range)} else NULL 
  if (graph.only==T){options("show.error.messages"=F);stop()}
  
  ###COMPUTATION#############################################################################################

	d<-Sys.time()
	f<-fisher() 

	p<-f$p 
	pp<-f$pp 
	se<-f$se 
	cv<-f$cv 
	sh<-f$sh
	mfisher<-f$somme 
	determinant<-f$det 
	crit<-f$crit 
	
	if (FIM=="B") {
      beta.estim<-!beta.fixed[diag(omega)!=0] 
      theta<-theta[diag(omega)!=0][beta.estim]
      beta<-beta[diag(omega)!=0][beta.estim]
      parameters.estim<-parameters[diag(omega)!=0][beta.estim]
  }else{
	beta.estim<-!beta.fixed 
  theta<-theta[beta.estim]
  beta<-beta[beta.estim]
  parameters.estim<-parameters[beta.estim]
	}
  
	#correlation matrix
	corr.mat<-cov2cor(solve(mfisher))
	
	#eigenvalues
	#for fixed effects
	EigenValues<-eigen(mfisher)$values
	EVmin1<-min(EigenValues[1:p] )
	EVmax1<-max(EigenValues[1:p] )
	EVratio1<-EVmax1/EVmin1
	
	#for variance components
  if (FIM!="B"){
	EVmin2<-min(EigenValues[(p+1):pp])
	EVmax2<-max(EigenValues[(p+1):pp])
	EVratio2<-EVmax2/EVmin2
  }
  else{
  EVmin2<-NA
	EVmax2<-NA
	EVratio2<-NA
  }
  
	EV.names<-c("min","max","max/min") 
	EV<-data.frame(c(EVmin1,EVmax1,EVratio1),c(EVmin2,EVmax2,EVratio2),row.names=EV.names)
	names(EV)<-c("FixedEffects","VarianceComponents")
  
  StdError<-se[1:p] 
	RSE<-cv[1:p] 
	Beta<-theta 

	a<-data.frame(Beta,StdError,RSE,row.names=parameters.estim) 
  .<-c(rep("%",p)) 	
	if (FIM=="B"){a<-cbind(a,.,Shrinkage=sh,.); names(a)[c(dim(a)[2]-2,dim(a)[2])]<-" "}  
  else {a<-cbind(a,.)
  names(a)[dim(a)[2]]<-" "} 
	
  
  b_condition<-F	
  if (FIM=="P"){
	Omega<-diag(as.matrix(omega1)) 
	StdError<-se[(p+1):(pp-lSIG)] 
	
	RSE<-cv[(p+1):(pp-lSIG)] 
	b_condition<-F
	if(length(which(vec==F))!=length(vec)) {
	b<-data.frame(omega2=Omega,StdError,RSE,row.names=parameters[vec]) 
	.<-c(rep("%",(pp-p-lSIG))) 
	b<-cbind(b,.)
  names(b)[dim(b)[2]]<-" "  
	b_condition<-T
	}
	}
	
	if (FIM!="B"){
	StdError<-se[(pp-(lSIG-1)):pp]
	RSE<-cv[(pp-(lSIG-1)):pp]
	.<-rep("%",lSIG) 
	#sig.names<-c('sig.inter','sig.slope','sig.interPD','sig.slopePD')[vec2] 
	sig.names<-unlist(c(f$xn))[vec2] 
	c<-cbind(data.frame(SIG,StdError,RSE,row.names=sig.names),.)
	names(c)[c(1,dim(c)[2])]<-c("Sigma"," ")
	}
	
	pm<-0
	ll<-lapply(prots,length)
	if(length(which(unlist(lapply(1:nr,function(nr) lapply(1:ll[[nr]], function(ll,prots,nr) length(prots[[nr]][[ll]]),prots=prots,nr=nr)))==0))==0){
  	e<-lapply(1:nr,function(nr,subjects,prots) data.frame(times=as.character(prots[[nr]][1:length(prots[[nr]])]),subjects), subjects<-subjects,prots<-prots)
  	pm<-1
  }

	sink((paste(directory,"\\",output,sep=""))) 
	cat("PFIM 4.0 Option 2",'\n',"\n")
	cat("Project: ", project) 
	cat("\n","\n")
	cat('Date: ', date()) 
	cat("\n","\n")
	cat("\n","\n")
	cat("**************************** INPUT SUMMARY ********************************") 
	cat("\n","\n") 
	cat('Differential Equations form of the model: ','\n','\n') 
	print(attributes(formED)$srcref) 
	cat("\n","\n")
	cat('Design: ','\n')
  cat("\n","\n")
for (i in 1:nr){
    if (pm==1){
        cat('Sample times for response:',LETTERS[i],'\n')
        print(e[[i]])
        cat("\n","\n")
    } else{
        cat('Sample times for response:',LETTERS[i],'                         Number of subjects per group \n')
        cat(paste(prots[[i]],"      ", t(t(subjects)),'\n'))
        cat("\n","\n")
    }
  }
  
	cat('Initial Conditions at time', time.condinit,':','\n','\n')
	ff<-sapply(1:lp, function(lp, form) cat(paste(form[[lp]][-1]),'\n'), form=condinit)
	cat("\n","\n") 

 
  if (FIM!="I"){
  cat("Random effect model: Trand = ",Trand) 
  cat("\n","\n")
  }
  
  ff1<-sapply(1:nr,function(nr,sigmainter,sigmaslope) cat(paste('Variance error model response',LETTERS[nr],': (',sigmainter[[nr]],'+',sigmaslope[[nr]],"*f)^2\n")),sigmainter<-sigmainter,sigmaslope<-sigmaslope)
	cat("\n","\n")
	
	cat("Error tolerance for solving differential equations system: RtolEQ =", RtolEQ,", AtolEQ =", AtolEQ,", Hmax = ",Hmax)
	cat("\n","\n") 

 if (FIM=="I") cat("Computation of the Individual Fisher information matrix") 
 if (FIM=="B") cat("Computation of the Bayesian Fisher information matrix") 
 if (FIM=="P") cat("Computation of the Population Fisher information matrix: option = ",option) 
 	cat("\n","\n","\n")
 	
 if (previous.FIM!=  "")

 {	

cat('Previous FIM from file',previous.FIM)
 	cat("\n","\n")  }

   if (outputFIM!=  "")

{
cat('FIM saved in',outputFIM)

cat("\n","\n","\n") }
  
  cat("******************* FISHER INFORMATION MATRIX ******************") 
  cat("\n","\n") 
	print(unclass(mfisher)) 
	cat("\n","\n") 
cat("************************** EXPECTED STANDARD ERRORS ************************") 
	cat("\n","\n") 
	cat("------------------------ Fixed Effects Parameters -------------------------") 
  	 cat("\n","\n") 
	print(a) 
	cat("\n","\n") 
	if (b_condition==T&&FIM=="P"){
	cat("------------------------- Variance of Inter-Subject Random Effects ----------------------") 
	cat("\n","\n") 
	print(b) 
	cat("\n","\n") 
  }
  if (FIM!="B"){
  cat("------------------------ Standard deviation of residual error ------------------------ ")
  cat("\n","\n")
	print(c) 
	cat("\n","\n")
  }  
	cat("******************************* DETERMINANT ********************************") 
	cat("\n","\n") 
	cat(determinant) 
	cat("\n","\n") 
	cat("******************************** CRITERION *********************************") 
	cat("\n","\n") 
	cat(crit) 
	cat("\n","\n") 	
	cat("\n","\n") 
	
		 
	cat("******************* EIGENVALUES OF THE FISHER INFORMATION MATRIX ******************")
	cat("\n","\n") 
	print(EV)
	cat("\n","\n") 
  cat("******************* CORRELATION MATRIX ******************")
	cat("\n","\n") 
	print(unclass(corr.mat))
	cat("\n","\n") 
	
	d1<-Sys.time()
  print(d1-d)
	sink() 

  
	cat('OK ODE EVAL option 2',"\n") 	
	return(list(condinit=condinit,prots=prots,subjects=subjects,mfisher=mfisher,determinant=determinant,crit=crit,se=se,cv=cv,Shrinkage=sh,EigenValues=EigenValues,corr.matrix=corr.mat))
	
}



