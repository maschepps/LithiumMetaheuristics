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

#if dose is missing
missdose<-tryCatch(get("dose"), error=function(e) T)
if(is.numeric(missdose)==F) dose<-doses<-0



###########################################################################
# 	Computation of the Population Fisher Information matrix       	      #     
#	 for PK Design with analytical form of the model			 									#							       
###########################################################################

			
#---------------------- PRE - COMPUTATION ------------------------------------


subjects.init<-subjects
subjects.initA<-subjects.init
theta<-beta
if(NUM==F){
lf<-length(form) #nombre de forme analytique
#nombre de form dans chaque reponse pour savoir quoi prendre en tf
#recover of the model
formg<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("form",LETTERS[1:nr],sep=""))  #sous forme de list
lformg<-lapply(formg,length)  #length(of each model tf or not)
c<-c(0,cumsum(lformg))
}
doses<-dose
p<-length(theta)
omega<-as.matrix(omega)

#Partie liee a l'implementation pour iov
#-------------------------------------------------------------------------------
if (n_occ<=1) 
{
  n_occ=1
  covariate_occ.model<-F
}
if (covariate_occ.model==F)
{
  covariate_occ.name<-list(NULL)
  covariate_occ.category<-list(NULL)
  covariate_occ.sequence<-list(list(NULL))
  covariate_occ.proportions<-list(NULL)
  parameter_occ.associated<-list(NULL)
  beta.covariate_occ<-list(NULL)
  locc<-0
  nb_seq<-1
}
sequence.proportions<-covariate_occ.proportions

 
pn<-n_occ
v<-rep(0,p)
param.asso.occ<-list()
occ.effect<-list()
beta.occ<-list()
nb_seq<-c()
  
for (kj in 1:length(covariate_occ.name)){
    nb_seq[kj]=length(covariate_occ.sequence[[kj]])
    param.asso.occ[[kj]]<-list()
    occ.effect[[kj]]<-list()
    beta.occ[[kj]]<-list()
    for (l in 1:nb_seq[kj]) {
      param.asso.occ[[kj]][[l]]<-list()
      occ.effect[[kj]][[l]]<-list()
      beta.occ[[kj]][[l]]<-list()
      for (i in 1:n_occ)  {param.asso.occ[[kj]][[l]][[i]]=v;occ.effect[[kj]][[l]][[i]]=v;beta.occ[[kj]][[l]][[i]]=v} 
      }
     
    if (n_occ>1){
      for (l in 1:nb_seq[kj]) {
  
        num<-c()
        if (covariate_occ.model==T){
          num[kj]=0
          for (i in which(covariate_occ.sequence[[kj]][[l]]!=covariate_occ.category[[kj]][1])){
            num[kj]=num[kj]+1
            param.asso.occ[[kj]][[l]][[i]]=parameter_occ.associated[[kj]]
            
            for (j  in 1:length(param.asso.occ[[kj]][[l]][[i]])){
              
              occ.effect[[kj]][[l]][[i]][which(param.asso.occ[[kj]][[l]][[i]][j]==parameters)]=1
              beta.occ[[kj]][[l]][[i]][which(param.asso.occ[[kj]][[l]][[i]][j]==parameters)]=beta.covariate_occ[[kj]][[num[kj]]][j]
            } 
         }
       }
     }
  }
}

#les combinaisons des sequences
sequence.reco<-list()
for (kj in 1:length(covariate_occ.name)){
  sequence.reco[[kj]]=c(1:nb_seq[kj])
}
comb.seq.poss<-unique(t(combn(unlist(sequence.reco),length(covariate_occ.name))))
ll.cov_occ<-length(covariate_occ.name)
for (i in 1: ll.cov_occ)
{
  nline<-which(comb.seq.poss[,i] > max(sequence.reco[[i]]) | comb.seq.poss[,i] < min(sequence.reco[[i]]))
  ref<-integer(length = 0)
  if (!identical(nline,ref)) comb.seq.poss<-comb.seq.poss[-nline,]
}
#nb de sequences possibles          
lposs<-dim(comb.seq.poss)[1]

#les propportions de chaque combinaison de sequence
comb.seq.prop<-c()
nbcov_occ<-length(covariate_occ.name)
if (covariate_occ.model==T){
for (i in 1: lposs){
comb.seq.prop[i]<-sequence.proportions[[1]][comb.seq.poss[i,1]]
if (nbcov_occ>1){
  for (j in 2:nbcov_occ){
    comb.seq.prop[i]=comb.seq.prop[i]*sequence.proportions[[j]][comb.seq.poss[i,j]]
  }
 }
}
}else{comb.seq.prop<-1}

#le vecteur des parametres fixes prenant compte des effets covariables changeant entre les periodes

theta1<-list()
Pt1<-list()
vecP1<-list()
lP1<-list() 
for (kj in 1:length(covariate_occ.name)){
  theta1[[kj]]<-list()
  Pt1[[kj]]<-list()
  vecP1[[kj]]<-list()
  lP1[[kj]]<-list()
  for (l in 1:nb_seq[kj]) {
  theta1[[kj]][[l]]<-list()
  Pt1[[kj]][[l]]<-list()
  vecP1[[kj]][[l]]<-list()
  lP1[[kj]][[l]]<-rep(0,pn)
  for (i in 1:pn)  { 
    Pt1[[kj]][[l]][[i]]<-c(occ.effect[[kj]][[l]][[i]])
    vecP1[[kj]][[l]][[i]]<-Pt1[[kj]][[l]][[i]]!=0
    Pt1[[kj]][[l]][[i]]<-Pt1[[kj]][[l]][[i]][vecP1[[kj]][[l]][[i]]]
    lP1[[kj]][[l]][i]<-length(Pt1[[kj]][[l]][[i]])
    if (lP1[[kj]][[l]][i] == 0) {
    theta1[[kj]][[l]][[i]] = beta
    }
    if (lP1[[kj]][[l]][i] > 0){
      if (Trand==1) theta1[[kj]][[l]][[i]]<-beta+beta.occ[[kj]][[l]][[i]]
      if (Trand==2) theta1[[kj]][[l]][[i]]<-beta*exp(beta.occ[[kj]][[l]][[i]])
    }
  }
 }
}

vecP<-list()
beta.occ2<-list()
lP<-list() 
for (i in 1: lposs){
 vecP[[i]]<-list()
 lP[[i]]<-rep(0,pn)
 beta.occ2[[i]]<-list()
 for (occ in 1: n_occ){
  vecP[[i]][[occ]]<-vecP1[[1]][[comb.seq.poss[i,1]]][[occ]]
  beta.occ2[[i]][[occ]]<-beta.occ[[1]][[comb.seq.poss[i,1]]][[occ]]
  lP[[i]]<-lP1[[1]][[comb.seq.poss[i,1]]]
  if (nbcov_occ>1){
   for (j in 2:nbcov_occ){
    vecP[[i]][[occ]]=vecP[[i]][[occ]]*vecP1[[j]][[comb.seq.poss[i,j]]][[occ]]
    beta.occ2[[i]][[occ]]=beta.occ2[[i]][[occ]]+beta.occ[[j]][[comb.seq.poss[i,j]]][[occ]]
    lP[[i]]<-lP[[i]]+lP1[[j]][[comb.seq.poss[i,j]]]
   }
  }
 }
}


theta2<-list()
  for (l in 1:lposs) {
  theta2[[l]]<-list()
  for (i in 1:n_occ)  { 
    if (lP[[l]][i] == 0) {
    theta2[[l]][[i]] = beta
    }
    if (lP[[l]][i] > 0){
      if (Trand==1) theta2[[l]][[i]]<-beta+beta.occ2[[l]][[i]]
      if (Trand==2) theta2[[l]][[i]]<-beta*exp(beta.occ2[[l]][[i]])
    }
  }
 }

#valeurs des effets covariables qui dependent de l'occasion  
beta.occ1<-unlist(beta.covariate_occ)    
locc<-length(beta.occ1)

#pour prendre en compte les variances des effets aleatoires lies a IOV
if (pn==1) gamma<-diag(NULL)
OMEGA<-diag(c(diag(omega),rep(c(diag(gamma)),pn))) #la matrice OMEGA qui regroupe omega et gamma

#----------- Fin de la partie precomputation pour IOV -------------------------------------------



#Dform<-lapply(1:lf,function(p,form,parameters,lf) lapply(1:p,function(p,form,parameters) D(form,parameters[p]), form=form[lf],parameters=parameters),p=p,form=form,parameters=parameters)
#recover the inital data for all the models
if(identical.times==T){
prots<-lapply(1:nr,function(nr,x) get(x),x=paste("prot",LETTERS[1],sep=""))
}else
prots<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("prot",LETTERS[1:nr],sep=""))
 
Q<-length(prots[[1]])
#l_prote<-lapply(1:nr,function(p,forma,parameters,lf) lapply(1:p,function(p,forma,parameters) D(forma,parameters[p]), forma=forma[lf],parameters=parameters),p=p,forma=forma,parameters=parameters)
#nombre de prelevement par prot elementaire
#l_prote<-lapply(1:nr,function(nr,Q,prots) lapply(1:Q, function(nr,prots,Q) length(prots[[Q]]),Q=Q,prots=prots[[nr]]),Q=Q,prots=prots)
#somme des prevement spa rpto selon les reponses
l_prote<-list()
l2<-list()
for (l in 1:Q){
  l_prote[[l]]<-0
 
  for (i in 1:nr){
   l_prote[[l]]<-sum(l_prote[[l]],length(prots[[i]][[l]]))
 #l2[[l]]<-sum(l2[[l]],l_prote[[i]][[l]])
  
   l2[[i]]<-lapply(1:Q, function(Q,prots) length(prots[[Q]]),prots=prots[[i]])
   
 }
}

#lapply(1:Q,function(nr,l_prote,Q) lapply(1:nr, function(Q,l_prote,nr) sum(l_prote[[nr]][[Q]]),nr=nr,l_prote=l_prote),nr=nr,l_prote=l_prote)

sigmainter<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("sig.inter",LETTERS[1:nr],sep=""))
sigmaslope<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("sig.slope",LETTERS[1:nr],sep=""))

#tests si toutes les bornes
pres_b<-lapply(1:nr,function(nr,x) exists(x[nr]),x=paste("bound",LETTERS[1:nr],sep="")) 
TF<-which(pres_b==F) 
if(length(TF)==0) bornes<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("bound",LETTERS[1:nr],sep="")) else {

if(length(TF)==nr)  bornes<-rep(list(list(c(0,Inf))),nr) else { 
    bornes<-rep(list(list(c(0,Inf))),nr)
    for (i in 1:nr){ 
          if(length(which(TF==i))==0 ) {
             x=paste("bound",LETTERS[i],sep="")
             bornes[[i]]<-get(x)
          }
    }
 }
}
          
                    


#algo simplex
if (algo.option=="SIMP")  {
if (identical.times==T){
upper<-lapply(1:nr,function(nr,x) get(x),x=paste("upper",LETTERS[1],sep=""))
lower<-lapply(1:nr,function(nr,x) get(x),x=paste("lower",LETTERS[1],sep="")) }
else{
upper<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("upper",LETTERS[1:nr],sep=""))
lower<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("lower",LETTERS[1:nr],sep=""))
}
}

if (algo.option=="FW")  {
if (identical.times==T){
nwind<-lapply(1:nr,function(nr,x) get(x),x=paste("nwind",LETTERS[1],sep=""))
sampwin<-lapply(1:nr,function(nr,x) get(x),x=paste("sampwin",LETTERS[1],sep=""))
#sampwin<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("sampwin",LETTERS[1:nr],sep=""))
nsamp<-lapply(1:nr,function(nr,x) get(x),x=paste("nsamp",LETTERS[1],sep=""))
nmaxpts<-lapply(1:nr,function(nr,x) get(x),x=paste("nmaxpts",LETTERS[1],sep=""))
nminpts<-lapply(1:nr,function(nr,x) get(x),x=paste("nminpts",LETTERS[1],sep=""))
}else{
nwind<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("nwind",LETTERS[1:nr],sep=""))
sampwin<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("sampwin",LETTERS[1:nr],sep=""))
#sampwin<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("sampwin",LETTERS[1:nr],sep=""))
nsamp<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("nsamp",LETTERS[1:nr],sep=""))
nmaxpts<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("nmaxpts",LETTERS[1:nr],sep=""))
nminpts<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("nminpts",LETTERS[1:nr],sep=""))
}
}
if(NUM==F){
#recover tf in order to compute the sensibility
tf<-vector(mode="list")
for (i in 1 : length(lformg)){

  if (lformg[i]!=1){ 
    tf[[i]]<-vector()
    tf[[i]][1]<-bornes[[i]][[1]][2]
    for (j in 2:lformg[[i]]){
      
      tf[[i]][j]<-c(bornes[[i]][[j]][2]) 
    }
  }
  else {
      tf[[i]]<-bornes[[i]][[1]][2]
  }     
} 
}

#attrape si il  ya une erreur#
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

#allow to build graph 
tt1<-list()
	  for (l in 1:nr){
          tt1[[l]]<-list()
          if (length(bornes[[l]])==1){
            tt1[[l]][[1]]<-seq(bornes[[l]][[1]][1],graphsup[[l]],(graphsup[[l]]-bornes[[l]][[1]][1])/1000)
           } else {
  		      for (i in 1: (length(bornes[[l]])-1)){
  		        tt1[[l]][[i]]<-seq(bornes[[l]][[i]][1],bornes[[l]][[i]][2],(bornes[[l]][[i]][2]-bornes[[l]][[i]][1])/1000)
  		        if(i==length(bornes[[l]])-1){
  		          tt1[[l]][[i+1]]<-seq(bornes[[l]][[i]][2],graphsup[[l]],(graphsup[[l]]-bornes[[l]][[i]][2])/1000)}
             }
          }
    }
tt1b<-lapply(1:nr,function(nr,tt1) unlist(tt1[[nr]]),tt1<-tt1) 
tt1<-lapply(1:nr,function(nr,tt1b) tt1b[[nr]][tt1b[[nr]]>=graphinf[[nr]][1] & tt1b[[nr]]<=graphsup[[nr]][1]],tt1b<-tt1b) 

#test if we have NULL in the prots in this case we can not optimise with same sampling time 
ll<-lapply(prots,length)
		if(identical.times==T && length(which(unlist(lapply(1:nr,function(nr) lapply(1:ll[[nr]], function(ll,prots,nr) length(prots[[nr]][[ll]]),prots=prots,nr=nr)))==0))!=0)
		       identical.times=F

if (subjects.input==1) 
	{
	Nt<-lapply(1:nr, function(nr,subjects.init,l2) sum(subjects.init*unlist(l2[[nr]])),subjects.init=subjects.init,l2=l2)
	Ntot<-sum(unlist(Nt))
	subjects.total<-sum(subjects.init)
	subjects.initN<-subjects.init
	subjects.init<-subjects.init/subjects.total
	subjects<-subjects.init
	} else    
	
subjects<-subjects.init


if (dose.identical==T) doses<-rep(doses,length(prots[[1]])) else NULL
condinit<-NULL

#-----------------------------------------------------------------------------
#------------------------ COMPUTE --------------------------------------------

#Function used in Population Fisher information matrix  function 
#------------------------------------
pts<-function(nr,i,prots){
protk<-c()
lprotk<-c()
for (k in 1:nr){
protk <-c(protk,prots[[k]][i])  
#lprotk<-c(protk,prots[[k]][i])
}
return(protk)
}



#Sensitivity matrix :
#--------------------
sensibilite<-function(nr,protk,k,lprotk,doses) #m=modele pk ou pd proti=prelevemet pk et protiPD=prelevementPD#
{	            #protk represente touts les temps de prelevement de toute les reponses pour le k ieme protocle elementaire
  #c<-length(proti) #nombre de prelevement par protocle elementaire
	#c1<-length(protiPD)
	#lprotk<-lapply(1:nr,function(nr,protk,lprotk) length(protk[[nr]]),protk=protk )
  dose<-doses[k]
	s1<-matrix(nrow=p,ncol=0)
	mod<-numeric()
	mod1<-numeric()
	df1<-list()
	df2<-list()
  s3d<-list()
  
	for (i in 1:p) {assign(parameters[i],theta[i]) }
  
	 for (l in 1 :nr){
	    #dat<-se(lprotk,protk,dose,form,tf,l)
	    #mod1<-c(mod1,dat[[1]])
	    #s1<-cbind(s1,dat[[2]])
	    #mod<-numeric()
      s<-matrix(nrow=p,ncol=lprotk[[l]])    
      s2d<-list() 
      s4d<-list()
      
  #symbolic derivatives of model expression
  if (NUM==F){
      vlf<-1:lformg[[l]]
      forma<-formg[[l]]
      lf<-length(forma)
      #derivees premieres par reponses
      Dforma<-lapply(1:lf,function(p,forma,parameters,lf) lapply(1:p,function(p,forma,parameters) D(forma,parameters[p]), forma=forma[lf],parameters=parameters),p=p,forma=forma,parameters=parameters)
      #derivees secondes par reponses et par parametres
      lf2<-length(unlist(Dforma))
      p1<-p
      D2forma<-lapply(1:lf,function(p,Dforma,parameters,lf) lapply(1:p,function(p,Dforma,parameters) lapply(1:p1, function(p1,Dforma,parameters) D(Dforma,parameters[p1]),Dforma=Dforma[[p]],parameters=parameters),Dforma=Dforma[[lf]],parameters=parameters),p=p,Dforma=Dforma,parameters=parameters)
			#print(D2forma)    
      if (lprotk[[l]]!=0) { 
      for (i in 1:p){
      s2d[[i]]<-matrix(nrow=p,ncol=lprotk[[l]])
	    s4d[[i]]<-matrix(nrow=p,ncol=lprotk[[l]]) 
	         for (j in 1:lprotk[[l]])
		       {   
               t<-protk[[l]][j]
			         ff1<-formg[[l]][t<=tf[[l]]][1]	     
			         u<-eval(Dforma[[vlf[t<=tf[[l]]][1]]][[i]])#*dose # 
               if (Trand ==2)  u<-theta[i]*u
			         mod[j]<-eval(ff1)#*dose #
               s[i,j]<-u   
               
               for ( j1 in 1:p)
               {
			          if (Trand==2)
			          {
			          int<-theta[j1]
			          s4d[[i]][j1,j]<-eval(D2forma[[vlf[t<=tf[[l]]][1]]][[i]][[j1]])
			          s2d[[i]][j1,j]<-eval(D2forma[[vlf[t<=tf[[l]]][1]]][[i]][[j1]])*int#*dose#
			          if(j1==i) 
                {s2d[[i]][j1,j]<-s2d[[i]][j1,j]+s[i,j]/theta[i]
			          #print(s2d[[i]][j1,j])
			          }
			          }
			          else
			          {
			          s2d[[i]][j1,j]<-eval(D2forma[[vlf[t<=tf[[l]]][1]]][[i]][[j1]])
			          }
			          }  
	         }  
             
	     } 
     
	   mod1<-c(mod1,mod)  
     s1<-cbind(s1,s)
     s3d[[l]]<-s2d
     }
    }
    
    
 #numerical derivatives of model function
 dd2<-list()
  if (NUM==T){
      #model prediction
      mod<-form(protk[[l]],theta,X=dose)
      if (nr>1){ if (l==1) mod<-mod[1:lprotk[[l]]] else  mod<-mod[(lprotk[[l-1]]+1):(lprotk[[l-1]]+lprotk[[l]])]}
      
      #derivees premieres par reponses
      df1[[l]]<-jacobian(form,x=theta,t=protk[[l]],X=dose)
      if (nr>1){if (l==1) df1[[l]]<-df1[[l]][1:lprotk[[l]],] 
      else  df1[[l]]<-df1[[l]][(lprotk[[l-1]]+1):(lprotk[[l-1]]+lprotk[[l]]),]}
          
      s<-t(df1[[l]])
      if (Trand==2) {s<-s*theta}
          
      #derivees secondes par reponses et par param?tre

      s2dp<-list()
      dd2[[l]]<-list()
      for (j in 1:lprotk[[l]])
      {
     
      dd2[[l]][[j]]<-hessian(form,x=theta,t=protk[[l]][j],X=dose)
      #if (nr>1){if (l==1) df2[[l]]<-df2[[l]][1:lprotk[[l]],] 
      #else  df2[[l]]<-df2[[l]][(lprotk[[l-1]]+1):(lprotk[[l-1]]+lprotk[[l]]),]}
      
      }
      df2[[l]]<-matrix(unlist(dd2[[l]]),nrow=p*p,ncol=lprotk[[l]])
   
      #print(df2[[l]])
      
      for (i in 1:p)
        {s2dp[[i]]<-df2[[l]][(p*(i-1)+1):(i*p),1:lprotk[[l]]] 
          if (Trand==2){            
            for (j in 1:lprotk[[l]]){
              for (i2 in 1:p){
                s2dp[[i]][i2,j]<-s2dp[[i]][i2,j]*theta[i2]+s[i,j]/theta[i]
                if (i2==i) s2dp[[i]][i2,j]<-s2dp[[i]][i2,j]+s[i2,j]/theta[i2]
              }
            }
          } 
        }
  

     mod1<-c(mod1,mod)  
     s1<-cbind(s1,s)
     s3d[[l]]<-s2dp
      
}
     
}    
     
     #joindre les morceaux de chaque reponse par parametres#
     s3dp<-list()
    
    for (j4 in 1:p)
     {  
     mo1<-NULL
       
     for (j3 in 1:length(s3d))
      {
      
     mo<-NULL  
     mo<-s3d[[j3]][[j4]] 
     mo1<-cbind(mo1,mo)
    }
    s3dp[[j4]]<-mo1
     }

     l<-list(mod1,s1,s3dp)
    # print(l)
	return(l)
}

bloc_PFIM<-function(protk,sensi,var,invvar,dvar,mod,lprotk,bout4)
{
	resB<-matrix(nrow=p+nr*2,ncol=p+nr*2)     #matrice aleatoire #
	resC<-matrix(rep(0,(p+nr*2)*p),nrow=p+nr*2,ncol=p)     #matrice incomplete option 1 covariance==0#
	resA2<-matrix(ncol=p,nrow=p)
	if (Trand==2) sensia<-sensi/theta 	else sensia<-sensi 
	resA1<-(sensia %*% invvar) %*% t(sensia)  #matrice parametres fixe#
	#print(resA1)
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
      #print(resB1)
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
	h1<-cbind(rbind(resA,resC),rbind(t(resC),resB))
	return(h1) 
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


bloc_BFIM<-function(sensi,invvar,dvar,omega)  
{

#Block A only
#-------
  resA1<-matrix(nrow=p,ncol=p)
  #taking into account a priori knowlegde on random effect distribution
 	omega1<-as.matrix(as.matrix(omega[,diag(omega)!=0])[diag(omega)!=0,])
 	if (Trand==2) sensia<-sensi/theta 	else sensia<-sensi
  if (Trand==2) resA1<-((sensia %*% invvar) %*% t(sensia))[diag(omega)!=0,diag(omega)!=0] +solve(diag(beta[diag(omega)!=0])%*%omega1%*%diag(beta[diag(omega)!=0]))
  else resA1<-(sensia %*% invvar) %*% t(sensia) +solve(omega1)
 	
  resA2<-matrix(ncol=p,nrow=p)
  for (i in 1:p){
		resA3<-invvar%*%dvar[[i]]%*%invvar
		resA2[i,]<-sapply(1:p,function(dvar,resA3,p) sum(diag(resA3%*%dvar[[p]])),dvar=dvar,resA3=resA3)
  } 


 	resA<-resA1+0.5*resA2[diag(omega)!=0,diag(omega)!=0]
 	
	return(resA) 
}



#prots,subjects,condinit,doses,Ntot,ind=0,l_prote,nq,categories,names.cov,proportions.cov,beta.covariate,tab_mod_prop,pfixe,p,theta,tab.cov.param
fisher<-function(prots,subjects,condinit,doses,Ntot,ind=0,l_prote,nq,categories,names.cov,proportions.cov,beta.covariate,tab_mod_prop,pfixe,p,theta2,tab.cov.param,n_occ,lposs)
{	
  
  #Nbsubjects<-Ntot/sum((nq+np)*subjects) 
  Nbsubjects<-Ntot/sum((unlist(l_prote)*subjects)) 
  #n<-np+nq 
  n<-sum(unlist(l_prote))              #nq=nombre de prelevements par protocoles elementaires#
	if(ind==2) Nbsubjects<-1
  pp<-2*p+2*nr #2 parameters for the residual variance by response#
	if (FIM=="I"){pp<-p+2*nr}
	if (FIM=="B"){pp<-p-length(which(diag(omega)==0))}
  somme<-matrix(c(rep(0,pp*pp)),ncol=pp) 
	se<-numeric() 
	cv<-numeric() 
	
	#for(i in 1: length(prots[[1]])) 
  for(i in 1: length(l_prote))  #pour chaque protocole meme nombre de protocole elementaire en pk et pd#
	{

      protk<-pts(nr,i,prots)
	    lprotk<-lapply(1:nr,function(nr,protk) length(protk[[nr]]),protk=protk )
	    lcumprotk<-c(1,cumsum(unlist(lprotk)))
	    if (lcumprotk[2]>0){
			sensibili<-sensibilite(nr,protk,i,lprotk,doses)  #i the elementar protocol
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
			eltdvar[[n]]<-lapply(1:p, function(p,sensi,t,t1,sigmaslope,n,theta) sigmaslope[[n]]*sensi[p,t:t1]/theta[p],sensi=sensi,t=t,t1=t1,sigmaslope=sigmaslope,n=n,theta=theta)
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
            for (j6 in 1:nr)
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
      
			somme<-somme+ subjects[i]*Nbsubjects*fish 

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
      #theta<-theta[diag(omega)!=0][beta.estim]
      p<-p-length(which((beta.fixed+(diag(omega)==0))>0))
  }
  else{
  beta.estim<-!beta.fixed 
  #theta<-theta[beta.estim]
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
	previous_information<-read.table(file.path(directory,previous.FIM), header=FALSE,sep="")    #On teste si la previous information matrix a la meme taille que la matrice actuelle. Si ce n'est pas le cas, on arrete l'execution et on affiche un message d'erreur : 
    if (nrow(previous_information) != nrow(somme) | ncol(previous_information) != ncol(somme) ) {stop ("The previous FIM has not the same size as the current FIM","\n") }
    somme<-as.matrix(previous_information)+somme
    }	
	
	if (ind==2) l<-list(somme) else
    {
	   pp<-dim(as.matrix(somme))[1] 
	   
   
	   det<-det(as.matrix(somme))

	   if (det==0) crit<-1e-99	else crit<-(det)^(1/pp)	  #dans le cas ou le protocole est degenere
     invcrit<-1/crit 
     bbbbb <<- c(bbbbb, crit)
	  
     if (ind==0) l<-invcrit else
		    {
		    f<-fisher.cv(somme,pp,SIG,omega1)
		    inv<-f[[1]]
		    se<-f[[2]]
		    cv<-f[[3]]
		    l<-list(somme=somme,doses=dose,prots=prots,inv=inv,se=se,cv=cv,det=det,crit=crit,p=p,pp=pp,subjects=subjects,xn=xn,pfixe=pfixe,sh=sh)
 		     }
 	}
 	
 return(l)
  
}

fisher.cv<-function(somme,pfixe,pp,SIG,OMEGA1,beta.covariate)
{ 
if (FIM=="B") {
      beta.estim<-!beta.fixed[diag(omega)!=0] 
      theta<-theta[diag(omega)!=0][beta.estim]
      #p<-p-length(which((beta.fixed+(diag(omega)==0))>0))
      #pfixe<-pfixe-length(which((beta.fixed+(diag(omega)==0))>0))
      }
else{
      beta.estim<-!beta.fixed 
      theta<-theta[beta.estim]
      #ncovfixe<-pfixe-p
      #p<-length(beta.estim)
      #pfixe<-pfixe-length(which(beta.estim==F)) 
}
lSIG<<-length(SIG)
      	
  inv<-try(solve(as.matrix(somme)))
	if(!is.null(attributes(inv)$class))
	{
		se<-rep(NA,pp)
		cv<-se
	}		
	else
		{
		se<-sqrt(diag(inv))   
    if (n_occ>1) {      
    cv1<-se[1:pfixe]/c(theta,unlist(beta.covariate))*100 
    cv1bis<-se[(pfixe+1):(pfixe+locc)]/beta.occ1*100 
		cv2<-se[(pfixe+locc+1):(pp-lSIG)]/diag(as.matrix(OMEGA1))*100
		cv3<-se[(pp-(lSIG-1)):pp]/SIG*100 
    cv<-abs(c(cv1,cv1bis,cv2,cv3))
		}
		else{
		se<-sqrt(diag(inv))
    
		cv1<-se[1:pfixe]/theta*100 
		if (FIM=="P") {
		cv2<-se[(pfixe+1):(pp-lSIG)]/diag(as.matrix(OMEGA1))*100 
		cv3<-se[(pp-(lSIG-1)):pp]/SIG*100 
		cv<-abs(c(cv1,cv2,cv3))}
    if (FIM=="I") {
      cv3<-se[(pp-(lSIG-1)):pp]/SIG*100 
      cv<-abs(c(cv1,cv3))
    }	
    if (FIM=="B") cv<-abs(c(cv1))
    
    }
    }
	l<-list(inv=inv,se=se,cv=cv)
	return(l)
}

#function to build graph
#############################
graphAFSopti<-function(prots,protopti,graphinf,graphsup,names.datax,names.datay,y.range) 
{

		for (i in 1:p) {assign(parameters[i],theta[i]) }

		if (dose.identical==T)
		{
		    dose<-doses[1] 
		    ff1<-list()
		    ff2<-list()
		    
        par(mfrow=c(1,nr))
        for (l in 1 : nr){
            
		        ff1[[l]]<-vector(mode="numeric", length(tt1[[l]]))
		        	        
				    for (j in 1:length(tt1[[l]]))
	          {
              #ff1[[l]]<-vector(mode="numeric",length=length(tt[[l]]))
				      t<-tt1[[l]][j]
				         ff1[[l]][j]<-eval(formg[[l]][t<=tf[[l]]][1])
    	      
  		      }
  		      
  		      
  		      if(log.logical!=F) plot(tt1[[l]],ff1[[l]],type='l',xlab=paste(names.datax[l]),ylab=paste(names.datay[l]),log=log.logical,ylim=c(y.range[[l]],range(ff1[[l]]))[1:2]) 
            else   plot(tt1[[l]],ff1[[l]],type='l',xlab=paste(names.datax[l]),ylab=paste(names.datay[l]),ylim=c(y.range[[l]],range(ff1[[l]]))[1:2])
        
            #points for the design on the graph
				        
				          for ( i in 1:length(protopti[[l]]))
				          {
				          ff2[[l]]<-vector(mode="numeric", length(protopti[[l]][[i]]))
				            for ( j in 1:length(protopti[[l]][[i]])) {
				    	          t<-protopti[[l]][[i]][j]
					              ff2[[l]][j]<-eval(formg[[l]][t<=tf[[l]]][1])#*dose[1] #
		                }
				                points(protopti[[l]][[i]],ff2[[l]],pch=paste(i),cex=2,col=paste(i))
	                } 
                  
             
		      }
      }
		else
		{
	    ff1<-list()
	    par(mfrow=c(nr,length(protopti[[1]])))
			for (l in 1 : nr){

			 ff1[[l]]<-vector(mode="numeric", length(tt1[[l]]))
			     for (kjk in 1:length(protopti[[l]]))
			     {
			#windows()
			#par(mfrow=c(1,2))
			       dose<-doses[kjk]
				      for (j in 1:length(tt1[[l]]))
				      {
				      	 t<-tt1[[l]][j]
					       ff1[[l]][j]<-eval(formg[[l]][t<=tf[[l]]][1])#*dose[k]  #
                
				      }
	         
				      if (log.logical!=F) 
				          plot(tt1[[l]],ff1[[l]],type='l',xlab=paste(names.datax[l]),ylab=paste(names.datay[l]),log=log.logical,ylim=c(y.range[[l]],range(ff1[[l]]))[1:2]) else plot(tt1[[l]],ff1[[l]],type='l',xlab=paste(names.datax[l]),ylab=paste(names.datay[l]),ylim=c(y.range[[l]],range(ff1[[l]]))[1:2])
			          	title(paste("Group",kjk,"- dose = ", doses[kjk],sep=" ")) 
            
              ff2<-list()
              ff2[[l]]<-vector(mode="numeric", length(protopti[[l]][[kjk]]))
				      for ( i in 1:length(protopti[[l]][[kjk]]))
				      {
	               dose<-doses[kjk]
					       t<-protopti[[l]][[kjk]][i]
					       ff2 [[l]][i]<-eval(formg[[l]][t<=tf[[l]]][1])#*dose[k] #				
		          }
           	    points(protopti[[l]][[kjk]],ff2[[l]],pch=paste(kjk),cex=2)
			          
			}
			}
      	
		}		
}

graphAFNopti<-function(prots,protopti,graphinf,graphsup,names.datax,names.datay,y.range)
{
		for (i in 1:p) {assign(parameters[i],theta[i]) }
		tt<-list()

    tt<-lapply(1:nr, function(nr,graphinf,graphsup) seq(graphinf[[nr]],graphsup[[nr]],(graphsup[[nr]]-graphinf[[nr]])/1000),graphinf=graphinf,graphsup=graphsup)
		if (dose.identical==T)
		{
				par(mfrow=c(1,nr))
        
				ff1<-list()
        for (li in 1 : nr){
				 
          ff1[[li]]<-form(tt[[li]],theta,doses[1])
          
           if (nr>1){if (li==1) ff1[[li]]<-ff1[[li]][1:length(tt[[li]])] 
           else  ff1[[li]]<-ff1[[li]][(length(tt[[li]])+1):(length(tt[[li-1]])+length(tt[[li]]))] }
        
        }
      
        ff2<-list()
				for (l in 1:nr){
				
          if (log.logical!=F) plot(tt[[l]],ff1[[l]],type='l',xlab=paste(names.datax[l]),ylab=paste(names.datay[l]),log=log.logical,ylim=c(y.range[[l]],range(ff1[[l]]))[1:2]) else  plot(tt[[l]],ff1[[l]],type='l',xlab=paste(names.datax[l]),ylab=paste(names.datay[l]),ylim=c(y.range[[l]],range(ff1[[l]]))[1:2])
 
   	

#points for the design on the graph
				      
				          for ( i in 1:length(protopti[[l]]))
				          {
				          ff2[[l]]<-vector(mode="numeric", length(protopti[[l]][[i]]))
				          ff2[[l]]<-form(protopti[[l]][[i]],theta,doses[1])
				          if (nr>1){if (l==1) ff2[[l]]<-ff2[[l]][1:length(protopti[[l]][[i]])] 
                  else  ff2[[l]]<-ff2[[l]][(length(protopti[[l]][[i]])+1):(length(protopti[[l-1]][[i]])+length(protopti[[l]][[i]]))] }
                  points(protopti[[l]][[i]],ff2[[l]],pch=paste(i),cex=2,col=paste(i))
	                } 
                  
        }

		}
		else     
		{ 
			par(mfrow=c(nr,length(protopti[[1]])))
      ff1<-list()
			ff2<-list()	
			for (l in 1:nr){
			 for (ki in 1:length(prots[[l]]))
			 {
				ff1<-list()
       
           ff1[[l]]<-form(tt[[l]],theta,doses[ki])
          
           if (nr>1){if (l==1) ff1[[l]]<-ff1[[l]][1:length(tt[[l]])] 
           else  ff1[[l]]<-ff1[[l]][(length(tt[[l]])+1):(length(tt[[l-1]])+length(tt[[l]]))] }
        
        if (log.logical!=F) plot(tt[[l]],ff1[[l]],type='l',xlab=paste(names.datax[l]),ylab=paste(names.datay[l]),log=log.logical,ylim=c(y.range[[l]],range(ff1[[l]]))[1:2]) else  plot(tt[[l]],ff1[[l]],type='l',xlab=paste(names.datax[l]),ylab=paste(names.datay[l]),ylim=c(y.range[[l]],range(ff1[[l]]))[1:2])
				 	title(paste("Group",ki,"- dose = ", doses[ki],sep=" ")) 
         if (length(prots[[l]][[ki]]) !=0){
        
        
        #design points on the curve
        ff2[[l]]<-form(protopti[[l]][[ki]],theta,doses[ki])
				          if (nr>1){if (l==1) ff2[[l]]<-ff2[[l]][1:length(protopti[[l]][[ki]])] 
                  else  ff2[[l]]<-ff2[[l]][(length(protopti[[l]][[ki]])+1):(length(protopti[[l-1]][[ki]])+length(protopti[[l]][[ki]]))]}
                  points(protopti[[l]][[ki]],ff2[[l]],pch=paste(ki),cex=2)
	       
         }
        
			 }
		  }
      } 
}






