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
  if (NUM==F){
lf<-length(form)#number of analytical form
#number of forms in each response
#recover of the model
#for multiple responses
  formg<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("form",LETTERS[1:nr],sep=""))  
  formg_init<-formg
  lformg<-lapply(formg,length)  #length(of each model tf or not)
  c<-c(0,cumsum(lformg))
}
 
  doses<-dose
  p<-length(theta)

#for IOV and covariates
#-------------------------------------------------------------------------------
#if (FIM=="I"||FIM=="B"){n_occ=1;covariate.model<-F}
if (n_occ<=1) 
{
  n_occ=1
  covariate_occ.model<-F
}

if(n_occ>1 && FIM!="P") {stop("Individual and bayesian matrices are only computed without inter-occasion variability (n_occ set to 1)")}

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
          for (i in which(covariate_occ.sequence[[kj]][[l]]!=covariate_occ.category[[kj]][1])){
            num[kj]=1
            param.asso.occ[[kj]][[l]][[i]]=parameter_occ.associated[[kj]]
            
            for (j  in 1:length(param.asso.occ[[kj]][[l]][[i]])){
              
              occ.effect[[kj]][[l]][[i]][which(param.asso.occ[[kj]][[l]][[i]][j]==parameters)]=1
              beta.occ[[kj]][[l]][[i]][which(param.asso.occ[[kj]][[l]][[i]][j]==parameters)]=beta.covariate_occ[[kj]][[j]][num[kj]]
   
            } 
            num[kj]=num[kj]+1
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
omega<-as.matrix(omega)
gamma<-as.matrix(gamma)
OMEGA<-diag(c(diag(omega),rep(c(diag(gamma)),pn))) #la matrice OMEGA qui regroupe omega et gamma
if (length(c(diag(omega),rep(c(diag(gamma)),pn)))==1) OMEGA<-as.matrix((c(diag(omega),rep(c(diag(gamma)),pn))))

#----------- IOV and covariates END-------------------------------------------



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

#recover tf in order to compute the sensibility
if(NUM==F){
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
if (length(which(unlist(erro)==1))>=1 || length(which(is.null(unlist(erro))==T))>=1)
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
for (kit in 1:nr){
protk <-c(protk,prots[[kit]][i])  
#lprotk<-c(protk,prots[[k]][i])
}
return(protk)
}


#Sensitivity matrix arguments:
#nr: number of responses
#protk: elementary design
#pfixe: number of fixed effects
#p: length of the vector beta in the stdin.r file
#kjk: number of elementray design
#lprotk: number of samples in the kith elementray design
#doses
#arguments if covariates
#categories.cov: categories for each covariate
#names.cov: name of each covariate
#beta.covaraite:
#theta: value of fixed effects in the stind.r file
#tab.cov.param
#--------------------
sensibilite<-function(nr,protk,pfixe,p,kjk,lprotk,doses,categories.cov,names.cov,beta.covariate,thetah,tab.cov.param,beta.occh)
{	
	#c<-length(proti)
  dose<-doses[kjk] 
  theta_init<-thetah
  s1<-matrix(nrow=p,ncol=0)
  s4<-matrix(nrow=p,ncol=0)
  #s1rand<-matrix(nrow=prand,ncol=0)
	
	if (NUM==T)
	{
	mod<-list()
	mod1<-list()
	dd<-list()
	df1<-list()
	}
	else{
  mod<-numeric()
	mod1<-numeric()
 }
	
	
	#IF COVARIATE  
	if (!is.null(categories.cov))
	{
        	parameters.cov.val<-rep(1,length(thetah))
        	n.theta<-c(theta_init,as.numeric(unlist(tab.cov.param[,5]))) 
          sum<-0
          for (i in 1:pfixe) {assign(c(parameters,tab.cov.param[,3])[i],n.theta[i])}
        	#attibuer les valeurs des modalites aux covariables
        
        	for (i in 1:length(names.cov)) {assign(names.cov[[i]],categories.cov[i])} #assign("trt",groupe)
     
          #boucle sur le nombre de parametres ayant une cov associee
          for (i in 1: length(unique(tab.cov.param[,1])))
        	{
          	   
            	 par.test<-unique(tab.cov.param[,1])[i]
            	 tab.cov.param1<-tab.cov.param[which(tab.cov.param[,1]==par.test),]
            	 tab.cov.param1<-matrix(tab.cov.param1,ncol=5)
             
             if (length(unique(tab.cov.param1[,2]))==1)
             { 
              # cas avec d'un parametre associe a une covariable  a une ou plusiers modalites
              # ex (Trand==2): Cl*exp(beta_CL*COV)
              #(Trand==1):Cl + beta_CL*COV
                test.bet<-tab.cov.param1[get(tab.cov.param1[1,2])==tab.cov.param1[,4],3]
                
                if (length(test.bet)!=0)  #si categories ==0
                {
                  if (Trand==2)   
                  {
                     
                      thetah[which(tab.cov.param1[1,1]==parameters)]<-thetah[which(tab.cov.param1[1,1]==parameters)]*exp(get(tab.cov.param1[get(tab.cov.param1[1,2])==tab.cov.param1[,4],3]))
                      parameters.cov.val[which(tab.cov.param1[1,1]==parameters)]<-exp(get(tab.cov.param1[get(tab.cov.param1[1,2])==tab.cov.param1[,4],3]))
                      
                  } else 
                  {
                  thetah[which(tab.cov.param1[1,1]==parameters)]<-thetah[which(tab.cov.param1[1,1]==parameters)]+ get(tab.cov.param1[get(tab.cov.param1[1,2])==tab.cov.param1[,4],3])  

                  }
                 }
               } else{  #cas d'un parametre associe a plusieurs covariables
                
                       p_val<-1
                       som<-0
                       for (i.par in 1:length(unique(tab.cov.param1[,2]))) #boucle sur les covaraibles pour le parametres
                       {
                             
                             par.test1<-unique(tab.cov.param1[,2])[i.par]
                             tab.cov.param2<-tab.cov.param1[which(tab.cov.param1[,2]==par.test1),]
                             tab.cov.param2<-matrix(tab.cov.param2,ncol=5)
                             
                             test.bet<-tab.cov.param2[get(tab.cov.param2[1,2])==tab.cov.param2[,4],3]
                           
                             if (length(test.bet)!=0)  #si categories ==0
                            {
                               if (Trand==2)   
                              {
                                 
                                  thetah[which(tab.cov.param2[1,1]==parameters)]<-thetah[which(tab.cov.param2[1,1]==parameters)]*exp(get(tab.cov.param2[get(tab.cov.param2[1,2])==tab.cov.param2[,4],3]))
                                  parameters.cov.val[which(tab.cov.param2[1,1]==parameters)]<-exp(get(tab.cov.param2[get(tab.cov.param2[1,2])==tab.cov.param2[,4],3]))
                                  
                              } else {
                              
                                  thetah[which(tab.cov.param[i,1]==parameters)]<-thetah[which(tab.cov.param[i,1]==parameters)]+ get(tab.cov.param2[get(tab.cov.param2[1,2])==tab.cov.param2[,4],3])  

                              }
                            }
                        }
                   }

                }
                    
            }
        		
  
	for (ip in 1:p) {assign(parameters[ip],thetah[ip])} 
	#for (i in 1:p) {assign(parameters[i],theta[i])}

  
  for (l in 1:nr){
  
  #symbolic derivatives of model expression
  if (NUM==F){
	   s<-matrix(nrow=p,ncol=lprotk[[l]]) 
	   s3<-matrix(nrow=p,ncol=lprotk[[l]])
     
     mod<-numeric() 
      for (i in 1:p){
	         for (j in 1:lprotk[[l]])
		       {    
                

                vlf<-1:lformg[[l]]
                forma<-formg[[l]]
                lf<-length(forma)
                Dforma<-lapply(1:lf,function(p,forma,parameters,lf) lapply(1:p,function(p,forma,parameters) D(forma,parameters[p]), forma=forma[lf],parameters=parameters),p=p,forma=forma,parameters=parameters)
              
               
                t<-protk[[l]][j]
                ff1<-formg[[l]][t<=tf[[l]]][1]
			          u<-eval(Dforma[[vlf[t<=tf[[l]]][1]]][[i]])
                
                             		           	          
			          
                #IF COVARIATE  only
			          if(!is.null(categories.cov)) s3[i,j]<-u*parameters.cov.val[i]
                
                if (Trand==2)  u<-beta[i]*u*exp(beta.occh[i])
		            mod[j]<-eval(ff1)
		            #IF COVARIATE  only
                if(!is.null(categories.cov)) s[i,j]<-u*parameters.cov.val[i] else s[i,j]<-u
           
		        }
		     
	     }                                       
	   MOD<-mod1<-c(mod1,mod)
	   s1<-cbind(s1,s)  
	   
	   #IF COVARIATE ONLY
     if(!is.null(categories.cov)) s4<-cbind(s4,s3) # pas de multiplication par valeur du parametre Trand==1
    
 }

#numerical derivatives of model function
 if (NUM==T)

 {
  if (lprotk[[l]]==0){s<-matrix(nrow=p,ncol=0); s3<-matrix(nrow=p,ncol=0) } else {	
    	 
          s<-matrix(nrow=p,ncol=lprotk[[l]]) 
    	    s3<-matrix(nrow=p,ncol=lprotk[[l]]) 
	        	        
          mod[[l]]<-form(protk[[l]],thetah,X=dose)
      
          if (nr>1){ if (l==1) mod[[l]]<-mod[[l]][1:lprotk[[l]]] else  mod[[l]]<-mod[[l]][(lprotk[[l-1]]+1):(lprotk[[l-1]]+lprotk[[l]])]}
         
                         
          df1[[l]]<-jacobian(form,x=thetah,t=protk[[l]],X=dose)

          if (nr>1){if (l==1) df1[[l]]<-df1[[l]][1:lprotk[[l]],] 
          else  df1[[l]]<-df1[[l]][(lprotk[[l-1]]+1):(lprotk[[l-1]]+lprotk[[l]]),]}
          
          s<-df1[[l]]
          s3<-df1[[l]] 

          
          if(lprotk[[l]]==1){s<-matrix(s,nrow=pfixe,ncol=1) ; s3<-s }  else {s<-t(s); s3<-t(s3)}
          if(!is.null(categories.cov))
          {
            s<-as.matrix(s)[1:length(theta),]*parameters.cov.val
            s3<-as.matrix(s3)[1:length(theta),]*parameters.cov.val
          }
          if (Trand==2) {
          s<-as.matrix(s)[1:length(theta),]*beta*exp(beta.occh)
          }
          s1<-cbind(s1,s)
          #print(s1)
          if(!is.null(categories.cov)) s4<-cbind(s4,s3)
          }
      MOD<-unlist(mod)
 }
 
     
  }

  
  #IF COVARIATE ONLY
  if(!is.null(categories.cov)) {
  #covariate parameters
     for (icat in 1: dim(tab.cov.param)[1])
     {
      
      if (!tab.cov.param[icat,4]==get(tab.cov.param[icat,2])) s2<-s1[which(parameters==tab.cov.param[[icat]]),]*0  else   s2<-s1[which(parameters==tab.cov.param[[icat]]),]
       
       
       s1<-rbind(s1,s2)
       s4<-rbind(s4,s2)
     }
      s1rand<-as.matrix(s1[1:length(theta),])
      l1<-list(MOD,s1,s1rand,n.theta,s4)
      thetah<-theta_init
  }else{  
    # (no covariate)                          
    if (dim(s1)[2]==0) l1<-list(MOD,s1,s1,theta_init,s1)
    if (dim(s1)[2]!=0) l1<-list(MOD,s1,s1,theta_init,s1/beta)
  }
  
 
	return(l1)

}


#Deriver par rapport aux effets fixes
sensifixe<-function(nr,protk,pfixe,p,kjk,lprotk,doses,categories.cov,names.cov,beta.covariate,theta2,tab.cov.param,n_occ,sequ)
{
  mod<-numeric()
  c<-sum(unlist(lprotk))
  s<-matrix(nrow=pfixe+locc,ncol=0)
  
  for(occ in 1:n_occ){
    
    sensibili<-sensibilite(nr,protk,pfixe,p,kjk,lprotk,doses,categories.cov,names.cov,beta.covariate,theta2[[sequ]][[occ]],tab.cov.param,beta.occ2[[sequ]][[occ]])
    if (Trand ==2) 
    {
      sensibilif<-sensibili[[5]]
      }else 
    { 
      sensibilif<-sensibili[[2]] 
    }
   
    
    mod_int<-sensibili[[1]]
    mod<-c(mod,mod_int)
    
    bloc0<-matrix(rep(0,locc*c),nrow=locc,ncol=c)
    
    if (lP[[sequ]][occ]==0) s_int<-rbind(sensibilif,bloc0)
    
    if (lP[[sequ]][occ]>0){ 
      block<-matrix(nrow=0,ncol=c)
      for (kj in 1:nbcov_occ){
        Seq<-comb.seq.poss[sequ,kj]
        l<-which(covariate_occ.sequence[[kj]][[Seq]][occ]==covariate_occ.category[[kj]])-1 
        block_int<-matrix(0,nrow=length(parameter_occ.associated[[kj]])*(length(covariate_occ.category[[kj]])-1),ncol=c)
        if (length(which(vecP1[[kj]][[Seq]][[occ]]==T))>0) 
        {
          for (j in 1:length(which(vecP1[[kj]][[Seq]][[occ]]==T))){
          i=l+(j-1)*(length(covariate_occ.category[[kj]])-1)
          if (Trand==2) s_int_ij<-as.matrix(sensibili[[2]][1:p,])[which(vecP1[[kj]][[Seq]][[occ]]==T)[j],]
          else  s_int_ij<-as.matrix(sensibili[[5]][1:p,])[which(vecP1[[kj]][[Seq]][[occ]]==T)[j],]
          block_int[i,]<-s_int_ij
        }
      }
      block<-rbind(block,block_int)
       if (Trand ==2 && covariate_occ.model==T) 
    {
     sensibilif[which(vecP1[[kj]][[Seq]][[occ]]==T),]<-sensibilif[which(vecP1[[kj]][[Seq]][[occ]]==T),]*exp(beta.occ2[[sequ]][[occ]][which(vecP1[[kj]][[Seq]][[occ]]==T)])
    }
      
    }
    s_int<-rbind(sensibilif,block)
   }                        
   s<-cbind(s,s_int)
  }
  l<-list(mod,s)
  return(l)
}


#Deriver par rapport aux effets aleatoires
sensialea<-function(nr,protk,pfixe,p,kjk,lprotk,doses,categories.cov,names.cov,beta.covariate,theta2,tab.cov.param,n_occ,sequ)
{
  c<-sum(unlist(lprotk))
  mod<-numeric()
  s1<-matrix(nrow=p,ncol=0)
  s2<-matrix(0,nrow=p*pn,ncol=c*pn)
  s<-matrix()
   
  if (n_occ==1) s<-sensibilite(nr,protk,pfixe,p,kjk,lprotk,doses,categories.cov,names.cov,beta.covariate,theta2[[sequ]][[1]],tab.cov.param,beta.occ2[[sequ]][[1]])[[3]]
  else{ 
    for(occ in 1:pn)
    { sensibili<-sensibilite(nr,protk,pfixe,p,kjk,lprotk,doses,categories.cov,names.cov,beta.covariate,theta2[[sequ]][[occ]],tab.cov.param,beta.occ2[[sequ]][[occ]])
      mod_int<-sensibili[[1]]
      s_int<-sensibili[[3]]
    
      mod<-c(mod,mod_int)
      s1<-cbind(s1,s_int)
      s2[(p*(occ-1)+1) : (p*occ),(c*(occ-1)+1) :(c*occ)]=s_int
    }
    s<-rbind(s1,s2)
  }
  l<-list(mod,s)
  return(l)
}

# Fisher matrix for 1 ind in 1 elementary design by block: block A (resA), block B(resB), block C (resC)
# matrix from a list of diagonal matrices
#-------------------------------------------------------------------------------------------------------

#fonction pour ecrire une matrice block diagonal a partir d'une liste de matrices
blockdiag <-function(x){
     if(!is.list(x)) stop("x not a list")
     n <- length(x)
     if(n==0) return(NULL)
     x <- lapply(x, function(y) if(length(y)) as.matrix(y) else
stop("Zero-length component in x"))
     d <- array(unlist(lapply(x, dim)), c(2, n))
     rr <- d[1,]
     cc <- d[2,]
     rsum <- sum(rr)
     csum <- sum(cc)
     out <- array(0, c(rsum, csum))
     ind <- array(0, c(4, n))
     rcum <- cumsum(rr)
     ccum <- cumsum(cc)
     ind[1,-1] <- rcum[-n]
     ind[2,] <- rcum
     ind[3,-1] <- ccum[-n]
     ind[4,] <- ccum
     imat <- array(1:(rsum * csum), c(rsum, csum))
     iuse <- apply(ind, 2, function(y, imat) imat[(y[1]+1):y[2],
(y[3]+1):y[4]], imat=imat)
     iuse <- as.vector(unlist(iuse))
     out[iuse] <- unlist(x)
     return(out)
} 

#population matrix
#-------------------------------------------------------------------------------------------------------
bloc_PFIM<-function(protk,sensif,sensia,var,invvar,pfixe,p,dvar,mod,lprotk,bout4,categories.cov,sequ)   
{

#Block A
#-------
  resA<-matrix(nrow=pfixe+locc,ncol=pfixe+locc)
  resA<-(sensif %*% invvar) %*% t(sensif)  #matrice parametres fixe#


#Block B
#-------
  n_occB<-n_occ
  if (length(which(gamma==0))==length(gamma)) n_occB<-1
  if (n_occB>1) l<-2*p else l<-p
	resB<-matrix(nrow=l+nr*2,ncol=l+nr*2)     #matrice aleatoire (block B) #
	resB_int<-matrix(nrow=l+nr*2,ncol=l+nr*2)     #calcul intermediaire
	
	l1<-list()
	li<-list()
	l2<-list()
	
	#premiere reponse
	c<-list()
	long<-c()
	a<-NULL
	for (occ in 1:n_occB){
	c[[occ]]<-cumsum(unlist(lprotk))+length(unlist(protk))*(occ-1)
	long[occ]<-c[[occ]][nr]-c[[occ]][1]
	a0<-c(bout4[[occ]][[1]],rep(0,long[occ]))
	a<-c(a,a0)
	}
	a1<-2*diag(a,length(unlist(protk))*n_occB)
	ai<-invvar%*%a1%*%invvar
	resB_int[l+1,l+1]<-sum(diag(ai%*%a1))

	a2<-a1*mod
	resB_int[l+2,l+2]<-sum(diag(a2%*%invvar%*%a2%*%invvar))
	resB_int[l+1,l+2]<-sum(diag(ai%*%a2))
	resB_int[l+2,l+1]<-resB_int[l+1,l+2]
	l1[[1]]<-a1
	li[[1]]<-ai
	l2[[1]]<-a2
	
	
	if (nr!=1){
	for (i in 1:(nr-1)){
	vec<-NULL
	 for (occ in 1:n_occB){
	   vec0<-rep(0,length(unlist(protk)))
	   if(length(protk[[i+1]])!=0){
	   
	     vec0<-c(rep(0,length(protk[[i]])),bout4[[occ]][[i+1]])
	   } else {vec0<-vec0}
	 
   vec<-c(vec,vec0)
   }
   
	a1<-2*diag(vec,length(unlist(protk))*n_occB)

	ai<-invvar%*%a1%*%invvar
	resB_int[l+i*2+1,l+i*2+1]<-sum(diag(ai%*%a1))
	#sigma slope #
	a2<-a1*mod
  resB_int[l+i*2+2,l+i*2+2]<-sum(diag(a2%*%invvar%*%a2%*%invvar))
	resB_int[l+(i*2)+1,l+i*2+2]<-sum(diag(ai%*%a2))
	resB_int[l+i*2+2,l+(i*2)+1]<-resB_int[l+(i*2)+1,l+i*2+2]
	l1[[i+1]]<-a1
	li[[i+1]]<-ai
	l2[[i+1]]<-a2
  }
 
 
  #sigma inter PD et sigma inter PK#
  for (i in 1: (nr)){
    for (j in  2 : (nr)){

	     resB_int[(l+j*2-1),l+i*2-1]<-sum(diag(li[[j]]%*%l1[[i]])) 
       resB_int[l+i*2-1,(l+j*2-1)]<-resB_int[(l+j*2-1),l+i*2-1] #resB[p+1,p+3]<-resB[p+3,p+1]
	     resB_int[l+j*2,l+i*2-1]<-sum(diag(li[[i]]%*%l2[[j]])) 	#resB[p+1,p+4]<-resB[p+4,p+1]
       resB_int[l+i*2-1,l+j*2]<-resB_int[l+j*2,l+i*2-1]
	    
       resB_int[l+j*2-1,l+i*2]<-sum(diag(li[[j]]%*%l2[[i]]))  #resB[p+2,p+3]<-resB[p+3,p+2]
	     resB_int[l+i*2,l+j*2-1]<-resB_int[l+j*2-1,l+i*2]
       resB_int[l+j*2,l+i*2]<-sum(diag(invvar%*%l2[[j]]%*%invvar%*%l2[[i]]))
       resB_int[l+i*2,l+j*2]<-resB_int[l+j*2,l+i*2] #resB[p+2,p+4]<-resB[p+4,p+2]
	}
	}
	
	for (i in 1:p)
	{
    resB1<-(sensia[i,])%*%t(sensia[i,])
    resB2<-invvar%*%resB1%*%invvar
    resB_int[i,1:p]<-sapply(1:p,function(resB2,sensia,p) sum(diag(resB2%*%(sensia[p,])%*%t(sensia[p,]))),resB2=resB2,sensia=sensia)
    resB_int[i,l+1]<-sum(diag(l1[[1]]%*%resB2))
		resB_int[i,l+2]<-sum(diag(l2[[1]]%*%resB2))
    for (j in 2:nr) {
		resB_int[i,l+j*2-1]<-sum(diag(l1[[j]]%*%resB2))
		resB_int[i,l+j*2]<-sum(diag(l2[[j]]%*%resB2)) 	
    }
		resB_int[(l+1):(l+(2*nr)),i]<-resB_int[i,(l+1):(l+(2*nr))]
		
    if (n_occB>1 && length(which(gamma==0))!=length(gamma))
    {
    resB3<-blockdiag(sapply(1:n_occB,function(mat,n_occB) list(((sensia[i,])%*%t(sensia[i,]))[1:length(unlist(protk)),1:length(unlist(protk))]),mat=((sensia[i,])%*%t(sensia[i,]))[1:length(unlist(protk)),1:length(unlist(protk))]))
    resB4<-invvar%*%resB3%*%invvar
    resB_int[p+i,1:p]<-sapply(1:p,function(resB4,sensia,p) sum(diag(resB4%*%(sensia[p,])%*%t(sensia[p,]))),resB4=resB4,sensia=sensia)
    resB_int[1:p,p+i]<-resB_int[p+i,1:p]
    resB_int[p+i,(p+1):l]<-sapply(1:p,function(resB4,sensia,p) sum(diag(resB4%*%blockdiag(sapply(1:n_occB,function(mat,n_occB) list(((sensia[p,])%*%t(sensia[p,]))[1:length(unlist(protk)),1:length(unlist(protk))]),mat=((sensia[p,])%*%t(sensia[p,]))[1:length(unlist(protk)),1:length(unlist(protk))])))),resB4=resB4,sensia=sensia)
		resB_int[p+i,l+1]<-sum(diag(l1[[1]]%*%resB4))
		resB_int[p+i,l+2]<-sum(diag(l2[[1]]%*%resB4))
    for (j in 2:nr) {
		resB_int[p+i,l+j*2-1]<-sum(diag(l1[[j]]%*%resB4))
		resB_int[p+i,l+j*2]<-sum(diag(l2[[j]]%*%resB4)) 	
    }
		resB_int[(l+1):(l+(2*nr)),p+i]<-resB_int[p+i,(l+1):(l+(2*nr))]
		}
	}
# Deduire B a partir de B_intermediare
 resB<-resB_int 	
 }
 

  
  else {
	 for (i in 1:p)
	{
    resB1<-(sensia[i,])%*%t(sensia[i,])
    resB2<-invvar%*%resB1%*%invvar
    resB_int[i,1:p]<-sapply(1:p,function(resB2,sensia,p) sum(diag(resB2%*%(sensia[p,])%*%t(sensia[p,]))),resB2=resB2,sensia=sensia)
    resB_int[i,l+1]<-sum(diag(l1[[1]]%*%resB2))
		resB_int[i,l+2]<-sum(diag(l2[[1]]%*%resB2))
		resB_int[(l+1):(l+2),i]<-resB_int[i,(l+1):(l+2)]
		
    if (n_occB>1 && length(which(gamma==0))!=length(gamma))
    {
    resB3<-blockdiag(sapply(1:n_occB,function(mat,n_occB) list(((sensia[i,])%*%t(sensia[i,]))[1:length(unlist(protk)),1:length(unlist(protk))]),mat=((sensia[i,])%*%t(sensia[i,]))[1:length(unlist(protk)),1:length(unlist(protk))]))
    resB4<-invvar%*%resB3%*%invvar
    resB_int[p+i,1:p]<-sapply(1:p,function(resB4,sensia,p) sum(diag(resB4%*%(sensia[p,])%*%t(sensia[p,]))),resB4=resB4,sensia=sensia)
    resB_int[1:p,p+i]<-resB_int[p+i,1:p]
    resB_int[p+i,(p+1):l]<-sapply(1:p,function(resB4,sensia,p) sum(diag(resB4%*%blockdiag(sapply(1:n_occB,function(mat,n_occB) list(((sensia[p,])%*%t(sensia[p,]))[1:length(unlist(protk)),1:length(unlist(protk))]),mat=((sensia[p,])%*%t(sensia[p,]))[1:length(unlist(protk)),1:length(unlist(protk))])))),resB4=resB4,sensia=sensia)
    resB_int[p+i,l+1]<-sum(diag(l1[[1]]%*%resB4))
		resB_int[p+i,l+2]<-sum(diag(l2[[1]]%*%resB4))
		resB_int[(l+1):(l+2),p+i]<-resB_int[p+i,(l+1):(l+2)]
		}
	}
  # Deduire B a partir de B_intermediare
 resB<-resB_int  	   
}

  
  # B final
  resB<-0.5*resB                               
  # Bloc C 
  resC<-matrix(rep(0,(l+nr*2)*(pfixe+locc)),nrow=l+nr*2,ncol=pfixe+locc)     #matrice incomplete option 1 covariance==0#
	# ou pfixe est 1 dim de A	
    
  h<-cbind(rbind(resA,resC),rbind(t(resC),resB))
	return(h) 
}



#individual matrix
#-------------------------------------------------------------------------------------------------------
bloc_IFIM<-function(protk,sensif,sensia,var,invvar,pfixe,p,dvar,mod,lprotk,bout4,categories.cov,sequ)   #proti,protiPD,sensi,var,invvar,dvar,mod
{

#Block A
#-------
  resA<-matrix(nrow=pfixe+locc,ncol=pfixe+locc)
  resA<-(sensif %*% invvar) %*% t(sensif)  #matrice parametres fixe#

#Block B
#-------
  n_occB<-n_occ
  if (length(which(gamma==0))==length(gamma)) n_occB<-1
  if (n_occB>1) l<-2*p else l<-p
  #no random effect
  l<-0
	resB<-matrix(nrow=l+nr*2,ncol=l+nr*2)     #matrice aleatoire (block B) #
	resB_int<-matrix(nrow=l+nr*2,ncol=l+nr*2)     #calcul intermediaire
	
	l1<-list()
	li<-list()
	l2<-list()
	
	#premiere reponse
	c<-list()
	long<-c()
	a<-NULL
	for (occ in 1:n_occB){
	c[[occ]]<-cumsum(unlist(lprotk))+length(unlist(protk))*(occ-1)
	long[occ]<-c[[occ]][nr]-c[[occ]][1]
	a0<-c(bout4[[occ]][[1]],rep(0,long[occ]))
	a<-c(a,a0)
	}
	a1<-2*diag(a,length(unlist(protk))*n_occB)
	ai<-invvar%*%a1%*%invvar
	resB_int[l+1,l+1]<-sum(diag(ai%*%a1))

	a2<-a1*mod
	resB_int[l+2,l+2]<-sum(diag(a2%*%invvar%*%a2%*%invvar))
	resB_int[l+1,l+2]<-sum(diag(ai%*%a2))
	resB_int[l+2,l+1]<-resB_int[l+1,l+2]
	l1[[1]]<-a1
	li[[1]]<-ai
	l2[[1]]<-a2
	
	
	if (nr!=1){
	for (i in 1:(nr-1)){
	vec<-NULL
	 for (occ in 1:n_occB){
	   vec0<-rep(0,length(unlist(protk)))
	   if(length(protk[[i+1]])!=0){
	   
	     vec0<-c(rep(0,length(protk[[i]])),bout4[[occ]][[i+1]])
	   } else {vec0<-vec0}
	 
   vec<-c(vec,vec0)
   }
   
	a1<-2*diag(vec,length(unlist(protk))*n_occB)

	ai<-invvar%*%a1%*%invvar
	resB_int[l+i*2+1,l+i*2+1]<-sum(diag(ai%*%a1))
	#sigma slope #
	a2<-a1*mod
  resB_int[l+i*2+2,l+i*2+2]<-sum(diag(a2%*%invvar%*%a2%*%invvar))
	resB_int[l+(i*2)+1,l+i*2+2]<-sum(diag(ai%*%a2))
	resB_int[l+i*2+2,l+(i*2)+1]<-resB_int[l+(i*2)+1,l+i*2+2]
	l1[[i+1]]<-a1
	li[[i+1]]<-ai
	l2[[i+1]]<-a2
  }
  #print(resB_int)
 
  #sigma inter PD et sigma inter PK#
  for (i in 1: (nr)){
    for (j in  2 : (nr)){

	     resB_int[(l+j*2-1),l+i*2-1]<-sum(diag(li[[j]]%*%l1[[i]])) 
       resB_int[l+i*2-1,(l+j*2-1)]<-resB_int[(l+j*2-1),l+i*2-1] #resB[p+1,p+3]<-resB[p+3,p+1]
	     resB_int[l+j*2,l+i*2-1]<-sum(diag(li[[i]]%*%l2[[j]])) 	#resB[p+1,p+4]<-resB[p+4,p+1]
       resB_int[l+i*2-1,l+j*2]<-resB_int[l+j*2,l+i*2-1]
	    
       resB_int[l+j*2-1,l+i*2]<-sum(diag(li[[j]]%*%l2[[i]]))  #resB[p+2,p+3]<-resB[p+3,p+2]
	     resB_int[l+i*2,l+j*2-1]<-resB_int[l+j*2-1,l+i*2]
       resB_int[l+j*2,l+i*2]<-sum(diag(invvar%*%l2[[j]]%*%invvar%*%l2[[i]]))
       resB_int[l+i*2,l+j*2]<-resB_int[l+j*2,l+i*2] #resB[p+2,p+4]<-resB[p+4,p+2]
	}
	}
	


 }
 
 # Deduire B a partir de B_intermediare
 resB<-resB_int 	
# B final
  resB<-0.5*resB                               
# Bloc C 
  resC<-matrix(rep(0,(l+nr*2)*(pfixe+locc)),nrow=l+nr*2,ncol=pfixe+locc)     #matrice incomplete option 1 covariance==0#
	# ou pfixe est 1 dim de A	
    
  h<-cbind(rbind(resA,resC),rbind(t(resC),resB))
  
	return(h) 
}


#bayesian matrix
#-------------------------------------------------------------------------------------------------------
bloc_BFIM<-function(sensif,invvar,pfixe,p,omega)   
{

#Block A only
#-------
  resA<-matrix(nrow=pfixe+locc,ncol=pfixe+locc)
  #taking into account a priori knowlegde on random effect distribution
 	omega1<-as.matrix(as.matrix(omega[,diag(omega)!=0])[diag(omega)!=0,])
  if (Trand==2) resA<-((sensif %*% invvar) %*% t(sensif))[diag(omega)!=0,diag(omega)!=0] +solve(as.matrix(diag((beta[diag(omega)!=0]),length(beta[diag(omega)!=0])))%*%as.matrix(omega1)%*%as.matrix(diag(((beta[diag(omega)!=0])),length(beta[diag(omega)!=0]))))
  else resA<-((sensif %*% invvar) %*% t(sensif))[diag(omega)!=0,diag(omega)!=0] +solve(omega1)
  
	return(resA) 
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

#Population Fisher information matrix
#------------------------------------
fisher<-function(prots,subjects,condinit,doses,Ntot,ind=0,l_prote,nq,categories,names.cov,proportions.cov,beta.covariate,tab_mod_prop,pfixe,p,theta2,tab.cov.param,n_occ,lposs)
{	 
  Nbsubjects<-Ntot/sum((unlist(l_prote)*subjects)) 
  n<-sum(unlist(l_prote))              #nq=nombre de prelevements par protocoles elementaires#
	if(ind==2) Nbsubjects<-1              
	 n_occB<-n_occ
  if (length(which(gamma==0))==length(gamma)) n_occB<-1
	if (n_occB>1) pp<-pfixe+length(beta.occ1)+p*2+2*nr  else pp<-pfixe+length(beta.occ1)+p+2*nr
	if (FIM=="I"){pp<-pfixe+2*nr}
	if (FIM=="B"){pp<-pfixe-length(which(diag(omega)==0))}

	somme<-matrix(c(rep(0,pp*pp)),ncol=pp) 
	   #test if covariate or not
	   #categories <-covariate.cat.poss
     #names.cov<-info.covariate[[1]]
  no.cov<-1
  
    if (is.null(categories)) 
        {
        categories<-matrix(0, nrow=1,ncol=1)
		    no.cov<-NULL
		    }
	#For each elementary design computation of the elementary fisher information matrix
for (sequ in 1:lposs){
    #somme1<-matrix(c(rep(0,pp*pp)),ncol=pp) 
	 for  (i in 1:length(l_prote))	#meme nombre de temps pour chaque reponse par groupe elementaire
		{  
		    
		    #categories:table with for each line a combination of the modalities for each covariates
		    #loop on all possible combination of modalities of all covariates 

		    for (ll in 1:dim(categories)[1])
        {
           
           if (is.null(no.cov)) categories.cov<-NULL else categories.cov<-c(categories[ll,])
           if (is.null(no.cov)) categories.prop<-1 else categories.prop<-prod(tab_mod_prop[[2]][ll,])  #si no covariate categories.prop=1

	     	   protk<-pts(nr,i,prots)
           lprotk<-lapply(1:nr,function(nr,protk) length(protk[[nr]]),protk=protk )
           sensi<-sensifixe(nr,protk,pfixe,p,i,lprotk,doses,categories.cov,names.cov,beta.covariate,theta2,tab.cov.param,n_occB,sequ) 
           mod<-sensi[[1]] 
           sensif<-sensi[[2]]   #matrice derivee premiere effets fixes (avec covariables si demande)
           
           sensia<-sensialea(nr,protk,pfixe,p,i,lprotk,doses,categories.cov,names.cov,beta.covariate,theta2,tab.cov.param,n_occB,sequ)[[2]]
           
           
           if (n_occB>1) bout1<-OMEGA%*%sensia else bout1<-omega%*%sensia
	       
			     bout4<-list()
           t1<-0
	      	 t<-0
	      	 for (occ in 1:n_occB){
      	   bout4[[occ]]<-list()
			     for (n in 1:(nr))
           {
		    	     t<-1+t1
			         t1<-t-1+lprotk[[n]]
			
			         if(lprotk[[n]]==0)
                {
			            bout4[[occ]][[n]]<-NULL 
                #if(n==nr){bout4[[occ]][[n+1]]=Inf}
                } 
                else 
                {
		        	     bout3<-sigmainter[[n]]+sigmaslope[[n]]*mod[t:t1]
		        	     bout4[[occ]][[n]]<-bout3
		            }
           }
           }
           
		    	 t<-NULL
			     t1<-NULL
           diagsig<-diag(unlist(bout4)^2,length(unlist(protk))*n_occB)

			     invdiagsig<-solve(as.matrix(diagsig))
			    
			     
           if (FIM=="I") {fish<-bloc_IFIM(protk,sensif,sensia,diagsig,invdiagsig,pfixe,p,dvar,mod,lprotk,bout4,categories.cov,sequ)}
           if (FIM=="B") {fish<-bloc_BFIM(sensif,invdiagsig,pfixe,p,omega); subjects<-c(1)} 
		                                                              
           var<-t(bout1) %*% sensia + diagsig 
			     if (FIM=="P") {
            if (length(which(dim(var)==0))!=2) {invvar<-solve(as.matrix(var))
            fish<-bloc_PFIM(protk,sensif,sensia,var,invvar,pfixe,p,dvar,mod,lprotk,bout4,categories.cov,sequ) 
           }
           else {fish<-matrix(c(rep(0,pp*pp)),ncol=pp)} }

          somme<-somme+subjects[i]*fish*Nbsubjects*categories.prop*comb.seq.prop[sequ]   
                
        }    
			 }
	   }
   

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
      pfixe<-pfixe-length(which((beta.fixed+(diag(omega)==0))>0))}
    else{
      beta.estim<-!beta.fixed 
      #theta<-theta[beta.estim]
      ncovfixe<-pfixe-p
      p<-length(beta.estim)
      pfixe<-pfixe-length(which(beta.estim==F)) }		

    OMEGABIS<-diag(c(diag(omega),diag(gamma)))
    OMEGABIS<-as.matrix(OMEGABIS)
    if (length(c(diag(omega),diag(gamma)))==1) OMEGABIS<-as.matrix(c(diag(omega),diag(gamma)))
		if (n_occB>1) vec1<<-diag(OMEGABIS)!=0 
		if (n_occB==1) vec1<<-diag(omega)!=0 
    OMEGA1<<-as.matrix(as.matrix(OMEGABIS[,vec1])[vec1,]) 
   	omega1<<-as.matrix(as.matrix(omega[,diag(omega)!=0])[diag(omega)!=0,])
   	gamma1<<-as.matrix(as.matrix(gamma[,diag(gamma)!=0])[diag(gamma)!=0,])
    if (n_occB==1) OMEGA1<-omega1
    
		#shrinkage
    sh<-c()
	  
 
    if (FIM=="P") {vec<-c(beta.estim,rep(T,ncovfixe+locc),vec1,vec2)} #remove the row i and the col i of somme if omega[i,i]=0  
 	  if (FIM=="I"){vec<-c(beta.estim,rep(T,locc),vec2)}
 	  if (FIM=="B"){
    #shrinkage matrix before removing the fixed parameters
    if (Trand==1) I_W<-solve(as.matrix(somme))*solve(as.matrix(OMEGA1))
    else I_W<-solve(as.matrix(somme))*solve(as.matrix(diag((beta[diag(omega)!=0]),length(beta[diag(omega)!=0])))%*%as.matrix(omega1)%*%as.matrix(diag(((beta[diag(omega)!=0])),length(beta[diag(omega)!=0]))))
    
    #sh<-diag(as.matrix(I_W))*100
    vec<-c(beta.estim,rep(T,locc))
    #shrinkage after removing fixed parameters
	  sh<-diag(as.matrix(I_W[vec,vec]))*100
	  sh<-sh[vec]
    }
    
    somme<-somme[vec,vec]
		
	
#taking into account previous information
    if (previous.FIM!=  ""){
    previous_information<-read.table(file.path(directory,previous.FIM),header=FALSE,sep="")
    #On teste si la previous information matrix a la meme taille que la matrice actuelle. Si ce n'est pas le cas, on arrete l'execution et on affiche un message d'erreur : 
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
	   #print(invcrit)
     if (ind==0) l<-invcrit else
		    { 
		    f<-fisher.cv(somme,pfixe,pp,SIG,OMEGA1,beta.covariate)
		    inv<-f[[1]]
		    se<-f[[2]]
		    cv<-f[[3]]
		    l<-list(somme=somme,doses=dose,prots=prots,inv=inv,se=se,cv=cv,det=det,crit=crit,p=p,pp=pp,subjects=subjects,xn=xn,pfixe=pfixe,sh=sh)
 		     }
 	}

 return(l)
  
}



#function to build graph


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
              if(graph.only==F){
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




