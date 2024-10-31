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




#--------------------------------------------------------------------------------
#To associate proportions of categories in a same tab
#--------------------------------------------------------------------------------
association_modality_covariate<-function(covariate.name,covariate.cat.reco,covariate.proportions,covariate.cat.poss)
{
pro.cov<-NULL
summary.pro.cov<-NULL

lcov.name<-length(covariate.name) 
lmod.poss<-dim(covariate.cat.poss)[1]
prop.comb<-matrix(nrow=lmod.poss,ncol=lcov.name)
for ( i in 1:lcov.name)
{
    pro.cov<-cbind(rep(covariate.name[[i]],length(covariate.cat.reco[[i]])),as.numeric(covariate.cat.reco[[i]]),as.numeric(covariate.proportions[[i]]) )
    #pro.cov<-cbind(covariate.category[[i]],covariate.proportions[[i]])
    summary.pro.cov<-rbind(summary.pro.cov,pro.cov)
}

for (i in 1:lcov.name){
                                          for (ij in 1: lmod.poss)
                                          {
                                          names.cov.ass<-covariate.name[[i]]
                                          mod.cov.ass<-covariate.cat.poss[ij,i]
                                          stab<-summary.pro.cov[summary.pro.cov[,1]==names.cov.ass,]
                                          prop.comb[ij,i]<-as.numeric(stab[stab[,2]==mod.cov.ass,3])
                                          }
                 }                                       

return(list(summary.pro.cov,prop.comb))
}



test<-function(i,protvec)
{
	res<-rep(F,nr)
	res2<-list()
	
	l3<-lapply(1:nr,function(nr,l2) sum(unlist(l2[[nr]])),l2=l2)
	p<-c(0,cumsum(unlist(l3)))
	l4<-length(p)
	protveck<-lapply(1:(l4-1), function (l4,protvec,p) c(protvec[(p[l4]+1):(p[l4+1])]),protvec=protvec,p=p)
	 
  lr<-lapply(1:(l4-1),function(l4,p,i) ifelse(i>(p[l4]) && i<=(p[l4+1]),(l4),0),p=p,i=i)	
  l<-max(unlist(lr))

          n<-length(lower[[l]])
          res2<-lapply(1:n,function(n,protvec,lower,upper,res,l) c(res[l]||(protvec[i]>=lower[[l]][[n]] && protvec[i]<=upper[[l]][n])),protvec=protvec,upper=upper, lower=lower,res=res,l=l)
  res3<-lapply(1:n,function(l,res2) length(which(unlist(res2)==T)),res2=res2)
  #print(protvec)
  #print(upper)
  #print(lower)
	ifelse((all(unlist(res3)==1)),0,1)

}

testsubjects<-function(subjects,Q)
{
	res<-0
	if ((sum(subjects[1:(Q-1)])>1)) res<-1
	else 
		for ( i in 1:Q)
	{
		if (subjects[i]<0) {
					res<-1
					break
					}
				else res<-0
	}
	return(res)
}


transveclist<-function(protvec,Q,l_prote,subjects)
{
	#protPK<-vector("list",Q)
	#protPD<-vector("list",Q)
	
  l3<-lapply(1:nr,function(nr,l2) sum(unlist(l2[[nr]])),l2=l2)
	p<-c(0,cumsum(unlist(l3)))
	l4<-length(p)
	protveck<-lapply(1:(l4-1), function (l4,protvec,p) c(protvec[(p[l4]+1):(p[l4+1])]),protvec=protvec,p=p)
	protveck1<-protveck
  prot<-list()	
	protvect<-list()
	
  #protvec1<-protvec[(sum(nq)+1 ):(length(protvec))]
  #protvec<-protvec[1:sum(nq)]
	for(l in 1:nr)
	{

	prot[[l]]<-list()
	for (i in 1:Q){
	if(l2[[l]][[i]]!=0) {
	prot[[l]][[i]]<-protveck[[l]][1:l2[[l]][[i]]]
	protveck[[l]]<-protveck[[l]][-(1:l2[[l]][[i]])]
	}
	if (l2[[l]][[Q]]==0){prot[[l]][[Q+1]]=Inf}
	}
	}
	p2<-c()
	if (length(unlist(protveck1))!=length(protvec))
	{
  p2<-protvec[-(1:length(unlist(protveck1)))]
  } 

	if (length(p2)==0) NULL else subjects<-c(unlist(p2),1-sum(p2))
   
	return(list(prot=prot,subjects=subjects))
		
}
		
		
constraints.simplex<-function(protvec,Q,l2,l_prote,subjects)
{
	#indic.res<-length(which((sapply(1:sum(nq+np),test,protvec))!=0))
	indic.res<-length(which((sapply(1:sum(unlist(l_prote)),test,protvec))!=0))
	if((indic.res!=0)) {
				cont<-1
				l<-list(cont)
				}
		else
		{	
		      trans<-transveclist(protvec,Q,l_prote,subjects)        #transveclist<-function(protvec,Q,nq,np)
					prots<-lapply(1:nr,function(nr,transprot) c(transprot[[nr]]),transprot=trans$prot)
          #prots[[1]]<-
					#prots[[2]]<-trans$prot[[2]]
					subjects<-trans$subjects
					res<-testsubjects(subjects,Q)                #testsubjects<-function(subjects,Q)
					#indic.delta1<-0
					#indic.delta2<-0
					
		  for (i in 1:Q)
		{

		indic.delta1<-lapply(1:nr,function(nr, prots,l2,delta.time) length(which(abs(prots[[nr]][[i]][2:l2[[nr]][[i]]]-prots[[nr]][[i]][1:(l2[[nr]][[i]]-1)])<delta.time)),prots=prots,l2=l2,delta.time=delta.time)

    }
			
		if((all(unlist(indic.delta1))!=0) ||(res==1))	cont<-1 else cont<-0
		l<-list(cont=cont,subjects=subjects,prots=prots)
		}
		
	return(l)				
}


fisher.simplex<-function(protvec,arg)

{
	Q<-arg[[2]]
	#nq<-arg[[3]]
	#np<-arg[[4]]
	l2<-arg[[3]]
	l_prote<-arg[[4]]
	subjects<-arg[[5]]
	Ntot<-arg[[6]]
	nq<-arg[[7]]
	categories<-arg[[8]]
	names.cov<-arg[[9]]
	proportions.cov<-arg[[10]]
	beta.covariate<-arg[[11]]
	tab_mod_prop<-arg[[12]]
	pfixe<-arg[[13]]
  p<-arg[[14]]
  theta2<-arg[[15]]
  tab.cov.param<-arg[[16]]
  n_occ<-arg[[17]]
  lposs<-arg[[18]]
  
	cont<-constraints.simplex(protvec,Q,l2,l_prote,subjects)   #constraints.simplex<-function(protvec,Q,nq,np,subjects)
	if (cont[[1]]==1) invcrit<-1e+99
		
	else
	{
	subjects<-cont[[2]]
	cat("subjects=",subjects,'\n')
	prots<-cont[[3]]
	count<<-count+1
	cat('count = ', count,'\n')	
	
  #prots<-lapply(1:nr,function(nr,Q,prots) lapply(1:Q,function(Q,prots) sort(prots[[nr]][[Q]]),prots=prots),Q=Q,prots=prots)

  invcrit<-fisher(prots,subjects,condinit,doses,Ntot,ind=0,l_prote,nq,categories,names.cov,proportions.cov,beta.covariate,tab_mod_prop,pfixe,p,theta2,tab.cov.param,n_occ,lposs)[[1]]
 #fisher<-function(prots,subjects,condinit,doses,Ntot,ind=0,l_prote,nq,categories,names.cov,proportions.cov,beta.covariate,tab_mod_prop,pfixe,p,theta2,tab.cov.param,n_occ,lposs)

	
	#cat('fisher ESSAI ')
	#print(fisher(prots,subjects,condinit,doses,Ntot,ind=0,l_prote,nq,categories,names.cov,proportions.cov,beta.covariate,tab_mod_prop,pfixe,p,theta2,tab.cov.param,n_occ,lposs))
	
	
	
	
	cat('crit=',1/invcrit,'\n')
	}
	return(invcrit)
}

bbbbb = c()
count<<-0

#fonction d'optimisation d'un protocole
#--------------------------------------
fisher.opti<-function(Q,l2,Ntot,nq,covariate.cat.poss,names.cov,proportions.cov,beta.covariate,tab_mod_prop,pfixe,p,theta2,tab.cov.param,n_occ,lposs)
{
	
	if ((subjects.opt==T) && (Q!=1)) 
    {#cas ou des NULL 
    #ll<-lapply(prots,length)
		#if(length(which(unlist(lapply(1:nr,function(nr) lapply(1:ll[[nr]], function(ll,prots,nr) length(prots[[nr]][[ll]]),prots=prots,nr=nr)))==0))==0)
    protkr<-c(unlist(prots),subjects.init[1:(Q-1)])
    } else {
    protkr<-unlist(prots) 
    }
         lpi<-length(protkr)
         ll<-lapply(prots,length)
  		   sv<-length(which(unlist(lapply(1:nr,function(nr) lapply(1:ll[[nr]], function(ll,prots,nr) length(prots[[nr]][[ll]]),prots=prots,nr=nr)))==0))==0
         if (identical.times==T && nr!=1 && sv==T){
            if(Q!=1){
            qq<-matrix(protkr,lpi,byrow=T,nrow=sum(unlist(l2[[1]]))) 
            diag(qq)<-diag(qq)-diag(qq)*simplex.parameter/100
            qqr<-qq[,1:sum(unlist(l2[[1]]))]
            qq1<-qqr
            for (i in 1:(nr-1)){
            qq1<-cbind(qq1,qqr)
            }
            qq<-cbind(qq1,matrix(protkr[(sum(unlist(l2))+1):lpi],ncol=length(protkr[(sum(unlist(l2))+1):lpi]),nrow=sum(unlist(l2[[1]]))))
            
            }else{
         qq<-matrix(protkr,l2[[1]][[1]],byrow=T,nrow=l2[[1]][[1]])
         diag(qq)<-diag(qq)-diag(qq)*simplex.parameter/100
         qq1<-qq
         for (i in 1:(nr-1)){
         qq1<-cbind(qq1,qq)
         }
         qq<-qq1
         qq1<-NULL
	       #qq<-matrix(rep(protPKPD,lpi),byrow=T,nrow=lpi)
	       }
         }else{
	       qq<-matrix(rep(protkr,lpi),byrow=T,nrow=lpi)	       
	       diag(qq)<-diag(qq)-diag(qq)*simplex.parameter/100
	       }
	       #matdeb<-rbind(protPKPD,qq)
         matdeb<-rbind(protkr,qq)
         #print(matdeb)
         # y<-matdeb - Single run
        y<-apply(matdeb,1,fisher.simplex,arg=list(matdeb,Q,l2,l_prote,subjects,Ntot,nq,covariate.cat.poss,names.cov,proportions.cov,beta.covariate,tab_mod_prop,pfixe,p,theta2,tab.cov.param,n_occ,lposs))     #fisher.simplex<-function(protvec,arg)
	degenere<-which(y==1e99)                   
	if (length(degenere)==(lpi+1)) 
	{
		cat("Error: the initial population design and all the vertices are very poor (criterion = 0); the Simplex algorithm cannot converge","\n")
		opti<-NULL
		eval<-NULL
		count<-NULL
		protopti<-NULL
		subjectsopti<-NULL
	}else
	{

	opti<-fun.amoeba(matdeb, y, Rctol, Max.iter,fisher.simplex,data=list(matdeb,Q,l2,l_prote,subjects,Ntot,nq,covariate.cat.poss,names.cov,proportions.cov,beta.covariate,tab_mod_prop,pfixe,p,theta2,tab.cov.param,n_occ,lposs),vrbs = iter.print) 
	
		
	if (Max.iter==0) protopti<-protkr else protopti<-(opti$p[which((opti$y==min(opti$y))==T)[1],])

	ind<<-1
	res<-transveclist(protopti,Q,l_prote,subjects)
	
	
	protopti<-res[[1]]
	subjectsopti<-res[[2]]
	#if (modelform=="DE") doses<-NULL
	eval<-fisher(protopti,subjectsopti,condinit,doses,Ntot,ind=1,l_prote,nq,covariate.cat.poss,names.cov,proportions.cov,beta.covariate,tab_mod_prop,pfixe,p,theta2,tab.cov.param,n_occ,lposs)
  # bbbbb <<- c(bbbbb, eval$crit)
	
	print(paste('JKDFKJDFKLSDFJLKDLKJFLKJ', eval$crit))
  }     
	return(list(opti=opti,eval=eval,count=count,protopti=protopti,subjectsopti=subjectsopti,degenere=degenere))
	
}
#Calculer l'intervalle de confiance
###################################

CI.computation<-function(alpha,beta.cov,beta.cov.se,nom.par,Trand)
{
norm.val<-abs(qnorm(alpha, mean=0, sd=1))

Beta<-beta.cov
lcov<-length(beta.cov)
#confidence interval for beta
CI_inf<-beta.cov-norm.val*beta.cov.se
CI_sup<-beta.cov+norm.val*beta.cov.se 
t<-unlist(lapply(1:lcov,function(lcov) paste("[",round(CI_inf[lcov],3),";",round(CI_sup[lcov],3),"]",sep="")))
lim<-paste(100*(1-2*alpha),"% CI")  
if (Trand==2)
{
#Confidence interval for exp(beta)   
Beta_exp<-exp(beta.cov)
vide<-rep(" ",lcov) 
t1<-unlist(lapply(1:lcov,function(lcov) paste("[",round(exp(CI_inf[lcov]),3),";",round(exp(CI_sup[lcov]),3),"]",sep=""))) 
ci<-data.frame(Beta,t,vide,Beta_exp,t1,row.names=paste(nom.par))
names(ci)<-c("Beta",lim," ","exp(Beta)",lim)
}
else
{ci<-data.frame(Beta,t,row.names=paste(nom.par))
names(ci)<-c("Beta",lim)
}

return(ci)

}


#Calculer la puissance et le nombre de sujet pour le test de comparaison
power_nni.computation<-function(alpha,given.power,init.num.subjects,beta.cov,beta.cov.se,nom.par,Trand)
{
norm.val<-abs(qnorm(alpha, mean=0, sd=1))

#Puissance avec cet effet la avec alpha =5%
power.prediction<-1-pnorm(norm.val-(beta.cov/beta.cov.se)) + pnorm(-norm.val-(beta.cov/beta.cov.se))

Beta<-beta.cov
Predicted_SE<-beta.cov.se
lcov<-length(beta.cov)  
Expected_power<-round(power.prediction,6)

#Calcul du nombre de sujet necessaires pour avoir une puissance donnee
#given.power
#calcul de la SE si je me met sous H1
NSE<-beta.cov/(norm.val-qnorm(1-given.power))

#calcul maintenant du nombre de sujet necessaires pour mettre en evidence cette effets
#avec une puissance donnee
NNI<-sum(init.num.subjects)*(beta.cov.se/NSE)^2

Number_subjects_needed<-NNI

#Affichage
vide<-rep(" ",lcov)
if (compute.power==T){
result<-data.frame(Expected_power,row.names=paste(nom.par))
names(result)<-c("Expected_power")
if (compute.nni==T){
result<-cbind(result,Number_subjects_needed,vide)
names(result)<-c("Expected_power","Number_subjects_needed",paste("(for a given power=",given.power,")",sep=""))
}
}
if (compute.power==F && compute.nni==T){
result<-data.frame(Number_subjects_needed,vide,row.names=paste(nom.par))
names(result)<-c("Number_subjects_needed",paste("(for a given power=",given.power,")",sep=""))
}
return(result)
}



#Calculer la puissance et le NSN pour le test d'equivalence
power_nni_eq.computation<-function(alpha,given.power,init.num.subjects,beta.cov,beta.cov.se,nom.par,Trand)
{
norm.val<-abs(qnorm(alpha, mean=0, sd=1))
Beta<-beta.cov
Predicted_SE<-beta.cov.se
lcov<-length(beta.cov)

#puissance
power.prediction<-c()
for (i in 1:length(beta.cov)){
if(beta.cov[i]==0){power.prediction[i]<-pnorm(-norm.val+(interval_eq[2]/beta.cov.se[i]))}
if(beta.cov[i] <0){power.prediction[i]<-1-pnorm(norm.val-(beta.cov[i]-interval_eq[1])/beta.cov.se[i])}
if(beta.cov[i] >0){power.prediction[i]<-pnorm(-norm.val-(beta.cov[i]+interval_eq[1])/beta.cov.se[i])}
#if (beta.cov[i] < interval_eq[1] || beta.cov[i] > interval_eq[2]) {power.prediction[i]<-c(".")
#test_eq<-F}
}




#Calcul du nombre de sujet necessaires pour avoir une puissance donnee%
#given.power
#calcul de la SE si je me met sous H1
NSE<-c()
for (i in 1:length(beta.cov)){
if (beta.cov[i]==0){NSE[i]<-interval_eq[2]/(norm.val+qnorm(given.power))}
if (beta.cov[i]<0){NSE[i]<-(-beta.cov[i]+interval_eq[1])/(-norm.val+qnorm(1-given.power))}
if (beta.cov[i]>0){NSE[i]<-(-beta.cov[i]-interval_eq[1])/(norm.val+qnorm(given.power))}
}
#calcul maintenant du nombre de sujet necessaires pour mettre en evidence cette effets
#avec une puissance donnee
NNI<-sum(init.num.subjects)*(beta.cov.se/NSE)^2



#si pas sous H1

for (i in 1:length(beta.cov)){
if (beta.cov[i] < interval_eq[1] || beta.cov[i] > interval_eq[2]) {
power.prediction[i]<-c(NA)
NNI[i]<-c(NA)
}
}

Expected_power<-power.prediction
Number_subjects_needed<-NNI
#Affichage
vide<-rep(" ",lcov)

if (compute.power_eq==T){
result<-data.frame(Expected_power,row.names=paste(nom.par))
names(result)<-c("Expected_power")
if (compute.nni_eq==T){
result<-cbind(result,Number_subjects_needed,vide)
names(result)<-c("Expected_power","Number_subjects_needed",paste("(for a given power=",given.power,")",sep=""))
}
}
if (compute.power_eq==F && compute.nni_eq==T){
result<-data.frame(Number_subjects_needed,vide,row.names=paste(nom.par))
names(result)<-c("Number_subjects_needed",paste("(for a given power=",given.power,")",sep=""))
}
return(result)
}


#-------------------------------------------------------------------------
#------------------------------ OUTPUT -----------------------------------

out.simplex<-function()
{

#graphs before optimisation
###########################

if (graph.only==T){ 
  if (modelform == "AF" && NUM==T){
    source(file.path(directory.program,"Pfim4.0op1.r"))
    if (graph.logical==T) {graphAFN(form,prots,graphinf,graphsup,names.datax,names.datay,y.range)}  else NULL
    if (graphsensi.logical==T) {graphsensiAFN(form,prots,graphinf,graphsup,names.datax,names.datay,y.range)}  else NULL
  }
  if (modelform == "AF" && NUM==F){
    source(file.path(directory.program,"Pfim4.0op1.r"))
    if (graph.logical==T) {graphAFS(formg_init,prots,graphinf,graphsup,names.datax,names.datay,y.range)} else NULL
    if (graphsensi.logical==T) {graphsensiAFS(formg_init,prots,graphinf,graphsup,names.datax,names.datay,y.range)} else NULL
  }
  if (modelform == "DE"){
    source(file.path(directory.program,"EQPfim4.0op1.r"))
    if (graph.logical==T)  {graph(formED,prots,graphinf,graphsup,names.datax,names.datay,y.range)} else NULL
    if (graphsensi.logical==T) {graphsensi(formED,prots,graphinf,graphsup,names.datax,names.datay,y.range)} else NULL
  }
  options("show.error.messages"=F);stop()
}
    

  
###COMPUTATION#############################################################################################

d<-Sys.time()
  g<-proc.time()
	stop<-1
  if (pn==1 || locc==0) no<-NULL
      if (n_occ>1 && locc>0){
      no<-NULL
      no1<-NULL
      lo<-NULL
  if (covariate_occ.model==T){
  
      #TEST COVARIATE_OCC MODEL
          if ((length(covariate_occ.name)!= length(covariate_occ.category))==T) 
          stop(" The list of the name of the covariate(s) depending also on the occasion has not the same length of the list of the category")
          
          if(length(which((sapply(covariate_occ.sequence,length) == sapply(sequence.proportions,length))==F))!=0 )
          stop(" The number of covariate categorie sequences is not equal to the number of proportions")
          
          if (length(covariate_occ.name)!= length(parameter_occ.associated))
          stop(" One or more covariate(s) has not at least one associated parameter ")
          
          if(sum(unlist(sapply(sequence.proportions,sum)))!=length(covariate_occ.name))
          stop(" The sum of the proportions of the categories is not equal to one")
      

      count_occ<-0
      for (k in 1:length(covariate_occ.name))
      {
          for (j in 1:length(parameter_occ.associated[[k]]))
          {
            
              lo<-length(covariate_occ.category[[k]])-1
              no1<-lapply(1:lo,function(lo) paste("beta_",parameter_occ.associated[[k]][j],"_",covariate_occ.name[[k]],"_",covariate_occ.category[[k]][lo+1],sep=""))
              no<-c(no,unlist(no1)) #liste des parametres
            }
             count_occ<-count_occ+lo
          }

      }
   }



	if((subjects.input==2)&&(sum(subjects.init)!=1)) 
	{
		cat("STOP: the sum of proportions of subjects must be equal to 1.")
		stop<-0  
	}
	
	
		
	if(((length(doses)!=Q)||(length(subjects.init)!=Q))&(modelform=="AF")) 
	{
		cat("STOP: the number of doses or the number of proportions of subjects must be equal to the number of elementary designs.")
		 stop<-0
	}

	
  for (i in 1:Q)
		{

		indic.delta1<-lapply(1:nr,function(nr, prots,l2,delta.time) length(which(abs(prots[[nr]][[i]][2:l2[[nr]][[i]]]-prots[[nr]][[i]][1:(l2[[nr]][[i]]-1)])<delta.time)), prots=prots,l2=l2,delta.time=delta.time)

    }
    
	if(sum(unlist(indic.delta1))!=0)
	{
		cat("STOP: the minimum delay between two sampling times should be of",delta.time,".")
		stop<-0
		   
	}

  ll<-lapply(prots,length)
		if(identical.times==T && length(which(unlist(lapply(1:nr,function(nr) lapply(1:ll[[nr]], function(ll,prots,nr) length(prots[[nr]][[ll]]),prots=prots,nr=nr)))==0))!=0)
		       identical.times=F
		       
    if (covariate.model==T && option==2)
        stop("The full expression of the Fisher information matrix (option 2) is not available with influence of covariates and interoccasions varabilities")
  

if (stop!=0)
{
	 #individual design
   b_condition<-T
   if (length(omega)!=0 && length(which(omega==0))==length(omega)) 
   {b_condition<-F;subjects.opt=F;Omega<-NULL;
    cat("All elements of the vector omega are equal to 0, no random effect are considered","\n", "You can not optimize subject proportions if you have several groups")
   
   }
   
   
	 if(covariate.model==T)
   { #TEST COVARIATE MODEL
          if(FIM!="P") stop("Individual and bayesian matrices are only computed without covariate")
          
          if ((length(covariate.name)!= length(covariate.category))==T || (length(covariate.category)!=length(covariate.proportions))==T) 
          stop(" The list of the name of the covariate(s) has not the same length of the list of the category or the proportions")
          
          if(length(which((sapply(covariate.category,length) == sapply(covariate.proportions,length))==F))!=0 )
          stop(" The number of categories is not equal to the number of proportions")
          
          if (length(covariate.name)!= length(parameter.associated))
          stop(" One or more covariate(s) has not at least one associated parameter ")
          
          if(sum(unlist(sapply(covariate.proportions,sum)))!=length(covariate.name))
          stop(" The sum of the proportions of the categories is not equal to one")
          
        
          covariate.cat.reco<-list()
          info.covariate<-list(covariate.name,covariate.category,covariate.proportions,parameter.associated,beta.covariate)
          info.power.nni<-list(alpha,compute.power,compute.nni,given.power)
          
          #recodage covariable si en character en vecteur numerique#METTRE UN lAPPLY ???
          for (i in 1: length(info.covariate[[1]]))
          {
          if(is.character(info.covariate[[2]][[i]]))
             covariate.cat.reco[[i]]<-c(1:length(info.covariate[[2]][[i]]))
          else 
             covariate.cat.reco[[i]]<-info.covariate[[2]][[i]]
          }
          
          #solution possibles pour les covariables
          covariate.cat.poss<-unique(t(combn(unlist(covariate.cat.reco),length(info.covariate[[1]]))))
          #verifier que les modalites presentes dans chque colonne sont bien des modalites possibles de la covariables
          #test si les modalites possibles pour cahque covaraibles ne sont pas abberantes
          ll.cov<-length(info.covariate[[1]])
          for (i in 1: ll.cov)
          {
            #nline<-which(covariate.cat.poss[,i] > max(info.covariate[[2]][[i]]))
            nline<-which(covariate.cat.poss[,i] > max(covariate.cat.reco[[i]]) | covariate.cat.poss[,i] < min(covariate.cat.reco[[i]]))
            ref<-integer(length = 0)
        
            if (!identical(nline,ref)) covariate.cat.poss<-covariate.cat.poss[-nline,]
          }
          
          #creation of a table with the parameters aand the covaraite associated 
          npar<-NULL
          npar1<-NULL
          lmodal<-NULL
          count<-0
          ilnp<-0
          s.cov<-matrix(nrow=length(unlist(info.covariate[[5]])),ncol=2)
          for (ijk in 1:length(info.covariate[[1]]))
          {
            for (ijk1 in 1:length(info.covariate[[4]][[ijk]]))
            {
              
              #s.cov[count+ijk1,2]<-info.covariate[[1]][[ijk]]
              #s.cov[count+ijk1,1]<-info.covariate[[4]][[ijk]][ijk1]
              s.cov[count+1,2]<-info.covariate[[1]][[ijk]]
              s.cov[count+1,1]<-info.covariate[[4]][[ijk]][ijk1]
              lnp<-length(info.covariate[[2]][[ijk]])-1
              npar1<-lapply(1:lnp,function(lnp) paste("beta_",info.covariate[[4]][[ijk]][ijk1],"_",info.covariate[[1]][[ijk]],"_",covariate.cat.reco[[ijk]][lnp]+1,sep=""))
              npar<-c(npar,unlist(npar1)) #liste des parametres
              
               #liste des modalites
              lmodal1<-lapply(1:lnp,function(lnp) c(covariate.cat.reco[[ijk]][lnp]+1))
              lmodal<-c(lmodal,unlist(lmodal1))
              
              
              if(lnp>1)
              {
                for (i in 2:(lnp))
                { 
                  ilnp<-ilnp+1
                  s.cov[count+1+ilnp,2]<-info.covariate[[1]][[ijk]]
                  s.cov[count+1+ilnp,1]<-info.covariate[[4]][[ijk]][ijk1] 
                  ilnp<-0
                }
              }
             count<-count+lnp
            }
            #count<-count+length(info.covariate[[4]][[ijk]])+ilnp 
            #ilnp<-0
          }
            tab.cov.param<-cbind(s.cov,npar,lmodal,unlist(info.covariate[[5]]))
            
            ###tranform in matrix
            #tab.cov.param1<-tab.cov.param[order(tab.cov.param[,1]),]
            #tab.cov.param<-matrix(tab.cov.param1,nrow=dim(tab.cov.param)[1],ncol=dim(tab.cov.param)[2])
            
            tab.cov.param<-matrix(tab.cov.param,nrow=dim(tab.cov.param)[1],ncol=dim(tab.cov.param)[2])
            
            #number of fixed effects with the covariate
            pfixe<-length(c(parameters,npar))
        
            #summary info covariate model
            #association des prportions aux differentes combianaisons
            tab_mod_prop<-list()
            tab_mod_prop<-association_modality_covariate(info.covariate[[1]],covariate.cat.reco,info.covariate[[3]],covariate.cat.poss)
            nq<-NULL
             
          	f<-fisher.opti(Q,l2,Ntot,nq,covariate.cat.poss,info.covariate[[1]],info.covariate[[3]],info.covariate[[5]],tab_mod_prop,pfixe,p,theta2,tab.cov.param,n_occ,lposs)     # fisher.opti<-function(Q,nq,np,Ntot)
          	}
            else 
            {
             covariate.cat.poss<-NULL
             info.covariate<-list(NULL,NULL,NULL,NULL,NULL)
             tab_mod_prop<-NULL
             pfixe<-length(parameters)
             tab.cov.param<-NULL
             nq<-NULL
             f<-fisher.opti(Q,l2,Ntot,nq,covariate.cat.poss,info.covariate[[1]],info.covariate[[3]],info.covariate[[5]],tab_mod_prop,pfixe,p,theta2,tab.cov.param,n_occ,lposs)     # fisher.opti<-function(Q,nq,np,Ntot) 
            }
  
  y<-f$opti$y 
	prot.opti<-f$protopti
	subjects.opti<-f$subjectsopti
	num.iter<-f$opti$iter
	conv<-f$opti$converge
	count<-f$count
	crit.init<-1/(fisher(prots,subjects,condinit,doses,Ntot,ind=0,l_prote,nq,covariate.cat.poss,info.covariate[[1]],info.covariate[[3]],info.covariate[[5]],tab_mod_prop,pfixe,p,theta2,tab.cov.param,n_occ,lposs)) 
	#fisher<-function(prots,subjects,condinit,doses,Q,Ntot,ind=0,nq,np)	
	
	if (subjects.input==1) subjects<-subjects.initN else subjects<-subjects.init
  
  subjects.prop<-subjects
	
	sink(file.path(directory,output)) 
	cat("PFIM 4.0",'\n',"\n")
	cat("Option",option)
	cat("\n","\n") 
	cat("Project:",project)
	cat("\n","\n") 
	cat('Date:',date()) 
	cat("\n","\n") 
	cat("\n","\n") 
	
	if(modelform=="DE")        #differential equation
	{
	  pm<-0
    ll<-lapply(prots,length)
		if(length(which(unlist(lapply(1:nr,function(nr) lapply(1:ll[[nr]], function(ll,prots,nr) length(prots[[nr]][[ll]]),prots=prots,nr=nr)))==0))==0)
    {
        e<-list()
        for (l in 1 : nr) 
        {	 
           e[[l]]<-data.frame(times=as.character(prots[[l]][1:Q]),subjects)
        }
        pm<-1
      }
  
  
	   cat("**************************** INPUT SUMMARY ********************************") 
		 cat("\n","\n") 
		 cat('Differential Equations form of the model: ','\n','\n') 
		 print(formED)
		 cat("\n","\n")
		 cat('Initial design: ','\n')
		 cat("\n","\n")
		  for (i in 1:nr){
            if (pm==1){
                cat('Sample times for response:',LETTERS[i],'\n')
                print(e[[i]])
            } 
            else{
              cat('Sample times for response:',LETTERS[i],'                         Number of subjects per group \n')
              cat(paste(prots[[i]],"      ", t(t(subjects)),'\n'))
              cat("\n","\n")
            }
      }
  cat("\n","\n")
  cat('Initial Conditions at time', time.condinit,':','\n','\n')
	   	ff<-sapply(1:Q, function(Q, form) cat(paste(form[[Q]][-1]),'\n'), form=condinit)
		  cat("\n","\n")
  
  if (b_condition==T && n_occ>1 && gamma !=rep(0,p)) 
  {
  cat("Number of occasions:",n_occ)
   cat("\n","\n")
  }
  
  

  if (FIM!="I"){
  cat("Random effect model: Trand = ",Trand) 
  cat("\n","\n")
  }
  
  ff1<-sapply(1:nr,function(nr,sigmainter,sigmaslope) cat(paste('Variance error model response',LETTERS[nr],': (',sigmainter[[nr]],'+',sigmaslope[[nr]],"*f)^2\n")),sigmainter<-sigmainter,sigmaslope<-sigmaslope)
	cat("\n","\n")
  cat("\n","\n")
  if (covariate.model==T || covariate_occ.model==T) {
    cat("Covariate model : ","\n","\n") 
    cat("\t")
    	     if (Trand==1) cat("NB: Covariates are additive on parameters") else cat("NB: Covariates are additive on log parameters")
          cat("\n","\n")
  }
  
		  
		  
  if (covariate.model==T)
	{
    	cat("Covariate model : ","\n","\n") 
    	     if (Trand==1) cat("NB: Covariates is additive on parameters") else cat("NB: Covariates is additive on log parameters")
          cat("\n","\n")
    	     e.cov<-list()
           lcov1<-length(info.covariate[[1]])
           #lcov<-length(info.covariate[[2]][[i]])
          for (i in 1:lcov1){
           #e.cov<-lapply(1:3,function(lcov1,Name,Categories,Proportions,Parameters_associated) data.frame(Name,Categories,Proportions,Parameters_associated),Name<-rep(paste(info.covariate[[1]][[lcov1]]),length(info.covariate[[2]][[lcov1]])), Categories<-info.covariate[[2]][[lcov1]],Proportions<-info.covariate[[3]][[lcov1]],Parameters_associated<-c(paste(info.covariate[[4]][[lcov1]]),rep(".",length(info.covariate[[2]][[lcov1]])-1)))
            Name<-rep(paste(info.covariate[[1]][[i]]),length(info.covariate[[2]][[i]]))
            lcov2<-length(info.covariate[[2]][[i]])
            prem_col<-lapply(1:lcov2,function(lcov2) paste("(",c(lcov2),")",sep=""))
            Categories<-info.covariate[[2]][[i]]
            Proportions<-info.covariate[[3]][[i]]
            References<-c("*",c(rep(" ",length(info.covariate[[2]][[i]])-1)))
            if(length(info.covariate[[4]][[i]])>1) {
                    charc1<-info.covariate[[4]][[i]][1]
                   for (ij in 2: length(info.covariate[[4]][[i]]))
                   {
                      
                      charc<-paste(charc1,info.covariate[[4]][[i]][ij],sep=",")
                      charc1<-charc
                   }
              Parameters_associated<-c(paste(charc1),rep(".",length(info.covariate[[2]][[i]])-1))
            } else 
            {
            Parameters_associated<-c(paste(info.covariate[[4]][[i]]),rep(".",length(info.covariate[[2]][[i]])-1))
            }
    
            e.cov[[i]]<-data.frame(Categories,References,Proportions,row.names=prem_col,check.names=F,check.rows=F)
          cat('Covariate',i,':',Name[1],'(',Parameters_associated[1],')',"\n")
                print(e.cov[[i]])
          cat("\n","\n")
          }
          
  }
  
  if (covariate_occ.model==T)
	{  
     
     cat("\n","\n")
     cat("\t") 
     cat("Covariates changing with occasion","\n","\n")
     e.cov1<-list()
     e.cov2<-list()
     
          for (k in 1: length(covariate_occ.name)){
          cat("\t")
          cat ("Covariate ",k,":",covariate_occ.name[[k]],"(",parameter_occ.associated[[k]],")","\n")
            
           lcov1<-length(covariate_occ.category[[k]])
           prem_col1<-lapply(1:lcov1,function(lcov1) paste("(",c(lcov1),")",sep=""))
           Categories<-covariate_occ.category[[k]]
           References<-c("*",c(rep(" ",length(covariate_occ.category[[k]])-1)))
           
           e.cov1[[k]]<-data.frame(Categories,References,row.names=prem_col1,check.names=F,check.rows=F)
           print(e.cov1[[k]])
           cat("\n","\n")
        
          lcov2<-length(sequence.proportions[[k]])
          prem_col2<-lapply(1:lcov2,function(lcov2) paste("(",c(lcov2),")",sep=""))
          Sequences<- list()

           for (sequ in 1:nb_seq[k]){
            Sequences[sequ]<-covariate_occ.sequence[[k]][[sequ]][1]
            for (occ in 2:n_occ){
              Sequences[sequ]<-paste(Sequences[sequ],covariate_occ.sequence[[k]][[sequ]][occ])
            }
           }

          #for (sequ in 1:nb_seq){
          #  Sequences[sequ]<-unlist(lapply(2:n_occ-1, function(n_occ) paste(covariate_occ.sequence[[k]][[sequ]][n_occ],covariate_occ.sequence[[k]][[sequ]][n_occ+1]))) 
          #}
           Proportions<-sequence.proportions[[k]]
          e.cov2[[k]]<-data.frame(unlist(Sequences),Proportions,row.names=prem_col2,check.names=F,check.rows=F)
          #. <- covariate_occ.sequence[[k]]
         
          #e.cov2[[k]]<-cbind(e.cov2[[k]],.,)
          names(e.cov2[[k]])<-c("Sequences","Proportions")
          
          print(e.cov2[[k]])
           cat("\n","\n")
          }
    }
  
		  cat("Computation of the Fisher information matrix: option = ",option) 
  		 cat("\n","\n") 
  	  cat('Total number of samples:',Ntot)
			cat("\n","\n") 
			cat('Associated criterion value:',round(crit.init,4))
			cat("\n","\n") 
			cat('Window of the allowed optimised sampling times for the ',nr,'responses:','\n')
      ff3<-sapply(1:nr,function(nr,lower,upper) cat(paste("Upper and lower admissible samples times for the response",nr,":",' [ ',lower[[nr]],":",upper[[nr]],"]",'\t ','\n'),sep=" "),lower=lower, upper=upper)
			cat("\n","\n") 
			cat('Minimum delay between two sampling times:',delta.time)
			cat("\n","\n") 
			cat('Optimisation of the proportions of subjects:',subjects.opt)
			cat("\n","\n") 
	    cat("\n","\n") 
		  cat("Error tolerance for solving differential equations system: RtolEQ =", RtolEQ,", AtolEQ =", AtolEQ,", Hmax = ",Hmax)
	    cat("\n","\n") 
	    
	}
	else            #IF ANALYTICAL FORM#
	{
  pm<-0
	  ll<-lapply(prots,length)
		if(length(which(unlist(lapply(1:nr,function(nr) lapply(1:ll[[nr]], function(ll,prots,nr) length(prots[[nr]][[ll]]),prots=prots,nr=nr)))==0))==0)
    {
        e<-list()
        for (l in 1 : nr) {	 
       e[[l]]<-data.frame(times=as.character(prots[[l]]),subjects.prop,doses)
         }
        pm<-1
    }
	
	   cat("**************************** INPUT SUMMARY ********************************") 
	   cat("\n","\n") 
	   cat('Analytical function model: ','\n','\n') 
	   if(NUM==T){print(form);cat("\n","\n")}
     else{ff<-sapply(1:lf, function(lf, form) cat(paste(form[lf]),'\n','\n'), form=form)}

     cat('Initial design: ','\n')
	   cat("\n","\n")
	   
	   
     for (i in 1:nr){
        cat('Sample times for response:',LETTERS[i],'\n')
        if (pm==1){
           print(e[[i]]) } 
        else{
           print(prots[[i]])
           cat('               Number of subjects per group \n',t(t(subjects.prop)),'\n')
           cat("\n","\n")
        }
           cat("\n","\n")
     }
  
     cat("\n","\n") 
     cat('Total number of samples (nr responses):',Ntot)
     cat("\n","\n") 
	   cat('Associated criterion value:', round(crit.init,4))
	   cat("\n","\n") 
	   if (nr!=1) cat('Identical sampling times for each response: ', identical.times) 
     cat("\n","\n")   
	   cat('Window of the allowed optimised sampling times:')
 	   cat("\n","\n") 
     ff3<-sapply(1:nr,function(nr,upper,lower) cat(paste('Upper and lower admissible samples times for the response',LETTERS[nr],': [' ,lower[[nr]],':',upper[[nr]],']','\n')),lower=lower, upper=upper)
	   cat("\n","\n")  
	   cat('Minimum delay between two sampling times:',delta.time)
	   cat("\n","\n") 
	   cat('Optimisation of the proportions of subjects:',subjects.opt)
	   cat("\n","\n") 
	   cat("\n","\n") 
	  
   if (b_condition==T && n_occ>1 && length(which(gamma==0))!=length(gamma)) 
  {
  cat("Number of occasions:",n_occ)
   cat("\n","\n")
  }
  
  if (FIM!="I"){
  cat("Random effect model: Trand = ",Trand) 
  cat("\n","\n")
  }
  
  ff1<-sapply(1:nr,function(nr,sigmainter,sigmaslope) cat(paste('Variance error model response',LETTERS[nr],': (',sigmainter[[nr]],'+',sigmaslope[[nr]],"*f)^2\n")),sigmainter<-sigmainter,sigmaslope<-sigmaslope)
	cat("\n","\n")
  
  if (covariate.model==T || covariate_occ.model==T) {
    cat("Covariate model : ","\n","\n") 
    cat("\t")
    	     if (Trand==1) cat("NB: Covariates are additive on parameters") else cat("NB: Covariates are additive on log parameters")
          cat("\n","\n")
  }
	 if (covariate.model==T)
	{
        cat("\t")
      	if (covariate_occ.model==T)     cat("Covariates not changing with occasion","\n","\n")
    	     e.cov<-list()
           lcov1<-length(info.covariate[[1]])
           #lcov<-length(info.covariate[[2]][[i]])
          for (i in 1:lcov1){
           #e.cov<-lapply(1:3,function(lcov1,Name,Categories,Proportions,Parameters_associated) data.frame(Name,Categories,Proportions,Parameters_associated),Name<-rep(paste(info.covariate[[1]][[lcov1]]),length(info.covariate[[2]][[lcov1]])), Categories<-info.covariate[[2]][[lcov1]],Proportions<-info.covariate[[3]][[lcov1]],Parameters_associated<-c(paste(info.covariate[[4]][[lcov1]]),rep(".",length(info.covariate[[2]][[lcov1]])-1)))
            Name<-rep(paste(info.covariate[[1]][[i]]),length(info.covariate[[2]][[i]]))
            lcov2<-length(info.covariate[[2]][[i]])
            prem_col<-lapply(1:lcov2,function(lcov2) paste("(",c(lcov2),")",sep=""))
            Categories<-info.covariate[[2]][[i]]
            Proportions<-info.covariate[[3]][[i]]
            References<-c("*",c(rep(" ",length(info.covariate[[2]][[i]])-1)))
            if(length(info.covariate[[4]][[i]])>1) {
                    charc1<-info.covariate[[4]][[i]][1]
                   for (ij in 2: length(info.covariate[[4]][[i]]))
                   {
                      
                      charc<-paste(charc1,info.covariate[[4]][[i]][ij],sep=",")
                      charc1<-charc
                   }
              Parameters_associated<-c(paste(charc1),rep(".",length(info.covariate[[2]][[i]])-1))
            } else 
            {
            Parameters_associated<-c(paste(info.covariate[[4]][[i]]),rep(".",length(info.covariate[[2]][[i]])-1))
            }
    
            e.cov[[i]]<-data.frame(Categories,References,Proportions,row.names=prem_col,check.names=F,check.rows=F)
          cat("\t",'Covariate',i,':',Name[1],'(',Parameters_associated[1],')',"\n")
                print(e.cov[[i]])
          cat("\n","\n")
          }
          
  }
  
  if (covariate_occ.model==T)
	{  
     
     cat("\n","\n")
     cat("\t") 
     cat("Covariates changing with occasion","\n","\n")
     e.cov1<-list()
     e.cov2<-list()
     
          for (k in 1: length(covariate_occ.name)){
          cat("\t")
          cat ("Covariate ",k,":",covariate_occ.name[[k]],"(",parameter_occ.associated[[k]],")","\n")
            
           lcov1<-length(covariate_occ.category[[k]])
           prem_col1<-lapply(1:lcov1,function(lcov1) paste("(",c(lcov1),")",sep=""))
           Categories<-covariate_occ.category[[k]]
           References<-c("*",c(rep(" ",length(covariate_occ.category[[k]])-1)))
           
           e.cov1[[k]]<-data.frame(Categories,References,row.names=prem_col1,check.names=F,check.rows=F)
           print(e.cov1[[k]])
           cat("\n","\n")
        
          lcov2<-length(sequence.proportions[[k]])
          prem_col2<-lapply(1:lcov2,function(lcov2) paste("(",c(lcov2),")",sep=""))
          Sequences<- list()

           for (sequ in 1:nb_seq[k]){
            Sequences[sequ]<-covariate_occ.sequence[[k]][[sequ]][1]
            for (occ in 2:n_occ){
              Sequences[sequ]<-paste(Sequences[sequ],covariate_occ.sequence[[k]][[sequ]][occ])
            }
           }

          #for (sequ in 1:nb_seq){
          #  Sequences[sequ]<-unlist(lapply(2:n_occ-1, function(n_occ) paste(covariate_occ.sequence[[k]][[sequ]][n_occ],covariate_occ.sequence[[k]][[sequ]][n_occ+1]))) 
          #}
           Proportions<-sequence.proportions[[k]]
          e.cov2[[k]]<-data.frame(unlist(Sequences),Proportions,row.names=prem_col2,check.names=F,check.rows=F)
          #. <- covariate_occ.sequence[[k]]
         
          #e.cov2[[k]]<-cbind(e.cov2[[k]],.,)
          names(e.cov2[[k]])<-c("Sequences","Proportions")
          
          print(e.cov2[[k]])
           cat("\n","\n")
          }
    }
  

	}
	
	prot.optiround<-list()
	if(length(f$degenere)!=(length(c(unlist(prots)+1))))
	{
	     for (l in 1 :nr){
	         if(length(which(l2[[l]]==0))==0){
                prot.optiround[[l]]<-lapply(prot.opti[[l]],round,3)
           } else 
			     {
                prot.optiround[[l]]<-list()
			           for ( v in 1: length(l2[[l]])){
			               if (is.numeric(prot.opti[[l]][[v]]))
			                  prot.optiround[[l]][[v]]<-unlist(lapply(prot.opti[[l]][[v]],round,3))
                     else { prot.optiround[[l]][[v]] <- prot.opti[[l]][[v]]
                     if( v==length(l2[[l]])) prot.optiround[[l]][[v+1]]<-Inf
                     } 
                }

            } 
      }   
			if (subjects.input==1) {
				subjects<-round((subjects.opti*Ntot/(unlist(l_prote))),4)
				subjects.prop<-round(subjects.opti,4)
			} else {
				subjects.prop<-round(subjects.opti,4)
				subjects<-round((subjects.opti*Ntot/(unlist(l_prote))),4)
			}



    	
		se<-f$eval[[5]]
		cv<-f$eval[[6]]
		pfixe<-f$eval[[13]]
    sh<-f$eval[[14]]		
		mfisher<-f$eval[[1]] 
		determinant<-f$eval[[7]] 
		crit<-f$eval[[8]]
		pp<-f$eval[[10]]

		StdError<-se[1:(pfixe+locc)] 
	  RSE<-cv[1:(pfixe+locc)] 
		
		if (FIM=="B") {
      beta.estim<-!beta.fixed[diag(omega)!=0] 
      theta<-theta[diag(omega)!=0][beta.estim]
      beta<-beta[diag(omega)!=0][beta.estim]
      parameters.estim<-parameters[diag(omega)!=0][beta.estim]
  
  }
  else{
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
	EVmin1<-min(EigenValues[1:(pfixe+locc)] )
	EVmax1<-max(EigenValues[1:(pfixe+locc)] )
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
		
		
		
	
	if (covariate.model==T) 
  {
  Beta<-c(theta,unlist(beta.covariate),beta.occ1) 
  if (n_occ==0) a<-data.frame(Beta,StdError,RSE,row.names=c(parameters.estim,npar)) 
  if (n_occ>0) a<-data.frame(Beta,StdError,RSE,row.names=c(parameters.estim,npar,no))  
  }else{
  Beta<-c(theta,beta.occ1)
    
  if (n_occ==0) a<-data.frame(Beta,StdError,RSE,row.names=c(parameters.estim)) 
  if (n_occ>0) a<-data.frame(Beta,StdError,RSE,row.names=c(parameters.estim,no)) 
	}
	
  .<-c(rep("%",pfixe+locc)) 	
	if (FIM=="B"){a<-cbind(a,.,Shrinkage=sh,.); names(a)[c(dim(a)[2]-2,dim(a)[2])]<-" "}  
  else {a<-cbind(a,.)
  names(a)[dim(a)[2]]<-" "} 
	
	
 
	b_condition<-F
	if (length(which(omega==0))!=length(omega)) {
	Omega<-diag(as.matrix(omega1)) 
	StdError<-se[(pfixe+locc+1):(pfixe+locc+length(Omega))] 
	RSE<-cv[(pfixe+locc+1):(pfixe+locc+length(Omega))]
	b<-data.frame(omega2=Omega,StdError,RSE,row.names=parameters[diag(omega)!=0]) 
	.<-c(rep("%",(length(Omega)))) 
	b<-cbind(b,.)
  names(b)[dim(b)[2]]<-" " 
	b_condition<-T
	}
	
	
	if (n_occ>1 && length(which(gamma==0))!=length(gamma)){
	Gamma<-diag(as.matrix(gamma1))
	StdError<-se[(pfixe+locc+length(Omega)+1):(pfixe+locc+length(Omega)+length(Gamma))] 
	RSE<-cv[(pfixe+locc+length(Omega)+1):(pfixe+locc+length(Omega)+length(Gamma))] 
	b2_condition<-F
	if(length(which(vec==F))!=length(vec)) {
	b2<-data.frame(gamma2=Gamma,StdError,RSE,row.names=parameters[diag(gamma)!=0]) 
	.<-c(rep("%",(length(Gamma)))) 
	b2<-cbind(b2,.) 
	names(b2)[dim(b2)[2]]<-" "
	b2_condition<-T
	}
	}
	if(n_occ==1 || length(which(gamma==0))!=length(gamma)){b2_condition<-F}
	
	if (FIM!="B"){
	StdError<-se[(pp-(lSIG-1)):pp]
	RSE<-cv[(pp-(lSIG-1)):pp]
	.<-rep("%",lSIG) 
	#sig.names<-c('sig.inter','sig.slope','sig.interPD','sig.slopePD')[vec2] 
	sig.names<-unlist(c(f$eval[[12]]))[vec2]
	c<-cbind(data.frame(SIG,StdError,RSE,row.names=sig.names),.)
	names(c)[c(1,dim(c)[2])]<-c("Sigma"," ")
	}
	
	
	
	
	
		if(modelform=="DE")  {
	    	  pm<-0
	         ll<-lapply(prot.optiround,length)
		      if(length(which(unlist(lapply(1:nr,function(nr) lapply(1:ll[[nr]], function(ll,prot.optiround,nr) length(prot.optiround[[nr]][[ll]]),prot.optiround=prot.optiround,nr=nr)))==0))==0)
        {
              f<-list()
              for (l in 1 : nr) {	 
                 f[[l]]<-data.frame(times=as.character(prot.optiround[[l]][1:Q]),subjects.prop,subjects)
              }
              pm<-1
        }
	  } else 
		{
     pm<-0
	         ll<-lapply(prot.optiround,length)
       if(length(which(unlist(lapply(1:nr,function(nr) lapply(1:ll[[nr]], function(ll,prot.optiround,nr) length(prot.optiround[[nr]][[ll]]),prot.optiround=prot.optiround,nr=nr)))==0))==0)
        {
        f<-list()
        for (l in 1 : nr) {	 
        f[[l]]<-data.frame(times=as.character(prot.optiround[[l]][1:Q]),subjects.prop,subjects)
         }
        pm<-1
        }
    }
		cat("**************************** OPTIMISED DESIGN *****************************") 
		cat("\n","\n") 
		cat("Number of iterations:",num.iter)
		cat("\n") 
		cat("Number of function evaluations:",count)
		cat("\n") 
		if (conv) cat("Convergence Achieved") else cat("False Convergence")
		cat("\n","\n") 
		cat("\n","\n") 
    cat('Design: ','\n')
		cat("\n","\n")
      for (i in 1:nr){
          cat('Sample times for response:',LETTERS[i],'\n')
          if (pm==1){
          print(f[[i]])
          cat("\n","\n") } else{
            print(prot.optiround[[i]][1:Q])
            cat('      Number of subjects per group \n')
            cat('     Proportions                 Numbers \n' )
            cat(t(t(subjects.prop)),'      ',  t(t(subjects)),'\n')
            cat("\n","\n")
          }
      }  
      cat("\n","\n")
		cat('Associated optimised criterion:', round(crit,4))  
		cat("\n","\n")
 	  if (nr!=1) cat('Identical sampling times for each response: ', identical.times)   
		cat("\n","\n") 
	
if (FIM=="I") cat("Computation of the Individual Fisher information matrix") 
if (FIM=="B") cat("Computation of the Bayesian Fisher information matrix") 
if (FIM=="P") cat("Computation of the Population Fisher information matrix: option = ",option)
cat("\n","\n") 
	
 #taking into account previous information
if (previous.FIM!=  "")
{	
#print(paste("Previous Fisher information matrix from file",previous.FIM)  )
cat('Previous FIM from file',previous.FIM)
cat("\n","\n")  }


    # Ecriture de la MF dans un fichier :	
    if (outputFIM !=  ""){
      if (ncol(mfisher) >= 4)
      {sommecomplete<-mfisher
      first.line<-c("# Project: ",project,                        '        Date: ',date(),  rep(" ",(ncol(mfisher)-4)))  }
      else
      {
        sommecomplete<-matrix("",nrow=nrow(mfisher),ncol=(4-ncol(mfisher)))
        sommecomplete<-cbind(mfisher,sommecomplete)
        first.line<-c("# Project: ",project,                        '        Date: ',date() ) }

      sommeoutput<-rbind(first.line,sommecomplete) 
      write.table(sommeoutput,file.path(directory,outputFIM),quote=F,append=F,row.names=F,col.names=F)
    }

	
if (outputFIM!=  "")
#print(paste('FIM saved in',outputFIM))
{
cat('FIM saved in',outputFIM)

cat("\n","\n","\n") 
}
		
		
  cat("******************* FISHER INFORMATION MATRIX ******************") 
  cat("\n","\n") 
	print(unclass(mfisher)) 
	cat("\n","\n")
		

	
		cat("\n","\n") 
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
  if (n_occ>1 && length(which(gamma==0))!=length(gamma)){
  cat("------------------------- Variance of Inter-Occasion Random Effects ----------------------") 
	cat("\n","\n") 
	print(b2) 
	cat("\n","\n")
	} 
  if (FIM!="B"){
  cat("------------------------ Standard deviation of residual error ------------------------ ")
  cat("\n","\n")
	print(c) 
	cat("\n","\n")}  
 
		cat("******************************* DETERMINANT ********************************") 
		cat("\n","\n") 
		cat(determinant) 
		cat("\n","\n") 
		cat("******************************** CRITERION *********************************") 
		cat("\n","\n") 
		cat(crit) 
		cat("\n","\n") 
		
	summary.ci<-NULL
	summary.power_nni<-NULL
	summary.ci_eq<-NULL
	summary.power_nni_eq<-NULL
	if (covariate.model==T || covariate_occ.model==T)
	{   	
  if (covariate.model==F)
	   {    
      beta.covariate<-NULL
      npar<-NULL
      }
    	if(compute.power==T || compute.nni==T)
    	{
    	cat("***************************** COMPARISON TEST ********************************") 
    	cat("\n","\n") 
    	summary.ci<-CI.computation(0.05/2,c(unlist(beta.covariate),beta.occ1),se[(p+1):(pfixe+locc)],c(npar,no),Trand)
    	summary.power_nni<-power_nni.computation(alpha/2,given.power,subjects,c(unlist(beta.covariate),beta.occ1),se[(p+1):(pfixe+locc)],c(npar,no),Trand)
    	print(summary.ci)
    	cat("\n","\n")
   	  cat("Type I error =",alpha)
      cat("\n","\n")
      print(summary.power_nni)                                         
    	}
    	cat("\n","\n") 
    	

    	cat("\n","\n") 
    	if(compute.power_eq==T || compute.nni_eq==T)
    	{
    	cat("***************************** EQUIVALENCE TEST ********************************") 
    	cat("\n","\n") 
    	summary.ci_eq<-CI.computation(0.05,c(unlist(beta.covariate),beta.occ1),se[(p+1):(pfixe+locc)],c(npar,no),Trand)
    	summary.power_nni_eq<-power_nni_eq.computation(alpha,given.power,subjects,c(unlist(beta.covariate),beta.occ1),se[(p+1):(pfixe+locc)],c(npar,no),Trand)
    	print(summary.ci_eq)
    	cat("\n","\n")
   	  cat("Type I error =",alpha)
   	  cat("\n")
 	  	cat("Equivalence interval = [log(",paste(exp(interval_eq[1])),"),log(",paste(exp(interval_eq[2])),")]",sep="")
      cat("\n","\n")
      print(summary.power_nni_eq)
      if (length(which(c(unlist(beta.covariate),beta.occ1) < interval_eq[1])) + length(which(c(unlist(beta.covariate),beta.occ1)> interval_eq[2])) > 0)
      cat("\n","NA: Unavailable because beta is not in the equivalence interval (not under H1)")                                          
    	}
    
    	cat("\n","\n") 
    	
	}
	
  cat("\n","\n") 
	cat("******************* EIGENVALUES OF THE FISHER INFORMATION MATRIX ******************")
	cat("\n","\n") 
	print(EV)
	cat("\n","\n") 
  cat("******************* CORRELATION MATRIX ******************")
	cat("\n","\n") 
	print(unclass(corr.mat))
	cat("\n") 
	
	  
	d1<-Sys.time()
	g1<-proc.time()
  print(d1-d) 
  print(g1[2]-g[2])
	sink() 
	#if(graph.logical==T) {graph(prots,prot.opti,graphinf,graphsup,names.datax,names.datay,y.range)} else NULL 


  
#graphs after optimisation
###########################

	if(graph.logical==T) {
    if (modelform == "AF" && NUM==T) {graphAFNopti(prots,prot.opti,graphinf,graphsup,names.datax,names.datay,y.range)}
    if (modelform == "AF" && NUM==F) {graphAFSopti(prots,prot.opti,graphinf,graphsup,names.datax,names.datay,y.range)}
    if (modelform == "DE" ) {graphDEopti(prots,prot.opti,graphinf,graphsup,names.datax,names.datay,y.range)}
  } else NULL 
	if(graphsensi.logical==T) {
    if (modelform == "AF" && NUM==T){source(file.path(directory.program,"Pfim4.0op1.r"));graphsensiAFN(form,prots,graphinf,graphsup,names.datax,names.datay,y.range)} 
    if (modelform == "AF" && NUM==F){source(file.path(directory.program,"Pfim4.0op1.r"));graphsensiAFS(formg_init,prots,graphinf,graphsup,names.datax,names.datay,y.range)}
    if (modelform == "DE"){source(file.path(directory.program,"EQPfim4.0op1.r")); graphsensi(formED,prots,graphinf,graphsup,names.datax,names.datay,y.range)}
  } else NULL 

		}
	   else
		{
		cat("Error: the initial population design and all the vertices are very poor (criterion = 0); the Simplex algorithm cannot converge.")
		sink()
	  }
	cat('OK',"\n")
  if (covariate.model==T)
	   return(list(prot.opti=prot.opti,subjects.opti=subjects.opti,mfisher=mfisher,determinant=determinant,crit=crit,se=se,cv=cv,summary.ci=summary.ci,summary.power_nni=summary.power_nni,summary.ci_eq=summary.ci_eq,summary.power_nni_eq=summary.power_nni_eq,EigenValues=EigenValues))
  else{
  	if (FIM!="B") return(list(prot.opti=prot.opti,subjects.opti=subjects.opti,mfisher=mfisher,determinant=determinant,crit=crit,se=se,cv=cv,EigenValues=EigenValues))
    else return(list(prot.opti=prot.opti,subjects.opti=subjects.opti,mfisher=mfisher,determinant=determinant,crit=crit,se=se,cv=cv,shrinkage=sh,EigenValues=EigenValues))
    }
	
}
}


ccccc <- list()

#add of a Splus function for the Simplex (taken on the Splus User Group) 
# --------------------------------------------------------------------

## ~/rtns/S/optfcn.q
## Splus functions useful in numerical optimization.
##   Daniel Heitjan, 13 September 1990
##   Revised, 10 March 1992
##   Style changes, 94.02.26
##
fun.amoeba<-function(p,y,ftol,itmax,funk,data,vrbs=F){
## Multidimensional minimization of the function funk(x,data) where x
## is an ndim-dimensional vector and data is a fixed set of data, by
## the downhill simplex method of Nelder and Mead.  The structure data
## is arbitrary and is evaluated only in funk.  Input is a matrix p
## whose ndim+1 rows are ndim-dimensional vectors which are the ## vertices of the starting simplex, and a data vector in a suitable
## format.  Also input is the vector y of length ndim+1, whose
## components must be pre-initialized to the values of funk evaluated
## at the ndim+1 vertices (rows) of p  and ftol the fractional
## convergence tolerance to be achieved in the function value (n.b.!).
## The output list will contain components p, ndim+1 new points all
## within ftol of a minimum function value, y, the function value,
## iter, the number of iterations taken, and converge, a convergence
## indicator.
##   Translated from *Numerical Recipes*, 13 March 1989
##   Revised for model fitting, 16 March 1989
##   Revised to handle univariate parameters, 29 July 1989
##   Comments and printing revised, 28 April 1990
##   Stopping criterion revised, 3 May 1990
##   Error in contraction routine fixed, 4 May 1990
##   Replaced 'order' with 'sort.list', 4 May 1990
##   Modified digits to print, 20 June 1990
##   Changed print to cat, 13 September 1990
##   Modify printing, 8 February 1991
##   Added verbose mode, 94.12.25
##   Daniel F. Heitjan
##
 if (vrbs) cat('entering fun.amoeba\n')
## Three parameters which define the expansions and contractions.
 alpha<-1.0
 beta<-0.5
 gamma<-2.0
## A parameter that governs stopping if the function is converging to
## zero.
 eps<-1.e-10
##
 mpts<-nrow(p)
 iter<-0
 contin<-T
 converge<-F
 bbbb <<- c()
 jjj <<- 1
 while (contin) {
## First we must determine which point is the highest (worst),
## next highest, and lowest (best).
  ord<-sort.list(y)
  ilo<-ord[1]
  ihi<-ord[mpts]
  inhi<-ord[mpts-1]
  bbbb <<- c(bbbb, y[ilo])
## Compute the fractional range from highest to lowest and return if
## satisfactory.
  rtol<-2.*abs(y[ihi]-y[ilo])/(abs(y[ihi])+abs(y[ilo])+eps)
  if ((rtol<ftol)||(iter==itmax)) { 
   contin<-F
   converge<-T
   if (iter==itmax) { converge<-F }
  } else { 
   if (iter==itmax) cat('Amoeba exceeding maximum iterations.\n')
   iter<-iter+1
   if (vrbs) cat('\n iteration',iter,'\n')
##
## Begin a new iteration.  Compute the vector average of all points
## except the highest, i.e. the center of the face of the simplex across
## from the high point.  We will subsequently explore along the ray from
## the high point through that center.
   pbar<-matrix(p[-ihi,],(mpts-1),(mpts-1))
   pbar<-apply(pbar,2,mean)
##
## Extrapolate by a factor alpha through the face, i.e. reflect the
## simplex from the high point.  Evaluate the function at the reflected
## point.
   pr<-(1+alpha)*pbar-alpha*p[ihi,]
   ypr<-funk(pr,data)
   print(ypr)    
if (ypr<=y[ilo]) { 
## Gives a result better than the best point, so try an additional
## extrapolation by a factor gamma, and check out the function there.
    prr<-gamma*pr+(1-gamma)*pbar
    yprr<-funk(prr,data)
    if (yprr<y[ilo]) { 
## The additional extrapolation succeeded, and replaces the highest point.
     p[ihi,]<-prr
     y[ihi]<-yprr
     if (vrbs) {
      cat('new p by extrapolated good flip\n')  
print(prr)
      cat('new y\n')  
print(yprr) }
    } else { 
## The additional extrapolation failed, but we can still use the
## reflected point.
     p[ihi,]<-pr
     y[ihi]<-ypr
     if (vrbs) {
      cat('new p by good flip\n')  
print(pr)
      cat('new y\n')  
print(ypr) }
    }
   } else { 
    if (ypr>=y[inhi]) { 
## The reflected point is worse than the second-highest.
     if (ypr<y[ihi]) { 
## If it's better than the highest, then replace the highest,
      p[ihi,]<-pr
      y[ihi]<-ypr
      if (vrbs) {
       cat('new p by bad flip\n')  
print(pr)
       cat('new y\n')  
print(ypr) }
     }
## but look for an intermediate lower point, in other words, perform a
## contraction of the simplex along one dimension.  Then evaluate the
## function.
     prr<-beta*p[ihi,]+(1-beta)*pbar
     yprr<-funk(prr,data)
     if (yprr<y[ihi]) { 
## Contraction gives an improvement, so accept it.
      p[ihi,]<-prr
      y[ihi]<-yprr
      if (vrbs) {
       cat('new p by contraction from p[ihi,]\n')  
print(prr)
       cat('new y\n')  
print(yprr) }
     } else { 
## Can't seem to get rid of that high point.  Better contract around the
## lowest (best) point.
      p<-0.5*(p+matrix(p[ilo,],nrow=nrow(p),ncol=ncol(p),byrow=T))
      for (j in 1:length(y)) {
       if (j!=ilo) {
        y[j]<-funk(p[j,],data)
       }
      }
      if (vrbs) {
       cat('new p by contraction around p[ilo,]\n')  
print(p)
       cat('new y\n')  
print(y) }
     }
    } else { 
## We arrive here if the original reflection gives a middling point.
## Replace the old high point and continue
     p[ihi,]<-pr
     y[ihi]<-ypr
     if (vrbs) {
      cat('new p by middling flip\n')  
print(pr)
      cat('new y\n')  
print(ypr) }
    }
   }
  }
  ccccc[[jjj]] <<-  data.frame(p=p[ihi,],
                             y=rep(1/y[ihi], length(p[ihi,])),
                             iter=rep(iter, length(p[ihi,])),
                             converge=rep(converge, length(p[ihi,])))
  jjj<<- jjj+1
 }
 if (vrbs) cat('  leaving fun.amoeba\n')
 list(p=p,y=y,iter=iter,converge=converge)
}
