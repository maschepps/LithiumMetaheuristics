#PFIM 4.0 
#April 2014
#Copyright (C) PFIM 4.0 - Universite Paris Diderot and INSERM.
#------------------------------------------------------------------------------------  

############################################################################
####                Data, model and first derivatives for all times     ####
############################################################################


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

#if fixed.times are missing
err13<-tryCatch(fixed.times_test<-lapply(1:nr,function(nr,x) get(x), x=paste("fixed.times",LETTERS[nr],sep="")), error=function(e) "T")
if(length(err13)==1 && err13=="T") for (i in 1:nr) {x=paste("fixed.times",LETTERS[i],sep="");assign(x,c())}  


###################################################################
#Function associated to covariate model
###################################################################
#--------------------------------------------------------------------------------
#association of the covariate categories to the corresponding proportions in a same table 
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


########################################
#Create table with times, doses, form and response for each time considered
create.timdos<-function(model.info,struc.info,sampl.info,init.info) {
#Needed to test which form should be used for each sampling time
  #comp.bnd<-function(t,boundaries) {
  #   max(which(unlist(lapply(1:length(boundaries),
  #   function(i) ((boundaries[[i]][1]-t)*(boundaries[[i]][2]-t)<=0)))))
  #}
  
  comp.bnd<-function(t,boundaries) {
    is.lower<-unlist(lapply(boundaries,function(vec) (vec[2]>=t)))
    is.higher<-unlist(lapply(boundaries,function(vec) (vec[1]<=t)))
    max(which((is.lower*is.higher)==1))
  } 
  
  
  nr<-model.info$nr
  tab<-NULL
  if(struc.info$modelform=="AF" ) {
    wlform<-0
    for(l in 1:nr) {
      wsamp<-sampl.info$sampwin[[l]]
      boundary<-struc.info$bound[[l]] #ECO TODO : v?rifier si ?a marche
      times<-sort(unique(unlist(wsamp)))
      tim.dos<-data.frame(times=times,doses=rep(init.info$doses[1],length(times)),
      condnum=rep(1,length(times)),irep=rep(l,length(times)))
      if (NUM==F){ifor<-unlist(lapply(tim.dos[,1],comp.bnd,boundaries=boundary))+wlform
      tim.dos<-cbind(tim.dos,iform=ifor)}
#eco TODO : sort this...
      if (NUM==F){
      wlform<-wlform+struc.info$lformg[[l]] }
      tab<-rbind(tab,tim.dos)
     
      }
    tim.dos<-tab
    if(!init.info$dose.identical & length(init.info$doses)>1) {
      tab1<-tim.dos;
      nl<-dim(tim.dos)[1]
      for(idos in 2:length(init.info$doses)) {
        tab1$doses<-rep(init.info$doses[idos],nl)
        tab1$condnum<-rep(idos,nl)
        #tab1$iform<-tab1$iform+(idos-1)*(struc.info$lf-1)
        if (NUM==F)#analytical form with standard way of fonction writing function() & numerical derivatives
        {tab1$iform<-tab1$iform+struc.info$lf}
        tim.dos<-rbind(tim.dos,tab1)
      }
    }
    tim.dos<-cbind(id=paste(tim.dos$times,tim.dos$irep,tim.dos$condnum,sep="x"),
      tim.dos)
  } else { #model is DE
    for(l in 1:nr) {
       wsamp<-sampl.info$sampwin[[l]]
       times<-sort(unique(unlist(wsamp)))
# eco : condnum contains the condinit group (1,..length(condinit))
       cond.number<-rep(1,length(times))
       
       tim.dos<-data.frame(times=times,condnum=cond.number,
       irep=rep(l,length(times))) #,time.init=time.init)         
    tab<-rbind(tab,tim.dos)
    }
    tim.dos<-tab
    if(!init.info$condinit.identical & length(init.info$condinit)>1) {
      tab1<-tim.dos;nl<-dim(tim.dos)[1]
      for(idos in 2:length(init.info$condinit)) {
        tab1$condnum<-rep(idos,nl)
        tim.dos<-rbind(tim.dos,tab1)
      }
    }
    tim.dos<-cbind(id=paste(tim.dos$times,tim.dos$irep,tim.dos$condnum,sep="x"),
      tim.dos)
  }
  sig.i<-unlist(lapply(1:dim(tim.dos)[1],function(i) model.info$sigmainter[[tim.dos$irep[i]]]))
  sig.s<-unlist(lapply(1:dim(tim.dos)[1],function(i) model.info$sigmaslope[[tim.dos$irep[i]]]))
  tim.dos<-cbind(tim.dos,sig.i=sig.i,sig.s=sig.s)

  return(tim.dos)
}

########################################
# Derivation symbolique des d?riv?es premi?res et secondes
#Speeden up the computation by only getting the symbolic derivatives of
# the forms we do need

symbolic.prem<-function(model.info,struc.info,kform,option) {
# Computation of 1st and 2nd order symbolic derivatives 

  p<-length(model.info$beta)
  for (i in 1:p) 
    assign(parameters[i],model.info$beta[i])
      #Deriv?e premi?re
      # liste de longueur Sum_l (n_l)
      # chaque ?l?ment de la liste est un vecteur de taille p (nb param?tres)
        df.forms<-vector(length(struc.info$form),mode="list")
        df.formsb<-vector(length(struc.info$form),mode="list")
        ddf.forms<-list()
        for(l in kform)
         {
          formi<-struc.info$form[[l]]
          df.formsb[[l]]<-lapply(1:p,function(i,formi) D(formi,parameters[i]),formi=formi)
          df.forms[[l]]<-unlist(df.formsb[[l]])
          #deriv?e seconde  lapply(1:p,function(p,Dforma,parameters) lapply(1:p1, function(p1,Dforma,parameters) D(Dforma,parameters[p1]),Dforma=Dforma[[p]],parameters=parameters),Dforma=Dforma[[lf]],parameters=parameters)
          if (option==2){
          p1<-p
          ddf.forms[[l]]<-lapply(1:p,function(p,df.formsb,parameters) lapply(1:p1, function (p1,df.formsb,parameters) D(df.formsb,parameters[p1]),df.formsb=df.formsb[[p]],parameters=parameters),parameters=parameters,df.formsb=df.formsb[[l]]) 
          }
      
       }
   
  return(list(df.forms=df.forms,ddf.forms=ddf.forms))
}



calc.modelprem<-function(tim.dos,model.info,struc.info,df.forms,ddf.forms,option,info.covariate,thetah,beta.occh) {
  p<-length(model.info$beta)
  theta_init<-thetah
  theta<-thetah
  for (i in 1:p) assign(parameters[i],thetah[i])
  if(info.covariate$covariate.model==T)  #calcul pour chacune des possiblit?s de covariables
  {
    f.p<-list()
    df.p<-list()
    df.pp<-list()
    ddf.p<-list()
    parameters.cov.val<-rep(1,p)
    pfixe<-length(c(parameters,info.covariate$tab.cov.param[,3]))
    n.theta<-c(thetah,as.numeric(unlist(info.covariate$tab.cov.param[,5]))) 
    sum<-0
    for (i in 1:pfixe) {assign(c(parameters,info.covariate$tab.cov.param[,3])[i],n.theta[i])} 
    for (ip in 1:p) {assign(parameters[ip],theta[ip])} 
    
    #calcul des derivees pour chaque combinaison de modelites des covariables
    for (ll in 1:dim(info.covariate$covariate.cat.poss)[1])
    {
        categories.cov<-c(info.covariate$covariate.cat.poss[ll,])
        categories.prop<-prod(info.covariate$tab_mod_prop[[2]][ll,])
        parameters.cov.val<-rep(1,p)
	   #attibuer les valeurs des modalit?s aux covariables
	   for (i in 1:length(info.covariate$covariate.name)) {assign(info.covariate$covariate.name[[i]],categories.cov[i])} #assign("trt",groupe)
     
          #calcul des exp(beta*cat) ou beta+cat selon les trand
          for (i in 1: length(unique(info.covariate$tab.cov.param[,1])))
      	   {
            	 #par.test<-tab.cov.param[i,1][
            	 par.test<-unique(info.covariate$tab.cov.param[,1])[i]
            	 tab.cov.param1<-info.covariate$tab.cov.param[which(info.covariate$tab.cov.param[,1]==par.test),]
            	 tab.cov.param1<-matrix(tab.cov.param1,ncol=5)
      	 
      	        if (length(unique(tab.cov.param1[,2]))==1)
               {
                                 test.bet<-tab.cov.param1[get(tab.cov.param1[1,2])==tab.cov.param1[,4],3]
                  if (length(test.bet)!=0)  #si categories ==0
                  {
                      if (Trand==2)   
                      {
                        
                          theta[which(tab.cov.param1[1,1]==parameters)]<-theta[which(tab.cov.param1[1,1]==parameters)]*exp(get(tab.cov.param1[get(tab.cov.param1[1,2])==tab.cov.param1[,4],3]))
                          parameters.cov.val[which(tab.cov.param1[1,1]==parameters)]<-exp(get(tab.cov.param1[get(tab.cov.param1[1,2])==tab.cov.param1[,4],3]))
                          
                      } else 
                      {
                      theta[which(tab.cov.param1[1,1]==parameters)]<-theta[which(tab.cov.param1[1,1]==parameters)]+ get(tab.cov.param1[get(tab.cov.param1[1,2])==tab.cov.param1[,4],3])  
    
                      }
                    }
                   }else {
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
                                   
                                    theta[which(tab.cov.param2[1,1]==parameters)]<-theta[which(tab.cov.param2[1,1]==parameters)]*exp(get(tab.cov.param2[get(tab.cov.param2[1,2])==tab.cov.param2[,4],3]))
                                    parameters.cov.val[which(tab.cov.param2[1,1]==parameters)]<-exp(get(tab.cov.param2[get(tab.cov.param2[1,2])==tab.cov.param2[,4],3]))
                                    
                                } else {
                                
                                    theta[which(tab.cov.param[i,1]==parameters)]<-theta[which(tab.cov.param[i,1]==parameters)]+ get(tab.cov.param2[get(tab.cov.param2[1,2])==tab.cov.param2[,4],3])  
  
                                }
                              }
                          }
                      }
              }   
        
  
         for (ip in 1:p) {assign(parameters[ip],theta[ip])}
         f.p[[ll]]<-unlist(lapply(1:dim(tim.dos)[1],function(i,tab,formi) {
         t<-tab$times[i]
         dose<-tab$doses[i]
         ifo<-tab$iform[i]
         return(eval(formi[[ifo]])) },tab=tim.dos,formi=struc.info$form))
         f.p[[ll]]<-matrix(f.p[[ll]],nrow=1)
         df.p[[ll]]<-matrix(nrow=p,ncol=dim(tim.dos)[1]) #without parameter covariate
         df.pp[[ll]]<-matrix(nrow=p,ncol=dim(tim.dos)[1]) #without parameter covariate
         ddf.p[[ll]]<-list()  
          
         for(i in 1:dim(tim.dos)[1]) {
            t<-tim.dos$times[i]
            dose<-tim.dos$doses[i]
            ifo<-tim.dos$iform[i]
            df.i<-df.forms[[ifo]]
            if (option==2) ddf.i<-ddf.forms[[ifo]]
            
              vec<-unlist(lapply(1:p,function(i) eval(df.i[[i]])))
              
              if(info.covariate$covariate.model==T) df.pp[[ll]][,i]<-vec*parameters.cov.val
              
              if(model.info$Trand==2) vec<-vec*model.info$beta*exp(beta.occh)
              
              if(info.covariate$covariate.model==T) df.p[[ll]][,i]<-vec*parameters.cov.val else df.p[,i]<-vec
               
            if(option==2) {
                p1<-p
                vec2<-lapply(1:p,function(p,ddf.i) lapply(1:p1,function(p1,ddf.i) eval(ddf.i[[p1]]),ddf.i=ddf.i[[p]]),ddf.i=ddf.i)
                for (j in 1:p){
                  if(i==1) ddf.p[[ll]][[j]]<-matrix(nrow=p,ncol=dim(tim.dos)[1])
                  vec3<-unlist(vec2[[j]])
                  if (model.info$Trand==2) {
                    ddfbeta.p<-vec3*thetah
                    ddfbeta.p[j]<-ddfbeta.p[j]+ df.p[j,i]/thetah[j] 
                    ddf.p[[ll]][[j]][,i]<-ddfbeta.p}
                  else {ddf.p[[ll]][[j]][,i]<-vec3}
                }     
            }
          }  
          #addition of derivatives for parameters covariates   
         for (icat in 1:dim(info.covariate$tab.cov.param)[1])
         {
                #if(model.info$Trand==2)
                 if(!info.covariate$tab.cov.param[icat,4]==get(info.covariate$tab.cov.param[icat,2])) s2<-df.p[[ll]][which(parameters==info.covariate$tab.cov.param[[icat]]),]*0 else  s2<-df.p[[ll]][which(parameters==info.covariate$tab.cov.param[[icat]]),]
                   df.p[[ll]]<-rbind(df.p[[ll]],s2)
                
         }  
      theta<-theta_init   
   }    
   
  }


 else {    
#model function
     f.p<-unlist(lapply(1:dim(tim.dos)[1],function(i,tab,formi) {
     t<-tab$times[i]
     dose<-tab$doses[i]
     ifo<-tab$iform[i]
     return(eval(formi[[ifo]])) },tab=tim.dos,formi=struc.info$form))
      f.p<-matrix(f.p,nrow=1)
#first derivatives
      df.p<-matrix(nrow=p,ncol=dim(tim.dos)[1])
      ddf.p<-list()
      for(i in 1:dim(tim.dos)[1]) {
        t<-tim.dos$times[i]
        dose<-tim.dos$doses[i]
        ifo<-tim.dos$iform[i]
        df.i<-df.forms[[ifo]]   #ifo=reponses
        if (option==2) ddf.i<-ddf.forms[[ifo]]
        vec<-unlist(lapply(1:p,function(i) eval(df.i[[i]])))
        if(model.info$Trand==2) vec<-vec*model.info$beta*exp(beta.occh)
        df.p[,i]<-vec
    #second derivatives
        if(option==2) {
        p1<-p
        vec2<-lapply(1:p,function(p,ddf.i) lapply(1:p1,function(p1,ddf.i) eval(ddf.i[[p1]]),ddf.i=ddf.i[[p]]),ddf.i=ddf.i)

      
        for (j in 1:p){                                                                                                                                         
          if(i==1) ddf.p[[j]]<-matrix(nrow=p,ncol=dim(tim.dos)[1])
          vec3<-unlist(vec2[[j]])
          if (model.info$Trand==2) {
           ddfbeta.p<-vec3*model.info$beta
           ddfbeta.p[j]<-ddfbeta.p[j]+ df.p[j,i]/thetah[j] 
           ddf.p[[j]][,i]<-ddfbeta.p}
          else {ddf.p[[j]][,i]<-vec3}
        }
        }
      }
 }

  return(list(mod=f.p,sensi=df.p,sensi2=ddf.p))
}


########################################
# Compute model predictions and first derivatives for all protocol times
# for a model described by differential equations

inter<-function(theta,model.info,kj,l,init.info,formED,condinit)  {
  p<-length(model.info$beta)
  for (i in 1:p) 
     assign(parameters[i],theta[i])
  cond<-eval(condinit[kj])
  res<-lsoda(cond,times=c(init.info$time.condinit,theta[p+1]),
    formED,theta,rtol=init.info$RtolEQ,
    atol=init.info$AtolEQ,hmax=init.info$Hmax)
  if (model.info$nr==1)
     res<-res[,dim(res)[2]][2]
  else 
     res<-res[,((dim(res)[2]-model.info$nr+1):(dim(res)[2]))][,l] [-1]
 return(res) 
}


calc.model.de<-function(tab,model.info,init.info,struc.info,option,info.covariate,thetah,beta.occh) {
  p<-length(model.info$beta)
  theta_init<-thetah
  theta<-thetah
  for (i in 1:p) assign(parameters[i],thetah[i])
  
   if(info.covariate$covariate.model==T)  #calcul pour chacune des possiblit?s de covariables
    {
       for (i in 1:p) assign(parameters[i],thetah[i])  
       df.p<-list()
       f.p<-NULL
       ddf.p<-list()
       df.pp<-list()
       l0<-0
       parameters_cov<-c(parameters,info.covariate$tab.cov.param[,3])
       theta_cov<-c(theta,as.numeric(unlist(info.covariate$tab.cov.param[,5])))
       p1<-length(theta_cov)
       pfixe<-length(c(parameters,info.covariate$tab.cov.param[,3]))
       n.theta<-c(thetah,as.numeric(unlist(info.covariate$tab.cov.param[,5]))) 
       sum<-0
       for (i in 1:pfixe) {assign(c(parameters,info.covariate$tab.cov.param[,3])[i],n.theta[i])} 
       for (ip in 1:p) {assign(parameters[ip],theta[ip])}
           
       for (ll in 1:dim(info.covariate$covariate.cat.poss)[1])
        {  f.p[[ll]]<-numeric()  
           df.p[[ll]]<-matrix(nrow=p,ncol=0)
           categories.cov<-c(info.covariate$covariate.cat.poss[ll,])
           categories.prop<-prod(info.covariate$tab_mod_prop[[2]][ll,])
           parameters.cov.val<-rep(1,p)
      	   #attibuer les valeurs des modalit?s aux covariables
      	   for (i in 1:length(info.covariate$covariate.name)) {assign(info.covariate$covariate.name[[i]],categories.cov[i])} #assign("trt",groupe)
           
          #calcul des exp(beta*cat) ou beta+cat selon les trand
            for (i in 1: length(unique(info.covariate$tab.cov.param[,1])))
      	   {
            	 #par.test<-tab.cov.param[i,1][
            	 par.test<-unique(info.covariate$tab.cov.param[,1])[i]
            	 tab.cov.param1<-info.covariate$tab.cov.param[which(info.covariate$tab.cov.param[,1]==par.test),]
            	 tab.cov.param1<-matrix(tab.cov.param1,ncol=5)
      	 
      	        if (length(unique(tab.cov.param1[,2]))==1)
               {
                                 test.bet<-tab.cov.param1[get(tab.cov.param1[1,2])==tab.cov.param1[,4],3]
                  if (length(test.bet)!=0)  #si categories ==0
                  {
                      if (Trand==2)   
                      {
                         
                          theta[which(tab.cov.param1[1,1]==parameters)]<-theta[which(tab.cov.param1[1,1]==parameters)]*exp(get(tab.cov.param1[get(tab.cov.param1[1,2])==tab.cov.param1[,4],3]))
                          parameters.cov.val[which(tab.cov.param1[1,1]==parameters)]<-exp(get(tab.cov.param1[get(tab.cov.param1[1,2])==tab.cov.param1[,4],3]))
                          
                      } else 
                      {
                      theta[which(tab.cov.param1[1,1]==parameters)]<-theta[which(tab.cov.param1[1,1]==parameters)]+ get(tab.cov.param1[get(tab.cov.param1[1,2])==tab.cov.param1[,4],3])    
                      }
                    }
                   }else {
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
                                   
                                    theta[which(tab.cov.param2[1,1]==parameters)]<-theta[which(tab.cov.param2[1,1]==parameters)]*exp(get(tab.cov.param2[get(tab.cov.param2[1,2])==tab.cov.param2[,4],3]))
                                    parameters.cov.val[which(tab.cov.param2[1,1]==parameters)]<-exp(get(tab.cov.param2[get(tab.cov.param2[1,2])==tab.cov.param2[,4],3]))
                                  
                                   
                                    
                                } else {
                                
                                    theta[which(tab.cov.param[i,1]==parameters)]<-theta[which(tab.cov.param[i,1]==parameters)]+ get(tab.cov.param2[get(tab.cov.param2[1,2])==tab.cov.param2[,4],3])  
  
                                }
                              }
                          }
                      }
              } 
            
            for (ip in 1:p) {assign(parameters[ip],theta[ip])} 
           
            for(kj in 1:length(unique(tab$condnum))) {
                
               cond<-eval(init.info$condinit[kj])
               fmod<-list()
               dd<-list()
               dd2<-list()
               s2dp<-list()
               for(l in 1:model.info$nr) {
                      
                      tim<-tab$times[tab$irep==l & tab$condnum==kj]
                      fmod[[l]]<-lsoda(cond,c(init.info$time.condinit,tim),struc.info$formED,
                       theta,rtol=init.info$RtolEQ,atol=init.info$AtolEQ,
                       hmax=init.info$Hmax)
                     if (model.info$nr==1) { #single response
                       fmod[[l]]<-fmod[[l]][,dim(fmod[[l]])[2]][-1] 
                     } else 
                       fmod[[l]]<-fmod[[l]][,
                         ((dim(fmod[[l]])[2]-model.info$nr+1):(dim(fmod[[l]])[2]))][,l][-1]
                       #if(l==1)
                       #{
                       #f.p[[ll]]<-fmod[[l]]
                       #}else{
                       #f.p[[ll]]<-c(f.p[[ll]],fmod[[l]])
                       #}
                       
                       #### 
                       f.p[[ll]]<-c(f.p[[ll]],fmod[[l]])
                       ####
                       
                       lfd<-length(tim)
                       dd[[l]]<-lapply(1:lfd,function(lfd,theta,tim,kj,t,l,model.info,init.info,
                       formED,condinit) fdHess(c(theta,tim[lfd]),inter,model.info=model.info,
                       kj=kj,l=l,init.info=init.info,formED=formED,condinit=condinit)$gradient,
                       theta=theta,tim=tim,kj=kj,l=l,model.info=model.info,
                       init.info=init.info,formED=struc.info$formED,
                       condinit=init.info$condinit)
                       df1<-t(matrix(unlist(dd[[l]]),ncol=length(tim)))[,1:p] 
                       if(length(tim)==1) df1<-matrix(df1,nrow=p,ncol=1)*parameters.cov.val else df1<-t(df1)*parameters.cov.val
                      
                      # if(l==1)
                       #{
                      #   df.pp[[ll]]<-df1
                      #   if (model.info$Trand==2) df1<-df1*model.info$beta*exp(beta.occh) 
                      #   df.p[[ll]]<-df1
                      # }else{
                      #   df.pp[[ll]]<-cbind(df.pp[[ll]],df1)
                      #   if (model.info$Trand==2) df1<-df1*model.info$beta*exp(beta.occh) 
                       #  df.p[[ll]]<-cbind(df.p[[ll]],df1)
                       #}
                       
                       #### 
                       if (model.info$Trand==2) df1<-df1*model.info$beta*exp(beta.occh) 
                       df.p[[ll]]<-cbind(df.p[[ll]],df1)
                       ####
                      
           
                       #matrice hessienne
                       if (option==2){
                         dd2[[l]]<-lapply(1:lfd,function(lfd,thetah,tim,kj,l,model.info,condinit,init.info,formED) hessian(inter,c(thetah,tim[lfd]),method="Richardson",kj=kj,l=l,formED=formED,model.info=model.info,condinit=condinit,init.info=init.info),thetah=thetah,tim=tim,kj=kj,l=l,model.info=model.info,condinit=init.info$condinit,init.info=init.info,formED=struc.info$formED)
                           
                           s2d<-NULL
                           s2dp[[l]]<-list()
                              for (j1 in 1:p)
                              {
                                   for (j in 1:length(tim))
                                    {
                                      if (model.info$Trand==2)
                                      {
                                        tran<-matrix(dd2[[l]][[j]][,j1][1:p]*thetah,ncol=1)
                                        tran[j1]<-tran[j1]+df.p[,which(tab$irep==l & tab$condnum==kj)][j1,j]/thetah[j1]
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
              }
          if (option==2){
              for (j4 in 1:p)
                 {
                    mo1<-NULL
                    for (j3 in 1:length(s2dp))
                    {
                      mo<-NULL
                      mo<-s2dp[[j3]][[j4]]
                       mo1<-cbind(mo1,mo)
                    }
                    ddf.p[[j4]]<-mo1
                 }    
          }
      
        for (icat in 1:dim(info.covariate$tab.cov.param)[1])
         {
                #if(model.info$Trand==2)
                 if(!info.covariate$tab.cov.param[icat,4]==get(info.covariate$tab.cov.param[icat,2])) s2<-df.p[[ll]][which(parameters==info.covariate$tab.cov.param[[icat]]),]*0 else  s2<-df.p[[ll]][which(parameters==info.covariate$tab.cov.param[[icat]]),]
                   df.p[[ll]]<-rbind(df.p[[ll]],s2)
                
         }  
      theta<-theta_init  
    }
    
   }else {
  
  #model without covaraite
     for (i in 1:p) 
      assign(parameters[i],thetah[i])
      #model function
       df.p<-NULL
       f.p<-NULL
       ddf.p<-list()
       l0<-0
       for(kj in 1:length(unique(tab$condnum))) {
         cond<-eval(init.info$condinit[kj])
         fmod<-list()
         dd<-list()
         dd2<-list()
         s2dp<-list()
         for(l in 1:model.info$nr) {
           tim<-tab$times[tab$irep==l & tab$condnum==kj]
           fmod[[l]]<-lsoda(cond,c(init.info$time.condinit,tim),struc.info$formED,
            thetah,rtol=init.info$RtolEQ,atol=init.info$AtolEQ,
             hmax=init.info$Hmax)
            
           if (model.info$nr==1) { #single response
             fmod[[l]]<-fmod[[l]][,dim(fmod[[l]])[2]][-1] 
           } else 
             fmod[[l]]<-fmod[[l]][,
               ((dim(fmod[[l]])[2]-model.info$nr+1):(dim(fmod[[l]])[2]))][,l][-1]
           f.p<-c(f.p,fmod[[l]])
           lfd<-length(tim)
           dd[[l]]<-lapply(1:lfd,function(lfd,thetah,tim,kj,t,l,model.info,init.info,
             formED,condinit) fdHess(c(thetah,tim[lfd]),inter,model.info=model.info,
             kj=kj,l=l,init.info=init.info,formED=formED,condinit=condinit)$gradient,
             thetah=thetah,tim=tim,kj=kj,l=l,model.info=model.info,
             init.info=init.info,formED=struc.info$formED,
             condinit=init.info$condinit)
           df1<-t(matrix(unlist(dd[[l]]),ncol=length(tim)))[,1:p] 
                  if(length(tim)==1) df1<-matrix(df1,nrow=p,ncol=1) else df1<-t(df1)
           if (model.info$Trand==2) df1<-df1*model.info$beta*exp(beta.occh) 
           df.p<-cbind(df.p,df1)
         
           #matrice hessienne
             if (option==2){
               dd2[[l]]<-lapply(1:lfd,function(lfd,thetah,tim,kj,l,model.info,condinit,init.info,formED) hessian(inter,c(thetah,tim[lfd]),method="Richardson",kj=kj,l=l,formED=formED,model.info=model.info,condinit=condinit,init.info=init.info),thetah=thetah,tim=tim,kj=kj,l=l,model.info=model.info,condinit=init.info$condinit,init.info=init.info,formED=struc.info$formED)
                 
                 s2d<-NULL
                 s2dp[[l]]<-list()
                    for (j1 in 1:p)
                    {
                         for (j in 1:length(tim))
                          {
                            if (model.info$Trand==2)
                            {
                              tran<-matrix(dd2[[l]][[j]][,j1][1:p]*thetah,ncol=1)
                              tran[j1]<-tran[j1]+df.p[,which(tab$irep==l & tab$condnum==kj)][j1,j]/thetah[j1]
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
            }
           if (option==2){
             for (j4 in 1:p)
             {
             mo1<-NULL
             for (j3 in 1:length(s2dp))
             {
             mo<-NULL
             mo<-s2dp[[j3]][[j4]]
             mo1<-cbind(mo1,mo)
             }
             ddf.p[[j4]]<-mo1
             }    
           }
   }
   
  return(list(mod=f.p,sensi=df.p,sensi2=ddf.p))
}


########################################
# Compute model predictions and numerical first derivatives for all protocol times
# for a user defined model described by analytical form
############################################################################

calc.model.num<-function(tim.dos,model.info,struc.info,option,info.covariate,thetah,beta.occh) {

  p<-length(model.info$beta)
  theta_init<-thetah
  theta<-thetah
  form<-struc.info$form
  for (i in 1:p) assign(parameters[i],thetah[i])
  
   if(info.covariate$covariate.model==T)  #calcul pour chacune des possiblit?s de covariables
    {
       for (i in 1:p) assign(parameters[i],thetah[i])  
       df.p<-list()
       f.p<-list()
       ddf.p<-list()
       df.pp<-list()
       l0<-0
       parameters_cov<-c(parameters,info.covariate$tab.cov.param[,3])
       theta_cov<-c(theta,as.numeric(unlist(info.covariate$tab.cov.param[,5])))
       p1<-length(theta_cov)
       pfixe<-length(c(parameters,info.covariate$tab.cov.param[,3]))
       n.theta<-c(thetah,as.numeric(unlist(info.covariate$tab.cov.param[,5]))) 
       sum<-0
       for (i in 1:pfixe) {assign(c(parameters,info.covariate$tab.cov.param[,3])[i],n.theta[i])} 
       for (ip in 1:p) {assign(parameters[ip],theta[ip])}
           
       for (ll in 1:dim(info.covariate$covariate.cat.poss)[1])
        {  f.p[[ll]]<-numeric()  
           df.p[[ll]]<-matrix(nrow=p,ncol=0)
           categories.cov<-c(info.covariate$covariate.cat.poss[ll,])
           categories.prop<-prod(info.covariate$tab_mod_prop[[2]][ll,])
           parameters.cov.val<-rep(1,p)
      	   #attibuer les valeurs des modalit?s aux covariables
      	   for (i in 1:length(info.covariate$covariate.name)) {assign(info.covariate$covariate.name[[i]],categories.cov[i])} #assign("trt",groupe)
           
          #calcul des exp(beta*cat) ou beta+cat selon les trand
            for (i in 1: length(unique(info.covariate$tab.cov.param[,1])))
      	   {
            	 #par.test<-tab.cov.param[i,1][
            	 par.test<-unique(info.covariate$tab.cov.param[,1])[i]
            	 tab.cov.param1<-info.covariate$tab.cov.param[which(info.covariate$tab.cov.param[,1]==par.test),]
            	 tab.cov.param1<-matrix(tab.cov.param1,ncol=5)
      	 
      	        if (length(unique(tab.cov.param1[,2]))==1)
               {
                                 test.bet<-tab.cov.param1[get(tab.cov.param1[1,2])==tab.cov.param1[,4],3]
                  if (length(test.bet)!=0)  #si categories ==0
                  {
                      if (Trand==2)   
                      {
                         
                          theta[which(tab.cov.param1[1,1]==parameters)]<-theta[which(tab.cov.param1[1,1]==parameters)]*exp(get(tab.cov.param1[get(tab.cov.param1[1,2])==tab.cov.param1[,4],3]))
                          parameters.cov.val[which(tab.cov.param1[1,1]==parameters)]<-exp(get(tab.cov.param1[get(tab.cov.param1[1,2])==tab.cov.param1[,4],3]))
                          
                      } else 
                      {
                      theta[which(tab.cov.param1[1,1]==parameters)]<-theta[which(tab.cov.param1[1,1]==parameters)]+ get(tab.cov.param1[get(tab.cov.param1[1,2])==tab.cov.param1[,4],3])    
                      }
                    }
                   }else {
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
                                   
                                    theta[which(tab.cov.param2[1,1]==parameters)]<-theta[which(tab.cov.param2[1,1]==parameters)]*exp(get(tab.cov.param2[get(tab.cov.param2[1,2])==tab.cov.param2[,4],3]))
                                    parameters.cov.val[which(tab.cov.param2[1,1]==parameters)]<-exp(get(tab.cov.param2[get(tab.cov.param2[1,2])==tab.cov.param2[,4],3]))
                                  
                                   
                                    
                                } else {
                                
                                    theta[which(tab.cov.param[i,1]==parameters)]<-theta[which(tab.cov.param[i,1]==parameters)]+ get(tab.cov.param2[get(tab.cov.param2[1,2])==tab.cov.param2[,4],3])  
  
                                }
                              }
                          }
                      }
              } 
            
            for (ip in 1:p) {assign(parameters[ip],theta[ip])} 
    
          
            for(kj in 1:length(unique(tim.dos$condnum))) {
                
              
               fmod<-list()
               dd<-list()
               dd2<-list()
               s2dp<-list()
               for(l in 1:model.info$nr) {
              
            tim<-tim.dos$times[tim.dos$irep==l & tim.dos$condnum==kj]
            dose<-tim.dos$doses[tim.dos$irep==l & tim.dos$condnum==kj]
           
          
          
            fmod[[l]]<-form(tim,theta,X=dose)  
           
           if (model.info$nr>1) { #single response
           if (l==1) fmod[[l]]<-fmod[[l]][1:length(tim)]
           else {
             tim0<-tim.dos$times[tim.dos$irep==(l-1) & tim.dos$condnum==kj]
             fmod[[l]]<-fmod[[l]][(length(tim0)+1):(length(tim0)+length(tim))]
             }
           }
           
           f.p[[ll]]<-c(f.p[[ll]],fmod[[l]])
          
          
           dd[[l]]<-jacobian(form,x=theta,t=tim,X=dose)

           if (nr>1){if (l==1) dd[[l]]<-dd[[l]][1:length(tim),] 
           else  {
            tim0<-tim.dos$times[tim.dos$irep==(l-1) & tim.dos$condnum==kj]
            dd[[l]]<-dd[[l]][(length(tim0)+1):(length(tim0)+length(tim)),]
            }
           }
           
           df1<-dd[[l]]
                        
                     
          if(length(tim)==1) df1<-matrix(df1,nrow=p,ncol=1)*parameters.cov.val else df1<-t(df1)*parameters.cov.val
                     
               
          #### 
          if (model.info$Trand==2) df1<-df1*model.info$beta*exp(beta.occh) 
          df.p[[ll]]<-cbind(df.p[[ll]],df1)
          ####
                      
           

                }
              }
              
              

     #################   
      for (icat in 1:dim(info.covariate$tab.cov.param)[1])
         {
                
                 if(!info.covariate$tab.cov.param[icat,4]==get(info.covariate$tab.cov.param[icat,2])) s2<-df.p[[ll]][which(parameters==info.covariate$tab.cov.param[[icat]]),]*0 else  s2<-df.p[[ll]][which(parameters==info.covariate$tab.cov.param[[icat]]),]
                   df.p[[ll]]<-rbind(df.p[[ll]],s2)
                
         }  
     theta<-theta_init  
      ######################
    }
    
   }else {
  
  #model without covariate
     for (i in 1:p) 
      assign(parameters[i],thetah[i])
      #model function
       df.p<-NULL
       f.p<-NULL
       ddf.p<-list()
       l0<-0
       
       
       for(kj in 1:length(unique(tim.dos$condnum))) {
         #cond<-eval(init.info$condinit[kj])
         fmod<-list()
         dd<-list()
         dd2<-list()
         s2dp<-list()
         for(l in 1:model.info$nr) {
          tim<-tim.dos$times[tim.dos$irep==l & tim.dos$condnum==kj]
          dose<-tim.dos$doses[tim.dos$irep==l & tim.dos$condnum==kj]
       
            
           fmod[[l]]<-form(tim,theta,X=dose)  
            
           if (model.info$nr>1) { #single response
           if (l==1) fmod[[l]]<-fmod[[l]][1:length(tim)]
           else {
             tim0<-tim.dos$times[tim.dos$irep==(l-1) & tim.dos$condnum==kj]
             fmod[[l]]<-fmod[[l]][(length(tim0)+1):(length(tim0)+length(tim))]
             }
           }
           
           f.p<-c(f.p,fmod[[l]])
          
          
           dd[[l]]<-jacobian(form,x=thetah,t=tim,X=dose)

           if (nr>1){if (l==1) dd[[l]]<-dd[[l]][1:length(tim),] 
           else  {
            tim0<-tim.dos$times[tim.dos$irep==(l-1) & tim.dos$condnum==kj]
            dd[[l]]<-dd[[l]][(length(tim0)+1):(length(tim0)+length(tim)),]
            }
           }
           
           df1<-dd[[l]]
           
           
           #dd[[l]]<-lapply(1:lfd,function(lfd,form,theta,tim,l) t(jacobian(form,x=theta,t=tim[lfd],X=dose)),form=form,theta=theta,tim=tim)
            
           #            ddl<-c()
           #             for (i in 1:lfd)
           #             {
           #             ddli<-dd[[l]][[i]][,l]
                        
            #            ddl<-c(ddl,ddli)
            #            }
            #            df1<-t(matrix(unlist(ddl),ncol=length(tim)))[,1:p]
                       
           
           if(length(tim)==1) df1<-matrix(df1,nrow=p,ncol=1) else df1<-t(df1)
           if (model.info$Trand==2) df1<-df1*model.info$beta*exp(beta.occh) 
           df.p<-cbind(df.p,df1)
          

              }
            }

   }
  
 
  return(list(mod=f.p,sensi=df.p))
}


######################################################################
#Compute model predictions and 1st and 2nd derivations including IOV #
######################################################################
                           
modelprem.cov.iov<-function(tim.dos,model.info,struc.info,sampl.info,init.info,option,info.covariate,sequ) {
  n_occ<-model.info$n_occ;
  if (length(which(model.info$gamma==0))==length(model.info$gamma)) n_occ<-1
  lP<-model.info$lP;lP1<-model.info$lP1;vecP1<-model.info$vecP1;comb.seq.poss<-model.info$comb.seq.poss;
  nbcov_occ<-model.info$nbcov_occ ;beta.occ<-model.info$beta.occ;beta.occ1<-model.info$beta.occ1;beta.occ2<-model.info$beta.occ2
  p<-length(model.info$beta);
  theta1<-model.info$theta1;theta2<-model.info$theta2
  if (info.covariate$covariate.model==T) r<-length(c(model.info$beta,unlist(info.covariate$beta.covariate)))+model.info$locc
  if (info.covariate$covariate.model==F) r<-length(model.info$beta)+model.info$locc
  co<-dim(tim.dos)[1]
  nr<-model.info$nr
  tim.dos<-create.timdos(model.info,struc.info,sampl.info,init.info)
  if(struc.info$modelform=="AF"&&NUM==F) { 
  df.forms<-symbolic.prem(model.info,struc.info,sort(unique(tim.dos$iform)),option)[[1]]
  ddf.forms<-symbolic.prem(model.info,struc.info,sort(unique(tim.dos$iform)),option)[[2]]
  }

  if (info.covariate$covariate.model==T){
    mod<-list()
    mod_int<-list()
    d1<-list()
    d1_int<-list()
    for (j in 1:dim(info.covariate$covariate.cat.poss)[1]){
      mod[[j]]<-numeric()
      mod_int[[j]]<-numeric()
      d1[[j]]<-matrix(nrow=r,ncol=0)
      for(occ in 1:n_occ){
        if(struc.info$modelform=="AF") {
          if (NUM==F){ 
          y<-calc.modelprem(tim.dos,model.info,struc.info,df.forms,ddf.forms,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[1]]
          dy<-calc.modelprem(tim.dos,model.info,struc.info,df.forms,ddf.forms,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[2]]
          }
          else{
          y<-calc.model.num(tim.dos,model.info,struc.info,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[1]]
          dy<-calc.model.num(tim.dos,model.info,struc.info,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[2]]
          } 
        }else 
        {
          y<-calc.model.de(tim.dos,model.info,init.info,struc.info,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[1]]
          dy<-calc.model.de(tim.dos,model.info,init.info,struc.info,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[2]]
        }
        if (Trand ==2) {
          dy[[j]]<-dy[[j]]
        }
        mod_int[[j]]<-y[[j]]
        mod[[j]]<-c(mod[[j]],mod_int[[j]])
        bloc0<-matrix(rep(0,model.info$locc*co),nrow=model.info$locc,ncol=co)
        if (lP[[sequ]][occ]==0) d1_int[[j]]<-rbind(dy[[j]],bloc0)
        if (lP[[sequ]][occ]>0){
          block<-matrix(nrow=0,ncol=co)
          d1_int0k<-matrix()
          for (kj in 1:nbcov_occ){
            Seq<-comb.seq.poss[sequ,kj]
            
            l<-which(covariate_occ.sequence[[kj]][[Seq]][occ]==covariate_occ.category[[kj]])-1 

    block_int<-matrix(0,nrow=length(parameter_occ.associated[[kj]])*(length(covariate_occ.category[[kj]])-1),ncol=co)
    if (length(which(vecP1[[kj]][[Seq]][[occ]]==T))>0) 
    {
    for (jj in 1:length(which(vecP1[[kj]][[Seq]][[occ]]==T))){
        i=l+(jj-1)*(length(covariate_occ.category[[kj]])-1)
        if (Trand==2) s_int_ij<-dy[[j]][1:p,][which(vecP1[[kj]][[Seq]][[occ]]==T)[jj],]
        else  s_int_ij<-dy[[j]][1:p,][which(vecP1[[kj]][[Seq]][[occ]]==T)[jj],]
        block_int[i,]<-s_int_ij
     }
    }
    block<-rbind(block,block_int)
  }

          d1_int[[j]]<-rbind(dy[[j]],block)
        }
        d1[[j]]<-cbind(d1[[j]],d1_int[[j]])
      }
    }
  }

  else
  {
    mod<-numeric()
    mod_int<-numeric()
    d1<-matrix(nrow=r,ncol=0)
    d1_int<-matrix(nrow=r,ncol=0)
    for(occ in 1:n_occ){ 
       if(struc.info$modelform=="AF") { 
          if (NUM==F){ 
          y<-calc.modelprem(tim.dos,model.info,struc.info,df.forms,ddf.forms,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[1]]
          dy<-calc.modelprem(tim.dos,model.info,struc.info,df.forms,ddf.forms,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[2]]
          }
          else{
          y<-calc.model.num(tim.dos,model.info,struc.info,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[1]]
          dy<-calc.model.num(tim.dos,model.info,struc.info,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[2]]
          } 
        }else 
        {
          y<-calc.model.de(tim.dos,model.info,init.info,struc.info,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[1]]
          dy<-calc.model.de(tim.dos,model.info,init.info,struc.info,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[2]]
        }
      
      mod_int<-y
      mod<-c(mod,mod_int)
      bloc0<-matrix(rep(0,model.info$locc*co),nrow=model.info$locc,ncol=co)
      if (lP[[sequ]][occ]==0) d1_int<-rbind(dy,bloc0)
      if (lP[[sequ]][occ]>0){
          block<-matrix(nrow=0,ncol=co)
          for (kj in 1:nbcov_occ){
            Seq<-comb.seq.poss[sequ,kj]
            
            l<-which(covariate_occ.sequence[[kj]][[Seq]][occ]==covariate_occ.category[[kj]])-1 

    block_int<-matrix(0,nrow=length(parameter_occ.associated[[kj]])*(length(covariate_occ.category[[kj]])-1),ncol=co)
    if (length(which(vecP1[[kj]][[Seq]][[occ]]==T))>0) 
    {
    for (jj in 1:length(which(vecP1[[kj]][[Seq]][[occ]]==T))){
        i=l+(jj-1)*(length(covariate_occ.category[[kj]])-1)
        if (Trand==2) s_int_ij<-dy[1:p,][which(vecP1[[kj]][[Seq]][[occ]]==T)[jj],]
        else  s_int_ij<-dy[1:p,][which(vecP1[[kj]][[Seq]][[occ]]==T)[jj],]
        block_int[i,]<-s_int_ij
     }
    }
    block<-rbind(block,block_int)
  }
          d1_int<-rbind(dy,block)
        }
      d1<-cbind(d1,d1_int)
    }
  } 
 l<-list(mod,d1)
 return(l)
}



#Pour les effets al?atoires

modelprem.cov.iov.rand<-function(tim.dos,model.info,struc.info,sampl.info,init.info,option,info.covariate,sequ) {
  n_occ<-model.info$n_occ
  if (length(which(model.info$gamma==0))==length(model.info$gamma)) n_occ<-1
  p<-length(model.info$beta)
  if (info.covariate$covariate.model==T) r<-length(c(model.info$beta,unlist(info.covariate$beta.covariate)))+model.info$locc
  else  r<-length(model.info$beta)+model.info$locc
  co<-dim(tim.dos)[1]
  theta2<-model.info$theta2
  beta.occ2<-model.info$beta.occ2
  tim.dos<-create.timdos(model.info,struc.info,sampl.info,init.info)
  #kform<-sort(unique(tim.dos$iform))
  if(struc.info$modelform=="AF"&&NUM==F) { 
  df.forms<-symbolic.prem(model.info,struc.info,sort(unique(tim.dos$iform)),option)[[1]]
  ddf.forms<-symbolic.prem(model.info,struc.info,sort(unique(tim.dos$iform)),option)[[2]]
  }
  
  if (info.covariate$covariate.model==T){
    mod<-list()
    mod_int<-list()
    d11<-list()
    d12<-list()
    d1<-list()
    d11_int<-list()
    for (j in 1:dim(info.covariate$covariate.cat.poss)[1]){
      mod[[j]]<-numeric()
      mod_int[[j]]<-numeric()
      d11[[j]]<-matrix(nrow=p,ncol=0)
      d12[[j]]<-matrix(0,nrow=p*n_occ,ncol=co*n_occ)
      d1[[j]]<-matrix()
      if (n_occ==1){
        if(struc.info$modelform=="AF") { 
                    if (NUM==F){ 
          y<-calc.modelprem(tim.dos,model.info,struc.info,df.forms,ddf.forms,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[1]]
          dy<-calc.modelprem(tim.dos,model.info,struc.info,df.forms,ddf.forms,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[2]]
          }
          else{
          y<-calc.model.num(tim.dos,model.info,struc.info,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[1]]
          dy<-calc.model.num(tim.dos,model.info,struc.info,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[2]]
          }
        }else 
        {
          y<-calc.model.de(tim.dos,model.info,init.info,struc.info,option,info.covariate,theta2[[sequ]][[1]],beta.occ2[[sequ]][[1]])[[1]]
          dy<-calc.model.de(tim.dos,model.info,init.info,struc.info,option,info.covariate,theta2[[sequ]][[1]],beta.occ2[[sequ]][[1]])[[2]]
        }
     
      mod[[j]]<-y[[j]]
      d1[[j]]<-dy[[j]][1:p,]}
      else{
        for(occ in 1:n_occ){
          if(struc.info$modelform=="AF") { 
                   if (NUM==F){ 
          y<-calc.modelprem(tim.dos,model.info,struc.info,df.forms,ddf.forms,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[1]]
          dy<-calc.modelprem(tim.dos,model.info,struc.info,df.forms,ddf.forms,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[2]]
          }
          else{
          y<-calc.model.num(tim.dos,model.info,struc.info,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[1]]
          dy<-calc.model.num(tim.dos,model.info,struc.info,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[2]]
          }
          }else 
        {
          y<-calc.model.de(tim.dos,model.info,init.info,struc.info,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[1]]
          dy<-calc.model.de(tim.dos,model.info,init.info,struc.info,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[2]]
        }
          mod_int[[j]]<-y[[j]]
          mod[[j]]<-c(mod[[j]],mod_int[[j]])

          d11_int[[j]]<-dy[[j]][1:p,]
          d11[[j]]<-cbind(d11[[j]],d11_int[[j]])
          d12[[j]][(p*(occ-1)+1) : (p*occ),(co*(occ-1)+1) :(co*occ)]<-d11_int[[j]]
        }
      d1[[j]]<-rbind(d11[[j]],d12[[j]])
      }
   }
  
 } else 
 {
   mod<-numeric()
   mod_int<-numeric()
   d11<-matrix(nrow=p,ncol=0)
   d12<-matrix(0,nrow=p*n_occ,ncol=co*n_occ)
   if (n_occ==1){ 
        if(struc.info$modelform=="AF") { 
                   if (NUM==F){ 
          y<-calc.modelprem(tim.dos,model.info,struc.info,df.forms,ddf.forms,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[1]]
          dy<-calc.modelprem(tim.dos,model.info,struc.info,df.forms,ddf.forms,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[2]]
          }
          else{
          y<-calc.model.num(tim.dos,model.info,struc.info,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[1]]
          dy<-calc.model.num(tim.dos,model.info,struc.info,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[2]]
          }
          }else 
        {
          y<-calc.model.de(tim.dos,model.info,init.info,struc.info,option,info.covariate,theta2[[sequ]][[1]],beta.occ2[[sequ]][[1]])[[1]]
          dy<-calc.model.de(tim.dos,model.info,init.info,struc.info,option,info.covariate,theta2[[sequ]][[1]],beta.occ2[[sequ]][[1]])[[2]]
        }   
        mod<-y
   d1<-dy[1:p,]}
   else{
    for(occ in 1:n_occ){
       if(struc.info$modelform=="AF") { 
         if (NUM==F){ 
          y<-calc.modelprem(tim.dos,model.info,struc.info,df.forms,ddf.forms,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[1]]
          dy<-calc.modelprem(tim.dos,model.info,struc.info,df.forms,ddf.forms,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[2]]
          }
          else{
          y<-calc.model.num(tim.dos,model.info,struc.info,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[1]]
          dy<-calc.model.num(tim.dos,model.info,struc.info,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[2]]
          }}else 
        {
          y<-calc.model.de(tim.dos,model.info,init.info,struc.info,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[1]]
          dy<-calc.model.de(tim.dos,model.info,init.info,struc.info,option,info.covariate,theta2[[sequ]][[occ]],beta.occ2[[sequ]][[occ]])[[2]]
        }
         mod_int<-y
     mod<-c(mod,mod_int)

     d11_int<-dy[1:p,]
     d11<-cbind(d11,d11_int)
     d12[(p*(occ-1)+1) : (p*occ),(co*(occ-1)+1) :(co*occ)]<-d11_int


   }
   d1<-rbind(d11,d12)

   }
 }
 l<-list(mod,d1)
 return(l)
}


 sensifixe<-function(tim.dos,model.info,struc.info,sampl.info,init.info,option,info.covariate)
{ 
  s<-list()
  if (model.info$lposs==1) s[[1]]<-modelprem.cov.iov(tim.dos,model.info,struc.info,sampl.info,init.info,option,info.covariate,sequ=1)
  if (model.info$lposs>1){
  for (sequ in 1:model.info$lposs)
    {
    s[[sequ]]<-modelprem.cov.iov(tim.dos,model.info,struc.info,sampl.info,init.info,option,info.covariate,sequ)
    }
  }
  return(s)
}
#sensi<-sensifixe(tim.dos,model.info,struc.info,sampl.info,init.info,option,info.covariate)

sensialea<-function(tim.dos,model.info,struc.info,sampl.info,init.info,option,info.covariate)
{
  s<-list()
  if (model.info$lposs==1) s[[1]]<-modelprem.cov.iov.rand(tim.dos,model.info,struc.info,sampl.info,init.info,option,info.covariate,sequ=1)
  if (model.info$lposs>1){
  for (sequ in 1:model.info$lposs)
    {
    s[[sequ]]<-modelprem.cov.iov.rand(tim.dos,model.info,struc.info,sampl.info,init.info,option,info.covariate,sequ)
    }
  }
  return(s)
}

#sensia<-sensialea(tim.dos,model.info,struc.info,sampl.info,init.info,option,info.covariate)



############################################################################
####                Fisher information matrix, method 1                 ####
############################################################################
#Only Fisher information matrix (for each elementary protocol)
#ecofisher.meth1(protk,tim.dos,f.p,df.p,model.info,option,ddf.p)
ecofisher.meth1<-function(ind1,tab,f.p,df.p,model.info,option,ddf.p) {
#No computation of C
#neglect all terms involving d2.f
#ind1 : indices des temps (rang?es de tab)
#tab : a table with times, doses, form and response for each time considered
#Computes only the Fisher information matrix (assumes elementary protocol)

  p<-length(model.info$beta)
  nr<-model.info$nr
#Calcul de V et V-1
   v1<-model.info$omega%*%df.p[,ind1]
   vvar<-t(v1)%*%df.p[,ind1]
   dia<-(tab$sig.i[ind1]+tab$sig.s[ind1]*f.p[ind1])**2
   vvar<-vvar+diag(dia,nrow=length(dia))
    invvar<-solve(vvar)     
    if (option==2) dvar<-lapply(1:p,function(p,ddf.p)
             t(ddf.p[[p]][,ind1])%*%v1+t(v1)%*%ddf.p[[p]][,ind1]+2*diag(sqrt(dia),length(dia))*diag(tab$sig.s[ind1]*df.p[p,ind1]/model.info$beta[p],length(tab$sig.s[ind1]*df.p[p,ind1]/model.info$beta[p])),ddf.p=ddf.p)
  
#Calcul des dV/dlambda_m.V-1 (m=1...p)
   dvdlam<-vector(p+2*nr,mode="list")
   dv2<-vector(p+2*nr,mode="list")
   for(m in 1:p) {
     dv1<- df.p[m,ind1]%*%t(df.p[m,ind1]) 
     dv2[[m]]<-dv1
     dvdlam[[m]]<-dv1%*%invvar
   }

#Calcul des dV/dlambda_m.V-1 (m>p)
   vecmi<-2*(tab$sig.i[ind1]+tab$sig.s[ind1]*f.p[ind1])
   vecms<-vecmi*f.p[ind1]
   for(l in 1:(2*nr))
      dvdlam[[p+l]]<-diag(rep(0,length(vecmi)),nrow=length(vecmi))
   for(l in 1:nr) { #sigma.inter
      if (model.info$sigmainter[[l]]>0) {
       diag(dvdlam[[p+(l-1)*2+1]])[tab$irep[ind1]==l]<-vecmi[tab$irep[ind1]==l]
       dvdlam[[p+(l-1)*2+1]]<-dvdlam[[p+(l-1)*2+1]]%*%invvar
       }
      }
   for(l in 1:nr) { #sigma.slope
      if (model.info$sigmaslope[[l]]>0) { 
       diag(dvdlam[[p+(l-1)*2+2]])[tab$irep[ind1]==l]<-vecms[tab$irep[ind1]==l]
       dvdlam[[p+(l-1)*2+2]]<-dvdlam[[p+(l-1)*2+2]]%*%invvar
       }
      }
   resB<-matrix(nrow=p+nr*2,ncol=p+nr*2)     #matrice al?atoire #
   resC<-matrix(rep(0,(p+nr*2)*p),nrow=p+nr*2,ncol=p)   #matrice incomplete option 1 covariance==0#
#Division par beta dns le cas Trand==2
   if(model.info$Trand==2) df1<-df.p[,ind1]/model.info$beta else df1<-df.p[,ind1]
   resA<-(df1 %*% invvar) %*% t(df1)  #matrice parametres fixe#
  
   for(i in 1:(p+2*nr)) {
      for(j in 1:(p+2*nr)) {
        resB[i,j]<-sum(diag(dvdlam[[i]]%*%dvdlam[[j]]))
        }
   }
   resB<-resB/2
#computation of block C for  matrice  complete
   if (option==2){
     #resA2
   resA2<-matrix(ncol=p,nrow=p)
   for (i in 1:p){
   resA3<-invvar%*%dvar[[i]]%*%invvar
   resA2[i,]<-sapply(1:p,function(dvar,resA3,p) sum(diag(resA3%*%dvar[[p]])),dvar=dvar,resA3=resA3)
   
   for(j in 1:p)
          {
            resC1<-invvar%*%dvar[[j]]%*%invvar
		        resC[i,j]<-sum(diag(resC1%*%dv2[[i]]))
            for ( j2 in 1:nr) {  
                vec_nulli<-rep(0,length(vecmi))
                vec_nulls<-rep(0,length(vecms))
                vec_nulli[which(tab$irep[ind1]==j2)]<-vecmi[tab$irep[ind1]==j2]
                vec_nulls[which(tab$irep[ind1]==j2)]<-vecms[tab$irep[ind1]==j2]
		            resC[p+j2*2-1,j]<-sum(diag(resC1%*%diag(vec_nulli,length(vec_nulli))))
		            resC[p+j2*2,j]<-sum(diag(resC1%*%diag(vec_nulls,length(vec_nulls))))
            }  
		      }
    
   }
   
  
   resA<-resA+resA2/2
   resC<-resC/2
   }
#Enlever les lignes de B et C pour lesquelles omega ou sigma=0
   if(length(model.info$nullpar)>0) {
     resB<-resB[-c(model.info$nullpar),-c(model.info$nullpar)]
     resC<-resC[-c(model.info$nullpar),]
   }

     fish<-cbind(rbind(resA,resC),rbind(t(resC),resB))
      
    #fixed parameters if beta.fixed	== T
    beta.estim<-!model.info$beta.fixed
    if (dim(fish)[1]>length(model.info$beta)) vec<-c(beta.estim,rep(T,dim(fish)[1]-length(model.info$beta)))  else vec<-beta.estim
    fish<-fish[vec,vec]
    
   return(fish)
}



####Individual Fisher matrix###############
indifisher.meth1<-function(ind1,tab,f.p,df.p,model.info) {
#No computation of C
#neglect all terms involving d2.f
#ind1 : indices des temps (rang?es de tab)
#tab : a table with times, doses, form and response for each time considered
#Computes only the Fisher information matrix (assumes elementary protocol)
  #print(ind1)
  p<-length(model.info$beta)
  nr<-model.info$nr
#Calcul de V et V-1
   #v1<-model.info$omega%*%df.p[,ind1]
   #vvar<-t(v1)%*%df.p[,ind1]
   dia<-(tab$sig.i[ind1]+tab$sig.s[ind1]*f.p[ind1])**2
   vvar<-diag(dia,nrow=length(dia))
   invvar<-solve(vvar)     
   

   dvdlam<-vector(2*nr,mode="list")
   
#Calcul des dV/dlambda_m.V-1 
   vecmi<-2*(tab$sig.i[ind1]+tab$sig.s[ind1]*f.p[ind1])
   vecms<-vecmi*f.p[ind1]
   for(l in 1:(2*nr))
      dvdlam[[l]]<-diag(rep(0,length(vecmi)),nrow=length(vecmi))
   for(l in 1:nr) { #sigma.inter
      if (model.info$sigmainter[[l]]>0) {
       diag(dvdlam[[(l-1)*2+1]])[tab$irep[ind1]==l]<-vecmi[tab$irep[ind1]==l]
       dvdlam[[(l-1)*2+1]]<-dvdlam[[(l-1)*2+1]]%*%invvar
       }
      }
   for(l in 1:nr) { #sigma.slope
      if (model.info$sigmaslope[[l]]>0) { 
       diag(dvdlam[[(l-1)*2+2]])[tab$irep[ind1]==l]<-vecms[tab$irep[ind1]==l]
       dvdlam[[(l-1)*2+2]]<-dvdlam[[(l-1)*2+2]]%*%invvar
       }
      }
   resB<-matrix(nrow=nr*2,ncol=nr*2)     #block for variance terms (no random effect)
   resC<-matrix(rep(0,(nr*2)*p),nrow=nr*2,ncol=p)   #matrice incomplete option 1 covariance==0#
#Division par beta dns le cas Trand==2
   if(model.info$Trand==2) df1<-df.p[,ind1]/model.info$beta else df1<-df.p[,ind1]
   resA<-(df1 %*% invvar) %*% t(df1)  #matrice parametres fixe#
  
   for(i in 1:(2*nr)) {
      for(j in 1:(2*nr)) {
        resB[i,j]<-sum(diag(dvdlam[[i]]%*%dvdlam[[j]]))
        }
   }
   resB<-resB/2

#Enlever les lignes de B et C pour lesquelles sigma=0
   if(length(model.info$nullsig)>0) {
     resB<-as.matrix(resB[-c(model.info$nullsig),-c(model.info$nullsig)])
     resC<-matrix(resC[-c(model.info$nullsig),],ncol=p)
   }

     fish<-cbind(rbind(resA,resC),rbind(t(resC),resB))
    
   #fixed parameters if beta.fixed	== T
    beta.estim<-!model.info$beta.fixed
    if (dim(fish)[1]>length(model.info$beta)) vec<-c(beta.estim,rep(T,dim(fish)[1]-length(model.info$beta)))  else vec<-beta.estim
    fish<-fish[vec,vec]

   return(fish)
}




####Bayesian Fisher matrix###############
bayfisher.meth1<-function(ind1,tab,f.p,df.p,model.info) {
#No computation of C
#neglect all terms involving d2.f
#ind1 : indices des temps (rang?es de tab)
#tab : a table with times, doses, form and response for each time considered
#Computes only the Fisher information matrix (assumes elementary protocol)


  p<-length(model.info$beta)
  nr<-model.info$nr
  omega1<-as.matrix(as.matrix(model.info$omega[,diag(model.info$omega)!=0])[diag(model.info$omega)!=0,])
  beta1<-model.info$beta[diag(model.info$omega)!=0]
  #solve(as.matrix(diag((beta[diag(omega)!=0]),length(beta[diag(omega)!=0])))%*%as.matrix(omega1)%*%as.matrix(diag(((beta[diag(omega)!=0])),length(beta[diag(omega)!=0]))))
    
#Calcul de V et V-1
   #v1<-model.info$omega%*%df.p[,ind1]
   #vvar<-t(v1)%*%df.p[,ind1]
   dia<-(tab$sig.i[ind1]+tab$sig.s[ind1]*f.p[ind1])**2
   vvar<-diag(dia,nrow=length(dia))
   invvar<-solve(vvar)     
   
#Division par beta dns le cas Trand==2
   if(model.info$Trand==2) df1<-df.p[,ind1]/model.info$beta else df1<-df.p[,ind1]
   #taking into account a priori knowlegde on random effect distribution
   if(model.info$Trand==2) solveomega<-solve(as.matrix(diag(beta1,length(beta1)))%*%as.matrix(omega1)%*%as.matrix(diag(beta1,length(beta1)))) else solveomega<-solve(omega1)
   #matrice parametres fixe#
   resA<-((df1 %*% invvar) %*% t(df1))[diag(model.info$omega)!=0,diag(model.info$omega)!=0] + solveomega
   
   
#fixed parameters	if beta.fixed	== T
    if (FIM=="B") {beta.estim<-!model.info$beta.fixed[diag(model.info$omega)!=0]}
    else{beta.estim<-!model.info$beta.fixed}
    if (dim(resA)[1]>length(model.info$beta)) vec<-c(beta.estim,rep(T,dim(resA)[1]-length(model.info$beta)))  else vec<-beta.estim
    resA<-resA[vec,vec]
    
  return(resA)
}


##with covariates##########
ecofisher.meth1.cov<-function(ind1,tab,f.p,df.p,model.info,option,ddf.p,info.covariate) {
#No computation of C
#neglect all terms involving d2.f
#ind1 : indices des temps (rang?es de tab)
#tab : a table with times, doses, form and response for each time considered
#Computes only the Fisher information matrix (assumes elementary protocol)
   p.cov<-length(c(model.info$beta,unlist(info.covariate$beta.covariate)))
   nr<-model.info$nr
   p<-length(model.info$beta)
   pp<-(p.cov+p)+2*nr-length(model.info$nullpar) #2 parameters for the residual variance by response#
   somme1<-matrix(c(rep(0,pp*pp)),ncol=pp) 
for (ll in 1: dim(info.covariate$covariate.cat.poss)[1])
    {
          categories.prop<-prod(info.covariate$tab_mod_prop[[2]][ll,])
          #Calcul de V et V-1
           v1<-model.info$omega%*%as.matrix(as.matrix(df.p[[ll]][,ind1])[1:p,])
           vvar<-t(v1)%*%as.matrix(as.matrix(df.p[[ll]][,ind1])[1:p,])
           dia<-(tab$sig.i[ind1]+tab$sig.s[ind1]*f.p[[ll]][ind1])**2
           vvar<-vvar+diag(dia,nrow=length(dia))
            invvar<-solve(vvar)     
            if (option==2) dvar<-lapply(1:p,function(p,ddf.p)
                     t(ddf.p[[p]][,ind1])%*%v1+t(v1)%*%ddf.p[[p]][,ind1]+2*diag(sqrt(dia),length(dia))*diag(tab$sig.s[ind1]*df.p[p,ind1]/model.info$beta[p],length(tab$sig.s[ind1]*df.p[p,ind1]/model.info$beta[p])),ddf.p=ddf.p)
          
        #Calcul des dV/dlambda_m.V-1 (m=1...p)
           dvdlam<-vector(p+2*nr,mode="list")
           dv2<-vector(p+2*nr,mode="list")
           for(m in 1:p) {
             dv1<- df.p[[ll]][m,ind1]%*%t(df.p[[ll]][m,ind1]) 
             dv2[[m]]<-dv1
             dvdlam[[m]]<-dv1%*%invvar
           }
        
        #Calcul des dV/dlambda_m.V-1 (m>p)
           vecmi<-2*(tab$sig.i[ind1]+tab$sig.s[ind1]*f.p[[ll]][ind1])
           vecms<-vecmi*f.p[[ll]][ind1]
           for(l in 1:(2*nr))
              dvdlam[[p+l]]<-diag(rep(0,length(vecmi)),nrow=length(vecmi))
           for(l in 1:nr) { #sigma.inter
              if (model.info$sigmainter[[l]]>0) {
               diag(dvdlam[[p+(l-1)*2+1]])[tab$irep[ind1]==l]<-vecmi[tab$irep[ind1]==l]
               dvdlam[[p+(l-1)*2+1]]<-dvdlam[[p+(l-1)*2+1]]%*%invvar
               }
              }
           for(l in 1:nr) { #sigma.slope
              if (model.info$sigmaslope[[l]]>0) { 
               diag(dvdlam[[p+(l-1)*2+2]])[tab$irep[ind1]==l]<-vecms[tab$irep[ind1]==l]
               dvdlam[[p+(l-1)*2+2]]<-dvdlam[[p+(l-1)*2+2]]%*%invvar
               }
              }
           resB<-matrix(nrow=p+nr*2,ncol=p+nr*2)     #matrice al?atoire #
        resC<-matrix(rep(0,(p+nr*2)*p.cov),nrow=p+nr*2,ncol=p.cov)     #matrice incomplete option 1 covariance==0#
        #Division par beta dns le cas Trand==2
           if(model.info$Trand==2)df1<-df.p[[ll]][,ind1]/c(model.info$beta,rep(1,length(unlist(info.covariate$beta.covariate)))) else df1<-df.p[[ll]][,ind1]
           resA<-(df1 %*% invvar) %*% t(df1)  #matrice parametres fixe#
          
           for(i in 1:(p+2*nr)) {
              for(j in 1:(p+2*nr)) {
                resB[i,j]<-sum(diag(dvdlam[[i]]%*%dvdlam[[j]]))
                }
           }
           resB<-resB/2
        #computation of block C for  matrice  complete
           if (option==2){
             #resA2
           resA2<-matrix(ncol=p,nrow=p)
           for (i in 1:p){
           resA3<-invvar%*%dvar[[i]]%*%invvar
           resA2[i,]<-sapply(1:p,function(dvar,resA3,p) sum(diag(resA3%*%dvar[[p]])),dvar=dvar,resA3=resA3)
           
           for(j in 1:p)
                  {
                    resC1<-invvar%*%dvar[[j]]%*%invvar
        		        resC[i,j]<-sum(diag(resC1%*%dv2[[i]]))
                    for ( j2 in 1:nr) {  
                        vec_nulli<-rep(0,length(vecmi))
                        vec_nulls<-rep(0,length(vecms))
                        vec_nulli[which(tab$irep[ind1]==j2)]<-vecmi[tab$irep[ind1]==j2]
                        vec_nulls[which(tab$irep[ind1]==j2)]<-vecms[tab$irep[ind1]==j2]
        		            resC[p+j2*2-1,j]<-sum(diag(resC1%*%diag(vec_nulli,length(vec_nulli))))
        		            resC[p+j2*2,j]<-sum(diag(resC1%*%diag(vec_nulls,length(vec_nulls))))
                    }  
        		      }
            
           }
           
          
           resA<-resA+resA2/2
           resC<-resC/2
           }
        #Enlever les lignes de B et C pour lesquelles omega ou sigma=0
           if(length(model.info$nullpar)>0){
             resB<-resB[-c(model.info$nullpar),-c(model.info$nullpar)]
             resC<-resC[-c(model.info$nullpar),]
           }
        
             fish<-cbind(rbind(resA,resC),rbind(t(resC),resB))
             somme1<-somme1+fish*categories.prop 
            
    } 
    #fixed parameters if beta.fixed	== T
    beta.estim<-!model.info$beta.fixed  
    if (dim(somme1)[1]>length(model.info$beta)) vec<-c(beta.estim,rep(T,dim(somme1)[1]-length(model.info$beta)))  else vec<-beta.estim
    somme1<-somme1[vec,vec] 
         
   return(somme1)
}


#fonction pour ?crire une matrice block diagonal ? partir d'une liste de matrices
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


#ind1,tab,f.p,df.p,model.info,option,ddf.p
ecofisher.meth1.iov<-function(ind1,tab,sensi,sensia,model.info,option,info.covariate) {
#No computation of C
#neglect all terms involving d2.f
#ind1 : indices des temps (rang?es de tab)
#tab : a table with times, doses, form and response for each time considered
#Computes only the Fisher information matrix (assumes elementary protocol)

  p<-length(model.info$beta)
  nr<-model.info$nr
  
  n_occ<-model.info$n_occ
  if (length(which(model.info$gamma==0))==length(model.info$gamma)) n_occ<-1
  if (n_occ==1) pr<-p else pr<-2*p
  pp<-(p+model.info$locc+pr)+2*nr-length(model.info$nullpar) #2 parameters for the residual variance by response#
  somme1<-matrix(c(rep(0,pp*pp)),ncol=pp) 
  
  for (sequ in 1:model.info$lposs){
       f.pr<-sensi[[sequ]][[1]]
       sensif<-sensi[[sequ]][[2]]
       df.pr<-sensia[[sequ]][[2]]
      
#Calcul de V et V-1
    if (length(which(model.info$gamma==0))!=length(model.info$gamma)) v1<-model.info$OMEGA%*%df.pr[,ind1]
    else  v1<-model.info$omega%*%df.pr[,ind1]
           vvar<-t(v1)%*%df.pr[,ind1]
           dia<-(rep(tab$sig.i[ind1[1:(length(ind1)/n_occ)]],n_occ)+rep(tab$sig.s[ind1[1:(length(ind1)/n_occ)]],n_occ)*f.pr[ind1])**2
           vvar<-vvar+diag(dia,nrow=length(dia))
           invvar<-solve(vvar) 
   # if (option==2) dvar<-lapply(1:pr,function(pr,ddf.pr)
    #         t(ddf.pr[[pr]][,ind1])%*%v1+t(v1)%*%ddf.pr[[pr]][,ind1]+2*diag(sqrt(dia),length(dia))*diag(tab$sig.s[ind1]*df.pr[pr,ind1]/model.info$beta[p],length(tab$sig.s[ind1]*df.p[p,ind1]/model.info$beta[p])),ddf.pr=ddf.pr)
  
#Calcul des dV/dlambda_m.V-1 (m=1...p)

   dvdlam<-vector(pr+2*nr,mode="list")
   dv2<-vector(pr+2*nr,mode="list")
   for(m in 1:pr) {
      dv1<- df.pr[m,ind1]%*%t(df.pr[m,ind1])
      dv2[[m]]<-dv1  
      dvdlam[[m]]<-dv1%*%invvar 
   }
   for(m in 1:p) {
      dv1<- df.pr[m,ind1]%*%t(df.pr[m,ind1])
      dv1bis<-blockdiag(sapply(1:n_occ,function(mat,n_occ) list((mat)%*%t(mat)),mat=(df.pr[m,ind1])[1:(length(ind1)/n_occ)])) 
      dv2[[m]]<-dv1  
      dvdlam[[m]]<-dv1%*%invvar 
      dv2[[p+m]]<-dv1bis  
      dvdlam[[p+m]]<-dv1bis%*%invvar 
   }

  
  

        #Calcul des dV/dlambda_m.V-1 (m>p)
           vecmi<-2*(rep(tab$sig.i[ind1[1:(length(ind1)/n_occ)]],n_occ)+rep(tab$sig.s[ind1[1:(length(ind1)/n_occ)]],n_occ)*f.pr[ind1])
           vecms<-vecmi*f.pr[ind1]
           for(l in 1:(2*nr))
              dvdlam[[pr+l]]<-diag(rep(0,length(vecmi)),nrow=length(vecmi))
           for(l in 1:nr) { #sigma.inter
              if (model.info$sigmainter[[l]]>0) {
               diag(dvdlam[[pr+(l-1)*2+1]])[rep(tab$irep[ind1[1:(length(ind1)/n_occ)]],n_occ)==l]<-vecmi[rep(tab$irep[ind1[1:(length(ind1)/n_occ)]],n_occ)==l]
               dvdlam[[pr+(l-1)*2+1]]<-dvdlam[[pr+(l-1)*2+1]]%*%invvar
               }
              }
           for(l in 1:nr) { #sigma.slope
              if (model.info$sigmaslope[[l]]>0) { 
               diag(dvdlam[[pr+(l-1)*2+2]])[rep(tab$irep[ind1[1:(length(ind1)/n_occ)]],n_occ)==l]<-vecms[rep(tab$irep[ind1[1:(length(ind1)/n_occ)]],n_occ)==l]
               dvdlam[[pr+(l-1)*2+2]]<-dvdlam[[pr+(l-1)*2+2]]%*%invvar
               }                                                
              }
              
            resB<-matrix(nrow=pr+nr*2,ncol=pr+nr*2)     #matrice al?atoire #
           resB2<-matrix(nrow=pr+nr*2,ncol=pr+nr*2) 
        resC<-matrix(rep(0,(pr+nr*2)*(p+model.info$locc)),nrow=pr+nr*2,ncol=(p+model.info$locc))     #matrice incomplete option 1 covariance==0#
 #Division par beta dans le cas Trand==2  
   if (Trand==2) df1<-sensif[,ind1]/c(model.info$beta,rep(1,model.info$locc)) 
     
   else df1<-sensif[,ind1]
      
   resA<-(df1 %*% invvar) %*% t(df1)  #matrice parametres fixe#
     
   for(i in 1:(pr+2*nr)) {
    for(j in 1:(pr+2*nr)) {
      resB[i,j]<-sum(diag(dvdlam[[i]]%*%dvdlam[[j]]))
    }
   }
           

  
        
        
#B final
        
           resBf<-resB/2
           
#Enlever les lignes de B et C pour lesquelles omega ou sigma=0
   if(length(model.info$nullpar)>0) {
     resBf<-resBf[-c(model.info$nullpar),-c(model.info$nullpar)]
     resC<-resC[-c(model.info$nullpar),]
   }

     fish<-cbind(rbind(resA,resC),rbind(t(resC),resBf))
    somme1<-somme1+fish*model.info$comb.seq.prop[sequ] 
}
    #fixed parameters if beta.fixed	== T
    beta.estim<-!model.info$beta.fixed
    if (dim(somme1)[1]>length(model.info$beta)) vec<-c(beta.estim,rep(T,dim(somme1)[1]-length(model.info$beta)))  else vec<-beta.estim
    somme1<-somme1[vec,vec] 
     
   return(somme1)
}


ecofisher.meth1.cov.iov<-function(ind1,tab,sensi,sensia,model.info,option,info.covariate) {
#df.pr and ddf.pr correspond to derivations / random     
#No computation of C
#neglect all terms involving d2.f
#ind1 : indices des temps (rang?es de tab)
#tab : a table with times, doses, form and response for each time considered
#Computes only the Fisher information matrix (assumes elementary protocol)
   p.cov<-length(c(model.info$beta,unlist(info.covariate$beta.covariate)))
   nr<-model.info$nr
   p<-length(model.info$beta)
   n_occ<-model.info$n_occ
   if (length(which(model.info$gamma==0))==length(model.info$gamma)) n_occ<-1
   if (n_occ==1) pr<-p else pr<-2*p
   pp<-(p.cov+model.info$locc+pr)+2*nr-length(model.info$nullpar) #2 parameters for the residual variance by response#
   somme1<-matrix(c(rep(0,pp*pp)),ncol=pp)                                                        
    #tab<-tim.dos ; ind1<-c(1:14); f.p<-modelprem[[1]]; df.p<-modelprem[[2]]
    
  
   
for (sequ in 1:model.info$lposs){     
for (ll in 1: dim(info.covariate$covariate.cat.poss)[1])
    {  f.pr<-sensi[[sequ]][[1]]
       sensif<-sensi[[sequ]][[2]]
       df.pr<-sensia[[sequ]][[2]]
       
          categories.prop<-prod(info.covariate$tab_mod_prop[[2]][ll,])
          #Calcul de V et V-1
    if (length(which(model.info$gamma==0))!=length(model.info$gamma)) v1<-model.info$OMEGA%*%df.pr[[ll]][,ind1]
    else  v1<-model.info$omega%*%df.pr[[ll]][,ind1]
           vvar<-t(v1)%*%df.pr[[ll]][,ind1]
           dia<-(rep(tab$sig.i[ind1[1:(length(ind1)/n_occ)]],n_occ)+rep(tab$sig.s[ind1[1:(length(ind1)/n_occ)]],n_occ)*f.pr[[ll]][ind1])**2
           vvar<-vvar+diag(dia,nrow=length(dia))
            invvar<-solve(vvar)     
            #if (option==2) dvar<-lapply(1:l,function(l,ddf.pr)
            #         t(ddf.pr[[pr]][,ind1])%*%v1+t(v1)%*%ddf.pr[[pr]][,ind1]+2*diag(sqrt(dia),length(dia))*diag(tab$sig.s[ind1]*df.p[pr,ind1]/c(model.info$beta[p],model.info$beta[p]),length(tab$sig.s[ind1]*df.pr[pr,ind1]/model.info$beta[p])),ddf.pr=ddf.pr)
          
        #Calcul des dV/dlambda_m.V-1 (m=1...p)
   dvdlam<-vector(pr+2*nr,mode="list")
   dv2<-vector(pr+2*nr,mode="list")
   for(m in 1:pr) {
      dv1<- df.pr[[ll]][m,ind1]%*%t(df.pr[[ll]][m,ind1])
      dv2[[m]]<-dv1  
      dvdlam[[m]]<-dv1%*%invvar 
   }
   for(m in 1:p) {
      dv1<- df.pr[[ll]][m,ind1]%*%t(df.pr[[ll]][m,ind1])
      dv1bis<-blockdiag(sapply(1:n_occ,function(mat,n_occ) list((mat)%*%t(mat)),mat=(df.pr[[ll]][m,ind1])[1:(length(ind1)/n_occ)])) 
      dv2[[m]]<-dv1  
      dvdlam[[m]]<-dv1%*%invvar 
      dv2[[p+m]]<-dv1bis  
      dvdlam[[p+m]]<-dv1bis%*%invvar 
   }
        
        #Calcul des dV/dlambda_m.V-1 (m>p)
           vecmi<-2*(rep(tab$sig.i[ind1[1:(length(ind1)/n_occ)]],n_occ)+rep(tab$sig.s[ind1[1:(length(ind1)/n_occ)]],n_occ)*f.pr[[ll]][ind1])
           vecms<-vecmi*f.pr[[ll]][ind1]
           for(l in 1:(2*nr))
              dvdlam[[pr+l]]<-diag(rep(0,length(vecmi)),nrow=length(vecmi))
           for(l in 1:nr) { #sigma.inter
              if (model.info$sigmainter[[l]]>0) {
               diag(dvdlam[[pr+(l-1)*2+1]])[rep(tab$irep[ind1[1:(length(ind1)/n_occ)]],n_occ)==l]<-vecmi[rep(tab$irep[ind1[1:(length(ind1)/n_occ)]],n_occ)==l]
               dvdlam[[pr+(l-1)*2+1]]<-dvdlam[[pr+(l-1)*2+1]]%*%invvar
               }
              }
           for(l in 1:nr) { #sigma.slope
              if (model.info$sigmaslope[[l]]>0) { 
               diag(dvdlam[[pr+(l-1)*2+2]])[rep(tab$irep[ind1[1:(length(ind1)/n_occ)]],n_occ)==l]<-vecms[rep(tab$irep[ind1[1:(length(ind1)/n_occ)]],n_occ)==l]
               dvdlam[[pr+(l-1)*2+2]]<-dvdlam[[pr+(l-1)*2+2]]%*%invvar
               }                                                
              }
           resB<-matrix(nrow=pr+nr*2,ncol=pr+nr*2)     #matrice al?atoire #
           resB2<-matrix(nrow=pr+nr*2,ncol=pr+nr*2) 
        resC<-matrix(rep(0,(pr+nr*2)*(p.cov+model.info$locc)),nrow=pr+nr*2,ncol=(p.cov+model.info$locc))     #matrice incomplete option 1 covariance==0#
        #Division par beta dns le cas Trand==2     
  if (Trand==2) df1<-sensif[[ll]][,ind1]/c(model.info$beta,rep(1,model.info$locc+length(unlist(info.covariate$beta.covariate)))) 
     else df1<-sensif[[ll]][,ind1] 
     
     resA<-(df1 %*% invvar) %*% t(df1)  #matrice parametres fixe#
          
           for(i in 1:(pr+2*nr)) {
              for(j in 1:(pr+2*nr)) {
                resB[i,j]<-sum(diag(dvdlam[[i]]%*%dvdlam[[j]]))
                }
           }
           
       
  
        
        
#B final
        
           resBf<-resB/2
      
        #Enlever les lignes de B et C pour lesquelles omega ou sigma=0
           if(length(model.info$nullpar)>0){
             resBf<-resBf[-c(model.info$nullpar),-c(model.info$nullpar)]
             resC<-resC[-c(model.info$nullpar),]
           }
        
             fish<-cbind(rbind(resA,resC),rbind(t(resC),resBf))
             somme1<-somme1+fish*categories.prop*model.info$comb.seq.prop[sequ] 
            
    } 
 } 
 #fixed parameters if beta.fixed	== T
    beta.estim<-!model.info$beta.fixed
    if (dim(somme1)[1]>length(model.info$beta)) vec<-c(beta.estim,rep(T,dim(somme1)[1]-length(model.info$beta)))  else vec<-beta.estim
    somme1<-somme1[vec,vec] 
    
   return(somme1)
}


#Population, individual or bayesian information matrix
fisher.cv<-function(somme,pp,SIG,omega1,beta) {
#somme=FIM for the population design
   
   #if (FIM!="B") {beta.estim<-!beta.fixed ; p<-length(beta)-length(which(beta.estim==F))}
   #else {beta.estim<-!beta.fixed[diag(omega)!=0]; p<-length(beta)-length(which((beta.fixed+(diag(omega)==0))>0))}
   
   p<-length(beta)
   
   inv<-try(solve(as.matrix(somme)))
   if(!is.null(attributes(inv)$class)) {
      se<-rep(NA,pp)
      cv<-se
   }      
   else {
      se<-sqrt(diag(inv)) 
      cv1<-se[1:p]/beta*100
     	if (FIM=="P") {cv2<-se[(p+1):(pp-length(SIG[SIG!=0]))]/diag(as.matrix(omega1))*100} 
      cv3<-se[(pp-(length(SIG[SIG!=0])-1)):pp]/SIG[SIG!=0]*100 
      
    if (FIM=="P") cv<-abs(c(cv1,cv2,cv3))	
    if (FIM=="I") cv<-abs(c(cv1,cv3))
    if (FIM=="B") cv<-abs(c(cv1))
      }
   l<-list(inv,se,cv)
   return(l)
}

#Compute the population information matrix for the final protocol (case ind=1)
ecofisher.meth1.final<-function(lprot,tab,f.p,df.p,subjects,model.info,
  cost,option,ddf.p,info.covariate){
   Nbsubjects<-cost/sum((unlist(lapply(lprot,length))*subjects)) 
   nr<-model.info$nr
   p<-length(model.info$beta)
   
   if (info.covariate$covariate.model==T){
   p.cov<-length(c(model.info$beta,unlist(info.covariate$beta.covariate)))
   pp<-(p+p.cov)+2*nr-length(model.info$nullpar)-length(which(beta.fixed==T)) #2 parameters for the residual variance by response#
   }else
   {
   if (FIM=="P") pp<-2*p+2*nr-length(model.info$nullpar)-length(which(beta.fixed==T)) #2 parameters for the residual variance by response#
   if (FIM=="I") pp<-p+2*nr-length(model.info$nullsig)-length(which(beta.fixed==T))
	 if (FIM=="B") pp<-p-length(which((beta.fixed+(diag(omega)==0))>0))#-length(which(diag(model.info$omega)==0))
   }
   somme<-matrix(c(rep(0,pp*pp)),ncol=pp) 
   se<-numeric() 
   cv<-numeric()
   if (info.covariate$covariate.model==T)
   {
   for(i in 1:length(lprot)) {   
      ind1<-lprot[[i]]    
      fish<-ecofisher.meth1.cov(ind1,tab,f.p,df.p,model.info,option,ddf.p,info.covariate)
      somme<-somme+subjects[i]*Nbsubjects*fish
     
   }
   }else
   {
   for(i in 1:length(lprot)) {
      ind1<-lprot[[i]]
      if (FIM=="P"){
      fish<-ecofisher.meth1(ind1,tab,f.p,df.p,model.info,option,ddf.p)
      
      if (previous.FIM == "") somme<-somme+subjects[i]*Nbsubjects*fish
      else {
       previous_information<-read.table(file.path(directory,previous.FIM),header=FALSE,sep="")
       somme<-somme+subjects[i]*Nbsubjects*fish + as.matrix(as.matrix(previous_information))
       } 
      } 
      if (FIM=="I"){
      fish<-indifisher.meth1(ind1,tab,f.p,df.p,model.info)
      somme<-somme+subjects[i]*Nbsubjects*fish
      } 
      if (FIM=="B"){
      somme<-bayfisher.meth1(ind1,tab,f.p,df.p,model.info)
      } 
      
   }
   }
   #in "somme": already removing omega and sigma = 0
   
   
   sigmaslopen<-lapply(1:nr,function(nr,x) x[nr],x=paste("sig.slope",
     LETTERS[1:nr],sep=""))
   sigmaintern<-lapply(1:nr,function(nr,x) x[nr],x=paste("sig.inter",
     LETTERS[1:nr],sep=""))
   xv<-rbind(model.info$sigmainter,model.info$sigmaslope)
   xn<-rbind(sigmaintern,sigmaslopen)
   SIG<-unlist(c(xv))     
   
   vec<-diag(model.info$omega)!=0 
   omega1<-as.matrix(model.info$omega[,vec])[vec,]
  
   
    #shrinkage
    sh<-c()
	  
 	  if (FIM=="B"){
    #shrinkage matrix before removing the fixed parameters
    if (Trand==1) I_W<-solve(as.matrix(somme))*solve(as.matrix(omega1))[!beta.fixed[diag(omega)!=0],!beta.fixed[diag(omega)!=0]]
    else I_W<-solve(as.matrix(somme))*solve(as.matrix(diag((beta[diag(omega)!=0]),length(beta[diag(omega)!=0])))%*%as.matrix(omega1)%*%as.matrix(diag(((beta[diag(omega)!=0])),length(beta[diag(omega)!=0]))))[!beta.fixed[diag(omega)!=0],!beta.fixed[diag(omega)!=0]]
    
    #shrinkage after removing fixed parameters
	  ##sh<-diag(as.matrix(I_W[vec,vec]))*100
	  sh<-diag(as.matrix(I_W))*100
    }

    ##somme<-somme[vec,vec]
#Matrice d?j? de la bonne dimension
   det<-det(as.matrix(somme))    #Fisher determinant
#det<-crit(criterion,somme,pp)
  pp<-dim(as.matrix(somme))[1] 
   if (det==0) crit<-1e-99   else crit<-(det)^(1/pp)     #dans le cas o? le protocole est degenere
   invcrit<-1/crit 
   
   if (FIM=="B") f<-fisher.cv(somme,pp,SIG[SIG!=0],omega1,model.info$beta[diag(model.info$omega)!=0][!model.info$beta.fixed[diag(model.info$omega)!=0]])
   else{  
   if (info.covariate$covariate.model==T)
   f<-fisher.cv(somme,pp,SIG[SIG!=0],omega1,c(model.info$beta[!model.info$beta.fixed],unlist(info.covariate$beta.covariate)))
   else
   f<-fisher.cv(somme,pp,SIG[SIG!=0],omega1,model.info$beta[!model.info$beta.fixed])
   }
   
   
  
   inv<-f[[1]]
   se<-f[[2]]
   cv<-f[[3]]
   l<-list(somme=somme,inv=inv,se=se,cv=cv,det=det,crit=crit,
     pp=pp,subjects=subjects,xn=xn,SIG=SIG,sh=sh)
   return(l)
}


fisher.cv.iov<-function(somme,pp,SIG,OMEGA1,beta) {
 
#somme=FIM for the population design
    p<-length(beta)
  
   inv<-try(solve(as.matrix(somme)))
   if(!is.null(attributes(inv)$class)) {
      se<-rep(NA,pp)
      cv<-se
   }      
   else {
      se<-sqrt((diag(inv))) 
      
      cv1<-se[1:p]/beta*100   
      cv2<-se[(p+1):(pp-length(SIG[SIG!=0]))]/diag(as.matrix(OMEGA1))*100 
      cv3<-se[(pp-(length(SIG[SIG!=0])-1)):pp]/SIG[SIG!=0]*100  
      cv<-abs(c(cv1,cv2,cv3))	
      }
   l<-list(inv,se,cv)
   return(l)
}



#Compute the population information matrix for the final protocol (case ind=1)
ecofisher.meth1.final.iov<-function(lprot,tab,sensi,sensia,subjects,model.info,
  cost,option,info.covariate){
  Nbsubjects<-cost/sum((unlist(lapply(lprot,length))*subjects))
   #Nbsubjects<-init.info$cost/sum((unlist(lapply(lprot,length))*subjects)) 
   nr<-model.info$nr
   p<-length(model.info$beta)
   #pour prendre en compte les variances des effets al?atoires li?s ? IOV
  if (model.info$n_occ==1) model.info$gamma<-diag(NULL)
  OMEGA<-model.info$OMEGA #la matrice OMEGA qui regroupe omega et gamma
   n_occ<-model.info$n_occ
   if (length(which(model.info$gamma==0))==length(model.info$gamma)) n_occ<-1
   if (n_occ==1) pr<-p else pr<-2*p
  

   if (info.covariate$covariate.model==T){
   p.cov<-length(c(model.info$beta,unlist(info.covariate$beta.covariate)))
   pp<-(p.cov+model.info$locc+pr)+2*nr-length(model.info$nullpar)-length(which(beta.fixed==T)) #2 parameters for the residual variance by response#
   }else
   {
    pp<-p+model.info$locc+pr+2*nr-length(model.info$nullpar)-length(which(beta.fixed==T))
    #2 parameters for the residual variance by response#
   }
   
   somme<-matrix(c(rep(0,pp*pp)),ncol=pp) 
   se<-numeric() 
   cv<-numeric()
   
   if (info.covariate$covariate.model==T)
    {
      for(i in 1:length(lprot)) {
        ind1<-lprot[[i]]
        if (model.info$n_occ>1 && length(which(gamma==0))!=length(gamma)){
        ind1<-lprot[[i]]
        for (occ in 2:n_occ) {ind1<-c(ind1,lprot[[i]]+(occ-1)*dim(tab)[1])}
      }
      
      
      fish<-ecofisher.meth1.cov.iov(ind1,tab,sensi,sensia,model.info,option,info.covariate)
        somme<-somme+subjects[i]*Nbsubjects*fish
     
      }
    }else
    {
      for(i in 1:length(lprot)) {
      ind1<-lprot[[i]]
      if (model.info$n_occ>1 && length(which(gamma==0))!=length(gamma)){
        ind1<-lprot[[i]]
        for (occ in 2:n_occ) {ind1<-c(ind1,lprot[[i]]+(occ-1)*dim(tab)[1])}
      }
      fish<-ecofisher.meth1.iov(ind1,tab,sensi,sensia,model.info,option,info.covariate)
      
      somme<-somme+subjects[i]*Nbsubjects*fish 
      }
    }
   sigmaslopen<-lapply(1:nr,function(nr,x) x[nr],x=paste("sig.slope",
     LETTERS[1:nr],sep=""))
   sigmaintern<-lapply(1:nr,function(nr,x) x[nr],x=paste("sig.inter",
     LETTERS[1:nr],sep=""))
   xv<-rbind(model.info$sigmainter,model.info$sigmaslope)
   xn<-rbind(sigmaintern,sigmaslopen)
   SIG<<-unlist(c(xv))
    
   OMEGABIS<-diag(c(diag(model.info$omega),diag(model.info$gamma)))
		if (n_occ>1) vec1<<-diag(OMEGABIS)!=0 
		if (n_occ==1) vec1<<-diag(omega)!=0 
    OMEGA1<<-as.matrix(OMEGABIS[,vec1])[vec1,] 
   	omega1<<-as.matrix(model.info$omega[,diag(model.info$omega)!=0])[diag(model.info$omega)!=0,]
   	gamma1<<-as.matrix(model.info$gamma[,diag(model.info$gamma)!=0])[diag(model.info$gamma)!=0,]
    if (n_occ==1) OMEGA1<-omega1
    #vec<-c(diag(model.info$omega),diag(model.info$gamma))!=0 
   #OMEGA1<-as.matrix(OMEGABIS[,vec])[vec,]
   
   
    
    pp<-dim(as.matrix(somme))[1] 
#Matrice d?j? de la bonne dimension
   det<-det(as.matrix(somme))    #Fisher determinant
#det<-crit(criterion,somme,pp)
   if (det==0) crit<-1e-99   else crit<-(det)^(1/pp)     #dans le cas o? le protocole est degenere
   invcrit<-1/crit 
     
   if (info.covariate$covariate.model==T)
   f<-fisher.cv.iov(somme,pp,SIG[SIG!=0],OMEGA1,c(model.info$beta[!beta.fixed],unlist(info.covariate$beta.covariate),model.info$beta.occ1))
   else   
   f<-fisher.cv.iov(somme,pp,SIG[SIG!=0],OMEGA1,c(model.info$beta[!beta.fixed],model.info$beta.occ1))
   inv<-f[[1]]
   se<-f[[2]]
   cv<-f[[3]]
   l<-list(somme=somme,inv=inv,se=se,cv=cv,det=det,crit=crit,
     pp=pp,subjects=subjects,xn=xn,SIG=SIG)
   return(l)
}


############################################################################
####                        Generate the elementary protocols                ####
############################################################################
# Objectif : passer d'une combinatoire sur les temps de protocole ? une
# combinatoire sur les indices de protocole

########################################
#Auxiliary functions to combine lists of protocols (several responses)
# Combine 2 protocols from different sampling windows
combn.protoc<-function(list1,list2) {
  list3<-list()
  ni<-1
  for(i in 1:length(list1)) {
    for(j in 1:length(list2)) {
       list3[[ni]]<-c(list1[[i]],list2[[j]])
       ni<-ni+1
    }
  }
  return(list3)
}
#Combine protocols for 2 responses, identical.times=F
combn.repo<-function(list1,list2,lirep,irep2) {
  list3<-list()
  lirep2<-list()
#  cat("length",length(list2),"\n")
  ni<-1
  for(i in 1:length(list1)) {
    for(j in 1:length(list2)) {
       list3[[ni]]<-c(list1[[i]],list2[[j]])
       lirep2[[ni]]<-c(lirep[[i]],rep(irep2,length(list2[[j]])))
       ni<-ni+1
    }
  }
  return(list(protoc=list3,lirep=lirep2))
}
#Combine protocols for 2 responses, identical.times=T
combn.repo.idT<-function(list1,list2,lirep,irep2) {
  list3<-list()
  lirep2<-list()
  for(i in 1:length(list1)) {
       list3[[i]]<-c(list1[[i]],list2[[i]])
       lirep2[[i]]<-c(lirep[[i]],rep(irep2,length(list2[[i]])))
  }
  return(list(protoc=list3,lirep=lirep2))
}

##################################################
#Collect information about sampling design
collect.sampling.info<-function(interactive=F,sampl.info,model.info) {        
   ipb<-0
   nr<-model.info$nr
   #length(unlist(lapply(1:nr, function(nr,nwind)!is.numeric(nwind[[nr]]),nwind=nwind))==T)==0
   if(length(unlist(lapply(1:nr,function(nr,nwind) !is.numeric(nwind[[nr]]),
      nwind=sampl.info$nwind))==T)==0) 
      {
      cat("Please give the number of separate","sampling windows available.\n")
      if (length(sampl.info$nwind)!=nr) 
        cat("nwind must be a list of the same size as the number of responses")
      ipb<-1
      nwind<-0
   } else { 
      vnsamp<-lapply(1:nr,function(nr,nsamp=nsamp,nwind=nwind) 
        c(nwind[[nr]]==length(nsamp[[nr]])),nsamp=sampl.info$nsamp,nwind=sampl.info$nwind)  
      if((length(unlist(lapply(1:nr,function(nr,nsamp)!is.numeric(nsamp[[nr]]),
        nsamp=sampl.info$nsamp))==T)==0)||(all(unlist(vnsamp))!=T)) {
          cat("Please give","the number of points to be taken in each sampling window.\n",
           "This must be a list with one vector for each window \n",
           "If the vector has several values, elementary protocols with different",
           "number of points will be generated. \n")
           ipb<-1
           }
      vsampwin<-lapply(1:nr,function(nr,sampwin,nwind) 
       c(nwind[[nr]]==length(sampwin[[nr]])),sampwin=sampl.info$sampwin,
       nwind=sampl.info$nwind)
      if((length(unlist(lapply(1:nr,function(nr,sampwin)
       !is.numeric(sampwin[[nr]]),sampwin=sampl.info$sampwin))==T)==0) ||
       (all(unlist(vsampwin))!=T)) {
        cat("Please give the sampling points for each sampling window.\n",
        "This must be a list with one vector of values for each window\n")
        ipb<-1
      }

      if((length(unlist(lapply(1:nr, function(nr,nmaxpts)
        !is.numeric(nmaxpts[[nr]]),nmaxpts=sampl.info$nmaxpts))==T)==0)) {
        cat("Please give the maximum total",
        "number of points per elementary protocol.\n")
        ipb<-1
      }

      if((length(unlist(lapply(1:nr, function(nr,nminpts)
        !is.numeric(nminpts[[nr]]),nminpts=sampl.info$nminpts))==T)==0)) {
         cat("No minimum total number of points per elementary protocol given,",
           "assuming for the first response",sampl.info$nminpts,".\n",
           "etc","\n")
            ipb<-1
      }
       
      t<-sapply(sampl.info$nsamp,length)
      mgf<-lapply(1:nr,function(nr,t,nsamp) lapply(1:(t[nr]),
        function(t,nsamp) max(nsamp[[nr]][[t]]),nsamp=nsamp),nsamp=sampl.info$nsamp,t=t)
#nmax par reponse
      nmgf<-lapply(1:nr,function(nr,mgf) sum(unlist(mgf[[nr]])),mgf=mgf)
      if (all(unlist(nmgf)==unlist(sampl.info$nmaxpts))!=T){
        cat("nmaxpts is not consistent with the nsamp for one response \n")
        ipb<-1
      }
#nminpts
      migf<-lapply(1:nr,function(nr,t,nsamp) lapply(1:(t[nr]),
        function(t,nsamp) min(nsamp[[nr]][[t]]),nsamp=nsamp),nsamp=sampl.info$nsamp,t=t)
      nmigf<-lapply(1:nr,function(nr,migf) sum(unlist(migf[[nr]])),migf=migf)
      if (all(unlist(nmigf)==unlist(sampl.info$nminpts))!=T) {
        cat("nminpts is not consistent with the nsamp for one response \n")
        ipb<-1}        
      }
      if(ipb==1) return(list(ipb=1))

#print the input
   cat("----------------------------------------------------------------\n")
   cat("            SUMMARY OF THE SAMPLING INFORMATION \n")
   cat("----------------------------------------------------------------\n")
   for ( l in 1: nr){
     cat("Sampling windows for the response:",LETTERS[l],"\n")
     for (i in 1:sampl.info$nwind[[l]]) {
       cat("Window",i,": t=",sampl.info$sampwin[[l]][[i]],"\n")
       cat("              Nb of sampling points to be taken in this window, n[",
       i,"]=",sampl.info$nsamp[[l]][[i]],"\n") }
       
       if(is.na(sampl.info$nmaxpts[[l]])) {
         sampl.info$nmaxpts[[l]]<-sum(unlist(lapply(sampl.info$nsamp[[l]],max)))
         cat("Computing the maximum number of points per protocol.\n")}
       cat("Maximum total number of points in one elementary protocol :",
         sampl.info$nmaxpts[[l]],"\n")
       if(is.na(sampl.info$nminpts[[l]])) {
         sampl.info$nminpts[[l]]<-sum(unlist(lapply(sampl.info$nsamp[[l]],min)))
         cat("Computing the minimum number of points per protocol.\n")}
       cat("Minimum total number of points in one elementary protocol :",sampl.info$nminpts[[l]],"\n")
       cat("----------------------------------------------------------------\n")
   }
   return(list(ipb=0,sampl.info=sampl.info))
}

##########################################
# Function to generate elementary protocols, simplified (no second derivatives)
generate.multprot.meth1<-function(filename="matelem.tmp",interactive=F,
 model.info,struc.info,sampl.info,init.info,option,info.covariate,categories.cov) {
   cat("Computing model and first derivative values for all possible times.\n")
   cat("This may take a while depending on the complexity of the model.\n")
#create table with all protocol times
   tim.dos<-create.timdos(model.info,struc.info,sampl.info,init.info)
# compute derivatives
  if (model.info$n_occ==1) {
   if(struc.info$modelform=="AF") {
    if(NUM==F){
     y<-symbolic.prem(model.info,struc.info,sort(unique(tim.dos$iform)),option)
     df.forms<-y[[1]]
     ddf.forms<-y[[2]]
# compute model and derivative values for all sampling times
     y2<-calc.modelprem(tim.dos,model.info,struc.info,df.forms,ddf.forms,option,info.covariate,model.info$beta,rep(0,length(model.info$beta)))
     f.p<-y2[[1]]
     df.p<-y2[[2]]
     ddf.p<-y2[[3]]
   }else{ #numerical derivatives
     y2<-calc.model.num(tim.dos,model.info,struc.info,option,info.covariate,model.info$beta,rep(0,length(model.info$beta)))
     f.p<-y2[[1]]
     df.p<-y2[[2]]
     ddf.p<-NULL
   }
   } else {
#Compute model predictions and first derivatives
# when model is described by differential equations
     y2<-calc.model.de(tim.dos,model.info,init.info,struc.info,option,info.covariate,model.info$beta,rep(0,length(model.info$beta)))
     f.p<-y2[[1]]
     df.p<-y2[[2]]
     ddf.p<-y2[[3]]
   }
 }  
#derivatives for models with iov 
  if (model.info$n_occ>1) { 
   sensi<-sensifixe(tim.dos,model.info,struc.info,sampl.info,init.info,option,info.covariate)
   sensia<-sensialea(tim.dos,model.info,struc.info,sampl.info,init.info,option,info.covariate)
  } 
#collect information 
  ipb<-collect.sampling.info(F,sampl.info,model.info)
  if (ipb[[1]]==1) return(list(ipb=1))
  sampl.info<-ipb$sampl.info
#Generate all possible protocols
  listprot<-list()
# creation d'un vecteur avec les indices au lieu des temps, idem pour temps contraints
  indwin<-sampl.info$sampwin
  iwin<-sampl.info$sampwin
  lwin<-iwin
  for(i in 1:length(indwin))
     lwin[[i]]<-lapply(iwin[[i]],function(vec,irep,idos) paste(vec,irep,idos,
     sep="x"),irep=i,idos=1)
  for(i in 1:length(indwin))
     indwin[[i]]<-lapply(lwin[[i]],function(vec) match(vec,tim.dos$id))
   nr<-model.info$nr
   
#liste des protocoles combin?s : list.cbn
#liste des protocoles possibles pour toutes les r?ponses
#liste de m?me taille avec la r?ponse correspondante : list.rep
#Generate the list for the first dose
  list.nr<-list() 
  if(sampl.info$identical.times) {#how many possible times
    nl<-dim(tim.dos[tim.dos$irep==1 & tim.dos$condnum==tim.dos$condnum[1],])[1]
  }
  for(l in 1:nr) {
    if(sampl.info$identical.times & l>1) {
      list.nr[[l]]<-lapply(list.nr[[1]],function(vec) vec+nl*(l-1))
    } else {
#generate all possible protocols for each sampling window
    list1<-list()
    for(ik in 1:sampl.info$nwind[[l]]) {      #number of windows of sampling times
#    print(sampl.info$nwind[[l]])
      np1<-0
      list1[[ik]]<-list()
      for(ij in 1:length(sampl.info$nsamp[[l]][[ik]]))  {  #number of sampling times for the design
#        print(sampl.info$nsamp[[l]][[ik]])
        if (length(indwin[[l]][[ik]])!=1)
        iprotoc<-t(combn(indwin[[l]][[ik]],sampl.info$nsamp[[l]][[ik]][ij]))
        else #cas ou on a seulement une fenetre de temps compos? de 1 seul temps
        iprotoc<- t(indwin[[l]][[ik]])
        
#        print(iprotoc)
        nprot<-dim(iprotoc)[1]
        l1<-lapply(1:nprot,function(i,vec) vec[i,],vec=iprotoc)
        for(i in 1:nprot) list1[[ik]][[i+np1]]<-l1[[i]]
        np1<-np1+nprot
      }#fin boucle sur le nombre de prelevements possibles dans chaque fenetre
    }#fin boucle sur le nombre de fenetres
#generate all possible protocols for each sampling window
    if(length(list1)>1) {
      list3<-list1[[1]]
      for(i in 2:length(list1)) {
        list3<-combn.protoc(list3,list1[[i]])
      }
      list1<-list3
    } else list1<-list1[[1]]
#prune protocols to respect the min/max nb of samples
    lprot<-unlist(lapply(list1,length))
    nmin<-sampl.info$nminpts[[l]];nmax<-sampl.info$nmaxpts[[l]];n1<-0
    list2<-list()
    for(i in 1:length(list1)) {
      if((lprot[i]-nmin)*(lprot[i]-nmax)<=0) {
        n1<-n1+1
        list2[[n1]]<-list1[[i]]
        }
     }
    list1<-list2
    # Keep only protocols with times in fixed.times
#    cat("List of protocols\n")
    if(length(sampl.info$fixed.times[[l]])>0) {
#    	cat("Fixed times for response",l,": ",unlist(sampl.info$fixed.times[[l]]),"\n")
    	nfixed.times<-length(sampl.info$fixed.times[[l]])
    	list.nr2<-list()
    	itrack<-1
    	for(iprot in 1:length(list1)) {
    		if(sum(sampl.info$fixed.times[[l]] %in% tim.dos$times[list1[[iprot]]]) == nfixed.times) {
    			list.nr2[[itrack]]<-list1[[iprot]]
    			itrack<-itrack+1
    		}   			 
    	}   
    	list1<-list.nr2
    }
    list.nr[[l]]<-list1
    } #end loop on identical.times
  } #end loop on number of responses
#   cat("list.nr:",length(list.nr)," ")
#   for(j in 1:nr) cat(length(list.nr[[j]])," ")
#   cat("\n")
#If several responses, combine the responses
  list.rep<-lapply(1:length(list.nr[[1]]),
     function(i,lis) rep(1,length(lis[[i]])),lis=list.nr[[1]])
  if(nr>1) {
    list1<-list.nr[[1]]
    if(sampl.info$identical.times) {
       for(j in 2:nr) {
          list.cb<-combn.repo.idT(list1,list.nr[[j]],list.rep,j)
          list1<-list.cb[[1]]
          list.rep<-list.cb[[2]]
       }
       list.cbn<-list1
    } else {
       for(j in 2:nr) {
          list.cb<-combn.repo(list1,list.nr[[j]],list.rep,j)
          list1<-list.cb[[1]]
          list.rep<-list.cb[[2]]
       }
       list.cbn<-list1
    }
  } else list.cbn<-list.nr[[1]]
  
#Duplicating the list for several doses
   if(struc.info$modelform=="AF") {
     if(!init.info$dose.identical & length(init.info$doses)>1) {
       len1<-length(list.cbn)
       l2<-list.cbn
       l3<-list.rep
       nl<-dim(tim.dos[tim.dos$doses==init.info$doses[1],])[1]
       for(idos in 2:length(init.info$doses)) {
         l1<-lapply(list.cbn,function(vec) vec+(idos-1)*nl)
         for(i in 1:length(l1)) {
           l2[[(i+len1)]]<-l1[[i]]
           l3[[(i+len1)]]<-list.rep[[i]]
           }
         len1<-len1+length(l1)
         }
       list.cbn<-l2
       list.rep<-l3
     }
   } else { #cas DE
     if(!init.info$condinit.identical & length(init.info$condinit)>1) {
       len1<-length(list.cbn)
       l2<-list.cbn
       l3<-list.rep
       nl<-dim(tim.dos[tim.dos$condnum==tim.dos$condnum[1],])[1]
       for(idos in 2:length(init.info$condinit)) {
         l1<-lapply(list.cbn,function(vec) vec+(idos-1)*nl)
         for(i in 1:length(l1)) {
           l2[[(i+len1)]]<-l1[[i]]
           l3[[(i+len1)]]<-list.rep[[i]]
           }
         len1<-len1+length(l1)
         }
       list.cbn<-l2
       list.rep<-l3
     }
   }
   
    cat("Now evaluating the Fisher Information Matrix for the",length(list.cbn),
    "protocols generated \n")
    cat(" This may take a while depending on the number of protocols.\n",
    " Please wait for the program to write the matrices in a file.\n")


#Compute the Fisher information matrix for each protocol
 if (length(which(model.info$gamma==0))==length(model.info$gamma))     lamb<-c(diag(model.info$omega))
 else  lamb<-c(diag(model.info$omega),diag(model.info$gamma))
 for(l in 1:nr) lamb<-c(lamb,model.info$sigmainter[[l]],model.info$sigmaslope[[l]])
    
 sig<-c(model.info$sigmainter[[1]],model.info$sigmaslope[[1]])
 if (nr>1){for(l in 2:nr) sig<-c(sig,model.info$sigmainter[[l]],model.info$sigmaslope[[l]])} 
 
  nullpar<-which(lamb==0)
  nullsig<-which(sig==0)
  model.info$nullpar<-nullpar
  model.info$nullsig<-nullsig
  lfish<-list()
  tprot<-list()
   ldet<-list()
  lk<-0
   ibest<-1
   fbest<-NA
   
   for(kprot in 1:length(list.cbn)) {
    protk<-list.cbn[[kprot]]
    repk<-list.rep[[kprot]]
    tprot[[kprot]]<-tim.dos$times[list.cbn[[kprot]]]
    
    if (model.info$n_occ>1 && length(which(gamma==0))!=length(gamma)){
    protk<-list.cbn[[kprot]]
    for (occ in 2:n_occ) {protk<-c(protk,list.cbn[[kprot]]+(occ-1)*dim(tim.dos)[1])}}
    
    if (info.covariate$covariate.mode==T){ 
       if (model.info$n_occ>1) y<-try(ecofisher.meth1.cov.iov(protk,tim.dos,sensi,sensia,model.info,option,info.covariate))
       else
       {  
       y<-try(ecofisher.meth1.cov(protk,tim.dos,f.p,df.p,model.info,option,ddf.p,info.covariate)) 
       }
    } 
     
    if (info.covariate$covariate.mode==F){
       if (model.info$n_occ>1) y<-try(ecofisher.meth1.iov(protk,tim.dos,sensi,sensia,model.info,option,info.covariate))
       else  {
       if (FIM=="P") y<-try(ecofisher.meth1(protk,tim.dos,f.p,df.p,model.info,option,ddf.p))
       if (FIM=="I") y<-try(indifisher.meth1(protk,tim.dos,f.p,df.p,model.info))
       if (FIM=="B") y<-try(bayfisher.meth1(protk,tim.dos,f.p,df.p,model.info))
       }
    }
    
    if(mode(y)=="character") {
       lfish[[kprot]]<-NULL 
       ldet[[kprot]]<-NULL 
       #cat("Pb in the computation of the Fisher information matrix for protocol",kprot,"\n")
       lk<-1
     } else {
       lfish[[kprot]]<-y
  
  #Cyrielle     
  #taking into account previous information
    if (previous.FIM!= ""){
    previous_information<-read.table(file.path(directory,previous.FIM),header=FALSE,sep="")
    #On teste si la previous information matrix a la même taille que la matrice actuelle. Si ce n'est pas le cas, on arrête l'exécution et on affiche un message d'erreur : 
    if (nrow(previous_information) != nrow(y) | ncol(previous_information) != ncol(y) ) {stop ("The previous FIM has not the same size as the current FIM","\n") }
   y<-as.matrix((previous_information)/subjects)+as.matrix(y)   
   lfish[[kprot]]<-y
   }
  
   
   ################################################################ 
       
       
       ldet[[kprot]]<-det(as.matrix(y))
       if(is.na(fbest) || fbest<ldet[[kprot]]) {
       	ibest<-kprot
       	fbest<-ldet[[kprot]]
       	
       }
     }
  }
#Save to files
#Cyrielle
 if(previous.FIM == "") { 
  if(lk>0) { #remove NULL protocols
    ik<-1
    lis1<-tp1<-lfi1<-lrep1<-ldet1<-list()
  
    for(k in 1:length(list.cbn)) {
      if(!is.null(lfish[[k]])) {
        lis1[[ik]]<-list.cbn[[k]]
        tp1[[ik]]<-tprot[[k]]
        lfi1[[ik]]<-lfish[[k]]
        lrep1[[ik]]<-list.rep[[k]]
        ldet1[[ik]]<-ldet[[k]]
        if(zapsmall(ldet[[k]])-fbest) ibest<-ik
        ik<-ik+1
      }
    }
   list.cbn<-lis1
   tprot<-tp1
   lfish<-lfi1
   list.rep<-lrep1
    ldet<-ldet1
  }
  for(kj in 1:length(list.cbn)) {
      write(kj,filename,append=kj>1)
      write(length(list.cbn[[kj]]),filename,append=T)
      write(tim.dos[list.cbn[[kj]],2],filename,ncol=length(list.cbn[[kj]]),append=T)
      write.table(lfish[[kj]],filename,col.names=F,row.names=F,quote=F,append=T)
  }
   } # end loop on use.priorinfo 
   #MODIFS 1 GROUP
   cat("----------------------------------------------------------------\n")
   cat("The best protocol with one group is:\n")
   nsuj<-init.info$cost/(length(list.cbn[[ibest]])) 
   cat("  Protocol      nb       Times       \n")

	cat("    ",ibest,"        ",nsuj,"         ",tim.dos[list.cbn[[ibest]],2],"\n")
   
   ffish<-lfish[[ibest]]*nsuj
   pp<-dim(ffish)[1] # Dimension
   xcri<-(det(as.matrix(ffish)))**(1/pp)
   cat("   determinant=",ldet[[ibest]],"  criterion=",xcri,"\n")
   
  #best 1 group protocol
  ######################
  #prot, dose/condinit, idbest, numbest
  
   bestprot1<-dose1<-list()
  
   for (i in 1:nr){
   if (modelform=="AF"){
    bestprot1[[i]]<-as.matrix(tim.dos[list.cbn[[ibest]],c(2,5)][tim.dos[list.cbn[[ibest]],c(2,5)]$irep==i,])[,1]
    dose1[[i]]<- as.matrix(tim.dos[list.cbn[[ibest]],c(3,5)][tim.dos[list.cbn[[ibest]],c(3,5)]$irep==i,])[1,1]
    }
    else{
    bestprot1[[i]]<-as.matrix(tim.dos[list.cbn[[ibest]],c(2,4)][tim.dos[list.cbn[[ibest]],c(2,4)]$irep==i,])[,1]
    dose1[[i]]<-as.character(condinit)
    }
   }

   idbest<-list(list.cbn[[ibest]])
   numbest<-ibest
   #######################
   
  if (model.info$n_occ==1) return(list(ipb=0,protoc=list.cbn,tprotoc=tprot,fisher=lfish,
    sampl.info=sampl.info,tab=tim.dos,f.p=f.p,df.p=df.p,ddf.p=ddf.p,nullpar=nullpar,nullsig=nullsig,
    list.rep=list.rep,model.info=model.info,bestprot1=bestprot1,bestcrit1=xcri,dose1=dose1,nsuj=nsuj,idbest=idbest,numbest=numbest))
  
  if (model.info$n_occ>1) return(list(ipb=0,protoc=list.cbn,tprotoc=tprot,fisher=lfish,
    sampl.info=sampl.info,tab=tim.dos,sensi=sensi,sensia=sensia,nullpar=nullpar,nullsig=nullsig,
    list.rep=list.rep,model.info=model.info,bestprot1=bestprot1,bestcrit1=xcri,dose1=dose1,nsuj=nsuj,idbest=idbest,numbest=numbest))  
  
}

#FW under steroids
fedorov.wynn.sg<-function(interactive=F,
 fileres=file.path(directory,"details.r"),
 filematelem=file.path(directory,"matelem.tmp"),
 model.info,struc.info,sampl.info,init.info,option,info.covariate) {
   nr<-model.info$nr
   prots<-init.info$prots
   
#Generate matelem.tmp or collect sampling information
   ipb<-generate.multprot.meth1(filematelem,interactive,model.info=model.info,
   struc.info=struc.info,sampl.info=sampl.info,init.info=init.info,option=option,info.covariate=info.covariate)
#Information to be passed on (including modified structures)
 
   tab=ipb$tab;
   if (model.info$n_occ==1) {f.p=ipb$f.p;df.p=ipb$df.p;ddf.p=ipb$ddf.p}
   if (model.info$n_occ>1) {sensi=ipb$sensi;sensia=ipb$sensia}
   nullpar=ipb$nullpar; nullsig<-ipb$nullsig
   model.info<-ipb$model.info;sampl.info<-ipb$sampl.info
   list.rep<-ipb$list.rep
   protoc<-ipb$protoc
    dimsig<-lapply(1:nr,function(nr,sigmainter,sigmaslope) 
     sum(as.integer(unlist(sigmainter[[nr]])>0)+as.integer(unlist(sigmaslope[[nr]])>0)),
     sigmainter=model.info$sigmainter,sigmaslope=model.info$sigmaslope)
  
  if (info.covariate$covariate.model==T){
   ndim<-length(c(model.info$beta, unlist(info.covariate$beta.covariate)))+model.info$locc+sum(as.integer(diag(model.info$omega)>0)) + sum(as.integer(diag(model.info$gamma)>0))+sum(unlist(dimsig))-length(which(model.info$beta.fixed==T))
   } else
   {
   ndim<-length(model.info$beta)+sum(as.integer(c(diag(model.info$omega))>0)) + sum(as.integer(diag(model.info$gamma)>0))+model.info$locc+sum(unlist(dimsig))-length(which(model.info$beta.fixed==T))
   }
   
   

   

   
#Transform initial protocols in indices on tab elements
   initnp<-length(prots[[1]])
   iprot<-prots
   if(struc.info$modelform=="AF") {
     dprot<-match(init.info$doses,unique(init.info$doses))
     if(length(dprot)<initnp) #recycle doses
        dprot<-rep(dprot,initnp)
   } else {
     dprot<-c(1:length(init.info$condinit))
     if(length(init.info$condinit)<initnp) #recycle doses
        dprot<-rep(dprot,initnp)
   }
   for(l in 1:nr) {
     for(i in 1:length(iprot[[l]])) {
      iprot[[l]][[i]]<-unlist(lapply(iprot[[l]][[i]],function(t,tab,icond,l) 
      which((tab$times==t & tab$condnum==icond & tab$irep==l)),tab=tab,
      icond=dprot[i],l=l))
     }
   }
#if several responses, combine the protocols to create the initial protocols
   if(nr>1) {
      initindx<-inittps<-list()
      for(i in 1:length(iprot[[1]])) {
        vec<-vec2<-c()
        for(l in 1:nr) {
          vec<-c(vec,iprot[[l]][[i]])
          vec2<-c(vec2,prots[[l]][[i]])
        }
        initindx[[i]]<-vec
        inittps[[i]]<-vec2
      }
   } else {
      initindx<-iprot[[1]]
      inittps<-prots[[1]]
   }


#Running Fedorov
   datprob<-list(filematelem,fileres,initindx,init.info$subjects,inittps)      
   yini<-fedorov.init(datprob,ndim,protoc,tab,init.info)
#print("federovinit")
  # if(yini$ipb==1) return(list(ipb=1)) #now affects 1st protocol if pb
   init.np<-yini$zeprot[1]
   init.freq<-yini$freq
   init.num<-yini$zeprot[2:(init.np+1)]
# garder les conditions d'initialisations
   init.fed<-list(init.np=init.np,init.freq=init.freq[1:init.np],init.num=init.num)
   if (model.info$n_occ==1) {data.fed<-list(tab=tab,f.p=f.p,df.p=df.p,ddf.p=ddf.p,nullpar=nullpar,nullsig=nullsig)}
   if (model.info$n_occ>1) {data.fed<-list(tab=tab,sensi=sensi,sensia=sensia,nullpar=nullpar)}
   protoc.fed<-list(protoc=protoc,tprotoc=ipb$tprotoc)

   yfed.inp<-list(namfich=yini$namfich,ndimen=yini$ndimen,zeprot=yini$zeprot,
   freq=yini$freq,datprob=datprob,protoc=protoc,tprotoc=ipb$tprotoc)
   
   if (FIM=="P" & previous.FIM==""){
   yfed<-fedorov.alg(yfed.inp)
#note : protopti, protfreq, nq returned by yfed never used...
   np.f<-yfed$y[[3]]
   opti.fed<-list(opti.np=np.f,opti.freq=yfed$protfreq[1:np.f],
      opti.num=yfed$y[[4]][1:np.f],opti.tps=yfed$protopti)
# print(list(ipb=0,yfed=yfed,protopti=yfed$protopti,protfreq=yfed$protfreq))
  }
###BEST 1 GROUP
bestprot1<-ipb$bestprot1
bestcrit1<-ipb$bestcrit1
dose1<-ipb$dose1
nsuj<-ipb$nsuj
idbest<-ipb$idbest
numbest<-ipb$numbest
##########


   
   if (FIM=="P" & previous.FIM=="") return(list(yfed=yfed,protopti=yfed$protopti,protfreq=yfed$protfreq,
   protoc.fed=protoc.fed,init.fed=init.fed,data.fed=data.fed,opti.fed=opti.fed,
   sampl.info=sampl.info,model.info=model.info,bestprot1=bestprot1,bestcrit1=bestcrit1,dose1=dose1,nsuj=nsuj,idbest=idbest,numbest=numbest))
  
   else return(list(protoc.fed=protoc.fed,init.fed=init.fed,data.fed=data.fed,
   sampl.info=sampl.info,model.info=model.info,bestprot1=bestprot1,bestcrit1=bestcrit1,dose1=dose1,nsuj=nsuj,idbest=idbest,numbest=numbest))
}

##############################################################################
#Fonction pour faire tourner Fedorov avec initialisation
fedorov.init<-function(datprob,ndim,protoc,tab,init.info) {
#datprob[[1]]--[[2]] :2 chaines de caract?res
# datprob[[3]]: index des tps pour chaque protocole (selon tab)
# datprob[[4]]: le vecteur de fr?quences
# datprob[[5]]: tps pour chaque protocole
   namfich<-c(datprob[[1]],datprob[[2]])
   vectps<-rep(0,ndim*4);
   fisher<-rep(0,ndim*(ndim+1)/2)
#Protocoles ?lementaires constituant le protocole de population initial
   initindx<-datprob[[3]] 
   zefreq<-datprob[[4]]
   freq<-rep(0,ndim*2)
   inittps<-datprob[[5]]
   ndimen<-c(length(protoc),ndim,init.info$cost)
# d?terminer quels sont les protocoles initiaux => idx.initprot
   idx.initprot<-c()
   for(i in 1:length(initindx)) {
     vec<-unlist(lapply(protoc,function(v1,vec) identical(as.integer(v1),vec),
        vec=as.integer(initindx[[i]])))
     idx.initprot<-c(idx.initprot,which(vec)[1])
   }
   
   kna<-which(is.na(idx.initprot))
   if(length(kna)>0) {
     for(i in kna)
        cat("Initial protocol",inittps[[i]],
        "does not correspond to the constraints, discarding it.\n")
   }
   kdup<-which(!is.na(idx.initprot) & duplicated(idx.initprot))
   if(length(kdup)>0) {
     for(i in kdup)
        cat("You have entered initial protocol",inittps[[i]],
        "more than one time, collapsing identical protocols in the input.\n")
   }
   kprot<-unique(idx.initprot[!is.na(idx.initprot)])
   if(length(kprot)==0) {
     cat("The initial protocols are not correctly specified.\n")
     cat("They must correspond to protocols possible given the constraints on time in the sampling windows.\n")
     cat("Replacing initial protocol by the first possible protocol\n")
     kprot<-c(1)
     zefreq<-c(1)
     cat("Protocol",kprot," times=",tab$times[protoc[[kprot]]],"\n")
     ipb<-1
   } else {
     ipb<-0
     idx.initprot[is.na(idx.initprot)]<-0
     for(ip in 1:length(kprot))
       freq[ip]<-sum(zefreq[kprot[ip]==idx.initprot])
     if(sum(freq)!=1) {
       cat("Standardising the frequencies of the initial protocols to a sum of 1.\n")
       xcal<-sum(freq)
       freq<-freq/xcal
     }
   }
   zeprot<-c(length(kprot),kprot,rep(0,(ndim*2-length(kprot))))
#   cat(zeprot[1],"initial elementary protocols:",zeprot[2:(1+length(kprot))],"\n")
#   cat("Frequencies:",freq[1:length(kprot)],"\n")
#   cat("Dimensions of the problem:",ndimen,"\n")
   yli<-list(ipb=ipb,namfich=namfich,ndimen=ndimen,zeprot=zeprot,freq=freq)
   return(yli)
}

fedorov.alg<-function(inputcheck) {
#inputcheck:liste de 2 chaines de caract?res, un entier, une liste de vecteurs
#de temps, et un vecteur de fr?quences
#   cat("Checking input\n")
   #print(inputcheck)
   namfich<-inputcheck[[1]] #fichiers matelem et r?sultats
   ndimen<-inputcheck[[2]] #dimensions du pb (nb prot elem, ndim, cout)
   zeprot<-inputcheck[[3]] #nb puis num?ros des protocoles ?l?mentaires ds prot init
   zefreq<-inputcheck[[4]] #fr?quences protocoles ?l?mentaires
   datprob<-inputcheck[[5]]
   protoc<-inputcheck[[6]]
   tprotoc<-inputcheck[[7]]
   np.init<-zeprot[1] #nb de protocoles ?l?mentaires ds protocole initial
   ndim<-ndimen[2]
#Dimensionnement des vecteurs
   numprot<-rep(0,ndim*2);
   nbdata<-numprot;
   freq<-numprot
   vectps<-rep(0,length(numprot)*max(unlist(lapply(protoc,length))));
   fisher<-rep(0,ndim*(ndim+1)/2)
   nok<-0
   cat("----------------------------------------------------------------\n")
   cat("Initial protocol\n")
   cat("----------------------------------------------------------------\n")
   cat("NUMBER OF GROUPS   ",np.init,"\n")
   cat("   Freq.   Protocol nb       Times       \n")
   for(i in 1:np.init) {
      cat(i,".   ",format(zefreq[i],width=5),"    ",format(zeprot[i+1],width=5),
        "           ")
      cat(tprotoc[[zeprot[[i+1]]]],"\n")
   }
   cat("\n")
   cat("----------------------------------------------------------------\n")

   y<-.C("FedorovInit_R",as.character(namfich),as.integer(ndimen),
   as.integer(np.init),as.integer(numprot),as.double(freq),as.integer(nbdata),
   as.double(vectps),as.double(fisher),as.integer(nok),as.integer(zeprot),
   as.double(zefreq))

   
   if(y[[9]]>=0){ #Test si l'optimisation s'est bien pass?e
     cat("----------------------------------------------------------------\n")
     cat("Optimisation results\n",y[[9]]," iterations \n")
     cat("----------------------------------------------------------------\n")
     cat("NUMBER OF GROUPS   ",y[[3]],"\n")
     cat("   Freq.   Protocol nb       Times       \n")
     ij<-1
     for(i in 1:y[[3]]) {
       cat(i,".   ",format(y[[5]][i],width=5),"    ",format(y[[4]][i],width=5),
       "           ")
       for(j in 1:y[[6]][i]) {
         cat(format(y[[7]][ij],width=5)," ")
         ij<-ij+1
       }
       cat("\n")
     }
     cat("----------------------------------------------------------------\n")
     nq<-c(y[[6]][1:y[[3]]])
     protopti<-vector(y[[3]],mode="list")
     protfreq<-y[[5]][1:y[[3]]]
     n1<-1
     for(i in 1:y[[3]]) {
       protopti[[i]]<-c(y[[7]][n1:(n1+y[[6]][i]-1)])
       n1<-n1+y[[6]][i]
     }
     if(y[[9]]==0) {
       cat("The program has returned after 0 iterations.\n")
       cat("It is possible that the initial population protocol was optimal,\n")
       cat("but if any other error message appeared it is more likely\n")
       cat("that it was too poor for the program to be able to start.\n")
       cat("We suggest running PFIM again after adding more protocols to\n")
       cat("the initial population protocol.\n")
     }
     return(list(y=y,protopti=protopti,protfreq=protfreq,nq=nq))
   } else 
     return(y=y) # Problem in optimisation
}


###########################################################################
###                          Function for graph                          ###                                                                                                  #                                                             
###########################################################################

#function to build graph, model in analytical form, expression
graph.fw.afs<-function(gr.tps,gr.list,tab,icond,graph.info,struc.info,model.info,
  init.info,info.covariate) {
   p<-length(model.info$beta)
   nr<-model.info$nr
   for (i in 1:p) 
   assign(parameters[i],model.info$beta[i])
   graphinf<-graph.info$graphinf
   graphsup<-graph.info$graphsup
   tt1<-list()
   #tt<-lapply(1:nr,function(nr,graphinf,graphsup) seq(graphinf,graphsup,0.1),graphinf<-graphinf[[nr]],graphsup<-graphsup[[nr]])   
     for (l in 1:nr){
       tt1[[l]]<-list()
       if (length(struc.info$bound[[l]])==1){
         if( struc.info$bound[[l]][[1]][1]> graphsup[[l]]) {
            warning(" The upper limit for the graph  ",l,
              "  is too small, cannot draw graph \n")
           return(0)
         }
         tt1[[l]][[1]]<-seq(struc.info$bound[[l]][[1]][1],graphsup[[l]],(graphsup[[l]]-struc.info$bound[[l]][[1]][1])/1000)
         
       } else {
       for (i in 1: (length(struc.info$bound[[l]])-1)){
          tt1[[l]][[i]]<-seq(struc.info$bound[[l]][[i]][1],struc.info$bound[[l]][[i]][2],(struc.info$bound[[l]][[i]][2]-struc.info$bound[[l]][[i]][1])/1000)
          if(i==length(struc.info$bound[[l]])-1){
             if (struc.info$bound[[l]][[i]][2]>graphsup[[l]]) {
                warning(" The upper limit (graph.sup)  for the response  ",
                 LETTERS[ l],"  is too small")
             return(0)
              }
             tt1[[l]][[i+1]]<-seq(struc.info$bound[[l]][[i]][2],graphsup[[l]],(graphsup[[l]]-struc.info$bound[[l]][[i]][2])/1000)
          }
       }
     }
    }
   tt1b<-lapply(1:nr,function(nr,tt1) unlist(tt1[[nr]]),tt1<-tt1)
   tt1<-lapply(1:nr,function(nr,tt1b) tt1b[[nr]][tt1b[[nr]]>=graphinf[[nr]][1] & tt1b[[1]]<=graphsup[[nr]][1]],tt1b<-tt1b) 
   ff1<-list()
   if ((init.info$dose.identical==T) | length(unique(icond))==1) {
     idos<-which(c(1:length(init.info$doses))==icond[1])
     dose<-init.info$doses[icond[1]] 
     par(mfrow=c(1,nr))
     for (l in 1:nr){
        ff1[[l]]<-vector(mode="numeric", length(tt1[[l]]))
        for (j in 1:length(tt1[[l]])) {
          t<-tt1[[l]][j]
       
        if(info.covariate$covariate.model==T)
          ff1[[l]][j]<-eval(struc.info$formg_init[[(l+nr*(idos-1))]][t<=struc.info$tf[[l]]][1]) 
           else
          ff1[[l]][j]<-eval(struc.info$formg[[(l+nr*(idos-1))]][t<=struc.info$tf[[l]]][1]) 
        }
        if(graph.info$log.logical!=F) 
          plot(tt1[[l]],ff1[[l]],type='l',xlab=paste(graph.info$names.datax[l]),ylab=paste(graph.info$names.datay[l]),
            log=graph.info$log.logical,ylim=c(graph.info$y.range[[l]],range(ff1[[l]]))[1:2]) else 
          plot(tt1[[l]],ff1[[l]],type='l',xlab=paste(graph.info$names.datax[l]),ylab=paste(graph.info$names.datay[l]),
            ylim=c(graph.info$y.range[[l]],range(ff1[[l]]))[1:2])
            
        #points for the design on the graph
        for (i in 1:length(gr.tps)) {
           tps<-gr.tps[[i]][gr.list[[i]]==l]
           tab1<-tab[tab$irep==l & tab$condnum==icond[1],]
           if (info.covariate$covariate.model==T || model.info$n_occ!=1)
           {      f.tps<-vector() 
                 for (ij in 1: length(tps))
                  {
                    t<-tps[ij]
                    if (model.info$n_occ==1) f.tps[ij]<-eval(struc.info$formg_init[[(l+nr*(idos-1))]][t<=struc.info$tf[[l]]][1])
                    else f.tps[ij]<-eval(struc.info$formg[[(l+nr*(idos-1))]][t<=struc.info$tf[[l]]][1])
                  }
            }else  f.tps<-tab1$f.p[match(tps,tab1$times)]
           
           points(tps,f.tps,pch=paste(i),cex=2,col=paste(i))
        }
      
        title(paste(graph.info$names.datay[l],"model ","\n Dose = ",
          dose,sep=" "))
      }
   } else {
#Regroup doses
     uni.cond<-sort(unique(icond))
     par(mfrow=c(length(uni.cond),nr))
     for(kd in 1:length(uni.cond)) {
#       dose<-doses[unique(icond)[kd]]
       idos<-which(c(1:length(init.info$doses))==uni.cond[kd])
       dose<-init.info$doses[uni.cond[kd]]
       for (l in 1:nr){
          ff1[[l]]<-vector(mode="numeric", length(tt1[[l]]))
          for (j in 1:length(tt1[[l]])) {
            t<-tt1[[l]][j]
            if(info.covariate$covariate.model==T)
            ff1[[l]][j]<-eval(struc.info$formg_init[[(l+nr*(idos-1))]][t<=struc.info$tf[[l]]][1]) else
            ff1[[l]][j]<-eval(struc.info$formg[[(l+nr*(idos-1))]][t<=struc.info$tf[[l]]][1])
          }
          if(graph.info$log.logical!=F) 
         plot(tt1[[l]],ff1[[l]],type='l',xlab=paste(graph.info$names.datax[l]),ylab=paste(graph.info$names.datay[l]),
            log=graph.info$log.logical,ylim=c(graph.info$y.range[[l]],range(ff1[[l]]))[1:2]) else 
         plot(tt1[[l]],ff1[[l]],type='l',xlab=paste(graph.info$names.datax[l]),ylab=paste(graph.info$names.datay[l]),
            ylim=c(graph.info$y.range[[l]],range(ff1[[l]]))[1:2])
            
      #points for the design on the graph
          for(kj in which(icond==uni.cond[kd])) {
            tps<-gr.tps[[kj]][gr.list[[kj]]==l]
            tab1<-tab[tab$irep==l & tab$condnum==icond[kj],]
             if (info.covariate$covariate.model==T || model.info$n_occ==1)
           {      f.tps<-vector() 
                 for (ij in 1: length(tps))
                  {
                    t<-tps[ij]
                    if (model.info$n_occ==1) f.tps[ij]<-eval(struc.info$formg_init[[(l+nr*(idos-1))]][t<=struc.info$tf[[l]]][1])
                    else f.tps[ij]<-eval(struc.info$formg[[(l+nr*(idos-1))]][t<=struc.info$tf[[l]]][1])
                  }
            }
            else f.tps<-tab1$f.p[match(tps,tab1$times)]
            points(tps,f.tps,pch=paste(kj),cex=2,col=paste(kj))
          }
        
          title(paste(graph.info$names.datay[l],"model ","\n Dose = ",
           init.info$doses[uni.cond[kd]],sep=" "))
       }
     }
   }
}



#function to build graph, model in analytical form, standard model function writing
graph.fw.afn<-function(gr.tps,gr.list,tab,icond,graph.info,struc.info,model.info,
  init.info,info.covariate) {
   p<-length(model.info$beta)
   form<-struc.info$form
   nr<-model.info$nr
  
   
   for (i in 1:p) 
   assign(parameters[i],model.info$beta[i])
   theta<-model.info$beta
   
   graphinf<-graph.info$graphinf
   graphsup<-graph.info$graphsup
   
   tt<-list()
   tt<-lapply(1:nr, function(nr,graphinf,graphsup) seq(graphinf[[nr]],graphsup[[nr]],(graphsup[[nr]]-graphinf[[nr]])/1000),graphinf=graphinf,graphsup=graphsup)   
   
   ff1<-list()
   
   if ((init.info$dose.identical==T) | length(unique(icond))==1) {
     idos<-which(c(1:length(init.info$doses))==icond[1])
     dose<-init.info$doses[icond[1]] 
     par(mfrow=c(1,nr))
     
     ff1<-list()
        for (li in 1 : nr){
				 
          ff1[[li]]<-form(tt[[li]],theta,dose)
          
           if (nr>1){if (li==1) ff1[[li]]<-ff1[[li]][1:length(tt[[li]])] 
           else  ff1[[li]]<-ff1[[li]][(length(tt[[li]])+1):(length(tt[[li-1]])+length(tt[[li]]))] }
        
        }
        
        for (l in 1:nr){    
        if(graph.info$log.logical!=F) 
          plot(tt[[l]],ff1[[l]],type='l',xlab=paste(graph.info$names.datax[l]),ylab=paste(graph.info$names.datay[l]),
            log=graph.info$log.logical,ylim=c(graph.info$y.range[[l]],range(ff1[[l]]))[1:2]) else 
          plot(tt[[l]],ff1[[l]],type='l',xlab=paste(graph.info$names.datax[l]),ylab=paste(graph.info$names.datay[l]),
            ylim=c(graph.info$y.range[[l]],range(ff1[[l]]))[1:2])
            
      #points for the design on the graph

        for (i in 1:length(gr.tps)) {
           tps<-gr.tps[[i]][gr.list[[i]]==l]
           if (nr>1) tps1<-gr.tps[[i]][gr.list[[i]]==l-1]
           tab1<-tab[tab$irep==l & tab$condnum==icond[1],]
           if (info.covariate$covariate.model==T || model.info$n_occ!=1)
           {      f.tps<-vector() 
                  f.tps<-form(tps,theta,dose)
				          if (nr>1){if (l==1) f.tps<-f.tps[1:length(tps)] 
                  else  f.tps<-f.tps[(length(tps1)+1):(length(tps1)+length(tps))] 
                  }
                  }else f.tps<-tab1$f.p[match(tps,tab1$times)]
           
           points(tps,f.tps,pch=paste(i),cex=2,col=paste(i))
        }
     
        title(paste(graph.info$names.datay[l],"model ","\n Dose = ",
          dose,sep=" "))
      }
   } else {
#Regroup doses
     uni.cond<-sort(unique(icond))
     par(mfrow=c(length(uni.cond),nr))
     for(kd in 1:length(uni.cond)) {
#       dose<-doses[unique(icond)[kd]]
       idos<-which(c(1:length(init.info$doses))==uni.cond[kd])
       dose<-init.info$doses[uni.cond[kd]]
       for (l in 1:nr){

           ff1[[l]]<-form(tt[[l]],theta,dose)
          
           if (nr>1){if (l==1) ff1[[l]]<-ff1[[l]][1:length(tt[[l]])] 
           else  ff1[[l]]<-ff1[[l]][(length(tt[[l]])+1):(length(tt[[l-1]])+length(tt[[l]]))] }
       	 	
          if(graph.info$log.logical!=F) 
         plot(tt[[l]],ff1[[l]],type='l',xlab=paste(graph.info$names.datax[l]),ylab=paste(graph.info$names.datay[l]),
            log=graph.info$log.logical,ylim=c(graph.info$y.range[[l]],range(ff1[[l]]))[1:2]) else 
         plot(tt[[l]],ff1[[l]],type='l',xlab=paste(graph.info$names.datax[l]),ylab=paste(graph.info$names.datay[l]),
            ylim=c(graph.info$y.range[[l]],range(ff1[[l]]))[1:2])
            
      #points for the design on the graph

          for(kj in which(icond==uni.cond[kd])) {
            tps<-gr.tps[[kj]][gr.list[[kj]]==l]
            tab1<-tab[tab$irep==l & tab$condnum==icond[kj],]
            if (info.covariate$covariate.model==T || model.info$n_occ!=1)
           {      f.tps<-vector() 
                  f.tps<-form(tps,theta,dose)
				          if (nr>1){if (l==1) f.tps<-f.tps[1:length(tps)] 
                  else  f.tps<-f.tps[(length(tps1)+1):(length(tps1)+length(tps))] 
                  }
                  }else f.tps<-tab1$f.p[match(tps,tab1$times)]

            points(tps,f.tps,pch=paste(kj),cex=2,col=paste(kj))
          }
        
          title(paste(graph.info$names.datay[l],"model ","\n Dose = ",
           init.info$doses[uni.cond[kd]],sep=" "))
       }
     }
   }
}


#function to build graph, model described by differential equations
graph.fw.de<-function(gr.tps,gr.list,tab,icond,graph.info,struc.info,model.info,
  init.info,info.covariate) {
   p<-length(model.info$beta)
   formED<-struc.info$formED
   nr<-model.info$nr
   for (i in 1:p) assign(parameters[i],model.info$beta[i])
   graphinf<-graph.info$graphinf
   graphsup<-graph.info$graphsup
   
   tt<-list()
   ff1<-list()
   e<-list()
   e1<-list()
   tt<-lapply(1:nr, function(nr,graphinf,graphsup) seq(graphinf[[nr]],
     graphsup[[nr]],(graphsup[[nr]]-graphinf[[nr]])/1000),graphinf=graphinf,graphsup=graphsup)
   if (init.info$condinit.identical==T | length(unique(icond))==1) {
     par(mfrow=c(1,nr))
     cond<-eval(init.info$condinit[icond[1]])
     e<-lapply(1:nr,function(nr,time.condinit,tt) c(time.condinit,tt[[nr]]),
       time.condinit=init.info$time.condinit,tt=tt)
     for (li in 1:nr){
        ff1[[li]]<-lsoda(cond,e[[li]],formED,beta,rtol=RtolEQ,atol=AtolEQ,hmax=Hmax)
        if (nr==1) ff1[[li]]<-ff1[[li]][,dim(ff1[[li]])[2]] [-1]
        else  ff1[[li]]<-ff1[[li]][,((dim(ff1[[li]])[2]-nr+1):(dim(ff1[[li]])[2]))][,li] [-1]
     }
     for (l in 1:nr){
        if (graph.info$log.logical!=F) 
          plot(tt[[l]],ff1[[l]],type='l',xlab=paste(graph.info$names.datax[l]),ylab=paste(graph.info$names.datay[l]),
            log=graph.info$log.logical,ylim=c(graph.info$y.range[[l]],range(ff1[[l]]))[1:2]) else 
          plot(tt[[l]],ff1[[l]],type='l',xlab=paste(graph.info$names.datax[l]),ylab=paste(graph.info$names.datay[l]),
             ylim=c(graph.info$y.range[[l]],range(ff1[[l]]))[1:2])
      
      #points for the design on the graph

        for (i in 1:length(gr.tps)) {
           tps<-gr.tps[[i]][gr.list[[i]]==l]
           tab1<-tab[tab$irep==l & tab$condnum==icond[1],]
           if (info.covariate$covariate.model==T || model.info$n_occ==1 )
           {      f.tps<-vector() 
                 #for (ij in 1: length(tps))
                  #{
                    #t<-tps[ij]
                    #f.tps[ij]<-eval(struc.info$formg_init[[(l+nr*(idos-1))]][t<=struc.info$tf[[l]]][1])
                    f.tps1<-lsoda(cond,c(init.info$time.condinit,tps),formED,beta,rtol=RtolEQ,atol=AtolEQ,hmax=Hmax)
                    if (nr==1) f.tps<-f.tps1[,dim(f.tps1)[2]][-1] else f.tps<-f.tps1[,((dim(f.tps1)[2]-nr+1):(dim(f.tps1)[2]))][,l][-1]
                       #}
            }
            else f.tps<-tab1$f.p[match(tps,tab1$times)]
            points(tps,f.tps,pch=paste(i),cex=2,col=paste(i))
        }
       
        title(paste(graph.info$names.datay[l],"model ","\n Initial Conditions = ",
          init.info$condinit[icond[1]],sep=" "))
     }
   } else {
#Regroup doses
     par(mfrow=c(length(unique(icond)),nr))
     for(kd in 1:length(unique(icond))) {
       cond<-eval(init.info$condinit[unique(icond)[kd]])
        e<-lapply(1:nr,function(nr,time.condinit,tt) c(time.condinit,
          tt[[nr]]),time.condinit=init.info$time.condinit,tt=tt)
        ff1<-list()
        for (li in 1 : nr) {
          ff1[[li]]<-lsoda(cond,e[[li]],formED,beta,rtol=RtolEQ,atol=AtolEQ,hmax=Hmax)
          if (nr==1) ff1[[li]]<-ff1[[li]][,dim(ff1[[li]])[2]] [-1]
          else 
            ff1[[li]]<-ff1[[li]][,((dim(ff1[[li]])[2]-nr+1):(dim(ff1[[li]])[2]))][,li] [-1]
        }
        for (l in 1:nr){
          if (graph.info$log.logical!=F) 
         plot(tt[[l]],ff1[[l]],type='l',xlab=paste(graph.info$names.datax[l]),ylab=paste(graph.info$names.datay[l]),
            log=graph.info$log.logical,ylim=c(graph.info$y.range[[l]],range(ff1[[l]]))[1:2]) else 
         plot(tt[[l]],ff1[[l]],type='l',xlab=paste(graph.info$names.datax[l]),ylab=paste(graph.info$names.datay[l]),
             ylim=c(graph.info$y.range[[l]],range(ff1[[l]]))[1:2])
         
      #points for the design on the graph

          for(kj in which(icond==unique(icond)[kd])) {
            tps<-gr.tps[[kj]][gr.list[[kj]]==l]
            tab1<-tab[tab$irep==l & tab$condnum==icond[kj],]
            if (info.covariate$covariate.model==T || model.info$n_occ==1 )
             {
                  f.tps<-vector() 
                  f.tps1<-lsoda(cond,c(init.info$time.condinit,tps),formED,beta,rtol=RtolEQ,atol=AtolEQ,hmax=Hmax)
                  if (nr==1) f.tps<-f.tps1[,dim(f.tps1)[2]][-1] else f.tps<-f.tps1[,((dim(f.tps1)[2]-nr+1):(dim(f.tps1)[2]))][,l][-1]
              }else  f.tps<-tab1$f.p[match(tps,tab1$times)]
                 points(tps,f.tps,pch=paste(kj),cex=2,col=paste(kj))
            }
      
          title(paste(graph.info$names.datay[l],"model ","\n Initial Conditions = ",
            init.info$condinit[unique(icond)[kd]],sep=" "))
        }
      }
   }
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
#########################################################################

power_nni.computation<-function(alpha,given.power,init.num.subjects,beta.cov,beta.cov.se,nom.par,Trand)
{
norm.val<-abs(qnorm(alpha, mean=0, sd=1))

#Puissance avec cet effet la avec alpha =5%
power.prediction<-1-pnorm(norm.val-(beta.cov/beta.cov.se)) + pnorm(-norm.val-(beta.cov/beta.cov.se))

Beta<-beta.cov
Predicted_SE<-beta.cov.se
lcov<-length(beta.cov)  
Expected_power<-round(power.prediction,6)

#Calcul du nombre de sujet necessaires pour avoir une puissance donn?e%
#given.power
#calcul de la SE si je me met sous H1
NSE<-beta.cov/(norm.val-qnorm(1-given.power))

#calcul maintenant du nombre de sujet necessaires pour mettre en evidence cette effets
#avec une puissance donn?e
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



#Calculer la puissance et le NSN pour le test d'?quivalence
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




#Calcul du nombre de sujet necessaires pour avoir une puissance donn?e%
#given.power
#calcul de la SE si je me met sous H1
NSE<-c()
for (i in 1:length(beta.cov)){
if (beta.cov[i]==0){NSE[i]<-interval_eq[2]/(norm.val+qnorm(given.power))}
if (beta.cov[i]<0){NSE[i]<-(-beta.cov[i]+interval_eq[1])/(-norm.val+qnorm(1-given.power))}
if (beta.cov[i]>0){NSE[i]<-(-beta.cov[i]-interval_eq[1])/(norm.val+qnorm(given.power))}
}
#calcul maintenant du nombre de sujet necessaires pour mettre en evidence cette effets
#avec une puissance donn?e
NNI<-sum(init.num.subjects)*(beta.cov.se/NSE)^2



#si pas sous H1

for (i in 1:length(beta.cov)){
if (beta.cov[i] < interval_eq[1] || beta.cov[i] > interval_eq[2]) {
power.prediction[i]<-c(NA)
NNI[i]<-c(NA)
}#else
#{#if (length(which(beta.cov[i] < interval_eq[1])) + length(which(beta.cov[i]> interval_eq[2])) < length(beta.cov))
#power.prediction[i]<-round(power.prediction[i], digits = 6)
#NNI[i]<-round(NNI[i], digits = 6)} 
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


###########################################################################
###                                 Fedorov output                      ###                                                                                                                                                            
###########################################################################
out.fedorov<-function(modelfile="Stdin.r",directory=".",
directory.program=".") {  
#directories  
#directory : control file ("Stdin.r"), models and results
#directory.program : FW library, model library, and this program
##############################


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
    if (graph.logical==T) {graphAFS(form,prots,graphinf,graphsup,names.datax,names.datay,y.range)} else NULL
    if (graphsensi.logical==T) {graphsensiAFS(form,prots,graphinf,graphsup,names.datax,names.datay,y.range)} else NULL
  }
  if (modelform == "DE"){
    source(file.path(directory.program,"EQPfim4.0op1.r"))
    if (graph.logical==T)  {graph(formED,prots,graphinf,graphsup,names.datax,names.datay,y.range)} else NULL
    if (graphsensi.logical==T) {graphsensi(formED,prots,graphinf,graphsup,names.datax,names.datay,y.range)} else NULL
  }
  options("show.error.messages"=F);stop()
}



   
#precomputation
   options(expressions=5000)  
   d<-Sys.time()
   g<-proc.time()
   under.unix<-!(version$os=='Microsoft Windows' ||
      version$os=='Win32' || version$os=='mingw32')
   #cat("Directory for PFIM3.1 programs:",directory.program,"\n")
   #cat("Directory for project:",directory,"\n")
#Global variables
   assign("directory",directory,env = .GlobalEnv)
   assign("directory.program",directory.program,env = .GlobalEnv)
   source(file.path(directory,modelfile))
   modfile<-file.path(directory,file.model)

#Loading libraries and Fedorov-Wynn
   require("combinat")
   require(deSolve)
   require(nlme)
   require(numDeriv)

   if(!is.loaded("FedorovInit_R")) {
      if(under.unix) {
         cmd<-paste("R CMD SHLIB -o ",file.path(directory.program,"libFED.so"),file.path(directory.program,"initfedoR.c"))
         system(cmd)
         dyn.load(file.path(directory.program,"libFED.so"))
      } else
      {
       if ( sessionInfo()[1]$R.version$arch == "x86_64" )  {dyn.load(file.path(directory.program,"libFED64.dll"))}
       else  {dyn.load(file.path(directory.program,"libFED.dll"))}
      }
   }
  
  
#Common for models in analytical form or differential equations
   doses<-NULL
   p<-length(beta)
   if (FIM=="B") {subjects<-c(1)}
   ##############################
   
#Partie li?e ? l'impl?mentation pour iov
#-------------------------------------------------------------------------------
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

#les combinaisons des s?quences
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

#les propportions de chaque combinaison de s?quence
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

#le vecteur des parametres fixes prenant compte des effets covariables changeant entre les p?riodes

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




   #############################
   sigmainter<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("sig.inter",LETTERS[1:nr],sep=""))
   sigmaslope<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("sig.slope",LETTERS[1:nr],sep=""))
#Default values
   if(!exists("identical.times")) identical.times<-F # or T???

#recover the inital protocols for all the models
   if(identical.times==T){
     prots<-lapply(1:nr,function(nr,x) get(x),x=paste("prot",LETTERS[1],sep=""))
   } else
     prots<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("prot",LETTERS[1:nr],sep=""))
     
#if several responses, check that the same nb of protocols are specified for
# each response
   if(nr>1) {
     if(length(unique(unlist(lapply(prots,length))))>1)
       stop("You must specify the same number of initial protocols for each response.\n")
   }
   ll<-lapply(prots,length)


#For models in analytical form
   if(modelform=="AF") {
#Problem with multiple doses
            prot.dose<-dose #doses associated to each initial protocol
            #doses<-sort(unique(dose)) # nb of possible doses
            doses<-unique(dose)
            assign("dose",doses[1],env=.GlobalEnv)
            dose<-doses[1]
            source(modfile,local=T)
            
            if (NUM==F){
            lf<-length(form) #nombre de forme analytique
            
      #nombre de form dans chaque r?ponse pour savoir quoi prendre en tf
            formg<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("form",LETTERS[1:nr],sep=""))  #sous forme de list
            lformg<-lapply(formg,length)  #length(of each model tf or not)
            #lf<-sum(unlist(lformg)) #vir? car inutile ?: c(0,cumsum(lformg))
            form<-form[1:lf]
            if(length(doses)>1) {
              eco.form1<-form
              for(idos in 2:length(doses)) {
                assign("dose",doses[idos],env=.GlobalEnv)
                 dose<-doses[idos] #seems to need both???
                source(modfile,local=T)
                for(ilf in 1:lf)
                  eco.form1[[(ilf+(idos-1)*lf)]]<-form[[ilf]]
                eco.formg1<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("form",
                  LETTERS[1:nr],sep=""))  #sous forme de list
                for(inr in 1:nr)
                  formg[[(inr+nr*(idos-1))]]<-eco.formg1[[inr]]
              }
              form<-eco.form1
            }
            nform<-lf*length(doses)
            }
#tests si tous les noms des X et Y pour toutes les reponses      
#if (length(names.datax)!=nr) names.datax<-rep("NULL",nr)
#if (length(names.datay)!=nr) names.datay<-rep("NULL",nr)      
      
#tests si toutes les bornes
      pres.b<-lapply(1:nr,function(nr,x) exists(x[nr]),
        x=paste("bound",LETTERS[1:nr],sep="")) 
      TF<-which(pres.b==F) 
      if(length(TF)==0) 
        bornes<-lapply(1:nr,function(nr,x) get(x[nr]),
          x=paste("bound",LETTERS[1:nr],sep="")) else {
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
      if (NUM==F){
#recover tf in order to compute the sensibility
      tf<-vector(mode="list")
      for (i in 1:length(lformg)){
        if (lformg[i]!=1){ 
          tf[[i]]<-vector()
          tf[[i]][1]<-bornes[[i]][[1]][2]
          for (j in 2:lformg[[i]])
             tf[[i]][j]<-c(bornes[[i]][[j]][2]) 
        } else tf[[i]]<-bornes[[i]][[1]][2]
      }}
       
      if (dose.identical==T) doses<-rep(prot.dose[1],length(prots[[1]])) else NULL
      condinit<-NULL
      
      #pas possible de faire du design individuel avec FW
       b_condition<-T
      if (length(omega)!=0 && length(which(omega==0))==length(omega)) 
     {
     #stop("You can not use the Fedorov-Wynn algorithm for optimization of individual design (i.e. no random effects).\n")
      FIM<-"I"
     }
          
   } else { #for models described by differential equations
      source(modfile)
      if (condinit.identical==T) 
        condinit<-rep(condinit[1],length(prots[[1]]))
 
          
      #tests si tous les noms des X et Y pour toutes les reponses      
        if (length(names.datax)!=nr) names.datax<-rep("NULL",nr)
        if (length(names.datay)!=nr) names.datay<-rep("NULL",nr)  

   }
   if (identical.times==T){
     nwind<-lapply(1:nr,function(nr,x) get(x),x=paste("nwind",LETTERS[1],sep=""))
     sampwin<-lapply(1:nr,function(nr,x) get(x),x=paste("sampwin",LETTERS[1],sep=""))
#sampwin<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("sampwin",LETTERS[1:nr],sep=""))
     fixed.times<-lapply(1:nr,function(nr,x) get(x),x=paste("fixed.times",LETTERS[1],sep=""))
     nsamp<-lapply(1:nr,function(nr,x) get(x),x=paste("nsamp",LETTERS[1],sep=""))
     nmaxpts<-lapply(1:nr,function(nr,x) get(x),x=paste("nmaxpts",LETTERS[1],sep=""))
     nminpts<-lapply(1:nr,function(nr,x) get(x),x=paste("nminpts",LETTERS[1],sep=""))
   } else {
     nwind<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("nwind",
       LETTERS[1:nr],sep=""))
     sampwin<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("sampwin",
       LETTERS[1:nr],sep=""))
     fixed.times<-lapply(1:nr,function(nr,x) get(x),x=paste("fixed.times",LETTERS[1:nr],sep=""))
     #sampwin<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("sampwin",LETTERS[1:nr],sep=""))
     nsamp<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("nsamp",
       LETTERS[1:nr],sep=""))
     nmaxpts<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("nmaxpts",
       LETTERS[1:nr],sep=""))
     nminpts<-lapply(1:nr,function(nr,x) get(x[nr]),x=paste("nminpts",
       LETTERS[1:nr],sep=""))
   }

   subjects.init<-subjects
   if (subjects.input==1) {
     l2<-lapply(prots,function(lis) lapply(lis,length))
     Nt<-lapply(1:nr, function(nr,subjects.init,l2) sum(subjects.init*unlist(l2[[nr]])),
       subjects.init=subjects.init,l2=l2)
     Ntot<-sum(unlist(Nt))
     subjects.total<-sum(subjects.init)
     subjects.init<-subjects.init/subjects.total
   }
   subjects<-subjects.init
 
 
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



#Creating structures to pass information along to subroutines
#All info pertaining to the model parameters
 
    
   if (n_occ<=1) gamma<-diag(NULL)
   model.info<-list(beta=beta,theta1=theta1,theta2=theta2,beta.occ=beta.occ,beta.occ1=beta.occ1,beta.occ2=beta.occ2,
   n_occ=n_occ,locc=locc,lposs=lposs,lP=lP,lP1=lP1,vecP=vecP,vecP1=vecP1,
   nbcov_occ=nbcov_occ,comb.seq.poss=comb.seq.poss,comb.seq.prop=comb.seq.prop,
   omega=omega,gamma=gamma,OMEGA=diag(c(diag(omega),rep(c(diag(gamma)),n_occ))),
   nullpar=which(c(diag(omega),diag(gamma))==0),
      sigmainter=sigmainter,sigmaslope=sigmaslope,nr=nr,Trand=Trand,
      parameters=parameters,names.datax=names.datax,names.datay=names.datay,beta.fixed=beta.fixed) 
   if(modelform=="AF"){
     if(NUM==F) struc.info<-list(modfile=modfile,modelform=modelform,form=form,lf=lf,
        formg=formg,lformg=lformg,tf=tf,bound=bornes,NUM=NUM) 
     else struc.info<-list(modelform=modelform,form=form,NUM=NUM)   
       } else
     struc.info<-list(modelform=modelform,formED=formED)
    
#All info pertaining to sampling
   sampl.info<-list(identical.times=identical.times,nwind=nwind,
     sampwin=sampwin,nsamp=nsamp,nmaxpts=nmaxpts,nminpts=nminpts,fixed.times=fixed.times)
#All info needed in graphs
   graph.info<-list(graph.logical=graph.logical,log.logical=log.logical,graphsensi.logical=graphsensi.logical,
      graphinf=graphinf,graphsup=graphsup,names.datax=names.datax,names.datay=names.datay,y.range=y.range)     
#All info pertaining to initial conditions, forms, etc...
   if(modelform=="AF")
     init.info<-list(dose.identical=dose.identical,doses=doses,
        prot.dose=prot.dose,prots=prots,subjects=subjects,cost=Ntot) else
     init.info<-list(condinit.identical=condinit.identical,condinit=condinit,
       prots=prots,subjects=subjects,cost=Ntot,
       time.condinit=time.condinit,RtolEQ=RtolEQ,AtolEQ=AtolEQ,Hmax=Hmax)
    info.covariate<-list(covariate.model=covariate.model)
    
#All info pertaining to covaraite model...
    if(covariate.model==T) {
                                
      #error messages
      if(option==2) stop(" You can not use option 2 with covariate currently")
      if(FIM!="P") stop("Individual and bayesian matrices are only computed without covariate")
          
       covariate.cat.reco<-list()
        info.covariate<-list(covariate.name=covariate.name,covariate.category=covariate.category,covariate.proportions=covariate.proportions,
        parameter.associated= parameter.associated,beta.covariate=beta.covariate,covariate.model=covariate.model)
        info.power.nni<-list(alpha=alpha,compute.power=compute.power,compute.nni=compute.nni,
        compute.power_eq=compute.power_eq,compute.nni_eq=compute.nni_eq,given.power=given.power,interval_eq=interval_eq)
     
          if ((length(covariate.name)!= length(covariate.category))==T || (length(covariate.category)!=length(covariate.proportions))==T) 
          stop(" The list of the name of the covariate(s) has not the same length of the list of the category or the proportions")
  
          if(length(which((sapply(covariate.category,length) == sapply(covariate.proportions,length))==F))!=0 )
          stop(" The number of categories is not equal to the number of proportions")
          
          if (length(covariate.name)!= length(parameter.associated))
          stop(" One or more covariate(s) has not at least one associated parameter ")
          
          if(sum(unlist(sapply(covariate.proportions,sum)))!=length(covariate.name))
          stop(" The sum of the proportions of the categories is not equal to one")
          for (i in 1: length(info.covariate[[1]]))
            {
                if(is.character(info.covariate[[2]][[i]]))
                   covariate.cat.reco[[i]]<-c(1:length(info.covariate[[2]][[i]]))
                else 
                   covariate.cat.reco[[i]]<-info.covariate[[2]][[i]]
            }
            
          for (i in 1: length(unlist(info.covariate$parameter.associated)))
          {
            if(length(which(unlist(info.covariate$parameter.associated)[i]==parameters))==0) 
            stop(" One of the parameter associated to one covariate does not correspond with the parameters of the model")
          
          }
            
          #solution possibles pour les covariables
          covariate.cat.poss<-unique(t(combn(unlist(covariate.cat.reco),length(info.covariate[[1]]))))
          #verifier que les modalit?s presentes dans chque colonne sont bien des modalit?s possibles de la covariables
           #test si les modalit?s possibles pour cahque covaraibles ne sont pas abb?rantes
          ll.cov<-length(info.covariate[[1]])
          for (i in 1: ll.cov)
          {
            nline<-which(covariate.cat.poss[,i] > max(covariate.cat.reco[[i]]))
            ref<-integer(length = 0)
        
            if (!identical(nline,ref)) covariate.cat.poss<-covariate.cat.poss[-nline,]
          }
          
          for (i in 1: ll.cov)
          {
            #nline<-which(covariate.cat.poss[,i] < min(covariate.cat.reco[[i]]))
                        nline<-which(covariate.cat.poss[,i] > max(covariate.cat.reco[[i]]) | covariate.cat.poss[,i] < min(covariate.cat.reco[[i]]))
            ref<-integer(length = 0)
        
            if (!identical(nline,ref)) covariate.cat.poss<-covariate.cat.poss[-nline,]
          }
          
           #creation de la formule du modele pour tenir compte des parametres des covariables
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
              
               #liste des modalit?s
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
            #tableau avec la liste des parametres sur lequels il y a une covariables ainsique la covariable associ?s et le nom du parametre de
           #la covariables
            ##tab.cov.param1<-tab.cov.param[order(tab.cov.param[,1]),]
            ##tab.cov.param<-matrix(tab.cov.param1,nrow=dim(tab.cov.param)[1],ncol=dim(tab.cov.param)[2])
            
            tab.cov.param<-matrix(tab.cov.param,nrow=dim(tab.cov.param)[1],ncol=dim(tab.cov.param)[2])
            
            p.parameters<-c(parameters,npar)
            pfixe<-length(c(parameters,npar))
            tab_mod_prop<-list()
            tab_mod_prop<-association_modality_covariate(info.covariate[[1]],covariate.cat.reco,info.covariate[[3]],covariate.cat.poss)
          
          
           #####
            if( modelform=="AF" && NUM==F)
            {
              formg_init<-formg
             # if (Trand==1)formg<-transf_model_trand1(nr,formg,lformg,parameters,tab.cov.param)  else  formg<-transf_model_trand2(nr,formg,lformg,parameters,tab.cov.param) 
                
              #summary info covariate model
              #association des prportions aux diff?rentes combianaisons
                
              struc.info$formg<-formg
              #lformg<-lapply(formg,length)  #length(of each model tf or not)
              #lf<-sum(unlist(lformg)) #vir? car inutile ?: c(0,cumsum(lformg))
              #form<-unlist(formg)[1:lf]
              #struc.info$form<-form   
              struc.info<-list(modfile=modfile,modelform=modelform,form=form,lf=lf,formg=formg,lformg=lformg,tf=tf,bound=bornes,formg_init=formg_init)
              } 
          #####
            
            info.covariate<-list(covariate.name=covariate.name,covariate.category=covariate.category,covariate.proportions=covariate.proportions,
                parameter.associated=parameter.associated,beta.covariate=beta.covariate,covariate.model=covariate.model,tab.cov.param=tab.cov.param,tab_mod_prop=tab_mod_prop,covariate.cat.poss=covariate.cat.poss,p.parameters=p.parameters)
           
             power.nni.info<-list(alpha=alpha,compute.power=compute.power,compute.nni=compute.nni,compute.power_eq=compute.power_eq,compute.nni_eq=compute.nni_eq,given.power=given.power,interval_eq=interval_eq)  
        }
         
#All info pertaining to model of covariates which depend on the subject and the occasion
  if (pn==1 || locc==0) no<-NULL
      if (n_occ>1 && locc>0){
      no<-NULL
      no1<-NULL
      lo<-NULL
  if (covariate_occ.model==T){
         #error messages
      if(option==2) stop(" You can not use option 2 with covariate currently")
      if(FIM!="P") stop("Individual and bayesian matrices are only computed without covariates")
      
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
      for (kj in 1:length(covariate_occ.name))
      {
          for (j in 1:length(parameter_occ.associated[[kj]]))
          {
            
              lo<-length(covariate_occ.category[[kj]])-1
              no1<-lapply(1:lo,function(lo) paste("beta_",parameter_occ.associated[[kj]][j],"_",covariate_occ.name[[kj]],"_",covariate_occ.category[[kj]][lo+1],sep=""))
              no<-c(no,unlist(no1)) #liste des parametres
            }
             count_occ<-count_occ+lo
          }
        power.nni.info<-list(alpha=alpha,compute.power=compute.power,compute.nni=compute.nni,compute.power_eq=compute.power_eq,compute.nni_eq=compute.nni_eq,given.power=given.power,interval_eq=interval_eq)  
       
      }
   }
   
 
 
       
##############################
#Running Fedorov-Wynn
   fed<-fedorov.wynn.sg(F,fileres=file.path(directory,"details.r"),
       filematelem=file.path(directory,"matelem.tmp"),
       model.info,struc.info,sampl.info,init.info,option,info.covariate)
       #if not 'P'FIM : return(list(protoc.fed=protoc.fed,init.fed=init.fed,data.fed=data.fed,
   #sampl.info=sampl.info,model.info=model.info,bestprot1=bestprot1,bestcrit1=bestcrit1,dose1=dose1,nsuj=nsuj))
#retrieve modified data structures
   model.info<-fed$model.info;
   sampl.info<-fed$sampl.info
#evaluation du protocle initial
   init.fed<-fed$init.fed
   initnp<-init.fed$init.np
   data.fed<-fed$data.fed
   tab<-data.fed$tab
   if (FIM=="P") {opti.fed<-fed$opti.fed}   
   protoc.fed<-fed$protoc.fed
   protoc<-protoc.fed$protoc;
   tprotoc<-protoc.fed$tprotoc
   iprote<-lapply(init.fed$init.num,function(i) protoc[[i]])
   lprote<-lapply(iprote,length)
   tprote<-lapply(init.fed$init.num,function(i) tprotoc[[i]])
#modif best 1 group protocol
bestprot1<-fed$bestprot1
bestcrit1<-fed$bestcrit1
dose1<-fed$dose1
nsuj<-fed$nsuj 
idbest<-fed$idbest
numbest<-fed$numbest


#changer pour refl?ter les fr?quences initiales (?ventuellement ajust?es en cas de pb)
   Nbsubjects.init<-init.info$cost/sum((unlist(lprote)*init.fed$init.freq))
   subjects.init<-init.fed$init.freq*Nbsubjects.init
   
   if (model.info$n_occ==1)  {call.fish<-ecofisher.meth1.final(iprote,data.fed$tab,data.fed$f.p,
   data.fed$df.p,init.fed$init.freq,model.info,
   init.info$cost,option,data.fed$ddf.p,info.covariate)}
   if (model.info$n_occ>1)  {call.fish<-ecofisher.meth1.final.iov(iprote,data.fed$tab,data.fed$sensi,
   data.fed$sensia,init.fed$init.freq,model.info,
   init.info$cost,option,info.covariate)}
   crit.init<-call.fish$crit
   #print("passe")
   #print(call.fish)
   
   SIG<-call.fish$SIG  
   lSIG<-length(SIG[SIG!=0])

#Eco : now saves to same directory as 'control' file
   sink(file.path(directory,output)) 
   
   cat("PFIM 4.0 ",'\n','\n') 
   cat("Option:",option,'\n')
   cat("\n","\n") 
   cat("Project:",project)
   cat("\n","\n") 
   cat('Date:',date()) 
   cat("\n","\n") 
   cat("\n","\n") 
 
   e<-list()
   if(modelform=="DE") {     
     for (l in 1:nr) {   
       Protocol<-c()
       for(iprot in 1:init.fed$init.np) {
         tab1<-tab[protoc[[init.fed$init.num[iprot]]],]
         Protocol<-c(Protocol,
           paste("c=(",paste(tab1$times[tab1$irep==l],collapse=", "),")",sep=""))
       }
        e[[l]]<-data.frame(Protocol,subjects=subjects.init,
          condinit=as.character(init.info$condinit))
     }
  
     cat("**************************** INPUT SUMMARY ********************************") 
     cat("\n","\n") 
     cat('Differential Equations form of the model: ','\n','\n') 
     print(struc.info$formED)
     cat("\n","\n")
    
      cat('Initial Conditions at time',init.info$time.condinit,':','\n','\n')
      ff<-sapply(1:initnp, function(initnp, form) cat(paste(form[[initnp]][-1]),'\n'), form=condinit)
      cat("\n")
      cat("Error tolerance for solving differential equations system: RtolEQ =", RtolEQ,", AtolEQ =", AtolEQ,", Hmax = ",Hmax)
      cat("\n","\n")
   } else { #modelform==AF
     for (l in 1:nr) {
       Protocol<-c()
       for(iprot in 1:init.fed$init.np) {
         tab1<-tab[protoc[[init.fed$init.num[iprot]]],]
         Protocol<-c(Protocol,
           paste("c=(",paste(tab1$times[tab1$irep==l],collapse=", "),")",sep=""))
       }
        if(length(prot.dose)<length(subjects.init)) 
          prot.dose<-rep(prot.dose,length(subjects.init))
       e[[l]]<-data.frame(Protocol,subjects=subjects.init,
          doses=prot.dose[1:length(subjects.init)])
     }
     cat("**************************** INPUT SUMMARY ********************************") 
     cat("\n","\n") 
     cat('Analytical function model: ','\n','\n') 
     if(NUM==T){print(attributes(form)$srcref);cat("\n","\n")}
     else{ff<-sapply(1:lf, function(lf, form) cat(paste(form[lf]),'\n','\n'), form=form)}
  
     #ff<-sapply(1:lf, function(lf, form) cat(paste(form[lf]),'\n','\n'), form=form)
     cat("\n","\n")
   }

cat('Initial design:','\n')   
     cat("\n","\n") 
     for (i in 1:nr){
        cat('Sample times for response:',LETTERS[i],'\n')
        print(e[[i]])
        cat("\n","\n")
     }
       
   cat('Total number of samples:',init.info$cost)
   cat("\n","\n") 
   cat('Associated criterion value:',round(crit.init,4))
   cat("\n","\n")
   cat('Identical sampling times for each response:',identical.times)
   cat("\n","\n")
   
   
if (n_occ>1 && length(which(gamma==0))!=length(gamma)) 
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
          cat("\t")
          cat('Covariate',i,':',Name[1],'(',Parameters_associated[1],')',"\n")
                print(e.cov[[i]])
          cat("\n")
          }
          
  }
   
  if (covariate_occ.model==T)
	{  
     
     cat("\n","\n")
     cat("\t") 
     cat("Covariates changing with occasion","\n","\n")
     e.cov1<-list()
     e.cov2<-list()
     
          for (kj in 1: length(covariate_occ.name)){
          cat("\t")
          cat ("Covariate ",kj,":",covariate_occ.name[[kj]],"(",parameter_occ.associated[[kj]],")","\n")
            
           lcov1<-length(covariate_occ.category[[kj]])
           prem_col1<-lapply(1:lcov1,function(lcov1) paste("(",c(lcov1),")",sep=""))
           Categories<-covariate_occ.category[[kj]]
           References<-c("*",c(rep(" ",length(covariate_occ.category[[kj]])-1)))
           
           e.cov1[[kj]]<-data.frame(Categories,References,row.names=prem_col1,check.names=F,check.rows=F)
           print(e.cov1[[kj]])
           cat("\n","\n")
        
          lcov2<-length(sequence.proportions[[kj]])
          prem_col2<-lapply(1:lcov2,function(lcov2) paste("(",c(lcov2),")",sep=""))
          Sequences<- list()

           for (sequ in 1:nb_seq[kj]){
            Sequences[sequ]<-covariate_occ.sequence[[kj]][[sequ]][1]
            for (occ in 2:n_occ){
              Sequences[sequ]<-paste(Sequences[sequ],covariate_occ.sequence[[kj]][[sequ]][occ])
            }
           }

          #for (sequ in 1:nb_seq){
          #  Sequences[sequ]<-unlist(lapply(2:n_occ-1, function(n_occ) paste(covariate_occ.sequence[[k]][[sequ]][n_occ],covariate_occ.sequence[[k]][[sequ]][n_occ+1]))) 
          #}
           Proportions<-sequence.proportions[[kj]]
          e.cov2[[kj]]<-data.frame(unlist(Sequences),Proportions,row.names=prem_col2,check.names=F,check.rows=F)
          #. <- covariate_occ.sequence[[k]]
         
          #e.cov2[[k]]<-cbind(e.cov2[[k]],.,)
          names(e.cov2[[kj]])<-c("Sequences","Proportions")
          
          print(e.cov2[[kj]])
           cat("\n","\n")
          }
    }
  
  
     
   cat("\n") 
   cat("Optimization step: ","\n","\n") 
   for ( l in 1: nr) {
     cat("Sampling windows for the response:",LETTERS[l],"\n")
     for (i in 1:sampl.info$nwind[[l]]) {
        cat("Window",i,": t=",sampl.info$sampwin[[l]][[i]],"\n")
        cat("    Nb of sampling points to be taken in this window, n[",i,"]=",
          sampl.info$nsamp[[l]][[i]],"\n") }
     if(is.na(sampl.info$nmaxpts[[l]])) {
       sampl.info$nmaxpts[[l]]<-sum(unlist(lapply(sampl.info$nsamp[[l]],max)))
       cat("Computing the maximum number of points per protocol.\n")}
     cat("Maximum total number of points in one elementary protocol :",
        sampl.info$nmaxpts[[l]],"\n")
     if(is.na(sampl.info$nminpts[[l]])) {
       sampl.info$nminpts[[l]]<-sum(unlist(lapply(sampl.info$nsamp[[l]],min)))
       cat("Computing the minimum number of points per protocol.\n")}
     cat("Minimum total number of points in one elementary protocol :",
       sampl.info$nminpts[[l]],"\n")
     cat("\n","\n") 
    }
    
   
  #Now do all this for the final protocol
  #Cyrielle
  if (FIM=="P" & previous.FIM==""){

   iprote.f<-lapply(opti.fed$opti.num,function(i) protoc[[i]])
   lprote.f<-lapply(iprote.f,length)
   tprote.f<-lapply(opti.fed$opti.num,function(i) tprotoc[[i]])
   subjects.opti<-opti.fed$opti.freq
   num.opti<-opti.fed$opti.num

#Recover initial condition (dose or condinit) for optimised protocols
   opti.cond<-unlist(lapply(opti.fed$opti.num,
      function(inum,protoc,tab) tab$condnum[protoc[[inum]][1]],protoc=protoc,
      tab=tab))
   opti.fed$cond<-opti.cond
   subjects1<-subjects.opti*init.info$cost/sum((unlist(lprote.f)*subjects.opti))
   if(modelform=="AF") doses.opti<-prot.dose[opti.cond] else 
      condinit.opti<-condinit[opti.cond]
   if (model.info$n_occ==1) {fse<-ecofisher.meth1.final(iprote.f,data.fed$tab,data.fed$f.p,
   data.fed$df.p,opti.fed$opti.freq,model.info,init.info$cost,option,data.fed$ddf.p,info.covariate) }
   if (model.info$n_occ>1) {fse<-ecofisher.meth1.final.iov(iprote.f,data.fed$tab,data.fed$sensi,
   data.fed$sensia,opti.fed$opti.freq,model.info,init.info$cost,option,info.covariate)}
   
   
   cat("\n")
   cat("Now evaluating the Fisher Information Matrix for the",fed$yfed$y[[2]][1],"protocols generated \n")

   }
   else{
   fse<-ecofisher.meth1.final(idbest,data.fed$tab,data.fed$f.p,
   data.fed$df.p,nsuj,model.info,init.info$cost,option,data.fed$ddf.p,info.covariate)
   opti.cond<-unlist(lapply(numbest,
      function(inum,protoc,tab) tab$condnum[protoc[[inum]][1]],protoc=protoc,
      tab=tab))
   subjects1<-nsuj 
   num.opti<-numbest
   subjects.opti<-nsuj
   if(modelform=="AF") doses.opti<-dose else condinit.opti<-as.character(condinit) 
   }
   
   
      
    cat("\n","\n")
    cat("BEST ONE GROUP PROTOCOL:","\n","\n") 
    #best1<-data.frame(times=as.matrix(bestprot1),freq=as.matrix(rep(1,nr)),Subjects=as.matrix(rep(nsuj,nr)),doses=as.matrix(dose1)) 
    for (i in 1:nr){  #col.names=c(times,freq,Subjects,doses)
      cat('Sample times for response:',LETTERS[i],'\n')
      if (modelform=="AF") best1<-data.frame(times=as.character(t(as.matrix(bestprot1)[i,])),freq=1,Subjects=nsuj,doses=dose1[[i]]) 
      else    best1<-data.frame(times=as.character(t(as.matrix(bestprot1)[i,])),freq=1,Subjects=nsuj,condinit=dose1[[i]]) 
      print(best1)
    }
    cat("\n","\n") 
    #print(typeof(bestcrit1))
    #print(bestcrit1)
    cat('Associated criterion:', round(bestcrit1,4))  
    cat("\n","\n") 
    
    
   se<-fse$se;cv<-fse$cv;sh<-fse$sh;mfisher<-fse$somme 
   determinant<-fse$det;crit<-fse$crit;pp<-fse$pp
   
   	if (FIM=="B") {
      beta.estim<-!model.info$beta.fixed[diag(model.info$omega)!=0] 
      #theta<-theta[diag(model.info$omega)!=0][beta.estim]
      beta<-model.info$beta[diag(model.info$omega)!=0][beta.estim]
      parameters.estim<-model.info$parameters[diag(model.info$omega)!=0][beta.estim]
  
  }
  else{
      beta.estim<-!model.info$beta.fixed 
      #theta<-theta[beta.estim]
      beta<-model.info$beta[beta.estim]
      parameters.estim<-model.info$parameters[beta.estim]
  }
	

   
  if (covariate.model==T) {
    	Beta<-c(beta,unlist(info.covariate$beta.covariate),model.info$beta.occ1) 
  	   p1<-length(Beta)
  	   P<-p1
       vec<-diag(omega)!=0 
      omega1<-as.matrix(omega[,vec])[vec,]
      vec2<-SIG!=0
      StdError<-se[1:p1] 
      RSE<-cv[1:p1] 
     	a.row<-data.frame(Beta,StdError,RSE,row.names=c(parameters.estim,npar,no)) 
  
      .<-c(rep("%",p1)) 	
	   if (FIM=="B"){a.row<-cbind(a.row,.,Shrinkage=sh,.); names(a.row)[c(dim(a.row)[2]-2,dim(a.row)[2])]<-" "}  
      else {a.row<-cbind(a.row,.)
      names(a.row)[dim(a.row)[2]]<-" "} 
      
      b_condition<-F
	   if (length(which(omega==0))!=length(omega)) {
      Omega<-diag(as.matrix(omega1)) 
      StdError<-se[(p1+1):(p1+length(Omega))] 
       RSE<-cv[(p1+1):(p1+length(Omega))]  
       b.row<-data.frame(omega2=Omega,StdError,RSE,row.names=parameters[diag(omega)!=0]) 
       .<-c(rep("%",(length(Omega)))) 
       b.row<-cbind(b.row,.)
       names(b.row)[dim(b.row)[2]]<-" " 
       b_condition<-T
	     }
       
       if (n_occ>1 && length(which(gamma==0))!=length(gamma)){
       Gamma<-diag(as.matrix(gamma1)) 
       StdError<-se[(p1+length(Omega)+1):(p1+length(Omega)+length(Gamma))] 
       RSE<-cv[(p1+length(Omega)+1):(p1+length(Omega)+length(Gamma))] 
       b2.row<-data.frame(gamma2=Gamma,StdError,RSE,row.names=parameters[diag(gamma)!=0]) 
	     .<-c(rep("%",(length(Gamma)))) 
	     b2.row<-cbind(b2.row,.)
       names(b2.row)[dim(b2.row)[2]]<-" " 
	     }
       
       if (FIM!="B"){
       StdError<-se[(pp-(lSIG-1)):pp]
       RSE<-cv[(pp-(lSIG-1)):pp]
       .<-rep("%",lSIG) 
       sig.names<-unlist(c(fse$xn))[SIG!=0] 
       c.row<-cbind(data.frame(SIG=SIG[SIG!=0],StdError,RSE,row.names=sig.names),.)
       names(c.row)[c(1,dim(c.row)[2])]<-c("Sigma"," ") 
       }
   }else{
   Beta<-c(beta,model.info$beta.occ1) 
   p2<-length(Beta)
   P<-p2
   vec<-diag(omega)!=0 
   omega1<-as.matrix(omega[,vec])[vec,]
   
   vec2<-SIG!=0
   StdError<-se[1:p2] 
   RSE<-cv[1:p2] 
   a.row<-data.frame(Beta,StdError,RSE,row.names=c(parameters.estim,no)) 
    .<-c(rep("%",p2)) 	
	   if (FIM=="B"){a.row<-cbind(a.row,.,Shrinkage=sh,.); names(a.row)[c(dim(a.row)[2]-2,dim(a.row)[2])]<-" "}  
      else {a.row<-cbind(a.row,.)
      names(a.row)[dim(a.row)[2]]<-" "} 

    b_condition<-F
    if (length(which(omega==0))!=length(omega)) {
    Omega<-diag(as.matrix(omega1)) 
      StdError<-se[(p2+1):(p2+length(Omega))] 
       RSE<-cv[(p2+1):(p2+length(Omega))]  
       b.row<-data.frame(omega2=Omega,StdError,RSE,row.names=parameters[diag(omega)!=0]) 
       .<-c(rep("%",(length(Omega)))) 
       b.row<-cbind(b.row,.)
       names(b.row)[dim(b.row)[2]]<-" "
     b_condition<-T
	   }  
       if (n_occ>1 && length(which(gamma==0))!=length(gamma)){
       Gamma<-diag(as.matrix(gamma1)) 
       StdError<-se[(p2+length(Omega)+1):(p2+length(Omega)+length(Gamma))] 
       RSE<-cv[(p2+length(Omega)+1):(p2+length(Omega)+length(Gamma))] 
       b2.row<-data.frame(gamma2=Gamma,StdError,RSE,row.names=parameters[diag(gamma)!=0]) 
	     .<-c(rep("%",(length(Gamma)))) 
       b2.row<-cbind(b2.row,.) 
       names(b2.row)[dim(b2.row)[2]]<-" "
	     }     
   if (FIM!="B"){
   StdError<-se[(pp-(lSIG-1)):pp]
   RSE<-cv[(pp-(lSIG-1)):pp]
   .<-rep("%",lSIG) 
   sig.names<-unlist(c(fse$xn))[SIG!=0] 
   c.row<-cbind(data.frame(SIG=SIG[SIG!=0],StdError,RSE,row.names=sig.names),.)
   names(c.row)[c(1,dim(c.row)[2])]<-c("Sigma"," ")
   }
   }
   
  #correlation matrix
	corr.mat<-cov2cor(solve(mfisher))
	
	#eigenvalues
	#for fixed effects
	EigenValues<-eigen(mfisher)$values
	EVmin1<-min(EigenValues[1:P] )
	EVmax1<-max(EigenValues[1:P] )
	EVratio1<-EVmax1/EVmin1
	
	#for variance components
  if (FIM!="B"){
	EVmin2<-min(EigenValues[P:pp])
	EVmax2<-max(EigenValues[P:pp])
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
   
   Subjects<-subjects1 
   f<-list()
   if(modelform=="DE")  {
#eco : v?rifier si ?a marche aussi pour DE et si oui, enlever le test!
      for (l in 1 : nr) {
           vec<-lapply(num.opti,function(i,protoc,tab,irep) {
              tab1<-tab[protoc[[i]],]
              tab1$times[tab1$irep==irep]},protoc=protoc,tab=tab,irep=l)
           f[[l]]<-data.frame(times=as.character(vec),freq=subjects.opti,Subjects,
             condinit=as.character(condinit.opti))
      }
     } else  {
        for (l in 1:nr) {
           vec<-lapply(num.opti,function(i,protoc,tab,irep) {
              tab1<-tab[protoc[[i]],]
              tab1$times[tab1$irep==irep]},protoc=protoc,tab=tab,irep=l)
           f[[l]]<-data.frame(times=as.character(vec),freq=subjects.opti,
             Subjects,doses=doses.opti)
        }
    }
    #Cyrielle
    if (FIM=="P" & previous.FIM=="" ){
    cat("\n","\n")
    cat("**************************** OPTIMISED DESIGN *****************************") 
    cat("\n","\n") 
    cat("\n","\n") 
    cat('Optimised design:','\n') 
    for (i in 1:nr){
      cat('Sample times for response:',LETTERS[i],'\n')
      print(f[[i]])
      cat("\n","\n")
    }
    cat("\n","\n") 
    cat('Associated optimised criterion:', round(crit,4))  
    cat("\n","\n") 
    }
    cat("\n","\n") 
    
  
# Ecriture de la MF dans un fichier :
if (outputFIM !=  "")
  {
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


#Cyrielle
  ###################

if (previous.FIM!="")
{
cat('Previous FIM from file',previous.FIM)
cat("\n","\n") }

##################

    
 if (outputFIM!=  "")

{
cat('FIM saved in',outputFIM)

cat("\n","\n","\n") }

  if (FIM=="I") cat("Computation of the Individual Fisher information matrix") 
  if (FIM=="B") cat("Computation of the Bayesian Fisher information matrix") 
  if (FIM=="P") cat("Computation of the Population Fisher information matrix: option = ",option) 
	
	cat("\n","\n") 

  
  cat("******************* FISHER INFORMATION MATRIX ******************") 
   cat("\n","\n") 
    print(unclass(mfisher)) 
     
    
    
    cat("\n","\n") 
    cat("\n","\n") 
    cat("************************** EXPECTED STANDARD ERRORS ************************") 
    cat("\n","\n") 
    cat("------------------------ Fixed Effects Parameters -------------------------") 
    cat("\n","\n") 
    print(a.row) 
    cat("\n","\n") 
    cat("\n","\n") 
  	if (FIM=="P"){
    cat("------------------------- Variance of Inter-Subject Random Effects ----------------------") 
	  cat("\n","\n") 
	  print(b.row) 
	  cat("\n","\n") 
    }
  
  if (n_occ>1 && length(which(gamma==0))!=length(gamma)){
  cat("------------------------- Variance of Inter-Occasion Random Effects ----------------------") 
	cat("\n","\n") 
	print(b2.row) 
	cat("\n","\n")
	}  
  if (FIM!="B"){
  cat("------------------------ Standard deviation of residual error ------------------------ ")
  cat("\n","\n")
	print(c.row) 
	cat("\n","\n")} 
    cat("******************************* DETERMINANT ********************************") 
    cat("\n","\n") 
    cat(determinant) 
    cat("\n","\n") 
    cat("******************************** CRITERION *********************************") 
    cat("\n","\n") 
    cat(crit) 
    cat("\n","\n") 
    cat("\n","\n") 
    
    
   
   if (FIM=="P" & previous.FIM=="") {
   opti.fed$nb.subj<-subjects1
   if(info.covariate$covariate.model==F && model.info$n_occ==1) tab<-cbind(tab,f.p=c(fed$data.fed$f.p))
   gr.list<-lapply(opti.fed$opti.num,function(j,tab,protoc) tab$irep[protoc[[j]]],tab=tab,protoc=protoc)
   gr.tps<-lapply(opti.fed$opti.num,function(j,protoc) protoc[[j]],protoc=tprotoc)
   } 
   else {
   tab<-cbind(tab,f.p=c(fed$data.fed$f.p))
   gr.list<-lapply(num.opti,function(j,tab,protoc) tab$irep[protoc[[j]]],tab=tab,protoc=protoc)
   gr.tps<-lapply(num.opti,function(j,protoc) protoc[[j]],protoc=tprotoc)
   }
    

  summary.ci<-NULL 
  summary.ci_eq<-NULL
	summary.power_nni<-NULL
	summary.power_nni_eq<-NULL
	 beta.covariate<-info.covariate$beta.covariate;beta.occ1<-model.info$beta.occ1;locc<-model.info$locc
if (covariate.model==T || covariate_occ.model==T) 
{compute.power<-power.nni.info$compute.power; compute.power_eq<-power.nni.info$compute.power_eq;
		compute.nni<-power.nni.info$compute.nni;	compute.nni_eq<-power.nni.info$compute.nni_eq;
    alpha<-power.nni.info$alpha;interval_eq<-power.nni.info$interval_eq;}
	if (covariate.model==T || covariate_occ.model==T)
	{   			if (covariate.model==F)
	   {    
  
      beta.covariate<-NULL
      npar<-NULL
      }
      
    	if(power.nni.info$compute.power==T || power.nni.info$compute.nni==T)
    	{
    	cat("***************************** COMPARISON TEST ********************************") 
    	cat("\n","\n")
      p1<-length(Beta)
      p<-length(model.info$beta) 
    	summary.ci<-CI.computation(0.05/2,c(unlist(beta.covariate),beta.occ1),se[(p+1):(p1)],c(npar,no),model.info$Trand)
    	summary.power_nni<-power_nni.computation(alpha/2,given.power,Subjects,c(unlist(beta.covariate),beta.occ1),se[(p+1):(p1)],c(npar,no),model.info$Trand)
    	print(summary.ci)
    	cat("\n","\n")
   	  cat("Type I error =",alpha)
      cat("\n","\n")
      print(summary.power_nni)                                         
    	}
    	cat("\n","\n") 
    

    	cat("\n","\n") 
    	if(power.nni.info$compute.power_eq==T || power.nni.info$compute.nni_eq==T)
    	{
    	cat("***************************** EQUIVALENCE TEST ********************************") 
    	cat("\n","\n") 
    	p1<-length(Beta)
      p<-length(model.info$beta)
      summary.ci_eq<-CI.computation(0.05,c(unlist(beta.covariate),beta.occ1),se[(p+1):(p1)],c(npar,no),model.info$Trand)
    	summary.power_nni_eq<-power_nni_eq.computation(alpha,given.power,Subjects,c(unlist(beta.covariate),beta.occ1),se[(p+1):(p1)],c(npar,no),model.info$Trand)
    	print(summary.ci_eq)
    	cat("\n","\n")
   	  cat("Type I error =",alpha)
   	  cat("\n")
 	  	cat("Equivalence interval = [log(",paste(exp(power.nni.info$interval_eq[1])),"),log(",paste(exp(power.nni.info$interval_eq[2])),")]",sep="")
      cat("\n","\n")
      print(summary.power_nni_eq)
      if (length(which(c(unlist(beta.covariate),beta.occ1) < power.nni.info$interval_eq[1])) + length(which(c(unlist(beta.covariate),beta.occ1)> power.nni.info$interval_eq[2])) > 0)
      cat("\n","NA: Unavailable because beta is not in the equivalence interval (not under H1)")                                          
    	}
    
    
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
    
    cat("\n","\n")
    d1<-Sys.time()
    print(d1-d) 
    g1<-proc.time()
    print(g1[2]-g[2])    
    


#graphs after optimisation
###########################

	if(graph.logical==T) {
    if (modelform == "AF" && NUM==T) {graph.fw.afn(gr.tps,gr.list,tab,opti.cond,graph.info,struc.info,
         model.info,init.info,info.covariate)}
    if (modelform == "AF" && NUM==F) {graph.fw.afs(gr.tps,gr.list,tab,opti.cond,graph.info,struc.info,
         model.info,init.info,info.covariate)}
    if (modelform == "DE" ) {graph.fw.de(gr.tps,gr.list,tab,opti.cond,graph.info,struc.info,
         model.info,init.info,info.covariate)}
  } else NULL 
	if(graphsensi.logical==T) {
    if (modelform == "AF" && NUM==T){source(file.path(directory.program,"Pfim4.0op1.r"));graphsensiAFN(form,prots,graphinf,graphsup,names.datax,names.datay,y.range)} 
    if (modelform == "AF" && NUM==F){source(file.path(directory.program,"Pfim4.0op1.r"));graphsensiAFS(formg_init,prots,graphinf,graphsup,names.datax,names.datay,y.range)}
    if (modelform == "DE"){source(file.path(directory.program,"EQPfim4.0op1.r")); graphsensi(formED,prots,graphinf,graphsup,names.datax,names.datay,y.range)}
  } else NULL 

  
   cat("\n","\n")
  
   sink()  
   if (covariate.model==T)  
   return(list(init.fed=fed$init.fed,data.fed=fed$data.fed,opti.fed=opti.fed,
   mfisher=mfisher,determinant=determinant,crit=crit,se=se,cv=cv,summary.ci=summary.ci,summary.power_nni=summary.power_nni,summary.ci_eq=summary.ci_eq,summary.power_nni_eq=summary.power_nni_eq,EigenValues=EigenValues))
   
   if (covariate.model==F)
   if (FIM=="P") {cat('\n','\n','Optimisation of the population design...','\n','\n') 
   if (previous.FIM=="") return(list(init.fed=fed$init.fed,data.fed=fed$data.fed,opti.fed=opti.fed,
   mfisher=mfisher,determinant=determinant,crit=crit,se=se,cv=cv,EigenValues=EigenValues))
   else return(list(init.fed=fed$init.fed,data.fed=fed$data.fed,
   mfisher=mfisher,determinant=determinant,crit=crit,se=se,cv=cv,EigenValues=EigenValues))
   }
   if (FIM=="I") {cat('\n','\n','Optimisation of the individual design...','\n','\n')
   return(list(init.fed=fed$init.fed,data.fed=fed$data.fed,
   mfisher=mfisher,determinant=determinant,crit=crit,se=se,cv=cv,EigenValues=EigenValues)) }
   if (FIM=="B") {cat('\n','\n','Optimisation of the bayesian design...','\n','\n')
   return(list(init.fed=fed$init.fed,data.fed=fed$data.fed,
   mfisher=mfisher,determinant=determinant,crit=crit,se=se,cv=cv,Shrinkage=sh,EigenValues=EigenValues)) }
}

