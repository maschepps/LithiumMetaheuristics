#PFIM 4.0 
#April 2014
#Copyright (C) PFIM 4.0 - Universite Paris Diderot and INSERM.
#------------------------------------------------------------------------------------  

#-----------------------------------------------------------------------------------
# dose: value of the dose 
# N: Number of doses
#	tau: interval between two doses
# Mutiple doses model (_md): model after N doses with a tau interval
# Steady-state model (_ss): model at steady state for doses given with a tau interval. 
#-----------------------------------------------------------------------------------

################################################################################
################################################################################
#######################                                  #######################
#######################       One-compartment models     #######################
#######################                                  #######################
################################################################################
################################################################################

#-------------------########################################--------------------
#-------------------########################################--------------------
#-------------------########                    ############--------------------
#-------------------######## Linear elimination ############--------------------
#-------------------########                    ############--------------------
#-------------------########################################--------------------
#-------------------########################################--------------------


#-------------------------
# Bolus
#-------------------------

#---------
# single dose
#---------
#$Model definition
#$V,k
#$
bolus_1cpt_Vk<-function()
{
  form1<-paste("dose/V*(exp(-k*t))")
  form1<-parse(text=form1,n=-1)
  tf<-Inf
  return(list(c(form1),tf))
}

#$Model definition
#$V,Cl
#$
bolus_1cpt_VCl<-function()
{
  form1<-paste("dose/V*(exp(-Cl/V*t))")
  form1<-parse(text=form1,n=-1)
  tf<-Inf
  return(list(c(form1),tf))
}

#---------
# multiples doses
#---------
#$Model definition
#$V,k
#$N,tau
bolus_1cpt_Vk_md<-function(N,tau)
{
  form1N<-paste("dose/(V)*(1-exp(-",N,"*k*",tau,"))/(1-exp(-k*",tau,"))*(exp(-k*(t-(",N,"-1)*",tau,")))")
  form1N<-parse(text=form1N,n=-1)
  tf<-Inf
  return(list(c(form1N),tf))
}

#$Model definition
#$V,Cl
#$N,tau
bolus_1cpt_VCl_md<-function(N,tau)
{
  form1N<-paste("dose/(V)*(1-exp(-",N,"*Cl/V*",tau,"))/(1-exp(-Cl/V*",tau,"))*(exp(-Cl/V*(t-(",N,"-1)*",tau,")))")
  form1N<-parse(text=form1N,n=-1)
  tf<-Inf
  return(list(c(form1N),tf))
}

#---------
# steady-state
#---------
#$Model definition
#$V,k
#$tau
bolus_1cpt_Vk_ss<-function(tau,TimeSS=0)
{
  formSS<-paste("dose/(V)/(1-exp(-k*",tau,"))*(exp(-k*(t-(",TimeSS,"))))")
  formSS<-parse(text=formSS,n=-1)	
  tf<-Inf
  return(list(c(formSS),tf))
}

#$Model definition
#$V,Cl
#$tau
bolus_1cpt_VCl_ss<-function(tau,TimeSS=0)
{
  formSS<-paste("dose/(V)/(1-exp(-Cl/V*",tau,"))*(exp(-Cl/V*(t-(",TimeSS,"))))")
  formSS<-parse(text=formSS,n=-1)	
  tf<-Inf
  return(list(c(formSS),tf))
}


#-------------------------
# Infusion
#-------------------------

#---------
# single dose
#---------
#$Model definition
#$V,k
#$TInf
infusion_1cpt_Vk<-function(TInf)
{
  R<-dose/TInf
  form1<-paste(R,"/(k*V)*(1-exp(-k*t))")
  form1<-parse(text=form1,n=-1)
  form2<-paste(R,"/(k*V)*(1-exp(-k*",TInf,"))*(exp(-k*(t-",TInf,")))")
  form2<-parse(text=form2,n=-1)
  tf<-c(TInf,Inf)
  return(list(c(form1,form2),tf))
}

#$Model definition
#$V,Cl
#$TInf
infusion_1cpt_VCl<-function(TInf)
{
  R<-dose/TInf
  form1<-paste(R,"/(Cl)*(1-exp(-(Cl/V)*t))")
  form1<-parse(text=form1,n=-1)
  form2<-paste(R,"/(Cl)*(1-exp(-(Cl/V)*",TInf,"))*(exp(-(Cl/V)*(t-",TInf,")))")
  form2<-parse(text=form2,n=-1)
  tf<-c(TInf,Inf)
  return(list(c(form1,form2),tf))
}

#---------
# multiples doses
#---------
#---------
#Needed intermediate function
#---------
infusion_1cpt_Ninter_Vk<-function(TInf,N,tau)
{
  R<-dose/TInf
  formN<-paste(R,"/(k*V)*(1-exp(-k*",TInf,"))*(1-exp(-",N,"*k*",tau,"))/(1-exp(-k*",tau,"))*(exp(-k*(t-(",N,"-1)*",tau,"-",TInf,")))")
  return(formN)
}
infusion_1cpt_Ninter_VCl<-function(TInf,N,tau)
{
  R<-dose/TInf
  formN<-paste(R,"/(Cl)*(1-exp(-(Cl/V)*",TInf,"))*(1-exp(-",N,"*(Cl/V)*",tau,"))/(1-exp(-(Cl/V)*",tau,"))*(exp(-(Cl/V)*(t-(",N,"-1)*",tau,"-",TInf,")))")
  return(formN)
}

#$Model definition
#$V,k
#$TInf,N,tau
infusion_1cpt_Vk_md<-function(TInf,N,tau)
{
  R<-dose/TInf	
  form1N<-paste(R,"/(k*V)*(1-exp(-k*(t-(",N,"-1)*",tau,")))+",infusion_1cpt_Ninter_Vk(TInf,N-1,tau))
  form2N<-infusion_1cpt_Ninter_Vk(TInf,N,tau)
  form1N<-parse(text=form1N,n=-1)
  form2N<-parse(text=form2N,n=-1)	
  tf<-c(TInf,Inf)
  return(list(c(form1N,form2N),tf))
}

#$Model definition
#$V,Cl
#$TInf,N,tau
infusion_1cpt_VCl_md<-function(TInf,N,tau)
{
  R<-dose/TInf	
  form1N<-paste(R,"/(Cl)*(1-exp(-(Cl/V)*(t-(",N,"-1)*",tau,")))+",infusion_1cpt_Ninter_VCl(TInf,N-1,tau))
  form2N<-infusion_1cpt_Ninter_VCl(TInf,N,tau)
  form1N<-parse(text=form1N,n=-1)
  form2N<-parse(text=form2N,n=-1)	
  tf<-c(TInf,Inf)
  return(list(c(form1N,form2N),tf))
}

#---------
# steady-state
#---------
#$Model definition
#$V,k
#$TInf,tau
infusion_1cpt_Vk_ss<-function(TInf,tau,TimeSS=0)
{
  R<-dose/TInf
  form1SS<-paste(R,"/(k*V)*((1-exp(-k*(t-",TimeSS,")))+exp(-k*",tau,")*(1-exp(-k*",TInf,"))*(exp(-k*(t-",TimeSS,"-",TInf,")))/(1-exp(-k*",tau,")))" )
  form2SS<-paste(R,"/(k*V)*(1-exp(-k*",TInf,"))*exp(-k*(t-",TimeSS,"-",TInf,"))/(1-exp(-k*",tau,"))" )
  form1SS<-parse(text=form1SS,n=-1)
  form2SS<-parse(text=form2SS,n=-1)	   
  tf<-c(TInf,Inf)
  return(list(c(form1SS,form2SS),tf))
}

#$Model definition
#$V,Cl
#$TInf,tau
infusion_1cpt_VCl_ss<-function(TInf,tau,TimeSS=0)
{
  R<-dose/TInf
  form1SS<-paste(R,"/(Cl)*((1-exp(-(Cl/V)*(t-",TimeSS,")))+exp(-(Cl/V)*",tau,")*(1-exp(-(Cl/V)*",TInf,"))*(exp(-(Cl/V)*(t-",TimeSS,"-",TInf,")))/(1-exp(-(Cl/V)*",tau,")))" )
  form2SS<-paste(R,"/(Cl)*(1-exp(-(Cl/V)*",TInf,"))*exp(-(Cl/V)*(t-",TimeSS,"-",TInf,"))/(1-exp(-(Cl/V)*",tau,"))" )
  form1SS<-parse(text=form1SS,n=-1)
  form2SS<-parse(text=form2SS,n=-1)	
  tf<-c(TInf,Inf)
  return(list(c(form1SS,form2SS),tf))
  
}

#-------------------------
# Oral1
#-------------------------

#---------
# single dose
#---------
#$Model definition                                                                                  
#$ka,V,k
#$
oral1_1cpt_kaVk<-function()
{
  form1<-paste("dose/V*ka/(ka-k)*(exp(-k*t)-exp(-ka*t))")
  form1<-parse(text=form1,n=-1)
  tf<-Inf
  return(list(c(form1),tf))
}

#$Model definition
#$ka,V,Cl
#$
oral1_1cpt_kaVCl<-function()
{
  form1<-paste("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*t)-exp(-ka*t))")
  form1<-parse(text=form1,n=-1)
  tf<-Inf
  return(list(c(form1),tf))
}

#---------
# multiples doses
#---------
#$Model definition
#$ka,V,k
#$N,tau 
oral1_1cpt_kaVk_md<-function(N,tau)
{
  formN<-paste("dose/V*ka/(ka-k)*(exp(-k*(t-(",N,"-1)*",tau,"))*(1-exp(-",N,"*k*",tau,"))/(1-exp(-k*",tau,"))-exp(-ka*(t-(",N,"-1)*",tau,"))*(1-exp(-",N,"*ka*",tau,"))/(1-exp(-ka*",tau,")))")
  formN<-parse(text=formN,n=-1)
  tf<-Inf
  return(list(c(formN),tf))
}

#$Model definition
#$ka,V,Cl
#$N,tau 
oral1_1cpt_kaVCl_md<-function(N,tau)
{
  formN<-paste("dose/V*ka/(ka-Cl/V)*(exp(-Cl/V*(t-(",N,"-1)*",tau,"))*(1-exp(-",N,"*(Cl/V)*",tau,"))/(1-exp(-(Cl/V)*",tau,"))-exp(-ka*(t-(",N,"-1)*",tau,"))*(1-exp(-",N,"*ka*",tau,"))/(1-exp(-ka*",tau,")))")
  formN<-parse(text=formN,n=-1)
  tf<-Inf
  return(list(c(formN),tf))
}

#---------
# steady-state
#---------
#$Model definition
#$ka,V,k
#$tau
oral1_1cpt_kaVk_ss<-function(tau,TimeSS=0)
{
  formSS<-paste("dose/V*ka/(ka-k)*(exp(-k*(t-",TimeSS,"))/(1-exp(-k*",tau,"))-exp(-ka*(t-",TimeSS,"))/(1-exp(-ka*",tau,")))")
  formSS<-parse(text=formSS,n=-1)
  tf<-Inf
  return(list(c(formSS),tf))
}

#$Model definition
#$ka,V,Cl
#$tau
oral1_1cpt_kaVCl_ss<-function(tau,TimeSS=0)
{
  formSS<-paste("dose/V*ka/(ka-(Cl/V))*(exp(-(Cl/V)*(t-",TimeSS,"))/(1-exp(-(Cl/V)*",tau,"))-exp(-ka*(t-",TimeSS,"))/(1-exp(-ka*",tau,")))")
  formSS<-parse(text=formSS,n=-1)
  tf<-Inf
  return(list(c(formSS),tf))
}


#-------------------########################################--------------------
#-------------------########################################--------------------
#-------------------##########                  ############--------------------
#-------------------########## Mickaelis Menten ############--------------------
#-------------------##########                  ############--------------------
#-------------------########################################--------------------
#-------------------########################################--------------------

#-------------------------
# Bolus
#-------------------------

#---------
# single dose
#---------
#$Model definition
#$V,Vm,km
bolus_1cpt_VVmkm<-function(){
  INT_bolus_1cpt_VVmkm<-function(t,y,p){
    V<-p[1]
    Vm <-p[2]
    km<-p[3]
    yd1<-(-Vm)*y[1]/(km+y[1])
    return(list(c(yd1),c(y[1]/V)))
  }
  return(INT_bolus_1cpt_VVmkm)
}

#-------------------------
# Infusion
#-------------------------

#---------
# single dose
#---------
#$Model definition
#$V,Vm,km
#$ dose, TInf
infusion_1cpt_VVmkm<-function(doseMM,TInf){
  INT_infusion_1cpt_VVmkm<-function(t,y, p){
    V<-p[1]
    Vm <-p[2]
    km<-p[3]
    if(t<=TInf){
      yd1<-(-Vm)*y[1]/(km+y[1])+(doseMM/(TInf*V))}
    else{
      yd1<-(-Vm)*y[1]/(km+y[1])}
    
    return(list(c(yd1),c(y[1])))
  }
  return(INT_infusion_1cpt_VVmkm)
}

#---------
# multiples doses
#---------
#$Model definition
#$V,Vm,km
#$ dose, TInf, tau
infusion_1cpt_VVmkm_md<-function(doseMM,TInf,tau){
  INT_infusion_1cpt_VVmkm_md<-function(t,y, p){
    V<-p[1]
    Vm <-p[2]
    km<-p[3]
    if(t%%tau<=TInf){
      yd1<-(-Vm)*y[1]/(km+y[1])+(doseMM/(TInf*V))}
    else{
      yd1<-(-Vm)*y[1]/(km+y[1])}
    return(list(c(yd1),c(y[1])))
  }
  return(INT_infusion_1cpt_VVmkm_md)
}


#-------------------------
# Oral1
#-------------------------

#---------
# single dose
#---------
#$Model definition
#$ka,V,Vm,km
#$ dose
oral1_1cpt_kaVVmkm<-function(doseMM){
  INT_oral1_1cpt_kaVVmkm<-function(t,y,p){
    ka<-p[1]
    V<-p[2]
    Vm <-p[3]
    km<-p[4]
    yd1<-(-Vm)*y[1]/(km+y[1])+(doseMM*ka/V)*exp(-ka*t)
    return(list(c(yd1),c(y[1])))
  }
  return(INT_oral1_1cpt_kaVVmkm)
}

#---------
# multiples doses
#---------
#---------
# Needed intermediate function
#---------
input_oral1<-function(ka,V,doseMM,n.libPK,tau,t){
  if(n<0){stop("error function input_oral1 n<=0")}
  else{ if(n==0){return(doseMM*ka/V*exp(-ka*t))}
    else{return(doseMM*ka/V*exp(-ka*(t-n.libPK*tau))+input_oral1(ka,V,doseMM,(n-1),tau,t))}
  } }

#$Model definition
#$ka,V,Vm,km
#$ dose, tau
oral1_1cpt_kaVVmkm_md<-function(doseMM,tau.lib){
  INT_oral1_1cpt_kaVVmkm_md<-function(t,y,p)
  {
    ka<-p[1]
    V<-p[2]
    Vm <-p[3]
    km<-p[4]
    n<-t%/%tau.lib
    input<-input_oral1(ka,V,doseMM,n,tau.lib,t)
    yd1<-(-Vm)*y[1]/(km+y[1])+input
    return(list(c(yd1),c(y[1])))
  }
  return(INT_oral1_1cpt_kaVVmkm_md)
}





################################################################################
################################################################################
#######################                                  #######################
#######################       Two-compartment models     #######################
#######################                                  #######################
################################################################################
################################################################################

#-------------------########################################--------------------
#-------------------########################################--------------------
#-------------------########                    ############--------------------
#-------------------######## Linear elimination ############--------------------
#-------------------########                    ############--------------------
#-------------------########################################--------------------
#-------------------########################################--------------------

#-------------------------
#Needed prefunction
#-------------------------
myprefunc<-function(k10,k12,k21,V)
{
  s<-paste("(",k10,"+",k21,"+",k12,")")
  rac2<-paste("sqrt(",s,"^2-4*",k10,"*",k21,")")
  alpha1<-paste("((",s,"+",rac2,")/2)")
  alpha2<-paste("((",s,"-",rac2,")/2)")
  
  A<-paste("1*(",alpha1,"-",k21,")/(",V,"*(",alpha1,"-",alpha2,"))")
  B<-paste("1*(",k21,"-",alpha2,")/(",V,"*(",alpha1,"-",alpha2,"))")
  return(list(alpha1,alpha2,A,B))
}
myprefuncoral<-function(ka,k10,k12,k21,V)
{
  s<-paste("(",k10,"+",k21,"+",k12,")")
  rac2<-paste("sqrt(",s,"^2-4*",k10,"*",k21,")")
  alpha1<-paste("((",s,"+",rac2,")/2)")
  alpha2<-paste("((",s,"-",rac2,")/2)")
  
  A<-paste(ka,"/",V,"*(",k21,"-",ka,")/((",alpha1,"-",ka,")*(",alpha2,"-",ka,"))")
  B<-paste(ka,"/",V,"*(",k21,"-",alpha1,")/((",ka,"-",alpha1,")*(",alpha2,"-",alpha1,"))")
  C<-paste(ka,"/",V,"*(",k21,"-",alpha2,")/((",ka,"-",alpha2,")*(",alpha1,"-",alpha2,"))")      
  return(list(alpha1,alpha2,A,B,C))
}

#-------------------------
# Bolus
#-------------------------

#---------
# single dose
#---------
#$Model definition
#$V,k,k12,k21
#$
bolus_2cpt_Vkk12k21<-function()
{
  prevalue<-myprefunc(k10="k",k12="k12",k21="k21",V="V")
  A<-prevalue[[3]]
  B<-prevalue[[4]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  form1<-paste("dose*(",A,"*exp(-",alpha1,"*t)+",B,"*exp(-",alpha2,"*t))")
  form1<-parse(text=form1,n=-1)
  tf<-Inf
  return(list(c(form1),tf))
}

#$Model definition
#$Cl,V1,Q,V2
#$
bolus_2cpt_ClV1QV2<-function()
{
  prevalue<-myprefunc(k10="Cl/V1",k12="Q/V1",k21="Q/V2",V="V1")
  A<-prevalue[[3]]
  B<-prevalue[[4]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  form1<-paste("dose*(",A,"*exp(-",alpha1,"*t)+",B,"*exp(-",alpha2,"*t))")
  form1<-parse(text=form1,n=-1)
  tf<-Inf
  return(list(c(form1),tf))
}

#---------
# multiples doses
#---------
#$Model definition
#$V,k,k12,k21
#$N,tau
bolus_2cpt_Vkk12k21_md<-function(N,tau)
{
  prevalue<-myprefunc(k10="k",k12="k12",k21="k21",V="V")
  A<-prevalue[[3]]
  B<-prevalue[[4]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  multiterm1<-paste("((1-exp(-",N,"*",alpha1,"*",tau,"))/(1-exp(-",alpha1,"*",tau,")))")
  multiterm2<-paste("((1-exp(-",N,"*",alpha2,"*",tau,"))/(1-exp(-",alpha2,"*",tau,")))")
  formN<-paste("dose*(",A,"*exp(-",alpha1,"*(t-(",N,"-1)*",tau,"))*",multiterm1,"+",B,"*exp(-",alpha2,"*(t-(",N,"-1)*",tau,"))*",multiterm2,")")	
  formN<-parse(text=formN,n=-1)
  tf<-Inf
  return(list(c(formN),tf))
}

#$Model definition
#$Cl,V1,Q,V2
#$N,tau
bolus_2cpt_ClV1QV2_md<-function(N,tau)
{
  prevalue<-myprefunc(k="Cl/V1",k12="Q/V1",k21="Q/V2",V="V1")
  A<-prevalue[[3]]
  B<-prevalue[[4]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  multiterm1<-paste("((1-exp(-",N,"*",alpha1,"*",tau,"))/(1-exp(-",alpha1,"*",tau,")))")
  multiterm2<-paste("((1-exp(-",N,"*",alpha2,"*",tau,"))/(1-exp(-",alpha2,"*",tau,")))")
  formN<-paste("dose*(",A,"*exp(-",alpha1,"*(t-(",N,"-1)*",tau,"))*",multiterm1,"+",B,"*exp(-",alpha2,"*(t-(",N,"-1)*",tau,"))*",multiterm2,")")
  formN<-parse(text=formN,n=-1)
  tf<-Inf
  return(list(c(formN),tf))
}

#---------
# steady-state
#---------
#$Model definition
#$V,k,k12,k21
#$tau
bolus_2cpt_Vkk12k21_ss<-function(tau,TimeSS=0)
{
  prevalue<-myprefunc(k10="k",k12="k12",k21="k21",V="V")
  A<-prevalue[[3]]
  B<-prevalue[[4]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  multiterm1<-paste("(1/(1-exp(-",alpha1,"*",tau,")))")
  multiterm2<-paste("(1/(1-exp(-",alpha2,"*",tau,")))")
  formSS<-paste("dose*(",A,"*exp(-",alpha1,"*(t-",TimeSS,"))*",multiterm1,"+",B,"*exp(-",alpha2,"*(t-",TimeSS,"))*",multiterm2,")")
  formSS<-parse(text=formSS,n=-1)	
  tf<-Inf
  return(list(c(formSS),tf))
}

#$Model definition
#$Cl,V1,Q,V2
#$tau
bolus_2cpt_ClV1QV2_ss<-function(tau,TimeSS=0)
{
  prevalue<-myprefunc(k="Cl/V1",k12="Q/V1",k21="Q/V2",V="V1")
  A<-prevalue[[3]]
  B<-prevalue[[4]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  multiterm1<-paste("(1/(1-exp(-",alpha1,"*",tau,")))")
  multiterm2<-paste("(1/(1-exp(-",alpha2,"*",tau,")))")
  formSS<-paste("dose*(",A,"*exp(-",alpha1,"*(t-",TimeSS,"))*",multiterm1,"+",B,"*exp(-",alpha2,"*(t-",TimeSS,"))*",multiterm2,")")
  formSS<-parse(text=formSS,n=-1)	
  tf<-Inf
  return(list(c(formSS),tf))
}

#-------------------------
# Infusion
#-------------------------

#---------
# single dose
#---------
#$Model definition
#$V,k,k12,k21
#$TInf
infusion_2cpt_Vkk12k21<-function(TInf)
{
  R<-dose/TInf
  prevalue<-myprefunc(k10="k",k12="k12",k21="k21",V="V")
  A<-prevalue[[3]]
  B<-prevalue[[4]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  form1<-paste(R,"*(",A,"/",alpha1,"*(1-exp(-",alpha1,"*t))+",B,"/",alpha2,"*(1-exp(-",alpha2,"*t)))")
  form1<-parse(text=form1,n=-1)
  form2<-paste(R,"*(",A,"/",alpha1,"*(",1,"-exp(-",alpha1,"*",TInf,"))*exp(-",alpha1,"*(t-",TInf,"))+",B,"/",alpha2,"*(",1,"-exp(-",alpha2,"*",TInf,"))*exp(-",alpha2,"*(t-",TInf,")))")
  form2<-parse(text=form2,n=-1)
  tf<-c(TInf,Inf)
  return(list(c(form1,form2),tf))
}

#$Model definition
#$Cl,V1,Q,V2
#$TInf
infusion_2cpt_ClV1QV2<-function(TInf)
{
  R<-dose/TInf
  prevalue<-myprefunc(k="Cl/V1",k12="Q/V1",k21="Q/V2",V="V1")      
  A<-prevalue[[3]]
  B<-prevalue[[4]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  form1<-paste(R,"*(",A,"/",alpha1,"*(1-exp(-",alpha1,"*t))+",B,"/",alpha2,"*(1-exp(-",alpha2,"*t)))")
  form1<-parse(text=form1,n=-1)
  form2<-paste(R,"*(",A,"/",alpha1,"*(",1,"-exp(-",alpha1,"*",TInf,"))*exp(-",alpha1,"*(t-",TInf,"))+",B,"/",alpha2,"*(",1,"-exp(-",alpha2,"*",TInf,"))*exp(-",alpha2,"*(t-",TInf,")))")
  form2<-parse(text=form2,n=-1)
  tf<-c(TInf,Inf)
  return(list(c(form1,form2),tf))
}

#---------
# multiples doses
#---------
#---------
#Needed intermediate function
#---------
infusion_2cpt_Ninter_Vkk12k21<-function(TInf,N,tau)
{
  R<-dose/TInf
  prevalue<-myprefunc(k10="k",k12="k12",k21="k21",V="V")
  A<-prevalue[[3]]
  B<-prevalue[[4]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  multiterm1<-paste("((1-exp(-",N,"*",alpha1,"*",tau,"))/(1-exp(-",alpha1,"*",tau,")))")
  multiterm2<-paste("((1-exp(-",N,"*",alpha2,"*",tau,"))/(1-exp(-",alpha2,"*",tau,")))")
  formN<-paste(R,"*(",A,"/",alpha1,"*(",1,"-exp(-",alpha1,"*",TInf,"))*",multiterm1,"*exp(-",alpha1,"*(t-(",N,"-",1,")*",tau,"-",TInf,"))+",B,"/",alpha2,"*(",1,"-exp(-",alpha2,"*",TInf,"))*",multiterm2,"*exp(-",alpha2,"*(t-(",N,"-",1,")*",tau,"-",TInf,")))")
  return(formN)
}
infusion_2cpt_Ninter_ClV1QV2<-function(TInf,N,tau)
{
  R<-dose/TInf
  prevalue<-myprefunc(k10="Cl/V1",k12="Q/V1",k21="Q/V2",V="V1")
  A<-prevalue[[3]]
  B<-prevalue[[4]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  multiterm1<-paste("((1-exp(-",N,"*",alpha1,"*",tau,"))/(1-exp(-",alpha1,"*",tau,")))")
  multiterm2<-paste("((1-exp(-",N,"*",alpha2,"*",tau,"))/(1-exp(-",alpha2,"*",tau,")))")
  formN<-paste(R,"*(",A,"/",alpha1,"*(",1,"-exp(-",alpha1,"*",TInf,"))*",multiterm1,"*exp(-",alpha1,"*(t-(",N,"-",1,")*",tau,"-",TInf,"))+",B,"/",alpha2,"*(",1,"-exp(-",alpha2,"*",TInf,"))*",multiterm2,"*exp(-",alpha2,"*(t-(",N,"-",1,")*",tau,"-",TInf,")))")
  return(formN)
}

#$Model definition
#$V,k,k12,k21
#$TInf,N,tau
infusion_2cpt_Vkk12k21_md<-function(TInf,N,tau)
{
  R<-dose/TInf
  prevalue<-myprefunc(k10="k",k12="k12",k21="k21",V="V")
  A<-prevalue[[3]]                                                                                                                                
  B<-prevalue[[4]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]	
  form1N<-paste(R,"*(",A,"/",alpha1,"*(1-exp(-",alpha1,"*(t-(",N,"-1)*",tau,")))+",B,"/",alpha2,"*(1-exp(-",alpha2,"*(t-(",N,"-1)*",tau,"))))+",infusion_2cpt_Ninter_Vkk12k21(TInf,N-1,tau))
  form2N<-infusion_2cpt_Ninter_Vkk12k21(TInf,N,tau)
  form1N<-parse(text=form1N,n=-1)
  form2N<-parse(text=form2N,n=-1)
  tf<-c(TInf+(N-1)*tau,Inf)
  return(list(c(form1N,form2N),tf))	
}

#$Model definition                                           
#$Cl,V1,Q,V2
#$TInf,N,tau
infusion_2cpt_ClV1QV2_md<-function(TInf,N,tau)
{
  R<-dose/TInf
  prevalue<-myprefunc(k10="Cl/V1",k12="Q/V1",k21="Q/V2",V="V1")
  A<-prevalue[[3]]
  B<-prevalue[[4]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]	
  form1N<-paste(R,"*(",A,"/",alpha1,"*(1-exp(-",alpha1,"*(t-(",N,"-1)*",tau,")))+",B,"/",alpha2,"*(1-exp(-",alpha2,"*(t-(",N,"-1)*",tau,"))))+",infusion_2cpt_Ninter_ClV1QV2(TInf,N-1,tau))
  form2N<-infusion_2cpt_Ninter_ClV1QV2(TInf,N,tau)
  form1N<-parse(text=form1N,n=-1)
  form2N<-parse(text=form2N,n=-1)	
  tf<-c(TInf+(N-1)*tau,Inf)
  return(list(c(form1N,form2N),tf))	
}

#---------
# steady-state
#---------
#$Model definition
#$V,k,k12,k21
#$TInf,tau
infusion_2cpt_Vkk12k21_ss<-function(TInf,tau)
{
  R<-dose/TInf
  prevalue<-myprefunc(k10="k",k12="k12",k21="k21",V="V")
  A<-prevalue[[3]]
  B<-prevalue[[4]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  
  form1N<-paste(R,"*(",A,"/",alpha1,"*((1-exp(-",alpha1,"*t))+exp(-",alpha1,"*",tau,")*
                (1-exp(-",alpha1,"*",TInf,"))*exp(-",alpha1,"*(t-",TInf,"))/(1-exp(-",alpha1,"*",tau,")))
                +",B,"/",alpha2,"*((1-exp(-",alpha2,"*t))+exp(-",alpha2,"*",tau,")*(1-exp(-",alpha2,"*",TInf,"))
                *exp(-",alpha2,"*(t-",TInf,"))/(1-exp(-",alpha2,"*",tau,"))))")
  form2N<-paste(R,"*(",A,"/",alpha1,"*((1-exp(-",alpha1,"*",TInf,"))*exp(-",alpha1,"*(t-",TInf,"))/(1-exp(-",alpha1,"*",tau,")))+",B,"/",alpha2,"*((1-exp(-",alpha2,"*",TInf,"))*exp(-",alpha2,"*(t-",TInf,"))/(1-exp(-",alpha2,"*",tau,"))))")
  
  form1N<-parse(text=form1N,n=-1)
  form2N<-parse(text=form2N,n=-1)
  tf<-c(TInf,Inf)
  return(list(c(form1N,form2N),tf))
  
}

#$Model definition
#$Cl,V1,Q,V2
#$TInf,tau
infusion_2cpt_ClV1QV2_ss<-function(TInf,tau)
{
  R<-dose/TInf                     
  prevalue<-myprefunc(k10="Cl/V1",k12="Q/V1",k21="Q/V2",V="V1")
  A<-prevalue[[3]]
  B<-prevalue[[4]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  
  form1N<-paste(R,"*(",A,"/",alpha1,"*((1-exp(-",alpha1,"*t))+exp(-",alpha1,"*",tau,")*
                (1-exp(-",alpha1,"*",TInf,"))*exp(-",alpha1,"*(t-",TInf,"))/(1-exp(-",alpha1,"*",tau,")))
                +",B,"/",alpha2,"*((1-exp(-",alpha2,"*t))+exp(-",alpha2,"*",tau,")*(1-exp(-",alpha2,"*",TInf,"))
                *exp(-",alpha2,"*(t-",TInf,"))/(1-exp(-",alpha2,"*",tau,"))))")
  form2N<-paste(R,"*(",A,"/",alpha1,"*((1-exp(-",alpha1,"*",TInf,"))*exp(-",alpha1,"*(t-",TInf,"))/(1-exp(-",alpha1,"*",tau,")))+",B,"/",alpha2,"*((1-exp(-",alpha2,"*",TInf,"))*exp(-",alpha2,"*(t-",TInf,"))/(1-exp(-",alpha2,"*",tau,"))))")
  
  form1N<-parse(text=form1N,n=-1)
  form2N<-parse(text=form2N,n=-1)
  tf<-c(TInf,Inf)
  return(list(c(form1N,form2N),tf))
  
}

#-------------------------
# Oral1
#-------------------------

#---------
# single dose
#---------
#$Model definition
#$ka,V,k,k12,k21
#$
oral1_2cpt_kaVkk12k21<-function()
{
  prevalue<-myprefuncoral(ka="ka",k10="k",k12="k12",k21="k21",V="V")
  A<-prevalue[[3]]
  B<-prevalue[[4]]
  C<-prevalue[[5]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  
  form1<-paste("dose*(",A,"*exp(-ka*t)+",B,"*exp(-",alpha1,"*t)+",C,"*exp(-",alpha2,"*t))")
  form1<-parse(text=form1,n=-1)
  tf<-Inf
  return(list(c(form1),tf))
}

#$Model definition
#$ka,Cl,V1,Q,V2
#$
oral1_2cpt_kaClV1QV2<-function()
{                                                 
  prevalue<-myprefuncoral(ka="ka",k10="Cl/V1",k12="Q/V1",k21="Q/V2",V="V1")
  A<-prevalue[[3]]
  B<-prevalue[[4]]
  C<-prevalue[[5]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  
  form1<-paste("dose*(",A,"*exp(-ka*t)+",B,"*exp(-",alpha1,"*t)+",C,"*exp(-",alpha2,"*t))")
  form1<-parse(text=form1,n=-1)
  tf<-Inf
  return(list(c(form1),tf))
}

#---------
# multiples doses
#---------
#$Model definition
#$ka,V,k,k12,k21
#$N,tau
oral1_2cpt_kaVkk12k21_md<-function(N,tau)
{
  prevalue<-myprefuncoral(ka="ka",k10="k",k12="k12",k21="k21",V="V")
  A<-prevalue[[3]]
  B<-prevalue[[4]]
  C<-prevalue[[5]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  formN<-paste("dose*(",A,"*exp(-ka*(t-(",N,"-1)*",tau,"))*(1-exp(-",N,"*ka*",tau,"))/(1-exp(-ka*",tau,"))+",B,"*exp(-",alpha1,"*(t-(",N,"-1)*",tau,"))*(1-exp(-",N,"*",alpha1,"*",tau,"))/(1-exp(-",alpha1,"*",tau,"))+",C,"*exp(-",alpha2,"*(t-(",N,"-1)*",tau,"))*(1-exp(-",N,"*",alpha2,"*",tau,"))/(1-exp(-",alpha2,"*",tau,")))")
  formN<-parse(text=formN,n=-1)
  tf<-Inf
  return(list(c(formN),tf))
}

#$Model definition
#$ka,Cl,V1,Q,V2
#$N,tau
oral1_2cpt_kaClV1QV2_md<-function(N,tau)
{
  prevalue<-myprefuncoral(ka="ka",k10="Cl/V1",k12="Q/V1",k21="Q/V2",V="V1")
  A<-prevalue[[3]]
  B<-prevalue[[4]]
  C<-prevalue[[5]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  formN<-paste("dose*(",A,"*exp(-ka*(t-(",N,"-1)*",tau,"))*(1-exp(-",N,"*ka*",tau,"))/(1-exp(-ka*",tau,"))+",B,"*exp(-",alpha1,"*(t-(",N,"-1)*",tau,"))*(1-exp(-",N,"*",alpha1,"*",tau,"))/(1-exp(-",alpha1,"*",tau,"))+",C,"*exp(-",alpha2,"*(t-(",N,"-1)*",tau,"))*(1-exp(-",N,"*",alpha2,"*",tau,"))/(1-exp(-",alpha2,"*",tau,")))")
  formN<-parse(text=formN,n=-1)
  tf<-Inf
  return(list(c(formN),tf))
}

#---------
# steady-state
#---------
#$Model definition
#$ka,V,k,k12,k21
#$tau                  
oral1_2cpt_kaVkk12k21_ss<-function(tau,TimeSS=0)
{
  prevalue<-myprefuncoral(ka="ka",k10="k",k12="k12",k21="k21",V="V")
  A<-prevalue[[3]]
  B<-prevalue[[4]]
  C<-prevalue[[5]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  formSS<-paste("dose*(",A,"*exp(-ka*(t-",TimeSS,"))/(1-exp(-ka*",tau,"))+",B,"*exp(-",alpha1,"*(t-",TimeSS,"))/(1-exp(-",alpha1,"*",tau,"))+",C,"*exp(-",alpha2,"*(t-",TimeSS,"))/(1-exp(-",alpha2,"*",tau,")))")
  formSS<-parse(text=formSS,n=-1)
  tf<-Inf
  return(list(c(formSS),tf))
}

#$Model definition
#$ka,Cl,V1,Q,V2
#$tau 
oral1_2cpt_kaClV1QV2_ss<-function(tau,TimeSS=0)
{                                                                    
  prevalue<-myprefuncoral(ka="ka",k10="Cl/V1",k12="Q/V1",k21="Q/V2",V="V1")  
  A<-prevalue[[3]]
  B<-prevalue[[4]]
  C<-prevalue[[5]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  formSS<-paste("dose*(",A,"*exp(-ka*(t-",TimeSS,"))/(1-exp(-ka*",tau,"))+",B,"*exp(-",alpha1,"*(t-",TimeSS,"))/(1-exp(-",alpha1,"*",tau,"))+",C,"*exp(-",alpha2,"*(t-",TimeSS,"))/(1-exp(-",alpha2,"*",tau,")))")
  formSS<-parse(text=formSS,n=-1)
  tf<-Inf
  return(list(c(formSS),tf))
}




#-------------------########################################--------------------
#-------------------########################################--------------------
#-------------------##########                  ############--------------------
#-------------------########## Mickaelis Menten ############--------------------
#-------------------##########                  ############--------------------
#-------------------########################################--------------------
#-------------------########################################--------------------

#-------------------------
# Bolus
#-------------------------

#---------
# single dose
#---------
#$Model definition
#$V,k12,k21,Vm,km
bolus_2cpt_Vk12k21Vmkm<-function(){
  INT_bolus_2cpt_Vk12k21Vmkm<-function(t,y,p){
    V<-p[1]
    k12<-p[2]
    k21<-p[3]
    Vm <-p[4]
    km<-p[5]
    yd1<-(-Vm)*y[1]/(km+y[1])-k12*y[1]+k21*y[2]
    yd2<-k12*y[1]-k21*y[2]
    return(list(c(yd1,yd2),c(y[1]/V)))
  }
  return(INT_bolus_2cpt_Vk12k21Vmkm)
}

#$Model definition
#$V1,Q,V2,Vm,km
bolus_2cpt_V1QV2Vmkm<-function(){
  INT_bolus_2cpt_V1QV2Vmkm<-function(t,y,p){
    V1<-p[1]
    Q<-p[2]
    V2<-p[3]
    Vm <-p[4]
    km<-p[5]
    k12<-Q/V1
    k21<-Q/V2
    yd1<-(-Vm)*y[1]/(km+y[1])-k12*y[1]+k21*y[2]
    yd2<-k12*y[1]-k21*y[2]
    return(list(c(yd1,yd2),c(y[1]/V1)))
  }
  return(INT_bolus_2cpt_V1QV2Vmkm)
}

#-------------------------
# Infusion
#-------------------------

#---------
# single dose
#---------
#$Model definition
#$V,k12,k21,Vm,km
#$ dose, TInf
infusion_2cpt_Vk12k21Vmkm<-function(doseMM,TInf){
  INT_infusion_2cpt_Vk12k21Vmkm<-function(t,y, p){
    V<-p[1]
    k12<-p[2]
    k21<-p[3]
    Vm <-p[4]
    km<-p[5]
    if(t<=TInf){
      yd1<-(-Vm)*y[1]/(km+y[1])+(doseMM/(TInf*V))-k12*y[1]+k21*y[2]}
    else{
      yd1<-(-Vm)*y[1]/(km+y[1])-k12*y[1]+k21*y[2]}
    yd2<-k12*y[1]-k21*y[2]
    return(list(c(yd1,yd2),c(y[1])))
  }
  return(INT_infusion_2cpt_Vk12k21Vmkm)
}

#$Model definition
#$V1,Q,V2,Vm,km
#$ dose, TInf
infusion_2cpt_V1QV2Vmkm<-function(doseMM,TInf){
  INT_infusion_2cpt_V1QV2Vmkm<-function(t,y, p){
    V1<-p[1]
    Q<-p[2]
    V2<-p[3]
    Vm <-p[4]
    km<-p[5]
    k12<-Q/V1
    k21<-Q/V2
    if(t<=TInf){
      yd1<-(-Vm)*y[1]/(km+y[1])+(dose/(TInf*V1))-k12*y[1]+k21*y[2]}
    else{
      yd1<-(-Vm)*y[1]/(km+y[1])-k12*y[1]+k21*y[2]}
    yd2<-k12*y[1]-k21*y[2]
    return(list(c(yd1,yd2),c(y[1])))
  }
  return(INT_infusion_2cpt_V1QV2Vmkm)
}

#---------
# multiples doses
#---------
#$Model definition
#$V,k12,k21,Vm,km
#$ dose, TInf, tau
infusion_2cpt_Vk12k21Vmkm_md<-function(doseMM,TInf,tau){
  INT_infusion_2cpt_Vk12k21Vmkm_md<-function(t,y, p){
    V<-p[1]
    k12<-p[2]
    k21<-p[3]
    Vm <-p[4]
    km<-p[5]
    if(t%%tau<=TInf){
      yd1<-(-Vm)*y[1]/(km+y[1])+(doseMM/(TInf*V))-k12*y[1]+k21*y[2]}
    else{
      yd1<-(-Vm)*y[1]/(km+y[1])-k12*y[1]+k21*y[2]}
    yd2<-k12*y[1]-k21*y[2]
    return(list(c(yd1,yd2),c(y[1])))
  }
  return(INT_infusion_2cpt_Vk12k21Vmkm_md)
}

#$Model definition
#$V1,Q,V2,Vm,km
#$ dose, TInf, tau
infusion_2cpt_V1QV2Vmkm_md<-function(dose,TInf,tau){
  INT_infusion_2cpt_V1QV2Vmkm_md<-function(t,y, p){
    V1<-p[1]
    Q<-p[2]
    V2<-p[3]
    Vm <-p[4]
    km<-p[5]
    k12<-Q/V1
    k21<-Q/V2
    if(t%%tau<=TInf){
      yd1<-(-Vm)*y[1]/(km+y[1])+(dose/(TInf*V1))-k12*y[1]+k21*y[2]}
    else{
      yd1<-(-Vm)*y[1]/(km+y[1])-k12*y[1]+k21*y[2]}
    yd2<-k12*y[1]-k21*y[2]
    return(list(c(yd1,yd2),c(y[1])))
  }
  return(INT_infusion_2cpt_V1QV2Vmkm_md)
}

#-------------------------
# Oral1
#-------------------------

#---------
# single dose
#---------
#$Model definition
#$ka,V,k12,k21,Vm,km
#$ dose
oral1_2cpt_kaVk12k21Vmkm<-function(doseMM){
  INT_oral1_2cpt_kaVk12k21Vmkm<-function(t,y,p){
    ka<-p[1]
    V<-p[2]
    k12<-p[3]
    k21<-p[4]
    Vm <-p[5]
    km<-p[6]
    yd1<-(-Vm)*y[1]/(km+y[1])+(doseMM*ka/V)*exp(-ka*t)-k12*y[1]+k21*y[2]
    yd2<-k12*y[1]-k21*y[2]
    return(list(c(yd1,yd2),c(y[1])))
  }
  return(INT_oral1_2cpt_kaVk12k21Vmkm)
}

#$Model definition
#$ka,V1,Q,V2,Vm,km
#$ dose
oral1_2cpt_kaV1QV2Vmkm<-function(doseMM){
  INT_oral1_2cpt_kaV1QV2Vmkm<-function(t,y,p){
    ka<-p[1]
    V1<-p[2]
    Q<-p[3]
    V2<-p[4]
    Vm <-p[5]
    km<-p[6]
    k12<-Q/V1
    k21<-Q/V2
    yd1<-(-Vm)*y[1]/(km+y[1])+(doseMM*ka/V1)*exp(-ka*t)-k12*y[1]+k21*y[2]
    yd2<-k12*y[1]-k21*y[2]
    return(list(c(yd1,yd2),c(y[1])))
  }
  return(INT_oral1_2cpt_kaV1QV2Vmkm)
}

#---------
# multiples doses
#---------
#---------
# Needed intermediate function
#---------
input_oral1<-function(ka,V,doseMM,n,tau,t){
  if(n<0){stop("error function input_oral1 n<=0")}
  else{ if(n==0){return(doseMM*ka/V*exp(-ka*t))}
    else{return(doseMM*ka/V*exp(-ka*(t-n*tau))+input_oral1(ka,V,doseMM,(n-1),tau,t))}
  } }

#$Model definition
#$ka,V,k12,k21,Vm,km
#$ dose, tau
oral1_2cpt_kaVk12k21Vmkm_md<-function(doseMM,tau){
  INT_oral1_2cpt_kaVk12k21Vmkm_md<-function(t,y,p)
  {
    ka<-p[1]
    V<-p[2]
    k12<-p[3]
    k21<-p[4]
    Vm <-p[5]
    km<-p[6]
    n<-t%/%tau
    input<-input_oral1(ka,V,doseMM,n,tau,t)
    yd1<-(-Vm)*y[1]/(km+y[1])+input-k12*y[1]+k21*y[2]
    yd2<-k12*y[1]-k21*y[2]
    return(list(c(yd1,yd2),c(y[1])))
  }
  return(INT_oral1_2cpt_kaVk12k21Vmkm_md)
}

#$Model definition
#$ka,V1,Q,V2,Vm,km
#$ dose,tau
oral1_2cpt_kaV1QV2Vmkm_md<-function(doseMM,tau){
  INT_oral1_2cpt_kaV1QV2Vmkm_md<-function(t,y,p)
  {
    ka<-p[1]
    V1<-p[2]
    Q<-p[3]
    V2<-p[4]
    Vm <-p[5]
    km<-p[6]
    k12<-Q/V1
    k21<-Q/V2
    n<-t%/%tau
    input<-input_oral1(ka,V1,doseMM,n,tau,t)
    yd1<-(-Vm)*y[1]/(km+y[1])+input-k12*y[1]+k21*y[2]
    yd2<-k12*y[1]-k21*y[2]
    return(list(c(yd1,yd2),c(y[1])))
  }
  return(INT_oral1_2cpt_kaV1QV2Vmkm_md)
}

################################################################################
################################################################################
#######################                                  #######################
#######################     Three-compartment models     #######################
#######################                                  #######################
################################################################################
################################################################################

#-------------------########################################--------------------
#-------------------########################################--------------------
#-------------------########                    ############--------------------
#-------------------######## Linear elimination ############--------------------
#-------------------########                    ############--------------------
#-------------------########################################--------------------
#-------------------########################################--------------------

#-------------------------
#Needed prefunction
#-------------------------
myprefunc_3cpt<-function(k10,k12,k21,k13,k31,V)
{
  a0<-paste("(",k10,"*",k21,"*",k31,")")
  a1<-paste("(",k10,"*",k31,"+",k21,"*",k31,"+",k21,"*",k13,"+",k10,"*",k21,"+",k31,"*",k12,")")
  a2<-paste("(",k10,"+",k12,"+",k13,"+",k21,"+",k31,")")
  
  b1<-paste("(",a1,"-",a2,"**",2,"/",3,")")
  b2<-paste("(",2,"*",a2,"**",3,"/",27,"-", a1,"*",a2,"/",3,"+",a0,")")
  
  r1<-paste("(","sqrt(-(",b1,"**",3,"/",27,")))")
  r2<-paste("(",2,"*",r1,"**(",1,"/",3,"))")
  
  alpha0 <- paste("(acos(-", b2,"/(", 2,"*", r1,"))/", 3, ")" )
  alpha1<- paste("(-(cos(",alpha0,")*",r2,"-", a2,"/",3,"))")
  alpha2<- paste("(-(cos(",alpha0,"+(",2,"*",pi,"/",3,"))*",r2,"-",a2,"/",3,"))")
  alpha3<- paste("(-(cos(",alpha0,"+(",4,"*",pi,"/",3,"))*",r2,"-",a2,"/",3,") )")
  
  A <- paste("((",1,"/",V,")*((",k21,"-",alpha1,")/(",alpha1,"-",alpha2,"))*((",k31,"-",alpha1,")/(",alpha1,"-",alpha3,")))")
  B <- paste("((",1,"/",V,")*((",k21,"-",alpha2,")/(",alpha2,"-",alpha1,"))*((",k31,"-",alpha2,")/(",alpha2,"-",alpha3,")))")
  C <- paste("((",1,"/",V,")*((",k21,"-",alpha3,")/(",alpha3,"-",alpha2,"))*((",k31,"-",alpha3,")/(",alpha3,"-",alpha1,")))")
  
  return(list(alpha1,alpha2,alpha3,A,B,C))
}
myprefuncoral_3cpt<-function(ka,k10,k12,k21,k13,k31,V)
{
  a0<-paste("(",k10,"*",k21,"*",k31,")")
  a1<-paste("(",k10,"*",k31,"+",k21,"*",k31,"+",k21,"*",k13,"+",k10,"*",k21,"+",k31,"*",k12,")")
  a2<-paste("(",k10,"+",k12,"+",k13,"+",k21,"+",k31,")")
  
  b1<-paste("(",a1,"-",a2,"**",2,"/",3,")")
  b2<-paste("(",2,"*",a2,"**",3,"/",27,"-", a1,"*",a2,"/",3,"+",a0,")")
  
  r1<-paste("(","sqrt(-(",b1,"**",3,"/",27,")))")
  r2<-paste("(",2,"*",r1,"**(",1,"/",3,"))")
  
  alpha0 <- paste("(acos(-", b2,"/(", 2,"*", r1,"))/", 3, ")" )
  alpha1<- paste("(-(cos(",alpha0,")*",r2,"-", a2,"/",3,"))")
  alpha2<- paste("(-(cos(",alpha0,"+(",2,"*",pi,"/",3,"))*",r2,"-",a2,"/",3,"))")
  alpha3<- paste("(-(cos(",alpha0,"+(",4,"*",pi,"/",3,"))*",r2,"-",a2,"/",3,") )")
  
  A <- paste("(",1,"/",V,")*(",ka,"/(",ka,"-",alpha1,"))*((",k21,"-",alpha1,")/(",alpha1,"-",alpha2,"))*((",k31,"-",alpha1,")/(",alpha1,"-",alpha3,"))")
  B <- paste("(",1,"/",V,")*(",ka,"/(",ka,"-",alpha2,"))*((",k21,"-",alpha2,")/(",alpha2,"-",alpha1,"))*((",k31,"-",alpha2,")/(",alpha2,"-",alpha3,"))")
  C<-  paste("(",1,"/",V,")*(",ka,"/(",ka,"-",alpha3,"))*((",k21,"-",alpha3,")/(",alpha3,"-",alpha2,"))*((",k31,"-",alpha3,")/(",alpha3,"-",alpha1,"))")
  
  return(list(alpha1,alpha2,alpha3,A,B,C))
}

#-------------------------
# Bolus
#-------------------------

#---------
# single dose
#---------
#$Model definition
#$V,k,k12,k21,k13,k31
#$
bolus_3cpt_Vkk12k21k13k31<-function()
{
  prevalue<-myprefunc_3cpt(k10="k",k12="k12",k21="k21",k13="k13",k31="k31",V="V")
  A<-prevalue[[4]]
  B<-prevalue[[5]]
  C<-prevalue[[6]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  alpha3<-prevalue[[3]]
  form1<-paste("dose*(",A,"*exp(-",alpha1,"*t)+",B,"*exp(-",alpha2,"*t)+",C,"*exp(-",alpha3,"*t))")
  form1<-parse(text=form1,n=-1)
  tf<-Inf
  return(list(c(form1),tf))
}

#$Model definition
#$Cl,V1,Q2,V2,Q3,V3
#$
bolus_3cpt_ClV1Q2V2Q3V3<-function()
{
  prevalue<-myprefunc_3cpt(k10="Cl/V1",k12="Q2/V1",k21="Q2/V2",k13="Q3/V1",k31="Q3/V3",V="V1")
  A<-prevalue[[4]]
  B<-prevalue[[5]]
  C<-prevalue[[6]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  alpha3<-prevalue[[3]]
  form1<-paste("dose*(",A,"*exp(-",alpha1,"*t)+",B,"*exp(-",alpha2,"*t)+",C,"*exp(-",alpha3,"*t))")
  form1<-parse(text=form1,n=-1)
  tf<-Inf
  return(list(c(form1),tf))
}

#---------
# multiples doses
#---------
#$Model definition
#$V,k,k12,k21,k13k31
#$N,tau
bolus_3cpt_Vkk12k21k13k31_md<-function(N,tau)
{
  prevalue<-myprefunc_3cpt(k10="k",k12="k12",k21="k21",k13="k13",k31="k31",V="V")
  A<-prevalue[[4]]
  B<-prevalue[[5]]
  C<-prevalue[[6]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  alpha3<-prevalue[[3]]
  multiterm1<-paste("((1-exp(-",N,"*",alpha1,"*",tau,"))/(1-exp(-",alpha1,"*",tau,")))")
  multiterm2<-paste("((1-exp(-",N,"*",alpha2,"*",tau,"))/(1-exp(-",alpha2,"*",tau,")))")
  multiterm3<-paste("((1-exp(-",N,"*",alpha3,"*",tau,"))/(1-exp(-",alpha3,"*",tau,")))")
  formN<-paste("dose*(",A,"*exp(-",alpha1,"*(t-(",N,"-1)*",tau,"))*",multiterm1,"+",B,"*exp(-",alpha2,"*(t-(",N,"-1)*",tau,"))*",multiterm2,"+",C,"*exp(-",alpha3,"*(t-(",N,"-1)*",tau,"))*",multiterm3,")")
  formN<-parse(text=formN,n=-1)
  tf<-Inf
  return(list(c(formN),tf))
}

#$Model definition
#$Cl,V1,Q2,V2,Q3,V3
#$N,tau
bolus_3cpt_ClV1Q2V2Q3V3_md<-function(N,tau)
{
  prevalue<-myprefunc_3cpt(k10="Cl/V1",k12="Q2/V1",k21="Q2/V2",k13="Q3/V1",k31="Q3/V3",V="V1")
  A<-prevalue[[4]]
  B<-prevalue[[5]]
  C<-prevalue[[6]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  alpha3<-prevalue[[3]]
  multiterm1<-paste("((1-exp(-",N,"*",alpha1,"*",tau,"))/(1-exp(-",alpha1,"*",tau,")))")
  multiterm2<-paste("((1-exp(-",N,"*",alpha2,"*",tau,"))/(1-exp(-",alpha2,"*",tau,")))")
  multiterm3<-paste("((1-exp(-",N,"*",alpha3,"*",tau,"))/(1-exp(-",alpha3,"*",tau,")))")
  formN<-paste("dose*(",A,"*exp(-",alpha1,"*(t-(",N,"-1)*",tau,"))*",multiterm1,"+",B,"*exp(-",alpha2,"*(t-(",N,"-1)*",tau,"))*",multiterm2,"+",C,"*exp(-",alpha3,"*(t-(",N,"-1)*",tau,"))*",multiterm3,")")
  formN<-parse(text=formN,n=-1)
  tf<-Inf
  return(list(c(formN),tf))
}

#---------
# steady-state
#---------
#$Model definition
#$V,k,k12,k21k13k31
#$tau
bolus_3cpt_Vkk12k21k13k31_ss<-function(tau,TimeSS=0)
{
  prevalue<-myprefunc_3cpt(k10="k",k12="k12",k21="k21",k13="k13",k31="k31",V="V")
  A<-prevalue[[4]]
  B<-prevalue[[5]]
  C<-prevalue[[6]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  alpha3<-prevalue[[3]]
  multiterm1<-paste("(1/(1-exp(-",alpha1,"*",tau,")))")
  multiterm2<-paste("(1/(1-exp(-",alpha2,"*",tau,")))")
  multiterm3<-paste("(1/(1-exp(-",alpha3,"*",tau,")))")
  formSS<-paste("dose*(",A,"*exp(-",alpha1,"*(t-",TimeSS,"))*",multiterm1,"+",B,"*exp(-",alpha2,"*(t-",TimeSS,"))*",multiterm2,"+",C,"*exp(-",alpha3,"*(t-",TimeSS,"))*",multiterm3,")")
  formSS<-parse(text=formSS,n=-1)
  tf<-Inf
  return(list(c(formSS),tf))
}

#$Model definition
#$Cl,V1,Q2,V2,Q3,V3
#$tau
bolus_3cpt_ClV1Q2V2Q3V3_ss<-function(tau,TimeSS=0)
{
  prevalue<-myprefunc_3cpt(k10="Cl/V1",k12="Q2/V1",k21="Q2/V2",k13="Q3/V1",k31="Q3/V3",V="V1")
  A<-prevalue[[4]]
  B<-prevalue[[5]]
  C<-prevalue[[6]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  alpha3<-prevalue[[3]]	
  multiterm1<-paste("(1/(1-exp(-",alpha1,"*",tau,")))")
  multiterm2<-paste("(1/(1-exp(-",alpha2,"*",tau,")))")
  multiterm3<-paste("(1/(1-exp(-",alpha3,"*",tau,")))")
  formSS<-paste("dose*(",A,"*exp(-",alpha1,"*(t-",TimeSS,"))*",multiterm1,"+",B,"*exp(-",alpha2,"*(t-",TimeSS,"))*",multiterm2,"+",C,"*exp(-",alpha3,"*(t-",TimeSS,"))*",multiterm3,")")
  formSS<-parse(text=formSS,n=-1)
  tf<-Inf
  return(list(c(formSS),tf))
}


#-------------------------
# Infusion
#-------------------------

#---------
# single dose
#---------
#$Model definition
#$V,k,k12,k21,k13,k31
#$TInf
infusion_3cpt_Vkk12k21k13k31<-function(TInf)
{
  R<-dose/TInf
  prevalue<-myprefunc_3cpt(k10="k",k12="k12",k21="k21",k13="k13",k31="k31",V="V")
  A<-prevalue[[4]]
  B<-prevalue[[5]]
  C<-prevalue[[6]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  alpha3<-prevalue[[3]]
  
  form1<-paste(R,"*(",A,"/",alpha1,"*(",1,"-exp(-",alpha1,"*t))+",B,"/",alpha2,"*(",1,"-exp(-",alpha2,"*t))+",C,"/",alpha3,"*(",1,"-exp(-",alpha3,"*t)))")
  form1<-parse(text=form1,n=-1)
  form2<-paste(R,"*(",A,"/",alpha1,"*(",1,"-exp(-",alpha1,"*",TInf,"))*exp(-",alpha1,"*(t-",TInf,"))+",B,"/",alpha2,"*(",1,"-exp(-",alpha2,"*",TInf,"))*exp(-",alpha2,"*(t-",TInf,"))+",C,"/",alpha3,"*(",1,"-exp(-",alpha3,"*",TInf,"))*exp(-",alpha3,"*(t-",TInf,")))")
  form2<-parse(text=form2,n=-1)
  tf<-c(TInf,Inf)
  return(list(c(form1,form2),tf))
}

#$Model definition
#$Cl,V1,Q2,V2,Q3,V3
#$TInf
infusion_3cpt_ClV1Q2V2Q3V3<-function(TInf)
{
  R<-dose/TInf
  prevalue<-myprefunc_3cpt(k10="Cl/V1",k12="Q2/V1",k21="Q2/V2",k13="Q3/V1",k31="Q3/V3",V="V1")
  A<-prevalue[[4]]
  B<-prevalue[[5]]
  C<-prevalue[[6]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  alpha3<-prevalue[[3]]
  
  form1<-paste(R,"*(",A,"/",alpha1,"*(",1,"-exp(-",alpha1,"*t))+",B,"/",alpha2,"*(",1,"-exp(-",alpha2,"*t))+",C,"/",alpha3,"*(",1,"-exp(-",alpha3,"*t)))")
  form1<-parse(text=form1,n=-1)
  form2<-paste(R,"*(",A,"/",alpha1,"*(",1,"-exp(-",alpha1,"*",TInf,"))*exp(-",alpha1,"*(t-",TInf,"))+",B,"/",alpha2,"*(",1,"-exp(-",alpha2,"*",TInf,"))*exp(-",alpha2,"*(t-",TInf,"))+",C,"/",alpha3,"*(",1,"-exp(-",alpha3,"*",TInf,"))*exp(-",alpha3,"*(t-",TInf,")))")
  form2<-parse(text=form2,n=-1)
  tf<-c(TInf,Inf)
  return(list(c(form1,form2),tf))
}

#---------
# multiples doses
#---------
#---------
# Needed intermediate functions
#---------
infusion_3cpt_Ninter_Vkk12k21k13k31<-function(TInf,N,tau)
{
  R<-dose/TInf
  prevalue<-myprefunc_3cpt(k10="k",k12="k12",k21="k21",k13="k13",k31="k31",V="V")
  A<-prevalue[[4]]
  B<-prevalue[[5]]
  C<-prevalue[[6]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  alpha3<-prevalue[[3]]	
  
  multiterm1<-paste("((1-exp(-",N,"*",alpha1,"*",tau,"))/(1-exp(-",alpha1,"*",tau,")))")
  multiterm2<-paste("((1-exp(-",N,"*",alpha2,"*",tau,"))/(1-exp(-",alpha2,"*",tau,")))")
  multiterm3<-paste("((1-exp(-",N,"*",alpha3,"*",tau,"))/(1-exp(-",alpha3,"*",tau,")))")
  
  formN<-paste(R,"*(",A,"/",alpha1,"*(",1,"-exp(-",alpha1,"*",TInf,"))*",multiterm1,"*exp(-",alpha1,"*(t-(",N,"-",1,")*",tau,"-",TInf,"))+",B,"/",alpha2,"*(",1,"-exp(-",alpha2,"*",TInf,"))*",multiterm2,"*exp(-",alpha2,"*(t-(",N,"-",1,")*",tau,"-",TInf,"))+",C,"/",alpha3,"*(",1,"-exp(-",alpha3,"*",TInf,"))*",multiterm3,"*exp(-",alpha3,"*(t-(",N,"-",1,")*",tau,"-",TInf,")))")
  return(formN)
}
infusion_3cpt_Ninter_ClV1Q2V2Q3V3<-function(TInf,N,tau)
{
  R<-dose/TInf
  prevalue<-myprefunc_3cpt(k10="Cl/V1",k12="Q2/V1",k21="Q2/V2",k13="Q3/V1",k31="Q3/V3",V="V1")
  A<-prevalue[[4]]
  B<-prevalue[[5]]
  C<-prevalue[[6]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  alpha3<-prevalue[[3]]	
  
  multiterm1<-paste("((1-exp(-",N,"*",alpha1,"*",tau,"))/(1-exp(-",alpha1,"*",tau,")))")
  multiterm2<-paste("((1-exp(-",N,"*",alpha2,"*",tau,"))/(1-exp(-",alpha2,"*",tau,")))")
  multiterm3<-paste("((1-exp(-",N,"*",alpha3,"*",tau,"))/(1-exp(-",alpha3,"*",tau,")))")
  
  formN<-paste(R,"*(",A,"/",alpha1,"*(",1,"-exp(-",alpha1,"*",TInf,"))*",multiterm1,"*exp(-",alpha1,"*(t-(",N,"-",1,")*",tau,"-",TInf,"))+",B,"/",alpha2,"*(",1,"-exp(-",alpha2,"*",TInf,"))*",multiterm2,"*exp(-",alpha2,"*(t-(",N,"-",1,")*",tau,"-",TInf,"))+",C,"/",alpha3,"*(",1,"-exp(-",alpha3,"*",TInf,"))*",multiterm3,"*exp(-",alpha3,"*(t-(",N,"-",1,")*",tau,"-",TInf,")))")
  return(formN)
}

#$Model definition
#$V,k,k12,k21,k13,k31
#$TInf,N,tau
infusion_3cpt_Vkk12k21k13k31_md<-function(TInf,N,tau)
{
  R<-dose/TInf
  prevalue<-myprefunc_3cpt(k10="k",k12="k12",k21="k21",k13="k13",k31="k31",V="V")
  A<-prevalue[[4]]
  B<-prevalue[[5]]
  C<-prevalue[[6]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  alpha3<-prevalue[[3]]	
  
  form1N<-paste(R,"*(",A,"/",alpha1,"*(",1,"-exp(-",alpha1,"*(t-(",N,"-",1,")*",tau,")))+",B,"/",alpha2,"*(",1,"-exp(-",alpha2,"*(t-(",N,"-",1,")*",tau,")))+",C,"/",alpha3,"*(",1,"-exp(-",alpha3,"*(t-(",N,"-",1,")*",tau,"))))+",infusion_3cpt_Ninter_Vkk12k21k13k31(TInf,N-1,tau))
  form2N<-infusion_3cpt_Ninter_Vkk12k21k13k31(TInf,N,tau)
  form1N<-parse(text=form1N,n=-1)
  form2N<-parse(text=form2N,n=-1)
  tf<-c(TInf+(N-1)*tau,Inf)
  return(list(c(form1N,form2N),tf))
}

#$Model definition
#$Cl,V1,Q2,V2,Q3,V3
#$TInf,N,tau
infusion_3cpt_ClV1Q2V2Q3V3_md<-function(TInf,N,tau)
{
  R<-dose/TInf
  prevalue<-myprefunc_3cpt(k10="Cl/V1",k12="Q2/V1",k21="Q2/V2",k13="Q3/V1",k31="Q3/V3",V="V1")
  A<-prevalue[[4]]
  B<-prevalue[[5]]
  C<-prevalue[[6]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  alpha3<-prevalue[[3]]	
  
  form1N<-paste(R,"*(",A,"/",alpha1,"*(",1,"-exp(-",alpha1,"*(t-(",N,"-",1,")*",tau,")))+",B,"/",alpha2,"*(",1,"-exp(-",alpha2,"*(t-(",N,"-",1,")*",tau,")))+",C,"/",alpha3,"*(",1,"-exp(-",alpha3,"*(t-(",N,"-",1,")*",tau,"))))+",infusion_3cpt_Ninter_ClV1Q2V2Q3V3(TInf,N-1,tau))
  form2N<-infusion_3cpt_Ninter_ClV1Q2V2Q3V3(TInf,N,tau)
  form1N<-parse(text=form1N,n=-1)
  form2N<-parse(text=form2N,n=-1)
  tf<-c(TInf+(N-1)*tau,Inf)
  return(list(c(form1N,form2N),tf))
}

#---------
# steady-state
#---------
#$Model definition
#$V,k,k12,k21,k13,k31
#$TInf,tau
infusion_3cpt_Vkk12k21k13k31_ss<-function(TInf,tau)
{
  R<-dose/TInf
  prevalue<-myprefunc_3cpt(k10="k",k12="k12",k21="k21",k13="k13",k31="k31",V="V")
  A<-prevalue[[4]]
  B<-prevalue[[5]]
  C<-prevalue[[6]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  alpha3<-prevalue[[3]]	
  
  form1N<-paste(R,"*(",A,"/",alpha1,"*((",1,"-exp(-",alpha1,"*t))+exp(-",alpha1,"*",tau,")*(",1,"-exp(-",alpha1,"*",TInf,"))*exp(-",alpha1,"*(t-",TInf,"))/(",1,"-exp(-",alpha1,"*",tau,")))+",
                B,"/",alpha2,"*((",1,"-exp(-",alpha2,"*t))+exp(-",alpha2,"*",tau,")*(",1,"-exp(-",alpha2,"*",TInf,"))*exp(-",alpha2,"*(t-",TInf,"))/(",1,"-exp(-",alpha2,"*",tau,")))+",
                C,"/",alpha3,"*((",1,"-exp(-",alpha3,"*t))+exp(-",alpha3,"*",tau,")*(",1,"-exp(-",alpha3,"*",TInf,"))*exp(-",alpha3,"*(t-",TInf,"))/(",1,"-exp(-",alpha3,"*",tau,"))))")
  form2N<-paste(R,"*(",A,"/",alpha1,"*((1-exp(-",alpha1,"*",TInf,"))*exp(-",alpha1,"*(t-",TInf,"))/(1-exp(-",alpha1,"*",tau,")))+",B,"/",alpha2,"*((1-exp(-",alpha2,"*",TInf,"))*exp(-",alpha2,"*(t-",TInf,"))/(1-exp(-",alpha2,"*",tau,")))+",C,"/",alpha3,"*((1-exp(-",alpha3,"*",TInf,"))*exp(-",alpha3,"*(t-",TInf,"))/(1-exp(-",alpha3,"*",tau,"))))")
  
  form1N<-parse(text=form1N,n=-1)
  form2N<-parse(text=form2N,n=-1)
  tf<-c(TInf,Inf)
  return(list(c(form1N,form2N),tf))
  
}

#$Model definition
#$Cl,V1,Q2,V2,Q3,V3
#$TInf,tau
infusion_3cpt_ClV1Q2V2Q3V3_ss<-function(TInf,tau)
{
  R<-dose/TInf
  prevalue<-myprefunc_3cpt(k10="Cl/V1",k12="Q2/V1",k21="Q2/V2",k13="Q3/V1",k31="Q3/V3",V="V1")
  A<-prevalue[[4]]
  B<-prevalue[[5]]
  C<-prevalue[[6]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  alpha3<-prevalue[[3]]	
  
  form1N<-paste(R,"*(",A,"/",alpha1,"*((",1,"-exp(-",alpha1,"*t))+exp(-",alpha1,"*",tau,")*(",1,"-exp(-",alpha1,"*",TInf,"))*exp(-",alpha1,"*(t-",TInf,"))/(",1,"-exp(-",alpha1,"*",tau,")))+",
                B,"/",alpha2,"*((",1,"-exp(-",alpha2,"*t))+exp(-",alpha2,"*",tau,")*(",1,"-exp(-",alpha2,"*",TInf,"))*exp(-",alpha2,"*(t-",TInf,"))/(",1,"-exp(-",alpha2,"*",tau,")))+",
                C,"/",alpha3,"*((",1,"-exp(-",alpha3,"*t))+exp(-",alpha3,"*",tau,")*(",1,"-exp(-",alpha3,"*",TInf,"))*exp(-",alpha3,"*(t-",TInf,"))/(",1,"-exp(-",alpha3,"*",tau,"))))")
  form2N<-paste(R,"*(",A,"/",alpha1,"*((1-exp(-",alpha1,"*",TInf,"))*exp(-",alpha1,"*(t-",TInf,"))/(1-exp(-",alpha1,"*",tau,")))+",B,"/",alpha2,"*((1-exp(-",alpha2,"*",TInf,"))*exp(-",alpha2,"*(t-",TInf,"))/(1-exp(-",alpha2,"*",tau,")))+",C,"/",alpha3,"*((1-exp(-",alpha3,"*",TInf,"))*exp(-",alpha3,"*(t-",TInf,"))/(1-exp(-",alpha3,"*",tau,"))))")
  
  form1N<-parse(text=form1N,n=-1)
  form2N<-parse(text=form2N,n=-1)
  tf<-c(TInf,Inf)
  return(list(c(form1N,form2N),tf))
  
}

#-------------------------
# Oral1
#-------------------------

#---------
# single dose
#---------
#$Model definition
#$ka,V,k,k12,k21,k13,k31
#$
oral1_3cpt_kaVkk12k21k13k31<-function()
{
  prevalue<-myprefuncoral_3cpt(ka="ka",k10="k",k12="k12",k21="k21",k13="k13",k31="k31",V="V")
  A<-prevalue[[4]]
  B<-prevalue[[5]]
  C<-prevalue[[6]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  alpha3<-prevalue[[3]]
  
  form1<-paste("dose*(",A,"*exp(-",alpha1,"*t)+",B,"*exp(-",alpha2,"*t)+",C,"*exp(-",alpha3,"*t)-(",A,"+",B,"+",C,")*exp(-ka*t))")
  form1<-parse(text=form1,n=-1)
  tf<-Inf
  return(list(c(form1),tf))
}

#$Model definition
#$ka,Cl,V1,Q2,V2,Q3,V3
#$
oral1_3cpt_kaClV1Q2V2Q3V3<-function()
{
  prevalue<-myprefuncoral_3cpt(ka="ka",k10="Cl/V1",k12="Q2/V1",k21="Q2/V2",k13="Q3/V1",k31="Q3/V3",V="V1")
  A<-prevalue[[4]]
  B<-prevalue[[5]]
  C<-prevalue[[6]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  alpha3<-prevalue[[3]]
  
  form1<-paste("dose*(",A,"*exp(-",alpha1,"*t)+",B,"*exp(-",alpha2,"*t)+",C,"*exp(-",alpha3,"*t)-(",A,"+",B,"+",C,")*exp(-ka*t))")
  form1<-parse(text=form1,n=-1)
  tf<-Inf
  return(list(c(form1),tf))
}


#---------
# multiples doses
#---------
#$Model definition
#$ka,V,k,k12,k21,k13,k31
#$N,tau
oral1_3cpt_kaVkk12k21k13k31_md<-function(N,tau)
{
  prevalue<-myprefuncoral_3cpt(ka="ka",k10="k",k12="k12",k21="k21",k13="k13",k31="k31",V="V")
  A<-prevalue[[4]]
  B<-prevalue[[5]]
  C<-prevalue[[6]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  alpha3<-prevalue[[3]]
  multiterm1<-paste("((1-exp(-",N,"*",alpha1,"*",tau,"))/(1-exp(-",alpha1,"*",tau,")))")
  multiterm2<-paste("((1-exp(-",N,"*",alpha2,"*",tau,"))/(1-exp(-",alpha2,"*",tau,")))")
  multiterm3<-paste("((1-exp(-",N,"*",alpha3,"*",tau,"))/(1-exp(-",alpha3,"*",tau,")))")
  multiterm4<-paste("((1-exp(-",N,"*ka*",tau,"))/(1-exp(-ka*",tau,")))")
  formN<-paste("dose*(",A,"*exp(-",alpha1,"*(t-(",N,"-1)*",tau,"))*",multiterm1,"+",B,"*exp(-",alpha2,"*(t-(",N,"-1)*",tau,"))*",multiterm2,"+",C,"*exp(-",alpha3,"*(t-(",N,"-1)*",tau,"))*",multiterm3,"-(",A,"+",B,"+",C,")*exp(-ka*(t-(",N,"-1)*",tau,"))*",multiterm4,")")
  formN<-parse(text=formN,n=-1)
  tf<-Inf
  return(list(c(formN),tf))
}

#$Model definition
#$ka,Cl,V1,Q2,V2,Q3,V3
#$
oral1_3cpt_kaClV1Q2V2Q3V3_md<-function(N,tau)
{
  prevalue<-myprefuncoral_3cpt(ka="ka",k10="Cl/V1",k12="Q2/V1",k21="Q2/V2",k13="Q3/V1",k31="Q3/V3",V="V1")
  A<-prevalue[[4]]
  B<-prevalue[[5]]
  C<-prevalue[[6]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  alpha3<-prevalue[[3]]
  multiterm1<-paste("((1-exp(-",N,"*",alpha1,"*",tau,"))/(1-exp(-",alpha1,"*",tau,")))")
  multiterm2<-paste("((1-exp(-",N,"*",alpha2,"*",tau,"))/(1-exp(-",alpha2,"*",tau,")))")
  multiterm3<-paste("((1-exp(-",N,"*",alpha3,"*",tau,"))/(1-exp(-",alpha3,"*",tau,")))")
  multiterm4<-paste("((1-exp(-",N,"*ka*",tau,"))/(1-exp(-ka*",tau,")))")
  formN<-paste("dose*(",A,"*exp(-",alpha1,"*(t-(",N,"-1)*",tau,"))*",multiterm1,"+",B,"*exp(-",alpha2,"*(t-(",N,"-1)*",tau,"))*",multiterm2,"+",C,"*exp(-",alpha3,"*(t-(",N,"-1)*",tau,"))*",multiterm3,"-(",A,"+",B,"+",C,")*exp(-ka*(t-(",N,"-1)*",tau,"))*",multiterm4,")")
  formN<-parse(text=formN,n=-1)
  tf<-Inf
  return(list(c(formN),tf))
}
#---------
# steady-state
#---------

#-----------------------

#$Model definition
#$ka,V,k,k12,k21,k13,k31
#$tau
oral1_3cpt_kaVkk12k21k13k31_ss<-function(tau,TimeSS=0)
{
  prevalue<-myprefuncoral_3cpt(ka="ka",k10="k",k12="k12",k21="k21",k13="k13",k31="k31",V="V")
  A<-prevalue[[4]]
  B<-prevalue[[5]]
  C<-prevalue[[6]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  alpha3<-prevalue[[3]]
  
  multiterm1<-paste("(1/(1-exp(-",alpha1,"*",tau,")))")
  multiterm2<-paste("(1/(1-exp(-",alpha2,"*",tau,")))")
  multiterm3<-paste("(1/(1-exp(-",alpha3,"*",tau,")))")
  multiterm4<-paste("(1/(1-exp(-ka*",tau,")))")
  
  formSS<-paste("dose*(",A,"*exp(-",alpha1,"*(t-",TimeSS,"))*",multiterm1,"+",B,"*exp(-",alpha2,"*(t-",TimeSS,"))*",multiterm2,"+",C,"*exp(-",alpha3,"*(t-",TimeSS,"))*",multiterm3,"-(",A,"+",B,"+",C,")*exp(-ka*(t-",TimeSS,"))*",multiterm4,")" )
  formSS<-parse(text=formSS,n=-1)
  tf<-Inf
  return(list(c(formSS),tf))
}

#$Model definition
#$ka,Cl,V1,Q2,V2,Q3,V3
#$tau
oral1_3cpt_kaClV1Q2V2Q3V3_ss<-function(tau,TimeSS=0)
{
  prevalue<-myprefuncoral_3cpt(ka="ka",k10="Cl/V1",k12="Q2/V1",k21="Q2/V2",k13="Q3/V1",k31="Q3/V3",V="V1")
  A<-prevalue[[4]]
  B<-prevalue[[5]]
  C<-prevalue[[6]]
  alpha1<-prevalue[[1]]
  alpha2<-prevalue[[2]]
  alpha3<-prevalue[[3]]
  
  multiterm1<-paste("(1/(1-exp(-",alpha1,"*",tau,")))")
  multiterm2<-paste("(1/(1-exp(-",alpha2,"*",tau,")))")
  multiterm3<-paste("(1/(1-exp(-",alpha3,"*",tau,")))")
  multiterm4<-paste("(1/(1-exp(-ka*",tau,")))")
  
  formSS<-paste("dose*(",A,"*exp(-",alpha1,"*(t-",TimeSS,"))*",multiterm1,"+",B,"*exp(-",alpha2,"*(t-",TimeSS,"))*",multiterm2,"+",C,"*exp(-",alpha3,"*(t-",TimeSS,"))*",multiterm3,"-(",A,"+",B,"+",C,")*exp(-ka*(t-",TimeSS,"))*",multiterm4,")" )
  formSS<-parse(text=formSS,n=-1)
  tf<-Inf
  return(list(c(formSS),tf))
}



#-------------------########################################--------------------
#-------------------########################################--------------------
#-------------------##########                  ############--------------------
#-------------------########## Mickaelis Menten ############--------------------
#-------------------##########                  ############--------------------
#-------------------########################################--------------------
#-------------------########################################--------------------

#-------------------------
# Bolus
#-------------------------

#---------
# single dose
#---------
#$Model definition
#$V,k12,k21,k13,k31,Vm,km
bolus_3cpt_Vk12k21k13k31Vmkm<-function(){
  INT_bolus_3cpt_Vk12k21k13k31Vmkm<-function(t,y,p){
    V<-p[1]
    k12<-p[2]
    k21<-p[3]
    k13<-p[4]
    k31<-p[5]
    Vm <-p[6]
    km<-p[7]
    yd1<-(-Vm)*y[1]/(km+y[1])-k12*y[1]+k21*y[2]-k13*y[1]+k31*y[3]
    yd2<-k12*y[1]-k21*y[2]
    yd3<-k13*y[1]-k31*y[3]
    return(list(c(yd1,yd2,yd3),c(y[1]/V)))
  }
  return(INT_bolus_3cpt_Vk12k21k13k31Vmkm)
}
#$Model definition
#$V1,Q2,V2,Q3,V3,Vm,km
bolus_3cpt_V1Q2V2Q3V3Vmkm<-function(){
  INT_bolus_3cpt_V1Q2V2Q3V3Vmkm<-function(t,y,p){
    V1<-p[1]
    Q2<-p[2]
    V2<-p[3]
    Q3<-p[4]
    V3<-p[5]
    Vm <-p[6]
    km<-p[7]
    k12<-Q2/V1
    k21<-Q2/V2
    k13<-Q3/V1
    k31<-Q3/V3
    yd1<-(-Vm)*y[1]/(km+y[1])-k12*y[1]+k21*y[2]-k13*y[1]+k31*y[3]
    yd2<-k12*y[1]-k21*y[2]
    yd3<-k13*y[1]-k31*y[3]
    return(list(c(yd1,yd2,yd3),c(y[1]/V1)))
  }
  return(INT_bolus_3cpt_V1Q2V2Q3V3Vmkm)
}

#-------------------------
# Infusion
#-------------------------

#---------
# single dose
#---------
#$Model definition
#$V,k12,k21,Vm,km
#$ dose, TInf
infusion_3cpt_Vk12k21k13k31Vmkm<-function(doseMM,TInf){
  INT_infusion_3cpt_Vk12k21k13k31Vmkm<-function(t,y, p){
    V<-p[1]
    k12<-p[2]
    k21<-p[3]
    k13<-p[4]
    k31<-p[5]
    Vm <-p[6]
    km<-p[7]
    if(t<=TInf){
      yd1<-(-Vm)*y[1]/(km+y[1])+(doseMM/(TInf*V))-k12*y[1]+k12*y[2]-k13*y[1]+k13*y[3]}
    else{
      yd1<-(-Vm)*y[1]/(km+y[1])-k12*y[1]+k12*y[2]-k13*y[1]+k13*y[3]}
    yd2<-k21*y[1]-k21*y[2]
    yd3<-k31*y[1]-k31*y[3]
    return(list(c(yd1,yd2,yd3),c(y[1])))
  }
  return(INT_infusion_3cpt_Vk12k21k13k31Vmkm)
}

#$Model definition
#$V1,Q,V2,Vm,km
#$ dose, TInf
infusion_3cpt_V1Q2V2Q3V3Vmkm<-function(doseMM,TInf){
  INT_infusion_3cpt_V1Q2V2Q3V3Vmkm<-function(t,y, p){
    V1<-p[1]
    Q2<-p[2]
    V2<-p[3]
    Q3<-p[4]
    V3<-p[5]
    Vm <-p[6]
    km<-p[7]
    k12<-Q2/V1
    k21<-Q2/V2
    k13<-Q3/V1
    k31<-Q3/V3
    if(t<=TInf){
      yd1<-(-Vm)*y[1]/(km+y[1])+(doseMM/(TInf*V1))-k12*y[1]+k12*y[2]-k13*y[1]+k13*y[3]}
    else{
      yd1<-(-Vm)*y[1]/(km+y[1])-k12*y[1]+k12*y[2]-k13*y[1]+k13*y[3]}
    yd2<-k21*y[1]-k21*y[2]
    yd3<-k31*y[1]-k31*y[3]
    return(list(c(yd1,yd2,yd3),c(y[1])))
  }
  return(INT_infusion_3cpt_V1Q2V2Q3V3Vmkm)
}

#---------
# multiples doses
#---------
#$Model definition
#$V,k12,k21,Vm,km
#$ dose, TInf, tau
infusion_3cpt_Vk12k21k13k31Vmkm_md<-function(doseMM,TInf,tau){
  INT_infusion_3cpt_Vk12k21k13k31Vmkm_md<-function(t,y, p){
    V<-p[1]
    k12<-p[2]
    k21<-p[3]
    k13<-p[4]
    k31<-p[5]
    Vm <-p[6]
    km<-p[7]
    if(t%%tau<=TInf){
      yd1<-(-Vm)*y[1]/(km+y[1])+(doseMM/(TInf*V))-k12*y[1]+k12*y[2]-k13*y[1]+k13*y[3]}
    else{
      yd1<-(-Vm)*y[1]/(km+y[1])-k12*y[1]+k12*y[2]-k13*y[1]+k13*y[3]}
    yd2<-k21*y[1]-k21*y[2]
    yd3<-k31*y[1]-k31*y[3]
    return(list(c(yd1,yd2,yd3),c(y[1])))
  }
  return(INT_infusion_3cpt_Vk12k21k13k31Vmkm_md)
}

#$Model definition
#$V1,Q,V2,Vm,km
#$ dose, TInf, tau
infusion_3cpt_V1Q2V2Q3V3Vmkm_md<-function(doseMM,TInf,tau){
  INT_infusion_3cpt_V1Q2V2Q3V3Vmkm_md<-function(t,y, p){
    V1<-p[1]
    Q2<-p[2]
    V2<-p[3]
    Q3<-p[4]
    V3<-p[5]
    Vm <-p[6]
    km<-p[7]
    k12<-Q2/V1
    k21<-Q2/V2
    k13<-Q3/V1
    k31<-Q3/V3
    if(t%%tau<=TInf){
      yd1<-(-Vm)*y[1]/(km+y[1])+(doseMM/(TInf*V1))-k12*y[1]+k12*y[2]-k13*y[1]+k13*y[3]}
    else{
      yd1<-(-Vm)*y[1]/(km+y[1])-k12*y[1]+k12*y[2]-k13*y[1]+k13*y[3]}
    yd2<-k21*y[1]-k21*y[2]
    yd3<-k31*y[1]-k31*y[3]
    return(list(c(yd1,yd2,yd3),c(y[1])))
  }
  return(INT_infusion_3cpt_V1Q2V2Q3V3Vmkm_md)
}


#-------------------------
# Oral1
#-------------------------

#---------
# single dose
#---------
#$Model definition
#$ka,V,k12,k21,Vm,km
#$ dose
oral1_3cpt_kaVk12k21k13k31Vmkm<-function(doseMM){
  INT_oral1_3cpt_kaVk12k21k13k31Vmkm<-function(t,y,p){
    ka<-p[1]
    V<-p[2]
    k12<-p[3]
    k21<-p[4]
    k13<-p[5]
    k31<-p[6]
    Vm <-p[7]
    km<-p[8]
    yd1<-(-Vm)*y[1]/(km+y[1])+(doseMM*ka/V)*exp(-ka*t)-k12*y[1]+k12*y[2]-k13*y[1]+k13*y[3]
    yd2<-k21*y[1]-k21*y[2]
    yd3<-k31*y[1]-k31*y[3]
    return(list(c(yd1,yd2,yd3),c(y[1])))
  }
  return(INT_oral1_3cpt_kaVk12k21k13k31Vmkm)
}

#$Model definition
#$ka,V1,Q,V2,Vm,km
#$ dose
oral1_3cpt_kaV1Q2V2Q3V3Vmkm<-function(doseMM){
  INT_oral1_3cpt_kaV1Q2V2Q3V3Vmkm<-function(t,y,p){
    ka<-p[1]
    V1<-p[2]
    Q2<-p[3]
    V2<-p[4]
    Q3<-p[5]
    V3<-p[6]
    Vm <-p[7]
    km<-p[8]
    k12<-Q2/V1
    k21<-Q2/V2
    k13<-Q3/V1
    k31<-Q3/V3
    yd1<-(-Vm)*y[1]/(km+y[1])+(doseMM*ka/V1)*exp(-ka*t)-k12*y[1]+k12*y[2]-k13*y[1]+k13*y[3]
    yd2<-k21*y[1]-k21*y[2]
    yd3<-k31*y[1]-k31*y[3]
    return(list(c(yd1,yd2,yd3),c(y[1])))
  }
  return(INT_oral1_3cpt_kaV1Q2V2Q3V3Vmkm)
}

#---------
# multiples doses
#---------
#---------
# Needed intermediate function
#---------
input_oral1<-function(ka,V,doseMM,n,tau,t){
  if(n<0){stop("error function input_oral1 n<=0")}
  else{ if(n==0){return(doseMM*ka/V*exp(-ka*t))}
    else{return(doseMM*ka/V*exp(-ka*(t-n*tau))+input_oral1(ka,V,doseMM,(n-1),tau,t))}
  } }

#$Model definition
#$ka,V,k12,k21,Vm,km
#$ dose, tau
oral1_3cpt_kaVk12k21k13k31Vmkm_md<-function(doseMM,tau){
  INT_oral1_3cpt_kaVk12k21k13k31Vmkm_md<-function(t,y,p)
  {
    ka<-p[1]
    V<-p[2]
    k12<-p[3]
    k21<-p[4]
    k13<-p[5]
    k31<-p[6]
    Vm <-p[7]
    km<-p[8]
    n<-t%/%tau
    input<-input_oral1(ka,V,doseMM,n,tau,t)
    yd1<-(-Vm)*y[1]/(km+y[1])+input-k12*y[1]+k12*y[2]-k13*y[1]+k13*y[3]
    yd2<-k21*y[1]-k21*y[2]
    yd3<-k31*y[1]-k31*y[3]
    return(list(c(yd1,yd2,yd3),c(y[1])))
  }
  return(INT_oral1_3cpt_kaVk12k21k13k31Vmkm_md)
}

#$Model definition
#$ka,V1,Q,V2,Vm,km
#$ dose,tau
oral1_3cpt_kaV1Q2V2Q3V3Vmkm_md<-function(doseMM,tau){
  INT_oral1_3cpt_kaV1Q2V2Q3V3Vmkm_md<-function(t,y,p)
  {
    ka<-p[1]
    V1<-p[2]
    Q2<-p[3]
    V2<-p[4]
    Q3<-p[5]
    V3<-p[6]
    Vm <-p[7]
    km<-p[8]
    k12<-Q2/V1
    k21<-Q2/V2
    k13<-Q3/V1
    k31<-Q3/V3
    n<-t%/%tau
    input<-input_oral1(ka,V1,doseMM,n,tau,t)
    yd1<-(-Vm)*y[1]/(km+y[1])+input-k12*y[1]+k12*y[2]-k13*y[1]+k13*y[3]
    yd2<-k21*y[1]-k21*y[2]
    yd3<-k31*y[1]-k31*y[3]
    return(list(c(yd1,yd2,yd3),c(y[1])))
  }
  return(INT_oral1_3cpt_kaV1Q2V2Q3V3Vmkm_md)
}


