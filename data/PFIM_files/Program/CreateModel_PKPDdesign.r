#PFIM 4.0 
#April 2014
#Copyright (C) PFIM 4.0 - Universite Paris Diderot and INSERM.
#------------------------------------------------------------------------------------  


#-----------------------------------------------------------------------------------
# function to create the file model_created.r containing an ODE PK/PD model 
#-----------------------------------------------------------------------------------
create_formED<-function(pk.name,pd.list,dose=NA,tau=NA, TInf=NA){
# beginning of the function formED
  fun_formED<-c("formED<-function(t,y,p){")

  pk.list<-get(paste(pk.name,".list",sep=""))
  #pd.list<-get(pd.name)
  pk.fun<-get(paste(pk.name,".fun",sep=""))


# adding PK and then PD parameters
  for(i in 1:length(pk.list$param)){
    fun_formED<-c(fun_formED,paste(pk.list$param[i],"<-p[",i,"]",sep="")) }
  
  for(i in (length(pk.list$param)+1):(length(pk.list$param)+length(pd.list$param))){
    fun_formED<-c(fun_formED,paste(pd.list$param[i-length(pk.list$param)],"<-p[",i,"]",sep="")) }

# intermediate parameters  
  fun_formED<-c(fun_formED,paste("pk<-y[1:",pk.list$node,"]",sep=""))
  fun_formED<-c(fun_formED,paste("pd<-y[",pk.list$node+1,":",pk.list$node+pd.list$node,"]",sep=""))

  if(as.numeric(regexpr("bolus",pk.name))==-1){
    fun_formED<-c(fun_formED,"conc<-y[1]")}
  else{
       if(as.numeric(regexpr("1cpt",pk.name))!=-1){
        fun_formED<-c(fun_formED,"conc<-y[1]/V")}
       else{           
         if(as.numeric(regexpr("V1",pk.name))!=-1){
          fun_formED<-c(fun_formED,"conc<-y[1]/V1")}
         else{
          fun_formED<-c(fun_formED,"conc<-y[1]/V")}
        }
  }
# adding the PK and then PD ODE  
  fun_formED<-c(fun_formED,pk.fun(dose,tau,TInf))
  fun_formED<-c(fun_formED,pd.list$fun)

# end of the function formED
  fun_formED<-c(fun_formED,paste("return(list(c(",pk.list$ret1,",",pd.list$ret1,")",",c(",pk.list$ret2,",",pd.list$ret2,")))",sep=""))
  fun_formED<-c(fun_formED,"}")
# write of the created function
  cat(fun_formED,file=file.path(directory,"model_created.r"),sep="\n")
}


#-----------------------------------------------------------------------------------
# PK model to use in the function create_formED: LINEAR ELIMINATION
#-----------------------------------------------------------------------------------
#-------------------------
# Bolus
#-------------------------
bolus_1cpt_Vk.list<-list(param=c("V","k"),ret1="dpk1",ret2="pk[1]",node=1)
bolus_1cpt_VCl.list<-list(param=c("V","Cl"),ret1="dpk1",ret2="pk[1]",node=1)
bolus_2cpt_Vkk12k21.list<-list(param=c("V","k","k12","k21"),ret1="dpk1,dpk2",ret2="pk[1]/V",node=2)
bolus_2cpt_ClV1QV2.list<-list(param=c("Cl","V1","Q","V2"),ret1="dpk1,dpk2",ret2="pk[1]/V1",node=2)
bolus_3cpt_Vkk12k21k13k31.list<-list(param=c("V","k","k12","k21","k13","k31"),ret1="dpk1,dpk2,dpk3",ret2="pk[1]/V",node=3)
bolus_3cpt_ClV1Q2V2Q3V3.list<-list(param=c("Cl","V1","Q2","V2","Q3","V3"),ret1="dpk1,dpk2,dpk3",ret2="pk[1]/V1",node=3)

bolus_1cpt_Vk.fun<-function(dose,tau,TInf){return(c("dpk1<--k*pk[1]"))}
bolus_1cpt_VCl.fun<-function(dose,tau,TInf){return(c("dpk1<--(Cl/V)*pk[1]"))}
bolus_2cpt_Vkk12k21.fun<-function(dose,tau,TInf){return(c("dpk1<--k*pk[1]-k12*pk[1]+k21*pk[2]","dpk2<-k12*pk[1]-k21*pk[2]"))}
bolus_2cpt_ClV1QV2.fun<-function(dose,tau,TInf){return(c("k12<-Q/V1","k21<-Q/V2","dpk1<--(Cl/V1)*pk[1]-k12*pk[1]+k21*pk[2]","dpk2<-k12*pk[1]-k21*pk[2]"))}
bolus_3cpt_Vkk12k21k13k31.fun<-function(dose,tau,TInf){return(c("dpk1<--k*pk[1]-k12*pk[1]+k21*pk[2]-k13*pk[1]+k31*pk[3]","dpk2<-k12*pk[1]-k21*pk[2]","dpk3<-k13*pk[1]-k31*pk[3]"))}
bolus_3cpt_ClV1Q2V2Q3V3.fun<-function(dose,tau,TInf){return(c("k12<-Q/V1","k21<-Q/V2","dpk1<--(Cl/V1)*pk[1]-k12*pk[1]+k21*pk[2]-k13*pk[1]+k31*pk[3]","dpk2<-k12*pk[1]-k21*pk[2]","dpk3<-k13*pk[1]-k31*pk[3]"))}

#-------------------------
# infusion
#-------------------------
infusion_1cpt_Vk.list<-list(param=c("V","k"),ret1="dpk1",ret2="pk[1]",node=1)
infusion_1cpt_VCl.list<-list(param=c("V","Cl"),ret1="dpk1",ret2="pk[1]",node=1)
infusion_1cpt_Vk_md.list<-list(param=c("V","k"),ret1="dpk1",ret2="pk[1]",node=1)
infusion_1cpt_VCl_md.list<-list(param=c("V","Cl"),ret1="dpk1",ret2="pk[1]",node=1)
infusion_2cpt_Vkk12k21.list<-list(param=c("V","k","k12","k21"),ret1="dpk1,dpk2",ret2="pk[1]",node=2)
infusion_2cpt_ClV1QV2.list<-list(param=c("Cl","V1","Q","V2"),ret1="dpk1,dpk2",ret2="pk[1]",node=2)
infusion_2cpt_Vkk12k21_md.list<-list(param=c("V","k","k12","k21"),ret1="dpk1,dpk2",ret2="pk[1]",node=2)
infusion_2cpt_ClV1QV2_md.list<-list(param=c("Cl","V1","Q","V2"),ret1="dpk1,dpk2",ret2="pk[1]",node=2)
infusion_3cpt_Vkk12k21k13k31.list<-list(param=c("V","k","k12","k21","k13","k31"),ret1="dpk1,dpk2,dpk3",ret2="pk[1]",node=3)
infusion_3cpt_ClV1Q2V2Q3V3.list<-list(param=c("Cl","V1","Q","V2","Q3","V3"),ret1="dpk1,dpk2,dpk3",ret2="pk[1]",node=3)
infusion_3cpt_Vkk12k21k13k31_md.list<-list(param=c("V","k","k12","k21","k13","k31"),ret1="dpk1,dpk2,dpk3",ret2="pk[1]",node=3)
infusion_3cpt_ClV1Q2V2Q3V3_md.list<-list(param=c("Cl","V1","Q","V2","Q3","V3"),ret1="dpk1,dpk2,dpk3",ret2="pk[1]",node=3)

infusion_1cpt_Vk.fun<-function(dose,tau,TInf){return(c(paste("if(t<=",TInf,"){",sep=""),paste("dpk1<-(",dose,"/(",TInf,"*V))-k*pk[1]}",sep=""),"else{","dpk1<--k*pk[1]}"))}
infusion_1cpt_VCl.fun<-function(dose,tau,TInf){return(c("k<-Cl/V",paste("if(t<=",TInf,"){",sep=""),paste("dpk1<-(",dose,"/(",TInf,"*V))-k*pk[1]}",sep=""),"else{","dpk1<--k*pk[1]}"))}
infusion_1cpt_Vk_md.fun<-function(dose,tau,TInf){return(c(paste("if(t%%",tau,"<=",TInf,"){",sep=""),paste("dpk1<-(",dose,"/(",TInf,"*V))-k*pk[1]}",sep=""),"else{","dpk1<--k*pk[1]}"))}
infusion_1cpt_VCl_md.fun<-function(dose,tau,TInf){return(c("k<-Cl/V",paste("if(t%%",tau,"<=",TInf,"){",sep=""),paste("dpk1<-(",dose,"/(",TInf,"*V))-k*pk[1]}",sep=""),"else{","dpk1<--k*pk[1]}"))}
infusion_2cpt_Vkk12k21.fun<-function(dose,tau,TInf){return(c(paste("if(t<=",TInf,"){",sep=""),paste("dpk1<-(",dose,"/(",TInf,"*V))-k*pk[1]-k12*pk[1]+k12*pk[2]}",sep=""),"else{","dpk1<--k*pk[1]-k12*pk[1]+k12*pk[2]}","dpk2<-k21*pk[1]-k21*pk[2]"))}
infusion_2cpt_ClV1QV2.fun<-function(dose,tau,TInf){return(c("k<-Cl/V1","k12<-Q/V1","k21<-Q/V2",paste("if(t<=",TInf,"){",sep=""),paste("dpk1<-(",dose,"/(",TInf,"*V))-k*pk[1]-k12*pk[1]+k12*pk[2]}",sep=""),"else{","dpk1<--k*pk[1]-k12*pk[1]+k12*pk[2]}","dpk2<-k21*pk[1]-k21*pk[2]"))}
infusion_2cpt_Vkk12k21_md.fun<-function(dose,tau,TInf){return(c(paste("if(t%%",tau,"<=",TInf,"){",sep=""),paste("dpk1<-(",dose,"/(",TInf,"*V))-k*pk[1]-k12*pk[1]+k12*pk[2]}",sep=""),"else{","dpk1<--k*pk[1]-k12*pk[1]+k12*pk[2]}","dpk2<-k21*pk[1]-k21*pk[2]"))}
infusion_2cpt_ClV1QV2_md.fun<-function(dose,tau,TInf){return(c("k<-Cl/V1","k12<-Q/V1","k21<-Q/V2",paste("if(t%%",tau,"<=",TInf,"){",sep=""),paste("dpk1<-(",dose,"/(",TInf,"*V))-k*pk[1]-k12*pk[1]+k12*pk[2]}",sep=""),"else{","dpk1<--k*pk[1]-k12*pk[1]+k12*pk[2]}","dpk2<-k21*pk[1]-k21*pk[2]"))}
infusion_3cpt_Vkk12k21k13k31.fun<-function(dose,tau,TInf){return(c(paste("if(t<=",TInf,"){",sep=""),paste("dpk1<-(",dose,"/(",TInf,"*V))-k*pk[1]-k12*pk[1]+k12*pk[2]-k13*pk[1]+k13*pk[3]}",sep=""),"else{","dpk1<--k*pk[1]-k12*pk[1]+k12*pk[2]-k13*pk[1]+k13*pk[3]}","dpk2<-k21*pk[1]-k21*pk[2]","dpk3<-k31*pk[1]-k31*pk[3]"))}
infusion_3cpt_ClV1Q2V2Q3V3.fun<-function(dose,tau,TInf){return(c("k<-Cl/V1","k12<-Q/V1","k21<-Q/V2","k13<-Q3/V1","k31<-Q3/V3",paste("if(t<=",TInf,"){",sep=""),paste("dpk1<-(",dose,"/(",TInf,"*V))-k*pk[1]-k12*pk[1]+k12*pk[2]-k13*pk[1]+k13*pk[3]}",sep=""),"else{","dpk1<--k*pk[1]-k12*pk[1]+k12*pk[2]-k13*pk[1]+k13*pk[3]}","dpk2<-k21*pk[1]-k21*pk[2]","dpk3<-k31*pk[1]-k31*pk[3]"))}
infusion_3cpt_Vkk12k21k13k31_md.fun<-function(dose,tau,TInf){return(c(paste("if(t%%",tau,"<=",TInf,"){",sep=""),paste("dpk1<-(",dose,"/(",TInf,"*V))-k*pk[1]-k12*pk[1]+k12*pk[2]-k13*pk[1]+k13*pk[3]}",sep=""),"else{","dpk1<--k*pk[1]-k12*pk[1]+k12*pk[2]-k13*pk[1]+k13*pk[3]}","dpk2<-k21*pk[1]-k21*pk[2]","dpk3<-k31*pk[1]-k31*pk[3]"))}
infusion_3cpt_ClV1Q2V2Q3V3_md.fun<-function(dose,tau,TInf){return(c("k<-Cl/V1","k12<-Q/V1","k21<-Q/V2","k13<-Q3/V1","k31<-Q3/V3",paste("if(t%%",tau,"<=",TInf,"){",sep=""),paste("dpk1<-(",dose,"/(",TInf,"*V))-k*pk[1]-k12*pk[1]+k12*pk[2]-k13*pk[1]+k13*pk[3]}",sep=""),"else{","dpk1<--k*pk[1]-k12*pk[1]+k12*pk[2]-k13*pk[1]+k13*pk[3]}","dpk2<-k21*pk[1]-k21*pk[2]","dpk3<-k31*pk[1]-k31*pk[3]"))}

#-------------------------
# Oral1
#-------------------------
oral1_1cpt_kaVk.list<-list(param=c("ka","V","k"),ret1="dpk1",ret2="pk[1]",node=1)
oral1_1cpt_kaVCl.list<-list(param=c("ka","V","Cl"),ret1="dpk1",ret2="pk[1]",node=1)
oral1_1cpt_kaVk_md.list<-list(param=c("ka","V","k"),ret1="dpk1",ret2="pk[1]",node=1)
oral1_1cpt_kaVCl_md.list<-list(param=c("ka","V","Cl"),ret1="dpk1",ret2="pk[1]",node=1)
oral1_2cpt_kaVkk12k21.list<-list(param=c("ka","V","k","k12","k21"),ret1="dpk1,dpk2",ret2="pk[1]",node=2)
oral1_2cpt_kaClV1QV2.list<-list(param=c("ka","Cl","V1","Q","V2"),ret1="dpk1,dpk2",ret2="pk[1]",node=2)
oral1_2cpt_kaVkk12k21_md.list<-list(param=c("ka","V","k","k12","k21"),ret1="dpk1,dpk2",ret2="pk[1]",node=2)
oral1_2cpt_kaClV1QV2_md.list<-list(param=c("ka","Cl","V1","Q","V2"),ret1="dpk1,dpk2",ret2="pk[1]",node=2)
oral1_3cpt_kaVkk12k21k13k31.list<-list(param=c("ka","V","k","k12","k21","k13","k31"),ret1="dpk1,dpk2,dpk3",ret2="pk[1]",node=3)
oral1_3cpt_kaClV1Q2V2Q3V3.list<-list(param=c("ka","Cl","V1","Q2","V2","Q3","V3"),ret1="dpk1,dpk2,dpk3",ret2="pk[1]",node=3)
oral1_3cpt_kaVkk12k21k13k31_md.list<-list(param=c("ka","V","k","k12","k21","k13","k31"),ret1="dpk1,dpk2,dpk3",ret2="pk[1]",node=3)
oral1_3cpt_kaClV1Q2V2Q3V3_md.list<-list(param=c("ka","Cl","V1","Q2","V2","Q3","V3"),ret1="dpk1,dpk2,dpk3",ret2="pk[1]",node=3)

oral1_1cpt_kaVk.fun<-function(dose,tau,TInf){return(c(paste("dpk1<--k*pk[1]+(",dose,"*ka/V)*exp(-ka*t)",sep="")))}
oral1_1cpt_kaVCl.fun<-function(dose,tau,TInf){return(c("k<-Cl/V",paste("dpk1<--k*pk[1]+(",dose,"*ka/V)*exp(-ka*t)",sep="")))}
oral1_1cpt_kaVk_md.fun<-function(dose,tau,TInf){return(c("input_oral1<-function(ka,V,dose,n,tau,t){","if(n<0){stop()} ","else{ if(n==0){return(dose*ka/V*exp(-ka*t))}","else{return(dose*ka/V*exp(-ka*(t-n*tau))+input_oral1(ka,V,dose,(n-1),tau,t))}}}",paste("n<-t%/%",tau,sep=""),paste("input<-input_oral1(ka,V,",dose,",n,",tau,",t)",sep=""),"dpk1<--k*pk[1]+input"))}
oral1_1cpt_kaVCl_md.fun<-function(dose,tau,TInf){return(c("k<-Cl/V","input_oral1<-function(ka,V,dose,n,tau,t){","if(n<0){stop()} ","else{ if(n==0){return(dose*ka/V*exp(-ka*t))}","else{return(dose*ka/V*exp(-ka*(t-n*tau))+input_oral1(ka,V,dose,(n-1),tau,t))}}}",paste("n<-t%/%",tau,sep=""),paste("input<-input_oral1(ka,V,",dose,",n,",tau,",t)",sep=""),"dpk1<--k*pk[1]+input"))}
oral1_2cpt_kaVkk12k21.fun<-function(dose,tau,TInf){return(c(paste("dpk1<--k*pk[1]+(",dose,"*ka/V)*exp(-ka*t)-k12*pk[1]+k12*pk[2]",sep=""),"dpk2<-k21*pk[1]-k21*pk[2]"))}
oral1_2cpt_kaClV1QV2.fun<-function(dose,tau,TInf){return(c("k<-Cl/V1","k12<-Q/V1","k21<-Q/V2",paste("dpk1<--k*pk[1]+(",dose,"*ka/V)*exp(-ka*t)-k12*pk[1]+k12*pk[2]",sep=""),"dpk2<-k21*pk[1]-k21*pk[2]")) }
oral1_2cpt_kaVkk12k21_md.fun<-function(dose,tau,TInf){return(c("input_oral1<-function(ka,V,dose,n,tau,t){","if(n<0){stop()} ","else{ if(n==0){return(dose*ka/V*exp(-ka*t))}","else{return(dose*ka/V*exp(-ka*(t-n*tau))+input_oral1(ka,V,dose,(n-1),tau,t))}}}",paste("n<-t%/%",tau,sep=""),paste("input<-input_oral1(ka,V,",dose,",n,",tau,",t)",sep=""),"dpk1<--k*pk[1]+input-k12*pk[1]+k12*pk[2]","dpk2<-k21*pk[1]-k21*pk[2]"))}
oral1_2cpt_kaClV1QV2_md.fun<-function(dose,tau,TInf){return(c("k<-Cl/V1","k12<-Q/V1","k21<-Q/V2","input_oral1<-function(ka,V1,dose,n,tau,t){","if(n<0){stop()} ","else{ if(n==0){return(dose*ka/V1*exp(-ka*t))}","else{return(dose*ka/V1*exp(-ka*(t-n*tau))+input_oral1(ka,V1,dose,(n-1),tau,t))}}}",paste("n<-t%/%",tau,sep=""),paste("input<-input_oral1(ka,V1,",dose,",n,",tau,",t)",sep=""),"dpk1<--k*pk[1]+input-k12*pk[1]+k12*pk[2]","dpk2<-k21*pk[1]-k21*pk[2]"))}
oral1_3cpt_kaVkk12k21k13k31.fun<-function(dose,tau,TInf){return(c(paste("dpk1<--k*pk[1]+(",dose,"*ka/V)*exp(-ka*t)-k12*pk[1]+k12*pk[2]-k13*pk[1]+k13*pk[3]",sep=""),"dpk2<-k21*pk[1]-k21*pk[2]","dpk3<-k31*pk[1]-k31*pk[3]"))}
oral1_3cpt_kaClV1Q2V2Q3V3.fun<-function(dose,tau,TInf){return(c("k<-Cl/V1","k12<-Q/V1","k21<-Q/V2","k13<-Q3/V1","k31<-Q3/V3",paste("dpk1<--k*pk[1]+(",dose,"*ka/V)*exp(-ka*t)-k12*pk[1]+k12*pk[2]-k13*pk[1]+k13*pk[3]",sep=""),"dpk2<-k21*pk[1]-k21*pk[2]","dpk3<-k31*pk[1]-k31*pk[3]"))}
oral1_3cpt_kaVkk12k21k13k31_md.fun<-function(dose,tau,TInf){return(c("input_oral1<-function(ka,V,dose,n,tau,t){","if(n<0){stop()} ","else{ if(n==0){return(dose*ka/V*exp(-ka*t))}","else{return(dose*ka/V*exp(-ka*(t-n*tau))+input_oral1(ka,V,dose,(n-1),tau,t))}}}",paste("n<-t%/%",tau,sep=""),paste("input<-input_oral1(ka,V,",dose,",n,",tau,",t)",sep=""),"dpk1<--k*pk[1]+input-k12*pk[1]+k12*pk[2]-k13*pk[1]+k13*pk[3]","dpk2<-k21*pk[1]-k21*pk[2]","dpk3<-k31*pk[1]-k31*pk[3]"))}
oral1_3cpt_kaClV1Q2V2Q3V3_md.fun<-function(dose,tau,TInf){return(c("k<-Cl/V1","k12<-Q2/V1","k21<-Q2/V2","k13<-Q3/V1","k31<-Q3/V3","input_oral1<-function(ka,V1,dose,n,tau,t){","if(n<0){stop()} ","else{ if(n==0){return(dose*ka/V1*exp(-ka*t))}","else{return(dose*ka/V1*exp(-ka*(t-n*tau))+input_oral1(ka,V1,dose,(n-1),tau,t))}}}",paste("n<-t%/%",tau,sep=""),paste("input<-input_oral1(ka,V1,",dose,",n,",tau,",t)",sep=""),"dpk1<--k*pk[1]+input-k12*pk[1]+k12*pk[2]-k13*pk[1]+k13*pk[3]","dpk2<-k21*pk[1]-k21*pk[2]","dpk3<-k31*pk[1]-k31*pk[3]"))}

#-----------------------------------------------------------------------------------
# PK model to use in the function create_formED: MICKAELIS MENTEN ELIMINATION
#-----------------------------------------------------------------------------------
#-------------------------
# Bolus
#-------------------------
bolus_1cpt_VVmkm.list<-list(param=c("V","Vm","km"),ret1="dpk1",ret2="pk[1]",node=1)
bolus_2cpt_Vk12k21Vmkm.list<-list(param=c("V","k12","k21","Vm","km"),ret1="dpk1,dpk2",ret2="pk[1]/V",node=2)
bolus_2cpt_V1QV2Vmkm.list<-list(param=c("V1","Q","V2","Vm","km"),ret1="dpk1,dpk2",ret2="pk[1]/V1",node=2)
bolus_3cpt_Vk12k21k13k31Vmkm.list<-list(param=c("V","k12","k21","k13","k31","Vm","km"),ret1="dpk1,dpk2,dpk3",ret2="pk[1]/V",node=3)
bolus_3cpt_V1Q2V2Q3V3Vmkm.list<-list(param=c("V1","Q2","V2","Q3","V3","Vm","km"),ret1="dpk1,dpk2,dpk3",ret2="pk[1]/V1",node=3)

bolus_1cpt_VVmkm.fun<-function(dose,tau,TInf){return(c("dpk1<-(-Vm)*pk[1]/(km+pk[1])"))}
bolus_2cpt_Vk12k21Vmkm.fun<-function(dose,tau,TInf){return(c("dpk1<-(-Vm)*pk[1]/(km+pk[1])-k12*pk[1]+k21*pk[2]","dpk2<-k12*pk[1]-k21*pk[2]"))}
bolus_2cpt_V1QV2Vmkm.fun<-function(dose,tau,TInf){return(c("k12<-Q/V1","k21<-Q/V2","dpk1<-(-Vm)*pk[1]/(km+pk[1])-k12*pk[1]+k21*pk[2]","dpk2<-k12*pk[1]-k21*pk[2]"))}
bolus_3cpt_Vk12k21k13k31Vmkm.fun<-function(dose,tau,TInf){return(c("dpk1<-(-Vm)*pk[1]/(km+pk[1])-k12*pk[1]+k21*pk[2]-k13*pk[1]+k31*pk[3]","dpk2<-k12*pk[1]-k21*pk[2]","dpk3<-k13*pk[1]-k31*pk[3]"))}
bolus_3cpt_V1Q2V2Q3V3Vmkm.fun<-function(dose,tau,TInf){return(c("k12<-Q/V1","k21<-Q/V2","dpk1<-(-Vm)*pk[1]/(km+pk[1])-k12*pk[1]+k21*pk[2]-k13*pk[1]+k31*pk[3]","dpk2<-k12*pk[1]-k21*pk[2]","dpk3<-k13*pk[1]-k31*pk[3]"))}

#-------------------------
# infusion
#-------------------------
infusion_1cpt_VVmkm.list<-list(param=c("V","Vm","km"),ret1="dpk1",ret2="pk[1]",node=1)
infusion_1cpt_VVmkm_md.list<-list(param=c("V","Vm","km"),ret1="dpk1",ret2="pk[1]",node=1)
infusion_2cpt_Vk12k21Vmkm.list<-list(param=c("V","k12","k21","Vm","km"),ret1="dpk1,dpk2",ret2="pk[1]",node=2)
infusion_2cpt_V1QV2Vmkm.list<-list(param=c("V1","Q","V2","Vm","km"),ret1="dpk1,dpk2",ret2="pk[1]",node=2)
infusion_2cpt_Vk12k21Vmkm_md.list<-list(param=c("V","k12","k21","Vm","km"),ret1="dpk1,dpk2",ret2="pk[1]",node=2)
infusion_2cpt_V1QV2Vmkm_md.list<-list(param=c("V1","Q","V2","Vm","km"),ret1="dpk1,dpk2",ret2="pk[1]",node=2)
infusion_3cpt_Vk12k21k13k31Vmkm.list<-list(param=c("V","k12","k21","k13","k31","Vm","km"),ret1="dpk1,dpk2,dpk3",ret2="pk[1]",node=3)
infusion_3cpt_V1Q2V2Q3V3Vmkm.list<-list(param=c("V1","Q","V2","Q3","V3","Vm","km"),ret1="dpk1,dpk2,dpk3",ret2="pk[1]",node=3)
infusion_3cpt_Vk12k21k13k31Vmkm_md.list<-list(param=c("V","k12","k21","k13","k31","Vm","km"),ret1="dpk1,dpk2,dpk3",ret2="pk[1]",node=3)
infusion_3cpt_V1Q2V2Q3V3Vmkm_md.list<-list(param=c("V1","Q","V2","Q3","V3","Vm","km"),ret1="dpk1,dpk2,dpk3",ret2="pk[1]",node=3)

infusion_1cpt_VVmkm.fun<-function(dose,tau,TInf){return(c(paste("if(t<=",TInf,"){",sep=""),paste("dpk1<-(",dose,"/(",TInf,"*V))+(-Vm)*pk[1]/(km+pk[1])}",sep=""),"else{","dpk1<-(-Vm)*pk[1]/(km+pk[1])}"))}
infusion_1cpt_VVmkm_md.fun<-function(dose,tau,TInf){return(c("k<-/V",paste("if(t%%",tau,"<=",TInf,"){",sep=""),paste("dpk1<-(",dose,"/(",TInf,"*V))+(-Vm)*pk[1]/(km+pk[1])}",sep=""),"else{","dpk1<-(-Vm)*pk[1]/(km+pk[1])}"))}
infusion_2cpt_Vk12k21Vmkm.fun<-function(dose,tau,TInf){return(c(paste("if(t<=",TInf,"){",sep=""),paste("dpk1<-(",dose,"/(",TInf,"*V))+(-Vm)*pk[1]/(km+pk[1])-k12*pk[1]+k12*pk[2]}",sep=""),"else{","dpk1<-(-Vm)*pk[1]/(km+pk[1])-k12*pk[1]+k12*pk[2]}","dpk2<-k21*pk[1]-k21*pk[2]"))}
infusion_2cpt_V1QV2Vmkm.fun<-function(dose,tau,TInf){return(c("k12<-Q/V1","k21<-Q/V2",paste("if(t<=",TInf,"){",sep=""),paste("dpk1<-(",dose,"/(",TInf,"*V))+(-Vm)*pk[1]/(km+pk[1])-k12*pk[1]+k12*pk[2]}",sep=""),"else{","dpk1<-(-Vm)*pk[1]/(km+pk[1])-k12*pk[1]+k12*pk[2]}","dpk2<-k21*pk[1]-k21*pk[2]"))}
infusion_2cpt_Vk12k21Vmkm_md.fun<-function(dose,tau,TInf){return(c(paste("if(t%%",tau,"<=",TInf,"){",sep=""),paste("dpk1<-(",dose,"/(",TInf,"*V))+(-Vm)*pk[1]/(km+pk[1])-k12*pk[1]+k12*pk[2]}",sep=""),"else{","dpk1<-(-Vm)*pk[1]/(km+pk[1])-k12*pk[1]+k12*pk[2]}","dpk2<-k21*pk[1]-k21*pk[2]"))}
infusion_2cpt_V1QV2Vmkm_md.fun<-function(dose,tau,TInf){return(c("k12<-Q/V1","k21<-Q/V2",paste("if(t%%",tau,"<=",TInf,"){",sep=""),paste("dpk1<-(",dose,"/(",TInf,"*V))+(-Vm)*pk[1]/(km+pk[1])-k12*pk[1]+k12*pk[2]}",sep=""),"else{","dpk1<-(-Vm)*pk[1]/(km+pk[1])-k12*pk[1]+k12*pk[2]}","dpk2<-k21*pk[1]-k21*pk[2]"))}
infusion_3cpt_Vk12k21k13k31Vmkm.fun<-function(dose,tau,TInf){return(c(paste("if(t<=",TInf,"){",sep=""),paste("dpk1<-(",dose,"/(",TInf,"*V))+(-Vm)*pk[1]/(km+pk[1])-k12*pk[1]+k12*pk[2]-k13*pk[1]+k13*pk[3]}",sep=""),"else{","dpk1<-(-Vm)*pk[1]/(km+pk[1])-k12*pk[1]+k12*pk[2]-k13*pk[1]+k13*pk[3]}","dpk2<-k21*pk[1]-k21*pk[2]","dpk3<-k31*pk[1]-k31*pk[3]"))}
infusion_3cpt_V1Q2V2Q3V3Vmkm.fun<-function(dose,tau,TInf){return(c("k12<-Q/V1","k21<-Q/V2","k13<-Q3/V1","k31<-Q3/V3",paste("if(t<=",TInf,"){",sep=""),paste("dpk1<-(",dose,"/(",TInf,"*V))+(-Vm)*pk[1]/(km+pk[1])-k12*pk[1]+k12*pk[2]-k13*pk[1]+k13*pk[3]}",sep=""),"else{","dpk1<-(-Vm)*pk[1]/(km+pk[1])-k12*pk[1]+k12*pk[2]-k13*pk[1]+k13*pk[3]}","dpk2<-k21*pk[1]-k21*pk[2]","dpk3<-k31*pk[1]-k31*pk[3]"))}
infusion_3cpt_Vk12k21k13k31Vmkm_md.fun<-function(dose,tau,TInf){return(c(paste("if(t%%",tau,"<=",TInf,"){",sep=""),paste("dpk1<-(",dose,"/(",TInf,"*V))+(-Vm)*pk[1]/(km+pk[1])-k12*pk[1]+k12*pk[2]-k13*pk[1]+k13*pk[3]}",sep=""),"else{","dpk1<-(-Vm)*pk[1]/(km+pk[1])-k12*pk[1]+k12*pk[2]-k13*pk[1]+k13*pk[3]}","dpk2<-k21*pk[1]-k21*pk[2]","dpk3<-k31*pk[1]-k31*pk[3]"))}
infusion_3cpt_V1Q2V2Q3V3Vmkm_md.fun<-function(dose,tau,TInf){return(c("k12<-Q/V1","k21<-Q/V2","k13<-Q3/V1","k31<-Q3/V3",paste("if(t%%",tau,"<=",TInf,"){",sep=""),paste("dpk1<-(",dose,"/(",TInf,"*V))+(-Vm)*pk[1]/(km+pk[1])-k12*pk[1]+k12*pk[2]-k13*pk[1]+k13*pk[3]}",sep=""),"else{","dpk1<-(-Vm)*pk[1]/(km+pk[1])-k12*pk[1]+k12*pk[2]-k13*pk[1]+k13*pk[3]}","dpk2<-k21*pk[1]-k21*pk[2]","dpk3<-k31*pk[1]-k31*pk[3]"))}

#-------------------------
# Oral1
#-------------------------
oral1_1cpt_kaVVmkm.list<-list(param=c("ka","V","Vm","km"),ret1="dpk1",ret2="pk[1]",node=1)
oral1_1cpt_kaVVmkm_md.list<-list(param=c("ka","V","Vm","km"),ret1="dpk1",ret2="pk[1]",node=1)
oral1_2cpt_kaVk12k21Vmkm.list<-list(param=c("ka","V","k12","k21","Vm","km"),ret1="dpk1,dpk2",ret2="pk[1]",node=2)
oral1_2cpt_kaV1QV2Vmkm.list<-list(param=c("ka","V1","Q","V2","Vm","km"),ret1="dpk1,dpk2",ret2="pk[1]",node=2)
oral1_2cpt_kaVk12k21Vmkm_md.list<-list(param=c("ka","V","k12","k21","Vm","km"),ret1="dpk1,dpk2",ret2="pk[1]",node=2)
oral1_2cpt_kaV1QV2Vmkm_md.list<-list(param=c("ka","V1","Q","V2","Vm","km"),ret1="dpk1,dpk2",ret2="pk[1]",node=2)
oral1_3cpt_kaVk12k21k13k31Vmkm.list<-list(param=c("ka","V","k12","k21","k13","k31","Vm","km"),ret1="dpk1,dpk2,dpk3",ret2="pk[1]",node=3)
oral1_3cpt_kaV1Q2V2Q3V3Vmkm.list<-list(param=c("ka","V1","Q2","V2","Q3","V3","Vm","km"),ret1="dpk1,dpk2,dpk3",ret2="pk[1]",node=3)
oral1_3cpt_kaVk12k21k13k31Vmkm_md.list<-list(param=c("ka","V","k12","k21","k13","k31","Vm","km"),ret1="dpk1,dpk2,dpk3",ret2="pk[1]",node=3)
oral1_3cpt_kaV1Q2V2Q3V3Vmkm_md.list<-list(param=c("ka","V1","Q2","V2","Q3","V3","Vm","km"),ret1="dpk1,dpk2,dpk3",ret2="pk[1]",node=3)

oral1_1cpt_kaVVmkm.fun<-function(dose,tau,TInf){return(c(paste("dpk1<-(-Vm)*pk[1]/(km+pk[1])+(",dose,"*ka/V)*exp(-ka*t)",sep="")))}
oral1_1cpt_kaVVmkm_md.fun<-function(dose,tau,TInf){return(c("input_oral1<-function(ka,V,dose,n,tau,t){","if(n<0){stop()} ","else{ if(n==0){return(dose*ka/V*exp(-ka*t))}","else{return(dose*ka/V*exp(-ka*(t-n*tau))+input_oral1(ka,V,dose,(n-1),tau,t))}}}",paste("n<-t%/%",tau,sep=""),paste("input<-input_oral1(ka,V,",dose,",n,",tau,",t)",sep=""),"dpk1<-(-Vm)*pk[1]/(km+pk[1])+input"))}
oral1_2cpt_kaVk12k21Vmkm.fun<-function(dose,tau,TInf){return(c(paste("dpk1<-(-Vm)*pk[1]/(km+pk[1])+(",dose,"*ka/V)*exp(-ka*t)-k12*pk[1]+k12*pk[2]",sep=""),"dpk2<-k21*pk[1]-k21*pk[2]"))}
oral1_2cpt_kaV1QV2Vmkm.fun<-function(dose,tau,TInf){return(c("k12<-Q/V1","k21<-Q/V2",paste("dpk1<-(-Vm)*pk[1]/(km+pk[1])+(",dose,"*ka/V)*exp(-ka*t)-k12*pk[1]+k12*pk[2]",sep=""),"dpk2<-k21*pk[1]-k21*pk[2]"))}
oral1_2cpt_kaVk12k21Vmkm_md.fun<-function(dose,tau,TInf){return(c("input_oral1<-function(ka,V,dose,n,tau,t){","if(n<0){stop()} ","else{ if(n==0){return(dose*ka/V*exp(-ka*t))}","else{return(dose*ka/V*exp(-ka*(t-n*tau))+input_oral1(ka,V,dose,(n-1),tau,t))}}}",paste("n<-t%/%",tau,sep=""),paste("input<-input_oral1(ka,V,",dose,",n,",tau,",t)",sep=""),"dpk1<-(-Vm)*pk[1]/(km+pk[1])+input-k12*pk[1]+k12*pk[2]","dpk2<-k21*pk[1]-k21*pk[2]"))}
oral1_2cpt_kaV1QV2Vmkm_md.fun<-function(dose,tau,TInf){return(c("k12<-Q/V1","k21<-Q/V2","input_oral1<-function(ka,V1,dose,n,tau,t){","if(n<0){stop()} ","else{ if(n==0){return(dose*ka/V1*exp(-ka*t))}","else{return(dose*ka/V1*exp(-ka*(t-n*tau))+input_oral1(ka,V1,dose,(n-1),tau,t))}}}",paste("n<-t%/%",tau,sep=""),paste("input<-input_oral1(ka,V1,",dose,",n,",tau,",t)",sep=""),"dpk1<-(-Vm)*pk[1]/(km+pk[1])+input-k12*pk[1]+k12*pk[2]","dpk2<-k21*pk[1]-k21*pk[2]"))}
oral1_3cpt_kaVk12k21k13k31Vmkm.fun<-function(dose,tau,TInf){return(c(paste("dpk1<-(-Vm)*pk[1]/(km+pk[1])+(",dose,"*ka/V)*exp(-ka*t)-k12*pk[1]+k12*pk[2]-k13*pk[1]+k13*pk[3]",sep=""),"dpk2<-k21*pk[1]-k21*pk[2]","dpk3<-k31*pk[1]-k31*pk[3]"))}
oral1_3cpt_kaV1Q2V2Q3V3Vmkm.fun<-function(dose,tau,TInf){return(c("k12<-Q/V1","k21<-Q/V2","k13<-Q3/V1","k31<-Q3/V3",paste("dpk1<-(-Vm)*pk[1]/(km+pk[1])+(",dose,"*ka/V)*exp(-ka*t)-k12*pk[1]+k12*pk[2]-k13*pk[1]+k13*pk[3]",sep=""),"dpk2<-k21*pk[1]-k21*pk[2]","dpk3<-k31*pk[1]-k31*pk[3]"))}
oral1_3cpt_kaVk12k21k13k31Vmkm_md.fun<-function(dose,tau,TInf){return(c("input_oral1<-function(ka,V,dose,n,tau,t){","if(n<0){stop()} ","else{ if(n==0){return(dose*ka/V*exp(-ka*t))}","else{return(dose*ka/V*exp(-ka*(t-n*tau))+input_oral1(ka,V,dose,(n-1),tau,t))}}}",paste("n<-t%/%",tau,sep=""),paste("input<-input_oral1(ka,V,",dose,",n,",tau,",t)",sep=""),"dpk1<-(-Vm)*pk[1]/(km+pk[1])+input-k12*pk[1]+k12*pk[2]-k13*pk[1]+k13*pk[3]","dpk2<-k21*pk[1]-k21*pk[2]","dpk3<-k31*pk[1]-k31*pk[3]"))}
oral1_3cpt_kaV1Q2V2Q3V3Vmkm_md.fun<-function(dose,tau,TInf){return(c("k12<-Q2/V1","k21<-Q2/V2","k13<-Q3/V1","k31<-Q3/V3","input_oral1<-function(ka,V1,dose,n,tau,t){","if(n<0){stop()} ","else{ if(n==0){return(dose*ka/V1*exp(-ka*t))}","else{return(dose*ka/V1*exp(-ka*(t-n*tau))+input_oral1(ka,V1,dose,(n-1),tau,t))}}}",paste("n<-t%/%",tau,sep=""),paste("input<-input_oral1(ka,V1,",dose,",n,",tau,",t)",sep=""),"dpk1<-(-Vm)*pk[1]/(km+pk[1])+input-k12*pk[1]+k12*pk[2]-k13*pk[1]+k13*pk[3]","dpk2<-k21*pk[1]-k21*pk[2]","dpk3<-k31*pk[1]-k31*pk[3]"))}



#-----------------------------------------------------------------------------------
# PD model to use in the function create_formED
#-----------------------------------------------------------------------------------
turn_input_Emax<-list(param=c("Rin","kout","Emax","C50"),ret1="dpd1",ret2="pd[1]",node=1,fun=c("dpd1<-Rin*(1+(Emax*conc)/(conc+C50))-kout*pd[1]"))
turn_input_gammaEmax<-list(param=c("Rin","kout","Emax","C50","gamma"),ret1="dpd1",ret2="pd[1]",node=1,fun=c("dpd1<-Rin*(1+(Emax*conc**gamma)/(conc**gamma+C50**gamma))-kout*pd[1]"))
turn_input_Imax<-list(param=c("Rin","kout","Imax","C50"),ret1="dpd1",ret2="pd[1]",node=1,fun=c("dpd1<-Rin*(1-(Imax*conc)/(conc+C50))-kout*pd[1]"))
turn_input_Imaxfull<-list(param=c("Rin","kout","C50"),ret1="dpd1",ret2="pd[1]",node=1,fun=c("dpd1<-Rin*(1-(conc)/(conc+C50))-kout*pd[1]"))
turn_input_gammaImax<-list(param=c("Rin","kout","Imax","C50","gamma"),ret1="dpd1",ret2="pd[1]",node=1,fun=c("dpd1<-Rin*(1-(Imax*conc**gamma)/(conc**gamma+C50**gamma))-kout*pd[1]"))
turn_input_gammaImaxfull<-list(param=c("Rin","kout","C50","gamma"),ret1="dpd1",ret2="pd[1]",node=1,fun=c("dpd1<-Rin*(1-(conc**gamma)/(conc**gamma+C50**gamma))-kout*pd[1]"))
turn_output_Emax<-list(param=c("Rin","kout","Emax","C50"),ret1="dpd1",ret2="pd[1]",node=1,fun=c("dpd1<-Rin-kout*pd[1]*(1+(Emax*conc)/(conc+C50))"))
turn_output_gammaEmax<-list(param=c("Rin","kout","Emax","C50","gamma"),ret1="dpd1",ret2="pd[1]",node=1,fun=c("dpd1<-Rin-kout*pd[1]*(1+(Emax*conc**gamma)/(conc**gamma+C50**gamma))"))
turn_output_Imax<-list(param=c("Rin","kout","Imax","C50"),ret1="dpd1",ret2="pd[1]",node=1,fun=c("dpd1<-Rin-kout*pd[1]*(1-(Imax*conc)/(conc+C50))"))
turn_output_Imaxfull<-list(param=c("Rin","kout","C50"),ret1="dpd1",ret2="pd[1]",node=1,fun=c("dpd1<-Rin-kout*pd[1]*(1-(conc)/(conc+C50))"))
turn_output_gammaImax<-list(param=c("Rin","kout","Imax","C50","gamma"),ret1="dpd1",ret2="pd[1]",node=1,fun=c("dpd1<-Rin-kout*pd[1]*(1-(Imax*conc**gamma)/(conc**gamma+C50**gamma))"))
turn_output_gammaImaxfull<-list(param=c("Rin","kout","C50","gamma"),ret1="dpd1",ret2="pd[1]",node=1,fun=c("dpd1<-Rin-kout*pd[1]*(1-(conc**gamma)/(conc**gamma+C50**gamma))"))

immed_lin_null<-list(param=c("Alin"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-Alin*conc"))
immed_lin_const<-list(param=c("Alin","S0"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0+Alin*conc"))
immed_lin_lin<-list(param=c("Alin","S0","kprog"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0+kprog*t+Alin*conc"))
immed_lin_exp<-list(param=c("Alin","S0","kprog"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0*exp(-kprog*t)+Alin*conc"))
immed_lin_dexp<-list(param=c("Alin","S0","kprog"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0*(1-exp(-kprog*t))+Alin*conc"))

immed_quad_null<-list(param=c("Alin","Aquad"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-Alin*conc+Aquad*conc**2"))
immed_quad_const<-list(param=c("Alin","Aquad","S0"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0+Alin*conc+Aquad*conc**2"))
immed_quad_lin<-list(param=c("Alin","Aquad","S0","kprog"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0+kprog*t+Alin*conc+Aquad*conc**2"))
immed_quad_exp<-list(param=c("Alin","Aquad","S0","kprog"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0*exp(-kprog*t)+Alin*conc+Aquad*conc**2"))
immed_quad_dexp<-list(param=c("Alin","Aquad","S0","kprog"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0*(1-exp(-kprog*t))+Alin*conc+Aquad*conc**2"))

immed_log_null<-list(param=c("Alog"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-Alog*log(conc)"))
immed_log_const<-list(param=c("Alog","S0"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0+Alog*log(conc)"))
immed_log_lin<-list(param=c("Alog","S0","kprog"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0+kprog*t+Alog*log(conc)"))
immed_log_exp<-list(param=c("Alog","S0","kprog"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0*exp(-kprog*t)+Alog*log(conc)"))
immed_log_dexp<-list(param=c("Alog","S0","kprog"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0*(1-exp(-kprog*t))+Alog*log(conc)"))


immed_Emax_null<-list(param=c("Emax","C50"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-(Emax*conc)/(conc+C50)"))
immed_Emax_const<-list(param=c("Emax","C50","S0"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0+((Emax*conc)/(conc+C50))"))
immed_Emax_lin<-list(param=c("Emax","C50","S0","kprog"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0+kprog*t+((Emax*conc)/(conc+C50))"))
immed_Emax_exp<-list(param=c("Emax","C50","S0","kprog"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0*exp(-kprog*t)+((Emax*conc)/(conc+C50))"))
immed_Emax_dexp<-list(param=c("Emax","C50","S0","kprog"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0*(1-exp(-kprog*t))+((Emax*conc)/(conc+C50))"))

immed_gammaEmax_null<-list(param=c("Emax","C50","gamma"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-(Emax*conc**gamma)/(conc**gamma+C50**gamma)"))
immed_gammaEmax_const<-list(param=c("Emax","C50","gamma","S0"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0+((Emax*conc**gamma)/(conc**gamma+C50**gamma))"))
immed_gammaEmax_lin<-list(param=c("Emax","C50","gamma","S0","kprog"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0+kprog*t+((Emax*conc**gamma)/(conc**gamma+C50**gamma))"))
immed_gammaEmax_exp<-list(param=c("Emax","C50","gamma","S0","kprog"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0*exp(-kprog*t)+((Emax*conc**gamma)/(conc**gamma+C50**gamma))"))
immed_gammaEmax_dexp<-list(param=c("Emax","C50","gamma","S0","kprog"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0*(1-exp(-kprog*t))+((Emax*conc**gamma)/(conc**gamma+C50**gamma))"))

immed_Imax_null<-list(param=c("Imax","C50"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-1-((Imax*conc)/(conc+C50))"))
immed_Imax_const<-list(param=c("Imax","C50","S0"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0*(1-((Imax*conc)/(conc+C50)))"))
immed_Imax_lin<-list(param=c("Imax","C50","S0","kprog"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-(kprog*t+S0)*(1-((Imax*conc)/(conc+C50)))"))
immed_Imax_exp<-list(param=c("Imax","C50","S0","kprog"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0*(exp(-kprog*t)-((Imax*conc)/(conc+C50)))"))
immed_Imax_dexp<-list(param=c("Imax","C50","S0","kprog"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0*(1-exp(-kprog*t)-((Imax*conc)/(conc+C50)))"))


immed_gammaImax_null<-list(param=c("Imax","C50","gamma"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-1-((Imax*conc**gamma)/(conc**gamma+C50**gamma))"))
immed_gammaImax_const<-list(param=c("Imax","C50","gamma","S0"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0*(1-((Imax*conc**gamma)/(conc**gamma+C50**gamma)))"))
immed_gammaImax_lin<-list(param=c("Imax","C50","gamma","S0","kprog"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-(kprog*t+S0)*(1-((Imax*conc**gamma)/(conc**gamma+C50**gamma)))"))
immed_gammaImax_exp<-list(param=c("Imax","C50","gamma","S0","kprog"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0*(exp(-kprog*t)-((Imax*conc**gamma)/(conc**gamma+C50**gamma)))"))
immed_gammaImax_dexp<-list(param=c("Imax","C50","gamma","S0","kprog"),ret1="dpd1",ret2="pdIm",node=1,fun=c("dpd1<-0","pdIm<-S0*(1-exp(-kprog*t)-((Imax*conc**gamma)/(conc**gamma+C50**gamma)))"))


#-----------------------------------------------------------------------------------
# character 
#-----------------------------------------------------------------------------------
bolus_1cpt_Vk<-"bolus_1cpt_Vk"
bolus_1cpt_VCl<-"bolus_1cpt_VCl"
bolus_2cpt_Vkk12k21<-"bolus_2cpt_Vkk12k21"
bolus_2cpt_ClV1QV2<-"bolus_2cpt_ClV1QV2"
bolus_3cpt_Vkk12k21k13k31<-"bolus_3cpt_Vkk12k21k13k31"
bolus_3cpt_ClV1Q2V2Q3V3<-"bolus_3cpt_ClV1Q2V2Q3V3"
bolus_1cpt_VVmkm<-"bolus_1cpt_VVmkm"
bolus_2cpt_Vk12k21Vmkm<-"bolus_2cpt_Vk12k21Vmkm"
bolus_2cpt_V1QV2Vmkm<-"bolus_2cpt_V1QV2Vmkm"
bolus_3cpt_Vk12k21k13k31Vmkm<-"bolus_3cpt_Vk12k21k13k31Vmkm"
bolus_3cpt_V1Q2V2Q3V3Vmkm<-"bolus_3cpt_V1Q2V2Q3V3Vmkm"
oral1_1cpt_kaVk<-"oral1_1cpt_kaVk"
oral1_1cpt_kaVCl<-"oral1_1cpt_kaVCl"
oral1_1cpt_kaVk_md<-"oral1_1cpt_kaVk_md"
oral1_1cpt_kaVCl_md<-"oral1_1cpt_kaVCl_md"
oral1_2cpt_kaVkk12k21<-"oral1_2cpt_kaVkk12k21"
oral1_2cpt_kaClV1QV2<-"oral1_2cpt_kaClV1QV2"
oral1_2cpt_kaVkk12k21_md<-"oral1_2cpt_kaVkk12k21_md"
oral1_2cpt_kaClV1QV2_md<-"oral1_2cpt_kaClV1QV2_md"
oral1_3cpt_kaVkk12k21k13k31<-"oral1_3cpt_kaVkk12k21k13k31"
oral1_3cpt_kaClV1Q2V2Q3V3<-"oral1_3cpt_kaClV1Q2V2Q3V3"
oral1_3cpt_kaVkk12k21k13k31_md<-"oral1_3cpt_kaVkk12k21k13k31_md"
oral1_3cpt_kaClV1Q2V2Q3V3_md<-"oral1_3cpt_kaClV1Q2V2Q3V3_md"
oral1_1cpt_kaVVmkm<-"oral1_1cpt_kaVVmkm"
oral1_1cpt_kaVVmkm_md<-"oral1_1cpt_kaVVmkm_md"
oral1_2cpt_kaVk12k21Vmkm<-"oral1_2cpt_kaVk12k21Vmkm"
oral1_2cpt_kaV1QV2Vmkm<-"oral1_2cpt_kaV1QV2Vmkm"
oral1_2cpt_kaVk12k21Vmkm_md<-"oral1_2cpt_kaVk12k21Vmkm_md"
oral1_2cpt_kaV1QV2Vmkm_md<-"oral1_2cpt_kaV1QV2Vmkm_md"
oral1_3cpt_kaVk12k21k13k31Vmkm<-"oral1_3cpt_kaVk12k21k13k31Vmkm"
oral1_3cpt_kaV1Q2V2Q3V3Vmkm<-"oral1_3cpt_kaV1Q2V2Q3V3Vmkm"
oral1_3cpt_kaVk12k21k13k31Vmkm_md<-"oral1_3cpt_kaVk12k21k13k31Vmkm_md"
oral1_3cpt_kaV1Q2V2Q3V3Vmkm_md<-"oral1_3cpt_kaV1Q2V2Q3V3Vmkm_md"
infusion_1cpt_Vk<-"infusion_1cpt_Vk"
infusion_1cpt_VCl<-"infusion_1cpt_VCl"
infusion_1cpt_Vk_md<-"infusion_1cpt_Vk_md"
infusion_1cpt_VCl_md<-"infusion_1cpt_VCl_md"
infusion_2cpt_Vkk12k21<-"infusion_2cpt_Vkk12k21"
infusion_2cpt_ClV1QV2<-"infusion_2cpt_ClV1QV2"
infusion_2cpt_Vkk12k21_md<-"infusion_2cpt_Vkk12k21_md"
infusion_2cpt_ClV1QV2_md<-"infusion_2cpt_ClV1QV2_md"
infusion_3cpt_Vkk12k21k13k31<-"infusion_3cpt_Vkk12k21k13k31"
infusion_3cpt_ClV1Q2V2Q3V3<-"infusion_3cpt_ClV1Q2V2Q3V3"
infusion_3cpt_Vkk12k21k13k31_md<-"infusion_3cpt_Vkk12k21k13k31_md"
infusion_3cpt_ClV1Q2V2Q3V3_md<-"infusion_3cpt_ClV1Q2V2Q3V3_md"
infusion_1cpt_VVmkm<-"infusion_1cpt_VVmkm"
infusion_1cpt_VVmkm_md<-"infusion_1cpt_VVmkm_md"
infusion_2cpt_Vk12k21Vmkm<-"infusion_2cpt_Vk12k21Vmkm"
infusion_2cpt_V1QV2Vmkm<-"infusion_2cpt_V1QV2Vmkm"
infusion_2cpt_Vk12k21Vmkm_md<-"infusion_2cpt_Vk12k21Vmkm_md"
infusion_2cpt_V1QV2Vmkm_md<-"infusion_2cpt_V1QV2Vmkm_md"
infusion_3cpt_Vk12k21k13k31Vmkm<-"infusion_3cpt_Vk12k21k13k31Vmkm"
infusion_3cpt_V1Q2V2Q3V3Vmkm<-"infusion_3cpt_V1Q2V2Q3V3Vmkm"
infusion_3cpt_Vk12k21k13k31Vmkm_md<-"infusion_3cpt_Vk12k21k13k31Vmkm_md"
infusion_3cpt_V1Q2V2Q3V3Vmkm_md<-"infusion_3cpt_V1Q2V2Q3V3Vmkm_md"
 


