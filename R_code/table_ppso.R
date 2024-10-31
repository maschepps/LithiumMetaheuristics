#try PSO lithium

#install.packages("devtools")
#devtools::install_github("TillF/ppso")
library(ppso)
#Load in 
setwd("~/Optimal_Lithium/src")
setwd("C:/Users/Admin/Desktop/Github/JRSS_Lithium2/Publish/src")
file_list = list.files()
sapply(file_list, source)

#####################################################
#######################Fix doses VA example#####################
#####################################################

###Lithium model

####################################

###two groups, 4 sampling times, reparametrized as Camille (macro)
#KSE devient CLSE, KES devient CLES, k devient ClS
#KSE= CLSE/V ; k = ClS/V ; KES = CLES/Ves = CLES/57.5

formSS3<-"dose * (KA/V * ((CLSE/V) - KA)/((((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) - KA) * (((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) - KA)) * exp(-KA * (t - 0))/(1 - exp(-KA * 24)) + KA/V * ((CLSE/V) - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2))/((KA - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2)) * (((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2))) * exp(-((((ClS/V) + 
(CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) * (t - 0))/(1 - exp(-((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) * 24)) + KA/V * ((CLSE/V) - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2))/((KA - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2)) * (((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2))) * exp(-((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) * (t - 
0))/(1 - exp(-((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) * 24)))" 
###one group 5 samp times
Dopt_lith_1grp_SS3<-function(time){
  res<-funFIMem(formSS3,c("KA","V","ClS", "CLES", "CLSE"),
                c(0.93,22.30,1.24,11.1,4.15),c(0.72,0.30,0.20,0.27,0),c(0, 0.137),
                list(time[1:5]),2,c(36),c(time[6]/sum(time[6])),Number_of_subjects)
  return(-res)
}

###evaluation of a design only (-D-criterion of the FIM)
time_test<-c(0,2,4,6,8,Number_of_subjects)
Dopt_lith_1grp_SS3(time = time_test)
#####

number_of_parameters<-6 
parameter_bounds = cbind(rep(0, number_of_parameters),rep(8, number_of_parameters))
ptm <- proc.time()
result6_cov_1grp<-optim_pso(objective_function = D1cov, number_of_parameters = number_of_parameters, number_of_particles = 40, 
                            max_number_of_iterations = 20,max_number_function_calls= 1000, 
                            parameter_bounds = parameter_bounds, 
                            initial_estimates=NULL, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                            lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
proc.time() - ptm

t_group1<-sort(round(result6_cov_1grp$par[1:5],2))

allocation1<-round(result6_cov_1grp$par[6]/sum(result6_cov_1grp$par[6])*Number_of_subjects*1, 2)

design.matrix<-matrix(c(t_group1, allocation1), byrow=TRUE, ncol=6)
colnames(design.matrix)<-c("time1", "time2", "time3", "time4",'time5', "N patients")

result6_cov_1grp$value
design.matrix

###two groups
Number_of_subjects<-100
Dopt_lith_2grp_SS3<-function(time){
  res<-funFIMem(formSS3,c("KA","V","ClS", "CLES", "CLSE"),
                c(0.93,22.30,1.24,11.1,4.15),c(0.72,0.30,0.20,0.27,0),c(0, 0.137),
                list(time[1:4],time[5:8]),2,c(36,36),c(time[9]/sum(time[9:10]),time[10]/sum(time[9:10])),Number_of_subjects)
  return(-res)
}


number_of_parameters<-10 
parameter_bounds = cbind(rep(0, number_of_parameters),rep(8, number_of_parameters))
ptm <- proc.time()
result7<-optim_pso(objective_function = D2cov, number_of_parameters = number_of_parameters, number_of_particles = 40, 
                   max_number_of_iterations = 20,max_number_function_calls= 1000, 
                   parameter_bounds = parameter_bounds, 
                   initial_estimates=NULL, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                   lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
proc.time() - ptm

t_group1<-round(result7$par[1:4],2)
t_group2<-round(result7$par[5:8],2)

allocation1<-round(result7$par[9]/sum(result7$par[9:10])*Number_of_subjects, 0)
allocation2<-round(result7$par[10]/sum(result7$par[9:10])*Number_of_subjects, 0)

design.matrix<-matrix(c(t_group1, allocation1, t_group2, allocation2), byrow=TRUE, ncol=5)
colnames(design.matrix)<-c("time1", "time2", "time3", "time4", "N patients")

result7$value
design.matrix




####idem 5groups

formSS3<-"dose * (KA/V * ((CLSE/V) - KA)/((((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) - KA) * (((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) - KA)) * exp(-KA * (t - 0))/(1 - exp(-KA * 24)) + KA/V * ((CLSE/V) - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2))/((KA - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2)) * (((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2))) * exp(-((((ClS/V) + 
(CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) * (t - 0))/(1 - exp(-((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) * 24)) + KA/V * ((CLSE/V) - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2))/((KA - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2)) * (((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2))) * exp(-((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) * (t - 
0))/(1 - exp(-((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) * 24)))" 

Number_of_subjects<-100
Dopt_lith_2grp_SS3<-function(time){
  res<-funFIMem(formSS3,c("KA","V","ClS", "CLES", "CLSE"),
                c(0.93,22.30,1.24,11.1,4.15),c(0.72,0.30,0.20,0.27,0),c(0, 0.137),
                list(time[1:4],time[5:8],time[9:12],time[13:16],time[17:20]),2,c(36,36,36,36,36),c(time[21]/sum(time[21:25]),time[22]/sum(time[21:25]),time[23]/sum(time[21:25]),time[24]/sum(time[21:25]),time[25]/sum(time[21:25])),Number_of_subjects)
  return(-res)
}


number_of_parameters<-25 
parameter_bounds = cbind(rep(0, number_of_parameters),rep(8, number_of_parameters))
ptm <- proc.time()
result7<-optim_pso(objective_function = Dopt_lith_2grp_SS3, number_of_parameters = number_of_parameters, number_of_particles = 40, 
                   max_number_of_iterations = 20,max_number_function_calls= 1000, 
                   parameter_bounds = parameter_bounds, 
                   initial_estimates=NULL, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                   lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
proc.time() - ptm

t_group1<-round(result7$par[1:4],2)
t_group2<-round(result7$par[5:8],2)
t_group3<-round(result7$par[9:12],2)
t_group4<-round(result7$par[13:16],2)
t_group5<-round(result7$par[17:20],2)


allocation1<-round(result7$par[21]/sum(result7$par[21:25])*Number_of_subjects, 2)
allocation2<-round(result7$par[22]/sum(result7$par[21:25])*Number_of_subjects, 2)
allocation3<-round(result7$par[23]/sum(result7$par[21:25])*Number_of_subjects, 2)
allocation4<-round(result7$par[24]/sum(result7$par[21:25])*Number_of_subjects, 2)
allocation5<-round(result7$par[25]/sum(result7$par[21:25])*Number_of_subjects, 2)


design.matrix<-matrix(c(t_group1, allocation1, t_group2, allocation2,t_group3, allocation3, t_group4, allocation4,t_group5, allocation5), byrow=TRUE, ncol=5)
colnames(design.matrix)<-c("time1", "time2", "time3", "time4", "N patients")

result7$value
design.matrix





#####################################################
#######################incorporate covariates#########
#####################################################
# use the dose as covariate



Number_of_subjects<-100
formSS3_cov<-"36 * (KA/V * ((CLSE/V) - KA)/((((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2) - KA) * (((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2) - KA)) * exp(-KA * (t - 0))/(1 - exp(-KA * 24)) + KA/V * ((CLSE/V) - ((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2))/((KA - ((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2)) * (((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2) - ((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2))) * exp(-((((ClS*exp(BETA*dose)/V) + 
(CLSE/V) + (CLES/57.5)) + sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2) * (t - 0))/(1 - exp(-((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2) * 24)) + KA/V * ((CLSE/V) - ((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2))/((KA - ((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2)) * (((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2) - ((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2))) * exp(-((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2) * (t - 
0))/(1 - exp(-((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2) * 24)))" 


###1grp 5 sampling times
Dopt_lith_1grp_SS3_cov<-function(time){
  res<-funFIMem(formSS3_cov,c("KA","V","ClS", "CLES", "CLSE","BETA"),
                c(0.93,22.30,1.24,11.1,4.15,0.32),c(0.72,0.30,0.20,0.27,0,0),c(0, 0.137),
                list(time[1:5],time[1:5]),2,c(0,1),c((time[6]/sum(time[6]))*0.8,(time[6]/sum(time[6]))*0.2),Number_of_subjects)
  return(-res)
}

number_of_parameters<-6 #nb de time (dans chaque bras) + nb de bras?
parameter_bounds = cbind(rep(0, number_of_parameters),rep(8, number_of_parameters))
ptm <- proc.time()
result6_cov_1grp<-optim_pso(objective_function = Dopt_lith_1grp_SS3_cov, number_of_parameters = number_of_parameters, number_of_particles = 40, 
                            max_number_of_iterations = 20,max_number_function_calls= 1000, 
                            parameter_bounds = parameter_bounds, 
                            initial_estimates=NULL, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                            lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
proc.time() - ptm

t_group1<-round(result6_cov_1grp$par[1:5],2)
t_group2<-round(result6_cov_1grp$par[1:5],2)


allocation1<-round(result6_cov_1grp$par[6]/sum(result6_cov_1grp$par[6])*Number_of_subjects*0.8, 2)
allocation2<-round(result6_cov_1grp$par[6]/sum(result6_cov_1grp$par[6])*Number_of_subjects*0.2, 2)

design.matrix<-matrix(c(t_group1, allocation1, t_group2, allocation2), byrow=TRUE, ncol=6)
colnames(design.matrix)<-c("time1", "time2", "time3", "time4",'time5', "N patients")

result6_cov_1grp$value
design.matrix
#



############2 grp
Dopt_lith_2grp_SS3_cov<-function(time){
  res<-funFIMem(formSS3_cov,c("KA","V","ClS", "CLES", "CLSE","BETA"),
                c(0.93,22.30,1.24,11.1,4.15,0.32),c(0.72,0.30,0.20,0.27,0,0),c(0, 0.137),
                list(time[1:4],time[1:4],time[5:8],time[5:8]),2,c(0,1,0,1),c((time[9]/sum(time[9:10]))*0.8,(time[9]/sum(time[9:10]))*0.2,
                                                                             (time[10]/sum(time[9:10]))*0.8,(time[10]/sum(time[9:10]))*0.2),Number_of_subjects)
  return(-res)
}

number_of_parameters<-10 #nb de time (dans chaque bras) + nb de bras?
parameter_bounds = cbind(rep(0, number_of_parameters),rep(8, number_of_parameters))
ptm <- proc.time()
result6_cov<-optim_pso(objective_function = Dopt_lith_2grp_SS3_cov, number_of_parameters = number_of_parameters, number_of_particles = 40, 
                       max_number_of_iterations = 20,max_number_function_calls= 1000, 
                       parameter_bounds = parameter_bounds, 
                       initial_estimates=NULL, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                       lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
proc.time() - ptm

t_group1<-round(result6_cov$par[1:4],2)
t_group2<-round(result6_cov$par[1:4],2)
t_group3<-round(result6_cov$par[5:8],2)
t_group4<-round(result6_cov$par[5:8],2)


allocation1<-round(result6_cov$par[9]/sum(result6_cov$par[9:10])*Number_of_subjects*0.8, 2)
allocation2<-round(result6_cov$par[9]/sum(result6_cov$par[9:10])*Number_of_subjects*0.2, 2)
allocation3<-round(result6_cov$par[10]/sum(result6_cov$par[9:10])*Number_of_subjects*0.8, 2)
allocation4<-round(result6_cov$par[10]/sum(result6_cov$par[9:10])*Number_of_subjects*0.2, 2)

design.matrix<-matrix(c(t_group1, allocation1, t_group2, allocation2,t_group3, allocation3, t_group4, allocation4), byrow=TRUE, ncol=5)
colnames(design.matrix)<-c("time1", "time2", "time3", "time4", "N patients")

result6_cov$value
design.matrix
##############5 groupes
Dopt_lith_5grp_SS3_cov<-function(time){
  res<-funFIMem(formSS3_cov,c("KA","V","ClS", "CLES", "CLSE","BETA"),
                c(0.93,22.30,1.24,11.1,4.15,0.32),c(0.72,0.30,0.20,0.27,0,0),c(0, 0.137),
                list(time[1:4],time[1:4],time[5:8],time[5:8],time[9:12],time[9:12],time[13:16],time[13:16],time[17:20],time[17:20]),
                2,c(0,1,0,1,0,1,0,1,0,1),
                c((time[21]/sum(time[21:25]))*0.8,(time[21]/sum(time[21:25]))*0.2,
                  (time[22]/sum(time[21:25]))*0.8,(time[22]/sum(time[21:25]))*0.2,
                  (time[23]/sum(time[21:25]))*0.8,(time[23]/sum(time[21:25]))*0.2,
                  (time[24]/sum(time[21:25]))*0.8,(time[24]/sum(time[21:25]))*0.2,
                  (time[25]/sum(time[21:25]))*0.8,(time[25]/sum(time[21:25]))*0.2),
                Number_of_subjects)
  return(-res)
}

number_of_parameters<-25 #nb de time (dans chaque bras) + nb de bras?
parameter_bounds = cbind(rep(0, number_of_parameters),rep(8, number_of_parameters))
ptm <- proc.time()
result6_cov_5grp<-optim_pso(objective_function = Dopt_lith_5grp_SS3_cov, number_of_parameters = number_of_parameters, number_of_particles = 40, 
                            max_number_of_iterations = 20,max_number_function_calls= 1000, 
                            parameter_bounds = parameter_bounds, 
                            initial_estimates=NULL, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                            lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
proc.time() - ptm

t_group1<-round(result6_cov_5grp$par[1:4],2)
t_group2<-round(result6_cov_5grp$par[1:4],2)
t_group3<-round(result6_cov_5grp$par[5:8],2)
t_group4<-round(result6_cov_5grp$par[5:8],2)
t_group5<-round(result6_cov_5grp$par[9:12],2)
t_group6<-round(result6_cov_5grp$par[9:12],2)
t_group7<-round(result6_cov_5grp$par[13:16],2)
t_group8<-round(result6_cov_5grp$par[13:16],2)
t_group9<-round(result6_cov_5grp$par[17:20],2)
t_group10<-round(result6_cov_5grp$par[17:20],2)



allocation1<-round(result6_cov_5grp$par[21]/sum(result6_cov_5grp$par[21:25])*Number_of_subjects*0.8, 2)
allocation2<-round(result6_cov_5grp$par[21]/sum(result6_cov_5grp$par[21:25])*Number_of_subjects*0.2, 2)
allocation3<-round(result6_cov_5grp$par[22]/sum(result6_cov_5grp$par[21:25])*Number_of_subjects*0.8, 2)
allocation4<-round(result6_cov_5grp$par[22]/sum(result6_cov_5grp$par[21:25])*Number_of_subjects*0.2, 2)
allocation5<-round(result6_cov_5grp$par[23]/sum(result6_cov_5grp$par[21:25])*Number_of_subjects*0.8, 2)
allocation6<-round(result6_cov_5grp$par[23]/sum(result6_cov_5grp$par[21:25])*Number_of_subjects*0.2, 2)
allocation7<-round(result6_cov_5grp$par[24]/sum(result6_cov_5grp$par[21:25])*Number_of_subjects*0.8, 2)
allocation8<-round(result6_cov_5grp$par[24]/sum(result6_cov_5grp$par[21:25])*Number_of_subjects*0.2, 2)
allocation9<-round(result6_cov_5grp$par[25]/sum(result6_cov_5grp$par[21:25])*Number_of_subjects*0.8, 2)
allocation10<-round(result6_cov_5grp$par[25]/sum(result6_cov_5grp$par[21:25])*Number_of_subjects*0.2, 2)

design.matrix<-matrix(c(t_group1, allocation1, t_group2, allocation2,t_group3, allocation3, t_group4, allocation4,
                        t_group5, allocation5, t_group6, allocation6,t_group7, allocation7, t_group8, allocation8,
                        t_group9, allocation9, t_group10, allocation10), byrow=TRUE, ncol=5)
colnames(design.matrix)<-c("time1", "time2", "time3", "time4", "N patients")

result6_cov_5grp$value
design.matrix

#############
####
#### # D_s-optimal
###################
#Load in 
#Source files
setwd("~/Optimal_Lithium/src")
file_list = list.files()
sapply(file_list, source)

#####################################################
#######################Fix doses VA example#####################
#####################################################

###Lithium model

####################################

###two groups, 4 sampling times, reparametrized as Camille (macro)
#KSE devient CLSE, KES devient CLES, k devient ClS
#KSE= CLSE/V ; k = ClS/V ; KES = CLES/Ves = CLES/57.5

formSS3<-"dose * (KA/V * ((CLSE/V) - KA)/((((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) - KA) * (((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) - KA)) * exp(-KA * (t - 0))/(1 - exp(-KA * 24)) + KA/V * ((CLSE/V) - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2))/((KA - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2)) * (((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2))) * exp(-((((ClS/V) + 
(CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) * (t - 0))/(1 - exp(-((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) * 24)) + KA/V * ((CLSE/V) - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2))/((KA - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2)) * (((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2))) * exp(-((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) * (t - 
0))/(1 - exp(-((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) * 24)))" 
###one group 5 samp times
Number_of_subjects<-100
Dopt_lith_1grp_SS3<-function(time){
  res<-funFIMem_ds2(formSS3,c("KA","V","ClS", "CLES", "CLSE"),
                    c(0.93,22.30,1.24,11.1,4.15),c(0.72,0.30,0.20,0.27,0),c(0, 0.137),
                    list(time[1:5]),2,c(36),c(time[6]/sum(time[6])),Number_of_subjects)
  return(-res)
}

###evaluation of a design only (-D-criterion of the FIM)
time_test<-c(0,2,4,6,8,Number_of_subjects)
Dopt_lith_1grp_SS3(time = time_test)
#####

number_of_parameters<-6 
parameter_bounds = cbind(rep(0, number_of_parameters),rep(8, number_of_parameters))
ptm <- proc.time()
result6_cov_1grp<-optim_pso(objective_function = Dopt_lith_1grp_SS3, number_of_parameters = number_of_parameters, number_of_particles = 40, 
                            max_number_of_iterations = 20,max_number_function_calls= 1000, 
                            parameter_bounds = parameter_bounds, 
                            initial_estimates=NULL, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                            lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
proc.time() - ptm

t_group1<-round(result6_cov_1grp$par[1:5],2)



allocation1<-round(result6_cov_1grp$par[6]/sum(result6_cov_1grp$par[6])*Number_of_subjects*1, 2)

design.matrix<-matrix(c(t_group1, allocation1), byrow=TRUE, ncol=6)
colnames(design.matrix)<-c("time1", "time2", "time3", "time4",'time5', "N patients")

result6_cov_1grp$value
design.matrix

###two groups
Number_of_subjects<-100
Dopt_lith_2grp_SS3<-function(time){
  res<-funFIMem_ds2(formSS3,c("KA","V","ClS", "CLES", "CLSE"),
                    c(0.93,22.30,1.24,11.1,4.15),c(0.72,0.30,0.20,0.27,0),c(0, 0.137),
                    list(time[1:4],time[5:8]),2,c(36,36),c(time[9]/sum(time[9:10]),time[10]/sum(time[9:10])),Number_of_subjects)
  return(-res)
}


number_of_parameters<-10 
parameter_bounds = cbind(rep(0, number_of_parameters),rep(8, number_of_parameters))
ptm <- proc.time()
result7<-optim_pso(objective_function = Dopt_lith_2grp_SS3, number_of_parameters = number_of_parameters, number_of_particles = 40, 
                   max_number_of_iterations = 20,max_number_function_calls= 1000, 
                   parameter_bounds = parameter_bounds, 
                   initial_estimates=NULL, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                   lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
proc.time() - ptm

t_group1<-round(result7$par[1:4],2)
t_group2<-round(result7$par[5:8],2)

allocation1<-round(result7$par[9]/sum(result7$par[9:10])*Number_of_subjects, 0)
allocation2<-round(result7$par[10]/sum(result7$par[9:10])*Number_of_subjects, 0)

design.matrix<-matrix(c(t_group1, allocation1, t_group2, allocation2), byrow=TRUE, ncol=5)
colnames(design.matrix)<-c("time1", "time2", "time3", "time4", "N patients")

result7$value
design.matrix




####idem 5groups

formSS3<-"dose * (KA/V * ((CLSE/V) - KA)/((((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) - KA) * (((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) - KA)) * exp(-KA * (t - 0))/(1 - exp(-KA * 24)) + KA/V * ((CLSE/V) - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2))/((KA - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2)) * (((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2))) * exp(-((((ClS/V) + 
(CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) * (t - 0))/(1 - exp(-((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) * 24)) + KA/V * ((CLSE/V) - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2))/((KA - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2)) * (((((ClS/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) - ((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2))) * exp(-((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) * (t - 
0))/(1 - exp(-((((ClS/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS/V) * (CLSE/V)))/2) * 24)))" 

Number_of_subjects<-100
Dopt_lith_2grp_SS3<-function(time){
  res<-funFIMem_ds2(formSS3,c("KA","V","ClS", "CLES", "CLSE"),
                    c(0.93,22.30,1.24,11.1,4.15),c(0.72,0.30,0.20,0.27,0),c(0, 0.137),
                    list(time[1:4],time[5:8],time[9:12],time[13:16],time[17:20]),2,c(36,36,36,36,36),c(time[21]/sum(time[21:25]),time[22]/sum(time[21:25]),time[23]/sum(time[21:25]),time[24]/sum(time[21:25]),time[25]/sum(time[21:25])),Number_of_subjects)
  return(-res)
}


number_of_parameters<-25 
parameter_bounds = cbind(rep(0, number_of_parameters),rep(8, number_of_parameters))
ptm <- proc.time()
result7<-optim_pso(objective_function = Dopt_lith_2grp_SS3, number_of_parameters = number_of_parameters, number_of_particles = 40, 
                   max_number_of_iterations = 20,max_number_function_calls= 1000, 
                   parameter_bounds = parameter_bounds, 
                   initial_estimates=NULL, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                   lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
proc.time() - ptm

t_group1<-round(result7$par[1:4],2)
t_group2<-round(result7$par[5:8],2)
t_group3<-round(result7$par[9:12],2)
t_group4<-round(result7$par[13:16],2)
t_group5<-round(result7$par[17:20],2)


allocation1<-round(result7$par[21]/sum(result7$par[21:25])*Number_of_subjects, 2)
allocation2<-round(result7$par[22]/sum(result7$par[21:25])*Number_of_subjects, 2)
allocation3<-round(result7$par[23]/sum(result7$par[21:25])*Number_of_subjects, 2)
allocation4<-round(result7$par[24]/sum(result7$par[21:25])*Number_of_subjects, 2)
allocation5<-round(result7$par[25]/sum(result7$par[21:25])*Number_of_subjects, 2)


design.matrix<-matrix(c(t_group1, allocation1, t_group2, allocation2,t_group3, allocation3, t_group4, allocation4,t_group5, allocation5), byrow=TRUE, ncol=5)
colnames(design.matrix)<-c("time1", "time2", "time3", "time4", "N patients")

result7$value
design.matrix





#####################################################
#######################incorporate covariates#########
#####################################################
# use the dose as covariate



Number_of_subjects<-100
formSS3_cov<-"36 * (KA/V * ((CLSE/V) - KA)/((((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2) - KA) * (((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2) - KA)) * exp(-KA * (t - 0))/(1 - exp(-KA * 24)) + KA/V * ((CLSE/V) - ((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2))/((KA - ((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2)) * (((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2) - ((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2))) * exp(-((((ClS*exp(BETA*dose)/V) + 
(CLSE/V) + (CLES/57.5)) + sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2) * (t - 0))/(1 - exp(-((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2) * 24)) + KA/V * ((CLSE/V) - ((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2))/((KA - ((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2)) * (((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) + sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2) - ((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2))) * exp(-((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2) * (t - 
0))/(1 - exp(-((((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5)) - sqrt(((ClS*exp(BETA*dose)/V) + (CLSE/V) + (CLES/57.5))**2 - 4 * (ClS*exp(BETA*dose)/V) * (CLSE/V)))/2) * 24)))" 


###1grp 5 sampling times
Dopt_lith_1grp_SS3_cov<-function(time){
  res<-funFIMem(formSS3_cov,c("KA","V","ClS", "CLES", "CLSE","BETA"),
                c(0.93,22.30,1.24,11.1,4.15,0.32),c(0.72,0.30,0.20,0.27,0,0),c(0, 0.137),
                list(time[1:5],time[1:5]),2,c(0,1),c((time[6]/sum(time[6]))*0.8,(time[6]/sum(time[6]))*0.2),Number_of_subjects)
  return(-res)
}

number_of_parameters<-6 #nb de time (dans chaque bras) + nb de bras?
parameter_bounds = cbind(rep(0, number_of_parameters),rep(8, number_of_parameters))
ptm <- proc.time()
result6_cov_1grp<-optim_pso(objective_function = Dopt_lith_1grp_SS3_cov, number_of_parameters = number_of_parameters, number_of_particles = 40, 
                            max_number_of_iterations = 20,max_number_function_calls= 1000, 
                            parameter_bounds = parameter_bounds, 
                            initial_estimates=NULL, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                            lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
proc.time() - ptm

t_group1<-round(result6_cov_1grp$par[1:5],2)
t_group2<-round(result6_cov_1grp$par[1:5],2)


allocation1<-round(result6_cov_1grp$par[6]/sum(result6_cov_1grp$par[6])*Number_of_subjects*0.8, 2)
allocation2<-round(result6_cov_1grp$par[6]/sum(result6_cov_1grp$par[6])*Number_of_subjects*0.2, 2)

design.matrix<-matrix(c(t_group1, allocation1, t_group2, allocation2), byrow=TRUE, ncol=6)
colnames(design.matrix)<-c("time1", "time2", "time3", "time4",'time5', "N patients")

result6_cov_1grp$value
design.matrix
#



############2 grp
Dopt_lith_2grp_SS3_cov<-function(time){
  res<-funFIMem_ds3(formSS3_cov,c("KA","V","ClS", "CLES", "CLSE","BETA"),
                    c(0.93,22.30,1.24,11.1,4.15,0.32),c(0.72,0.30,0.20,0.27,0,0),c(0, 0.137),
                    list(time[1:4],time[1:4],time[5:8],time[5:8]),2,c(0,1,0,1),c((time[9]/sum(time[9:10]))*0.8,(time[9]/sum(time[9:10]))*0.2,
                                                                                 (time[10]/sum(time[9:10]))*0.8,(time[10]/sum(time[9:10]))*0.2),Number_of_subjects)
  return(-res)
}

number_of_parameters<-10 #nb de time (dans chaque bras) + nb de bras?
parameter_bounds = cbind(rep(0, number_of_parameters),rep(8, number_of_parameters))
ptm <- proc.time()
result6_cov<-optim_pso(objective_function = Dopt_lith_2grp_SS3_cov, number_of_parameters = number_of_parameters, number_of_particles = 40, 
                       max_number_of_iterations = 20,max_number_function_calls= 1000, 
                       parameter_bounds = parameter_bounds, 
                       initial_estimates=NULL, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                       lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
proc.time() - ptm

t_group1<-round(result6_cov$par[1:4],2)
t_group2<-round(result6_cov$par[1:4],2)
t_group3<-round(result6_cov$par[5:8],2)
t_group4<-round(result6_cov$par[5:8],2)


allocation1<-round(result6_cov$par[9]/sum(result6_cov$par[9:10])*Number_of_subjects*0.8, 2)
allocation2<-round(result6_cov$par[9]/sum(result6_cov$par[9:10])*Number_of_subjects*0.2, 2)
allocation3<-round(result6_cov$par[10]/sum(result6_cov$par[9:10])*Number_of_subjects*0.8, 2)
allocation4<-round(result6_cov$par[10]/sum(result6_cov$par[9:10])*Number_of_subjects*0.2, 2)

design.matrix<-matrix(c(t_group1, allocation1, t_group2, allocation2,t_group3, allocation3, t_group4, allocation4), byrow=TRUE, ncol=5)
colnames(design.matrix)<-c("time1", "time2", "time3", "time4", "N patients")

result6_cov$value
design.matrix
##############5 groupes
Dopt_lith_5grp_SS3_cov<-function(time){
  res<-funFIMem_ds3(formSS3_cov,c("KA","V","ClS", "CLES", "CLSE","BETA"),
                    c(0.93,22.30,1.24,11.1,4.15,0.32),c(0.72,0.30,0.20,0.27,0,0),c(0, 0.137),
                    list(time[1:4],time[1:4],time[5:8],time[5:8],time[9:12],time[9:12],time[13:16],time[13:16],time[17:20],time[17:20]),
                    2,c(0,1,0,1,0,1,0,1,0,1),
                    c((time[21]/sum(time[21:25]))*0.8,(time[21]/sum(time[21:25]))*0.2,
                      (time[22]/sum(time[21:25]))*0.8,(time[22]/sum(time[21:25]))*0.2,
                      (time[23]/sum(time[21:25]))*0.8,(time[23]/sum(time[21:25]))*0.2,
                      (time[24]/sum(time[21:25]))*0.8,(time[24]/sum(time[21:25]))*0.2,
                      (time[25]/sum(time[21:25]))*0.8,(time[25]/sum(time[21:25]))*0.2),
                    Number_of_subjects)
  return(-res)
}

number_of_parameters<-25 #nb de time (dans chaque bras) + nb de bras?
parameter_bounds = cbind(rep(0, number_of_parameters),rep(8, number_of_parameters))
ptm <- proc.time()
result6_cov_5grp<-optim_pso(objective_function = Dopt_lith_5grp_SS3_cov, number_of_parameters = number_of_parameters, number_of_particles = 40, 
                            max_number_of_iterations = 20,max_number_function_calls= 1000, 
                            parameter_bounds = parameter_bounds, 
                            initial_estimates=NULL, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                            lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
proc.time() - ptm

t_group1<-round(result6_cov_5grp$par[1:4],2)
t_group2<-round(result6_cov_5grp$par[1:4],2)
t_group3<-round(result6_cov_5grp$par[5:8],2)
t_group4<-round(result6_cov_5grp$par[5:8],2)
t_group5<-round(result6_cov_5grp$par[9:12],2)
t_group6<-round(result6_cov_5grp$par[9:12],2)
t_group7<-round(result6_cov_5grp$par[13:16],2)
t_group8<-round(result6_cov_5grp$par[13:16],2)
t_group9<-round(result6_cov_5grp$par[17:20],2)
t_group10<-round(result6_cov_5grp$par[17:20],2)



allocation1<-round(result6_cov_5grp$par[21]/sum(result6_cov_5grp$par[21:25])*Number_of_subjects*0.8, 2)
allocation2<-round(result6_cov_5grp$par[21]/sum(result6_cov_5grp$par[21:25])*Number_of_subjects*0.2, 2)
allocation3<-round(result6_cov_5grp$par[22]/sum(result6_cov_5grp$par[21:25])*Number_of_subjects*0.8, 2)
allocation4<-round(result6_cov_5grp$par[22]/sum(result6_cov_5grp$par[21:25])*Number_of_subjects*0.2, 2)
allocation5<-round(result6_cov_5grp$par[23]/sum(result6_cov_5grp$par[21:25])*Number_of_subjects*0.8, 2)
allocation6<-round(result6_cov_5grp$par[23]/sum(result6_cov_5grp$par[21:25])*Number_of_subjects*0.2, 2)
allocation7<-round(result6_cov_5grp$par[24]/sum(result6_cov_5grp$par[21:25])*Number_of_subjects*0.8, 2)
allocation8<-round(result6_cov_5grp$par[24]/sum(result6_cov_5grp$par[21:25])*Number_of_subjects*0.2, 2)
allocation9<-round(result6_cov_5grp$par[25]/sum(result6_cov_5grp$par[21:25])*Number_of_subjects*0.8, 2)
allocation10<-round(result6_cov_5grp$par[25]/sum(result6_cov_5grp$par[21:25])*Number_of_subjects*0.2, 2)

design.matrix<-matrix(c(t_group1, allocation1, t_group2, allocation2,t_group3, allocation3, t_group4, allocation4,
                        t_group5, allocation5, t_group6, allocation6,t_group7, allocation7, t_group8, allocation8,
                        t_group9, allocation9, t_group10, allocation10), byrow=TRUE, ncol=5)
colnames(design.matrix)<-c("time1", "time2", "time3", "time4", "N patients")

result6_cov_5grp$value
design.matrix






