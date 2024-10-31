#########################################################################
##					        		                                               ##
##				            INPUT FILE FOR PFIM 4.0                          ##
#########################################################################


#Name of the project
#-------------------- 
project<-"w"

#Name of the file containing the PK or PD model
#----------------------------------------------
file.model<-"model_2cpss_macro_2.R"

#Name of the output file for the results and for the Fisher information matrix
#---------------------------------------
output<-"Stdout_optPSO_5times1grp_2cp_macro_optSimplex_Schepps_cov.r";
outputFIM<-"FIM_optPSO_5times1grp_2cp_macro_optSimplex_Schepps_cov.txt";
outputFIM<-"FIM_optPSO_5times1grp_2cp_macro_optSimplex_Schepps_cov.txt";

#FIM: Population (P) or Individual (I) or Bayesian (B) Fisher information matrix
#---------------------------------------
FIM<-"P"

#Previous information for population design (FIM<-"P") only:
#If previous information is available, please specify below the file name;
#otherwise leave it as the default
#-------------------------------------------------------
previous.FIM<-""

#RUN:  Evaluation (EVAL) or Optimisation (OPT) 
#-------------------------------------------------------
run<-"OPT"

#To display only  graphs of models and/or sensitivity functions before evaluating the Fisher Information matrix
graph.only<-F

#Block diagonal Fisher information matrix (option<-1) or complete Fisher information matrix (option<-2)
#----------------------------------------------------------
option<-1

#Number of responses
#--------------------------------------------------------------------
nr<-1



################### MODEL OPTION ###########################

#Model form: Differential equations (DE) or analytical form (AF)
#---------------------------------------------------------------

modelform<-"AF"

###### ANALYTICAL MODEL OPTION #############################
############################################################

#Identical dose in each elementary design (Yes=T, No=F)
#-------------------------------------------------------------
dose.identical<-F

# If 'Yes', enter the value of the dose, 
# else, enter the vector of the dose values for each elementary design
#--------------------------------------------------------------------
# dose <- c(36)
dose<-c(36)

#Vector of the times intervals of each expression  
#-----------------------------------------------------------
boundA<-list(c(0,Inf))
# boundB<-list(c(0,Inf))
# boundC<-list(c(0,Inf))

#Numerical derivatives  (Yes=T, No=F)
#If 'Yes', specify the model function "form" in the model file
#If 'No', specify the object "form" which is a vector of expressions in the model file
#-----------------------------------------------------------
NUM<-F 

###### END ANALYTICAL MODEL OPTION ########################



###### DIFFERENTIAL EQUATION OPTION ##########################
##############################################################

#Initial time for which initial conditions are given
#---------------------------------------------------
time.condinit<-0

#Identical initial conditions in each elementary design (Yes=T, No=F)
#-------------------------------------------------------------
condinit.identical<-T

# If 'Yes', enter once the expression of the initial values of the system at the initial time
# else, enter the vectors of the initial conditions for each elementary design
# If initial values depend on the parameters to be estimated, 
# enter this parameter into the expression without any quotation marks 
#---------------------------------------------------------
condinit<-expression(c(0,0))

# Error tolerance for solving differential equations
#----------------------------------------------------

RtolEQ<-1e-04
AtolEQ<-1e-04
Hmax<-Inf 

###### END DIFFERENTIAL EQUATION OPTION #################################



#Name of the fixed effects parameters
#-------------------------------------
parameters<-c("KA","V","ClS","CLSE","CLES")
parameters<-c("ka","V","ClS","CLSE","CLES")


#Fixed effects parameters values
#-------------------------------
beta<-c(0.93,22.3,1.24,4.15,11.1)


#Some parameters may not be estimated (not estimated = T, estimated = F)
#--------------------------------
beta.fixed<-c(F,F,F,F,F)

#Number of occasions
#--------------------------------------------------------------------------
n_occ<-1

#Random effect model (1) = additive  (2) = exponential 
#------------------------------------------------------------------
Trand<-2;

#Diagonal Matrix of variance for inter-subject random effects:
#---------------------------------------------------
omega<-diag(c(0.72,0.30,0.2,0,0.27))


#Diagonal Matrix of variance for inter-occasion random effects:
#---------------------------------------------------
gamma<-diag(c(0,0,0,0,0)) 

#Standard deviation of residual error (sig.inter+sig.slope*f)^2:
#------------------------------------------------------------------
sig.slopeA<-0.137
sig.interA<-0
# sig.interB<-0.02
# sig.slopeB<-0
# sig.interC<-1.85
# sig.slopeC<-0

#List of the vectors of sampling times for each elementary design 
#You can specify that a group has no sampling time by writing NULL 
#(ONLY if you have several response)
#-----------------------------------------------------------------
# protA<-list(c(0,0.4,2.3,4.9,8 ))
# protA<-list(c(0,2,5,8),c(0,1,4,8))
# protA<-list(c(0.00, 2,5, 8),c(0.00, 1,4, 8))

# protA<-list(c(0.00, 2,5, 8),c(0.00, 1,4, 8),c(0.00, 2,5, 8),c(0.00, 1,4, 8),c(0.00, 2,5, 8))
protA <- list(c(0,1,2, 4,6, 8))
# protA <- list(c(0,2,3,4,6, 8))
# protA <- list(c(0,2,4,8))

# protA<- list(c(sort(a1$result[3,1:5]), a1$result[6])) # 34.35098


#Vector of initial proportions or numbers of subjects for each elementary design 
#--------------------------------------------------------------
subjects<-c(100)

#Subjects input: (1) for number of subjects (2) for proportions of subjects
#---------------------------------------------------------------------------
subjects.input<-1

#If 'proportions of subjects' give the total number of samples
#-------------------------------------------------------------
# Ntot<-400




###################################################################
#                                                                 #
#                        Covariate model                          #
#                                                                 #
###################################################################

##########################################
# Covariates not changing with occasion  # 
##########################################

#Add covariate to the model  (Yes==T No==F)
#---------------------------------------------------------------------------
covariate.model<-T

#Vector of covariates
#---------------------------------------------------------------------
covariate.name<-list(c("gene"))

#Categories for each covariate (the first category is the reference)
#-----------------------------------------------------------------------
covariate.category<-list(gene=c("W","H"))

#Proportions of subjects in each category
#-------------------------------------------------------------------------
covariate.proportions<-list(gene=c(0.8,0.2))

#Parameter(s) associated with each covariate
#-------------------------------------------------------------------------
parameter.associated<-list(gene=c("ClS"))

# Values of covariate parameters in covariate model 
# (values of parameters for all other categories than the reference category (for which beta=0) 
# covariate is additive on parameter if additive random effect model (Trand=1)
# covariate is additive on log parameters if exponential random effect model (Trand=2)
#-----------------------------------------------------------------------
beta.covariate<-list(gene=list(c(0.32)))


#####################################
#Covariates changing with occasion  # 
#####################################


#Add covariate to the model   (Yes==T No==F)
#---------------------------------------------------------------------------
covariate_occ.model<-F

#Vector of covariates depending on the occasion
#---------------------------------------------------------------------
covariate_occ.name<-list(  c() )

#Categories for each covariate (the first category is the reference)
#-----------------------------------------------------------------------
#covariate_occ.category<-list(  admin=c() )

#Sequences of values of covariates at each occasion 
#Specify as many values in each sequence as number of occasions (n_occ) for each covariate
#-------------------------------------------------------------------------------------------------------
 
#covariate_occ.sequence<-list(  admin=list(c("A","B"),c("B","A"))  )
#covariate_occ.sequence<-list(  admin=list(c(),c())  )

#Proportions of elementary designs corresponding to each sequence of covariate values
#Specify as many values of proportion as number of sequences defined in covariate_occ.sequence for each covariate
#-----------------------------------------------------------------------------------------------------------------
covariate_occ.proportions<-list(  admin=c()  )

#Parameter(s) associated with each covariate
#-------------------------------------------------------------------------
parameter_occ.associated<-list(  admin=c()  )

# Values of covariate parameters in covariate model 
# (values of parameters for all other categories than the reference category (for which beta=0) 
# covariate is additive on parameter if additive random effect model (Trand=1)
# covariate is additive on log parameters if exponential random effect model (Trand=2)
#-----------------------------------------------------------------------
beta.covariate_occ<-list(  admin=list(c()  ))


#############################################
# Power and number of subjects              #
#############################################

#Type one error alpha 
#-----------------------------------------------------------------------------
alpha<-0.05

#Compute expected power for comparison test (Yes=T, No=F)
#---------------------------------------------------------------------------
compute.power<-T

#Compute the number of subjects needed for a given power for comparison test(Yes=T, No=F)
#----------------------------------------------------------------------------
compute.nni<-T

#Equivalence interval
interval_eq<-c(log(0.8),log(1.25))

#Compute expected power for equivalence test (Yes=T, No=F)
#---------------------------------------------------------------------------
compute.power_eq<-F

#Compute the number of subjects needed for a given power for equivalence test (Yes=T, No=F)
#----------------------------------------------------------------------------
compute.nni_eq<-F

#Set value the given power
#---------------------------------------------------------------------------
given.power<-0.9



############ONLY FOR OPTIMISATION ###############################

#Identical sampling times for each response
# (only if you do not have sampling times==NULL)
#-------------------------------------------------------------------------------------
identical.times<-T

######## OPTIMISATION ALGORITHM OPTION ###############

#Character string for choice of the optimisation algorithm: 
#	"FW" for the Fedorov-Wynn algorithm 
#	"SIMP" for the Simplex algorithm
#------------------------------------------

algo.option<-"SIMP"


########################
#SIMPLEX SPECIFICATION #
########################

#Optimisation of the proportions of subjects: (Yes=T, No=F)
#--------------------------------------------------------------

subjects.opt<-T

#Vector of lower and upper admissible sampling times
#---------------------------------------------------


lowerA<-c(0)
upperA<-c(8)

# lowerB<-c(0)
# upperB<-c(312)
# 
# lowerB<-c(0)
# upperB<-c(312)

#Minimum delay between two sampling times
#-------------------------------------------

delta.time<-0

#Print iteration step (Yes=T, No=F)
#---------------------------------

iter.print<-T


#Parameter for initial simplex building (%)
#------------------------------------------

simplex.parameter<-20


#Maximum iteration number
#------------------------

Max.iter<-5000


#Relative convergence tolerance
#------------------------------
Rctol<-1e-6



#############################
#FEDOROV-WYNN SPECIFICATION #
#############################


#Number of sampling windows
#--------------------------
nwindA<-1
nwindB<-1
# nwindC<-1


#List of vector of the allowed sampling times for each sampling window
#--------------------------------------------------------------------

sampwinA<-list(c(0,0.5,1,2,3,4,5,6,7,8))
# sampwinB<-list(c(288,289,290,291,293,295,311))
# sampwinC<-list(c(311))

#Fixed times (times which will be in all evaluated protocols, corresponding to fixed constraints)
#--------------------------------------------------------------------
fixed.timesA<-c()
# fixed.timesB<-c()
# fixed.timesC<-c()

#List of vector of allowed number of points to be taken from each sampling window
#------------------------------------------------------------------------------

nsampA<-list(c(4))
nsampB<-list(c(4))
# nsampC<-list(c(1))

#Maximum total number of sampling times per subject
#--------------------------------------------------

nmaxptsA<-4
nmaxptsB<-4
# nmaxptsC<-1

#Minimum total number of sampling times per subject
#--------------------------------------------------

nminptsA<-4
nminptsB<-4
# nminptsC<-1
############# END OF OPTIMISATION ALGORITHM OPTION ###############






############## GRAPH SPECIFICATION OPTION ###############

#graphical representation of the model (Yes=T, No=F)
#-------------------------------------
graph.logical<-T
#graphical representation of sensitivity functions (Yes=T, No=F)
#-------------------------------------
graphsensi.logical<-F


#Vector of Names on X axes for each response
#---------------------------------
# names.datax<-c("Time","Time","Time")
names.datax<-c("Time")

#Vector of Names on Y axes for each response
#---------------------------------
# names.datay<-c("Concentration","Concentration","Concentration")
names.datay<-c("Concentration")

#Controls logarithmic axes for the graphical representation of the model
#Values "xy", "x", or "y" produce log-log or log-x or log-y axes.
#(For standard graphic, log.logical<-F)
#--------------------------------------------------------------
log.logical<-F
#log.logical<-'y'

#Vector of lower and upper sampling times for the graphical representations
#-------------------------------------------------------------------------
graph.infA<-c(0)
graph.supA<-c(24)

# graph.infB<-c(288)
# graph.supB<-c(312)

# graph.infC<-c(0)
# graph.supC<-c(312)



#Vector of lower and upper concentration for the graphical representations
#------------------------------------------------------------------------
y.rangeA<-NULL # default range
# y.rangeB<-NULL
# y.rangeC<-NULL 

#y.rangeA<-c(1,10)

############# END OF GRAPH SPECIFICATION OPTION ###############













