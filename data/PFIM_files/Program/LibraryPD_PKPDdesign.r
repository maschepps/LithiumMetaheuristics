#PFIM 4.0 
#April 2014
#Copyright (C) PFIM 4.0 - Universite Paris Diderot and INSERM.
#------------------------------------------------------------------------------------  



#ANALYTICAL FORM

#-------------------------
# Linear drug action model
#-------------------------

#---------
# null baseline
#---------
#$Model definition
#$Alin
#$
immed_lin_null<-function(fconc)
{
  form1<-paste("Alin*",fconc)
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# constant baseline
#---------
#$Model definition
#$Alin,S0
#$
immed_lin_const<-function(fconc)
{
  form1<-paste("S0+Alin*",fconc)
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# linear progression
#---------
#$Model definition
#$Alin,S0,kprog
#$
immed_lin_lin<-function(fconc)
{
  form1<-paste("S0+kprog*t+Alin*",fconc)
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# exponential increase
#---------
#$Model definition
#$Alin,S0,kprog
#$
immed_lin_exp<-function(fconc)
{
  form1<-paste("S0*exp(-kprog*t)+Alin*",fconc)
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# exponential decrease
#---------
#$Model definition
#$Alin,S0,kprog
#$
immed_lin_dexp<-function(fconc)
{
  form1<-paste("S0*(1-exp(-kprog*t))+Alin*",fconc)
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}


#-------------------------
# Quadratic drug action model
#-------------------------

#---------
# null baseline
#---------
#$Model definition
#$Alin,Aquad
#$
immed_quad_null<-function(fconc)
{
  form1<-paste("Alin*",fconc,"+Aquad*",fconc,"**2")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# constant baseline
#---------
#$Model definition
#$Alin,Aquad,S0
#$
immed_quad_const<-function(fconc)
{
  form1<-paste("S0+Alin*",fconc,"+Aquad*",fconc,"**2")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# linear progression
#---------
#$Model definition
#$Alin,Aquad,S0,kprog
#$
immed_quad_lin<-function(fconc)
{
  form1<-paste("S0+kprog*t+Alin*",fconc,"+Aquad*",fconc,"**2")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# exponential increase
#---------
#$Model definition
#$Alin,Aquad,S0,kprog
#$
immed_quad_exp<-function(fconc)
{
  form1<-paste("S0*exp(-kprog*t)+Alin*",fconc,"+Aquad*",fconc,"**2")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# exponential decrease
#---------
#$Model definition
#$Alin,Aquad,S0,kprog
#$
immed_quad_dexp<-function(fconc)
{
  form1<-paste("S0*(1-exp(-kprog*t))+Alin*",fconc,"+Aquad*",fconc,"**2")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}


#-------------------------
# Logarithmic drug action model
#-------------------------

#---------
# null baseline
#---------
#$Model definition
#$Alog
#$
immed_log_null<-function(fconc)
{
  form1<-paste("Alog*log(",fconc,")")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# constant baseline
#---------
#$Model definition
#$Alog,S0
#$
immed_log_const<-function(fconc)
{
  form1<-paste("S0+Alog*log(",fconc,")")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# linear progression
#---------
#$Model definition
#$Alog,S0,krpog
#$
immed_log_lin<-function(fconc)
{
  form1<-paste("S0+kprog*t+Alog*log(",fconc,")")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# exponential increase
#---------
#$Model definition
#$Alog,S0,krpog
#$
immed_log_exp<-function(fconc)
{
  form1<-paste("S0*exp(-kprog*t)+Alog*log(",fconc,")")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# exponential decrease
#---------
#$Model definition
#$Alog,S0,krpog
#$
immed_log_dexp<-function(fconc)
{
  form1<-paste("S0*(1-exp(-kprog*t))+Alog*log(",fconc,")")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}


#-------------------------
# Emax drug action model
#-------------------------

#---------
# null baseline
#---------
#$Model definition
#$Emax,C50
#$
immed_Emax_null<-function(fconc)
{
  form1<-paste("Emax*",fconc,"/(",fconc,"+C50)")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# constant baseline
#---------
#$Model definition
#$Emax,C50,S0
#$
immed_Emax_const<-function(fconc)
{
  form1<-paste("Emax*",fconc,"/(",fconc,"+C50)+S0")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# linear progression
#---------
#$Model definition
#$Emax,C50,S0,kprog
#$
immed_Emax_lin<-function(fconc)
{
  form1<-paste("Emax*",fconc,"/(",fconc,"+C50)+S0+kprog*t")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# exponential increase
#---------
#$Model definition
#$Emax,C50,S0,kprog
#$
immed_Emax_exp<-function(fconc)
{
  form1<-paste("Emax*",fconc,"/(",fconc,"+C50)+S0*exp(-kprog*t)")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# exponential decrease
#---------
#$Model definition
#$Emax,C50,S0,kprog
#$
immed_Emax_dexp<-function(fconc)
{
  form1<-paste("Emax*",fconc,"/(",fconc,"+C50)+S0*(1-exp(-kprog*t))")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#-------------------------
# sigmoid Emax drug action model
#-------------------------

#---------
# null baseline
#---------
#$Model definition
#$Emax,gamma,C50
#$
immed_gammaEmax_null<-function(fconc)
{
  form1<-paste("Emax*",fconc,"**gamma/(",fconc,"**gamma+C50**gamma)")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# constant baseline
#---------
#$Model definition
#$Emax,gamma,C50,S0
#$
immed_gammaEmax_const<-function(fconc)
{
  form1<-paste("Emax*",fconc,"**gamma/(",fconc,"**gamma+C50**gamma)+S0")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# linear progression
#---------
#$Model definition
#$Emax,gamma,C50,S0,kprog
#$
immed_gammaEmax_lin<-function(fconc)
{
  form1<-paste("Emax*",fconc,"**gamma/(",fconc,"**gamma+C50**gamma)+S0+kprog*t")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# exponential increase
#---------
#$Model definition
#$Emax,gamma,C50,S0,kprog
#$
immed_gammaEmax_exp<-function(fconc)
{
  form1<-paste("Emax*",fconc,"**gamma/(",fconc,"**gamma+C50**gamma)+S0*exp(-kprog*t)")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# exponential decrease
#---------
#$Model definition
#$Emax,gamma,C50,S0,kprog
#$
immed_gammaEmax_dexp<-function(fconc)
{
  form1<-paste("Emax*",fconc,"**gamma/(",fconc,"**gamma+C50**gamma)+S0*(1-exp(-kprog*t))")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#-------------------------
# Imax drug action model
#-------------------------

#---------
# null baseline
#---------
#$Model definition
#$Imax,C50
#$
immed_Imax_null<-function(fconc)
{
  form1<-paste("1-(Imax*",fconc,"/(",fconc,"+C50))")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# constant baseline
#---------
#$Model definition
#$Imax,C50,S0
#$
immed_Imax_const<-function(fconc)
{
  form1<-paste("(1-(Imax*",fconc,"/(",fconc,"+C50)))*S0")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# linear progression
#---------
#$Model definition
#$Imax,C50,S0,krpog
#$
immed_Imax_lin<-function(fconc)
{
  form1<-paste("(1-(Imax*",fconc,"/(",fconc,"+C50)))*(S0+kprog*t)")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# exopnential increase
#---------
#$Model definition
#$Imax,C50,S0,krpog
#$
immed_Imax_exp<-function(fconc)
{
  form1<-paste("(exp(-kprog*t)-(Imax*",fconc,"/(",fconc,"+C50)))*S0")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# exopnential decrease
#---------
#$Model definition
#$Imax,C50,S0,krpog
#$
immed_Imax_dexp<-function(fconc)
{
  form1<-paste("(1-exp(-kprog*t)-(Imax*",fconc,"/(",fconc,"+C50)))*S0")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}


#-------------------------
# sigmoid Imax drug action model
#-------------------------

#---------
# null baseline
#---------
#$Model definition
#$Imax,C50,gamma
#$
immed_gammaImax_null<-function(fconc)
{
  form1<-paste("1-(Imax*",fconc,"**gamma/(",fconc,"**gamma+C50**gamma))")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# constant baseline
#---------
#$Model definition
#$Imax,C50,gamma,S0
#$
immed_gammaImax_const<-function(fconc)
{
  form1<-paste("(1-(Imax*",fconc,"**gamma/(",fconc,"**gamma+C50**gamma)))*S0")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# linear progression
#---------
#$Model definition
#$Imax,C50,gamma,S0,kprog
#$
immed_gammaImax_lin<-function(fconc)
{
  form1<-paste("(1-(Imax*",fconc,"**gamma/(",fconc,"**gamma+C50**gamma)))*(S0+kprog*t)")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# exponential increase
#---------
#$Model definition
#$Imax,C50,gamma,S0,kprog
#$
immed_gammaImax_exp<-function(fconc)
{
  form1<-paste("(exp(-kprog*t)-(Imax*",fconc,"**gamma/(",fconc,"**gamma+C50**gamma)))*S0")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}

#---------
# exponential decrease
#---------
#$Model definition
#$Imax,C50,gamma,S0,kprog
#$
immed_gammaImax_dexp<-function(fconc)
{
  form1<-paste("(1-exp(-kprog*t)-(Imax*",fconc,"**gamma/(",fconc,"**gamma+C50**gamma)))*S0")
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}


