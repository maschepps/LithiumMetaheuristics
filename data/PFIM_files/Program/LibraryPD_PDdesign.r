#PFIM 4.0 
#April 2014
#Copyright (C) PFIM 4.0 - Universite Paris Diderot and INSERM.
#------------------------------------------------------------------------------------  

#-------------------------
# Linear drug action model
#-------------------------

#---------
# null baseline
#---------
#$Model definition
#$Alin
#$
immed_lin_null<-function()
{
  form1<-"Alin*t"
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
immed_lin_const<-function()
{
  form1<-"S0+Alin*t"
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
immed_quad_null<-function()
{
  form1<-"Alin*t+Aquad*t**2"
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
immed_quad_const<-function()
{
  form1<-"S0+Alin*t+Aquad*t**2"
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
immed_log_null<-function()
{
  form1<-"Alog*log(t)"
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
immed_log_const<-function()
{
  form1<-"S0+Alog*log(t)"
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
immed_Emax_null<-function()
{
  form1<-"Emax*t/(t+C50)"
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
immed_Emax_const<-function()
{
  form1<-"Emax*t/(t+C50)+S0"
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
immed_gammaEmax_null<-function()
{
  form1<-"Emax*t**gamma/(t**gamma+C50**gamma)"
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
immed_gammaEmax_const<-function()
{
  form1<-"Emax*t**gamma/(t**gamma+C50**gamma)+S0"
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
immed_Imax_null<-function()
{
  form1<-"1-(Imax*t/(t+C50))"
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
immed_Imax_const<-function()
{
  form1<-"(1-(Imax*t/(t+C50)))*S0"
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
immed_gammaImax_null<-function()
{
  form1<-"1-(Imax*t**gamma/(t**gamma+C50**gamma))"
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
immed_gammaImax_const<-function()
{
  form1<-"(1-(Imax*t**gamma/(t**gamma+C50**gamma)))*S0"
	form1<-parse(text=form1,n=-1)
	tf<-Inf
	return(list(c(form1),tf))
}
