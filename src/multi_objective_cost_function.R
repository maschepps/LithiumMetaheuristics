#Multiobjective function
multi_obj_result <- function(time){
  funFIMem_nsga_ds3 <- function(equation,paramName,beta,o,sigma,t_group,Trand,d,PropSubjects,nbTot){
    #Name of the fixed effect parameters and sampling times
    paramF<-c(paramName,"t")
    
    #model equation
    form1<- equation
    
    PSI<-c(paramName,"sig.inter","sig.slope")
    lpsi<-length(PSI)
    
    
    #(Diagonal Matrix of) variance for inter-subject random effects:
    omega<-diag(o)
    
    
    #Random effect model (1) = additive  (2) = exponential 
    #------------------------------------------------------------------
    
    if ( Trand == 1 ) {
      form11 <- form1
      for(i in 1:length(paramName)){
        form11 <- gsub(paramName[i], paste0("(",paramName[i],"+b)"), form11)
        
      }
      
    } else {
      form11 <- form1
      for(i in 1:length(paramName)){
        form11 <- gsub(paramName[i], paste0("(",paramName[i],"*exp(b))"), form11)
      }
    }
    
    #gather all groups of protocol
    M_f<-list()
    M_F <- matrix(rep(0),nrow=lpsi+length(paramName),ncol=lpsi+length(paramName))
    for(q in 1:length(t_group)){
      t<-c()
      t<-c(t_group[[q]])
      
      #dose value
      
      dose<-d[q]
      
      #calculate matrix E for n individuals
      
      equatf <- parse(text = form11, n=-1)
      f<-function(paramF){eval(equatf[[1]])}
      
      #Fixed effects parameters values
      for(i in 1:length(paramName)){
        assign(paramName[i],beta[i])
      }
      
      # for(i in 1:length(PSI2)){
      #   assign(PSI2[i], test3[i])
      # }
      
      param <- c(beta,t)
      #calculate the predictions with fixed effects
      b <- 0
      fixed<-f(param)
      
      #Standard deviation of residual error (sig.inter+sig.slope*f)^2:
      sig.inter<-sigma[1]
      sig.slope<-sigma[2]
      
      var<-diag((sig.inter+sig.slope*fixed)^2)
      if(length(t) == 1){
        var = ((sig.inter+sig.slope*fixed)^2)
      }
      
      #calculate variance Var
      form2 <- paste0("( sig.inter + sig.slope * fixed )^2")
      Vmodel<- parse(text = form2)

      #get derivatives for fixed parameters
      df<-deriv(equatf[[1]],PSI)
      mdf<-attributes(eval(df))$gradient
      eval(df)[1]
      
      #delete the last two columns (correspond to sig.inter and sig.slope) 
      mdfi <- mdf[,-c(length(PSI)-1,length(PSI))]
      #complete derivative for exponential random effect model 
      if(Trand ==2 ){
        mdfie <- mdfi %*% diag(beta)
      }else {mdfie <- mdfi}
      
      
      #calculate variance Vi
      Vi <- mdfie %*% omega %*% t(mdfie) + var
      
      #get derivatives of sigma
      dv<-deriv(Vmodel[[1]],PSI)
      mdv<-attributes(eval(dv))$gradient
      
      
      #calculate matrix part A
      M_A <- t(mdfi) %*% solve(Vi) %*% mdfi
      
      
      #complete the rest of the matrix with 0
      for(i in 1:length(PSI)){
        M_A <- cbind(M_A,0)
        M_A <- rbind(M_A,0)
      }
      
      #calculate matrix part B
      #initialize the matrix with 0
      M_B <- matrix(rep(0),nrow=length(PSI)+length(paramName),ncol=length(PSI)+length(paramName))
      
      #prepare a list of matrix of derivatives of sigma to simplify usage
      m<-list()
      for(i in (length(paramName)+1):lpsi){
        if(length(t)==1){
          m[[i]] <- mdv[i]
        }else{
          m[[i]] <- diag(mdv[,i])
        }
        
      }
      #calculate first three rows of part B
      #Actually first 5 rows
      #Random effects f
      for(i in 1:length(paramName)){
        
        for(j in 1:length(paramName)){
          M_B[length(paramName)+i,length(paramName)+j] <- 1/2 * sum(diag(((mdfie[,i] %*% t(mdfie[,i])) %*% solve(Vi) %*% (mdfie[,j] %*% t(mdfie[,j])) %*% solve(Vi))))
        }
        for(j in (length(PSI)-1):length(PSI)){
          M_B[length(paramName)+i,length(paramName)+j] <- 1/2 * sum(diag(((mdfie[,i] %*% t(mdfie[,i])) %*% solve(Vi) %*% m[[j]] %*% solve(Vi))))
        }
        
      }
      #calculate the last two rows of partB
      for(i in (length(PSI)-1):length(PSI)){
        for(j in 1:length(paramName)){
          M_B[length(paramName)+i,length(paramName)+j] <- 1/2 * sum(diag( m[[i]] %*% solve(Vi) %*% (mdfie[,j] %*% t(mdfie[,j])) %*% solve(Vi)))
        }
        for(j in (length(PSI)-1):length(PSI)){
          M_B[length(paramName)+i,length(paramName)+j] <- 1/2 * sum(diag(m[[i]] %*% solve(Vi) %*% m[[j]] %*% solve(Vi)))
        }
      }
      
      M_f[[q]] <- (M_A+M_B)*PropSubjects[q]
      
      
      M_F <-M_F+M_f[[q]]
    }
    M_F <- M_F *nbTot
    
    
    #set names for vectors 
    fname<-c()
    for(n in 1:length(paramName)){
      fname<-c(fname,paste0("u_",paramName[n]))
    }
    for(n in 1:length(paramName)){
      fname<-c(fname,paste0("w_",paramName[n]))
    }
    fname<-c(fname,"sig.inter","sig.slope")
    rownames(M_F) <- fname
    colnames(M_F) <- fname
    
    
    if(0 %in% c(o,sigma)){
      
      M_F<-M_F[-c(length(paramName)+which(c(o,sigma)==0)),-c(length(paramName)+which(c(o,sigma)==0))]
      
      PSI<- PSI[-c(length(paramName)+which(c(o,sigma)==0))]
      
    }

    x = c(11,3,8,4,5,9,10,1,2,6,7)
    x = c(3,6,9,1,2,4,5,7,8,10,11)
    
    test = M_F[x, x]

    invmat = test[4:11,4:11]

    invmat[is.infinite(invmat)] = 10e-10
    invmat[is.na(invmat)] = 10e-10
   
    ans1 = test[1:3,1:3] - test[1:3,c(4:11)] %*% ginv(invmat) %*% t(test[1:3,c(4:11)])
    ans1 = test[1:3,1:3] - test[1:3,c(4:11)] %*% ginv(invmat) %*% t(test[1:3,c(4:11)])
    ans2 = invmat - test[c(4:11), c(1:3)] %*% ginv(test[1:3, 1:3]) %*% t(test[c(4:11), c(1:3)])
    
    # det(ans1)
    if(length(t)==1){
      return(M_F)
    }else{
      deter1 <- det(ans1) ^ (1/3)
      deter2 <- det(ans2) ^ (1/8)
      return(c(-deter1,-deter2))
    }
  }

  ans = funFIMem_nsga_ds3(formSS3_cov,c("KA","V","ClS", "CLES", "CLSE","BETA"),
                          c(0.93,22.30,1.24,11.1,4.15,0.32),c(0.72,0.30,0.20,0.27,0,0),c(0, 0.137),
                          list(time[1:5],time[1:5]),2,c(0,1),c((time[6]/sum(time[6]))*0.8,(time[6]/sum(time[6]))*0.2),Number_of_subjects)

  return(ans)
  
}

