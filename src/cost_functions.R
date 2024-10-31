# Universal Number of Subjects
Number_of_subjects<-100



###1grp 5 sampling times
Ds1cov<-function(time){
  res<-funFIMem_ds3(formSS3_cov,c("KA","V","ClS", "CLES", "CLSE","BETA"),
                    c(0.93,22.30,1.24,11.1,4.15,0.32),c(0.72,0.30,0.20,0.27,0,0),c(0, 0.137),
                    list(time[1:5],time[1:5]),2,c(0,1),c((time[6]/sum(time[6]))*0.8,(time[6]/sum(time[6]))*0.2),Number_of_subjects)
  return(-res)
}

D1cov<-function(time){
  res<-funFIMem(formSS3_cov,c("KA","V","ClS", "CLES", "CLSE","BETA"),
                c(0.93,22.30,1.24,11.1,4.15,0.32),c(0.72,0.30,0.20,0.27,0,0),c(0, 0.137),
                list(time[1:5],time[1:5]),2,c(0,1),c((time[6]/sum(time[6]))*0.8,(time[6]/sum(time[6]))*0.2),Number_of_subjects)
  return(-res)
}
# 2 group
############2 grp
D2cov<-function(time){
  res<-funFIMem(formSS3_cov,c("KA","V","ClS", "CLES", "CLSE","BETA"),
                c(0.93,22.30,1.24,11.1,4.15,0.32),c(0.72,0.30,0.20,0.27,0,0),c(0, 0.137),
                list(time[1:4],time[1:4],time[5:8],time[5:8]),2,c(0,1,0,1),c((time[9]/sum(time[9:10]))*0.8,(time[9]/sum(time[9:10]))*0.2,
                                                                             (time[10]/sum(time[9:10]))*0.8,(time[10]/sum(time[9:10]))*0.2),Number_of_subjects)
  return(-res)
}
Ds2cov<-function(time){
  res<-funFIMem_ds3(formSS3_cov,c("KA","V","ClS", "CLES", "CLSE","BETA"),
                    c(0.93,22.30,1.24,11.1,4.15,0.32),c(0.72,0.30,0.20,0.27,0,0),c(0, 0.137),
                    list(time[1:4],time[1:4],time[5:8],time[5:8]),2,c(0,1,0,1),c((time[9]/sum(time[9:10]))*0.8,(time[9]/sum(time[9:10]))*0.2,
                                                                                 (time[10]/sum(time[9:10]))*0.8,(time[10]/sum(time[9:10]))*0.2),Number_of_subjects)
  return(-res)
}

D1 <-function(time){
  res<-funFIMem(formSS3,c("KA","V","ClS", "CLES", "CLSE"),
                c(0.93,22.30,1.24,11.1,4.15),c(0.72,0.30,0.20,0.27,0),c(0, 0.137),
                list(time[1:5]),2,c(36),c(time[6]/sum(time[6])),Number_of_subjects)
  return(-res)
}
Ds1 <-function(time){
  res<-funFIMem_ds2(formSS3,c("KA","V","ClS", "CLES", "CLSE"),
                    c(0.93,22.30,1.24,11.1,4.15),c(0.72,0.30,0.20,0.27,0),c(0, 0.137),
                    list(time[1:5]),2,c(36),c(time[6]/sum(time[6])),Number_of_subjects)
  return(-res)
}

D2<-function(time){
  res<-funFIMem(formSS3,c("KA","V","ClS", "CLES", "CLSE"),
                c(0.93,22.30,1.24,11.1,4.15),c(0.72,0.30,0.20,0.27,0),c(0, 0.137),
                list(time[1:4],time[5:8]),2,c(36,36),c(time[9]/sum(time[9:10]),time[10]/sum(time[9:10])),Number_of_subjects)
  return(-res)
}

Ds2<-function(time){
  res<-funFIMem_ds2(formSS3,c("KA","V","ClS", "CLES", "CLSE"),
                    c(0.93,22.30,1.24,11.1,4.15),c(0.72,0.30,0.20,0.27,0),c(0, 0.137),
                    list(time[1:4],time[5:8]),2,c(36,36),c(time[9]/sum(time[9:10]),time[10]/sum(time[9:10])),Number_of_subjects)
  return(-res)
}

D5<-function(time){
  res<-funFIMem(formSS3,c("KA","V","ClS", "CLES", "CLSE"),
                c(0.93,22.30,1.24,11.1,4.15),c(0.72,0.30,0.20,0.27,0),c(0, 0.137),
                list(time[1:4],time[5:8],time[9:12],time[13:16],time[17:20]),2,c(36,36,36,36,36),c(time[21]/sum(time[21:25]),time[22]/sum(time[21:25]),time[23]/sum(time[21:25]),time[24]/sum(time[21:25]),time[25]/sum(time[21:25])),Number_of_subjects)
  return(-res)
}

Ds5<-function(time){
  res<-funFIMem_ds2(formSS3,c("KA","V","ClS", "CLES", "CLSE"),
                c(0.93,22.30,1.24,11.1,4.15),c(0.72,0.30,0.20,0.27,0),c(0, 0.137),
                list(time[1:4],time[5:8],time[9:12],time[13:16],time[17:20]),2,c(36,36,36,36,36),c(time[21]/sum(time[21:25]),time[22]/sum(time[21:25]),time[23]/sum(time[21:25]),time[24]/sum(time[21:25]),time[25]/sum(time[21:25])),Number_of_subjects)
  return(-res)
}

