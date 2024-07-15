#return the eigenvectors of the AMS inverse matrix for later unstrain
AMS_inv <- function(mat,type="v",prnt=TRUE, Shiny=FALSE){
  library(matlib)
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  AMS <- as.matrix(mat)
  #if type = v AMS is in vector form
  if (type=="v"){
    #AMS in the form of Vi, Vi_dec, Vi_inc
    A <- as.data.frame(matrix(c(AMS[2],AMS[3],
                                AMS[5],AMS[6],
                                AMS[8],AMS[9]),
                              ncol=2,byrow=T))
    colnames(A) <- c("dec","inc")
    A$x <- cos(d2r(A$dec))*cos(d2r(A$inc))
    A$y <- sin(d2r(A$dec))*cos(d2r(A$inc))
    A$z <- sin(d2r(A$inc))
    A_vec <- t(A[,3:5])
    A_val <- c(AMS[1],0,0,0,AMS[4],0,0,0,AMS[7])
    A_val <- matrix(A_val,nrow=3,byrow=T)
    AMS_M <- A_vec%*%A_val%*%inv(A_vec)
  }
  #if type = m AMS is in matrix form
  if (type=="m") {
    #AMS in the form of K11,K22,K33,K12,K23,K13
    AMS_M <- c(AMS[1],AMS[4],AMS[6],
               AMS[4],AMS[2],AMS[5],
               AMS[6],AMS[5],AMS[3])
    AMS_M <- matrix(AMS_M,nrow=3,byrow=T)
  }
  #take full 3x3 matrix
  if(type=="m3x3"){AMS_M <- mat}
  #inverts matrix
  invAMS <- inv(AMS_M)
  inv_AMSe <- eigen(invAMS,symmetric = TRUE)
  AMS_inv <- inv_AMSe$vectors
  AMS_inv_val <- inv_AMSe$values
  AMS_Me <- eigen(AMS_M,symmetric = TRUE)
  AMS_Mval <- AMS_Me$values
  #original anisotropy parameter
  L <- AMS_Mval[1]/AMS_Mval[2]
  F <- AMS_Mval[2]/AMS_Mval[3]
  P <- AMS_Mval[1]/AMS_Mval[3]
  #inverted anisotropy parameter
  Li <- AMS_inv_val[1]/AMS_inv_val[2]
  Fi <- AMS_inv_val[2]/AMS_inv_val[3]
  Pi <- AMS_inv_val[1]/AMS_inv_val[3]
  #calculate inverted direction of axese if shiny==true
  if(Shiny==TRUE){
    V1_inc <- r2d(asin(AMS_inv[3,1]/(sqrt((AMS_inv[1,1]^2)+(AMS_inv[2,1]^2)+(AMS_inv[3,1]^2)))))
    V1_dec <- (r2d(atan2(AMS_inv[2,1],AMS_inv[1,1])))%%360

    V2_inc <- r2d(asin(AMS_inv[3,2]/(sqrt((AMS_inv[1,2]^2)+(AMS_inv[2,2]^2)+(AMS_inv[3,2]^2)))))
    V2_dec <- (r2d(atan2(AMS_inv[2,2],AMS_inv[1,2])))%%360

    V3_inc <- r2d(asin(AMS_inv[3,3]/(sqrt((AMS_inv[1,3]^2)+(AMS_inv[2,3]^2)+(AMS_inv[3,3]^2)))))
    V3_dec <- (r2d(atan2(AMS_inv[2,3],AMS_inv[1,3])))%%360
    AMS_inv_eigen_tab <- round(matrix(c(V1_dec,V2_dec,V3_dec,Li,Fi,V1_inc,V2_inc,V3_inc),ncol = 5,byrow = T),digits=2)
    colnames(AMS_inv_eigen_tab) <- c("V1","V2","V3","L_inv","F_inv")
    rownames(AMS_inv_eigen_tab) <- c("Dec", "Inc")
    AMS_inv_eigen_tab[2,4:5] <- c("","")
  }

  if(prnt==TRUE){
    #print anisotropy parameters
    cat(paste("Anisotropy parameter:
L:",round(L,digits = 4),"
F:", round(F,digits = 4),"
P:", round(P,digits = 4),"
"))
  }

  #returns inverted anisotropy directions if Shiny is FALSE
  if(Shiny==FALSE){return(AMS_inv)}
  #return result list if Shiny is TRUE
  if(Shiny==TRUE){
    result <- list()
    result[[1]] <- AMS_inv
    result[[2]] <- AMS_inv_eigen_tab
    return(result)
  }
}

#Function that rotate geographic dec_inc pair(DI) into bedding coordinates
#if bedding ins not in the file as 3 and 4 column, can be added in function
bed_DI <- function(DI,in_file=TRUE, bed_az,bed_plunge,export=FALSE){
  #functions degree to radians and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  #functions spherical (Dec=x, Inc=y) to Cartesian
  s2cx <- function(x,y) {cos(d2r(x))*cos(d2r(y))}
  s2cy <- function(x,y) {sin(d2r(x))*cos(d2r(y))}
  s2cz <- function(y) {sin(d2r(y))}
  data <- DI
  data <- na.omit(data)
  #colnames change if BdecInc are in file or not
  if(in_file==FALSE){
    if(ncol(data)==3) {
      colnames(data) <- c("dec","inc","position")
      pos_exist=TRUE
      }
    if(ncol(data)==2) {
      colnames(data) <- c("dec","inc")
      pos_exist=FALSE
    }
  } else {
    if(ncol(data)==5) {
      colnames(data) <- c("dec","inc","baz","bplunge","position")
      pos_exist=TRUE
      }
    if(ncol(data)==4) {
      colnames(data) <- c("dec","inc","baz","bplunge")
      pos_exist=FALSE
      }
    }

  #sines and cosines if BdecInc are not in file
  if(in_file==FALSE){
    sbd <- -sin(d2r(bed_az))
    cbd <- cos(d2r(bed_az))
    sbi <- sin(d2r(bed_plunge))
    cbi <- cos(d2r(bed_plunge))
  }

  newDI <- as.data.frame(matrix(ncol=2,nrow=0))
  for(i in 1:nrow(data)){
    newDI_p <- as.data.frame(matrix(ncol=2,nrow=1))
    if(in_file==TRUE){
      sbd <- -sin(d2r(data[i,3]))
      cbd <- cos(d2r(data[i,3]))
      sbi <- sin(d2r(data[i,4]))
      cbi <- cos(d2r(data[i,4]))
    }
    x <- s2cx(data[i,1],data[i,2])
    y <- s2cy(data[i,1],data[i,2])
    z <- s2cz(data[i,2])
    xn <- x*(sbd^2+cbd^2*cbi)+
      y*(cbd*sbd*(1-cbi))+
      z*sbi*cbd
    yn <- x*cbd*sbd*(1-cbi)+
      y*(cbd^2+sbd*sbd*cbi)-
      z*sbd*sbi
    zn <- -(x*cbd*sbi-
              y*sbi*sbd-
              z*cbi)
    newdec <- r2d(atan2(yn,xn))
    newdec <- ifelse(newdec<0,newdec+360,newdec)
    newinc <- r2d(asin(zn))
    newDI_p[1,1:2] <- c(newdec,newinc)
    newDI <- rbind(newDI,newDI_p)
  }
  if(pos_exist==TRUE){
    newDI <- cbind(newDI,data$position)
    colnames(newDI) <- c("TCdec","TCinc","position")
  }else{colnames(newDI) <- c("TCdec","TCinc")}
  if(export==TRUE){write.csv(newDI,"tilt_corrected_directions.csv",row.names = FALSE)}
  return(newDI)
}

#check_bipolarity
bip_check <- function(DI){
  #fucnctions deg to rads and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  data <- DI[,1:2]
  data <- na.omit(data)
  colnames(data) <- c("dec", "inc")

  #directions in Cartesian coordinates
  data$x <- cos(d2r(data$dec))*cos(d2r(data$inc))
  data$y <- sin(d2r(data$dec))*cos(d2r(data$inc))
  data$z <- sin(d2r(data$inc))

  #elements of the distribution matrix
  T_elements <- c(sum((data$x)*(data$x)),sum(data$x*data$y),sum(data$x*data$z),
                  sum(data$y*data$x),sum(data$y*data$y),sum(data$y*data$z),
                  sum(data$z*data$x),sum(data$z*data$y),sum(data$z*data$z))

  #distribution matrix
  T <- matrix(T_elements,nrow=3, byrow=TRUE)

  #calculate and copy eigenvalues and vectors
  T_e <- eigen(T)
  T_vec <- T_e$vectors
  T_val <- T_e$values

  #calculate dec inc of max variance
  V1inc <- r2d(asin(T_vec[3,1]/(sqrt((T_vec[1,1]^2)+(T_vec[2,1]^2)+(T_vec[3,1]^2)))))
  V1dec <- r2d(atan2(T_vec[2,1],T_vec[1,1]))
  V1dec <- V1dec%%360

  #next  calculates difference between dec_inc and average
  data$Dec_aver <- rep(V1dec)
  data$Inc_aver <- rep(V1inc)
  data$delta <- abs(data$dec-data$Dec_aver)
  data$diff <- r2d(acos((sin(d2r(data$inc))*sin(d2r(data$Inc_aver)))+
                          (cos(d2r(data$inc))*cos(d2r(data$Inc_aver))*cos(d2r(data$delta)))))
  #Isolate modes
  m1ind <- as.numeric(which(data$diff<=90), arr.ind = TRUE)
  m2ind <- as.numeric(which(data$diff>90), arr.ind = TRUE)
  t <- ifelse(length(m2ind)>0 && length(m1ind)>0, TRUE, FALSE)
  return(t)
}

#function that generates resampled Data dec_inc
boots_DI <- function(DI) {
  library("tidyverse", warn.conflicts = FALSE)
  n <- nrow(DI)
  newDI <- DI[sample(n,n,replace = T),]
  return(newDI)
}


#flips all data toward common polarity
common_DI <- function(DI,down=TRUE, export=FALSE,name="common_dirs") {
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  data <- DI[,1:2]
  data <- na.omit(data)
  colnames(data) <- c("dec", "inc")
  #directions in Cartesian coordinates
  data$x <- cos(d2r(data$dec))*cos(d2r(data$inc))
  data$y <- sin(d2r(data$dec))*cos(d2r(data$inc))
  data$z <- sin(d2r(data$inc))
  #averaged Cartesian coordinates
  x_av <- mean(data$x)
  y_av <- mean(data$y)
  z_av <- mean(data$z)
  #elements of the distribution matrix
  T_elements <- c(sum((data$x)*(data$x)),sum(data$x*data$y),sum(data$x*data$z),
                  sum(data$y*data$x),sum(data$y*data$y),sum(data$y*data$z),
                  sum(data$z*data$x),sum(data$z*data$y),sum(data$z*data$z))
  #distribution matrix
  T <- matrix(T_elements,nrow=3, byrow=TRUE)
  #calculate and copy eigenvalues and vectors
  T_e <- eigen(T)
  T_vec <- T_e$vectors
  T_val <- T_e$values
  #calculate dec inc of max variance
  V1inc <- r2d(asin(T_vec[3,1]/(sqrt((T_vec[1,1]^2)+(T_vec[2,1]^2)+(T_vec[3,1]^2)))))
  V1dec <- r2d(atan2(T_vec[2,1],T_vec[1,1]))
  V1dec <- ifelse(V1dec<0,V1dec+360,V1dec)
  #force V1 down if down=TRUE or up if FALSE
  if(down==TRUE){
    if(V1inc<0){
      V1inc <- -V1inc
      V1dec <- ifelse((V1dec+180)>360,V1dec-180,V1dec+180)
    }
  }
  if(down==FALSE){
    if(V1inc>=0){
      V1inc <- -V1inc
      V1dec <- ifelse((V1dec+180)>360,V1dec-180,V1dec+180)
    }
  }
  #next  calculates difference between dec_inc and average
  data$Dec_aver <- rep(V1dec)
  data$Inc_aver <- rep(V1inc)
  data$delta <- abs(data$dec-data$Dec_aver)
  data$diff <- r2d(acos((sin(d2r(data$inc))*sin(d2r(data$Inc_aver)))+
                          (cos(d2r(data$inc))*cos(d2r(data$Inc_aver))*cos(d2r(data$delta)))))
  #flips all directions
  data$dec_N <- ifelse(data$diff>90,ifelse((data$dec+180)>360,data$dec-180,data$dec+180),data$dec)
  data$inc_N <- ifelse(data$diff>90,-data$inc,data$inc)
  new_dec_inc <- subset(data,select=c(dec_N,inc_N))
  colnames(new_dec_inc) <- c("dec", "inc")
  new_dec_inc
  #if export==TRUE export flipped data into file
  if(export==TRUE){
    write.csv(new_dec_inc,paste(name, ".csv"), row.names = FALSE)
    # cat(paste("File saved as",name,".csv"))
  }
  return(new_dec_inc)
}

#find crossing point of two dataset of xy coordinates
curve_cross <- function(a, b) {
  colnames(a) <- c("x","y")
  colnames(b) <- c("x","y")
  curve1 <- approxfun(a$x, a$y, rule = 2)
  curve2 <- approxfun(b$x, b$y, rule = 2)

  inters_x <- uniroot(function(x) curve1(x) - curve2(x),
                      c(min(a$x), max(a$x)))$root

  inters_y <- curve2(inters_x)
  cross <- cbind(inters_x,inters_y)
  return(cross)
}

#Dynamic VANDAMME or VGP(45) cutoff (EI before the filter)
#(Physics of the Earth and Planetary Interiors 85;1994)
cut_DI <- function(DI,VD=TRUE,lat=0,long=0,cutoff=40, geo=FALSE,inc_f=TRUE, export=FALSE, name="cut_dirs",Shiny=F){
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #if coordinates are geographic, takes the four columns to calculate TC directions
  if(geo==TRUE){
    DIAP <- DI
    data <- bed_DI((DIAP))
    data <- data[,1:2]
  }else{data <- DI[,1:2]}

  data <- na.omit(data)
  colnames(data) <- c("dec", "inc")
  #enter longitude an latitude and convert radians
  Slat_d <- lat
  Slong_d <- long
  Slat_r <- d2r(Slat_d)
  Slong_r <- d2r(Slong_d)
  #reiterate Vandamme filter applying incl-flattening correction first
  n <- 0
  repeat{
    #number of reiteration
    n <- n+1
    if(inc_f==TRUE){
      #calculate f factor of distribution and f=1 if it is not flattened
      I_E_Edec_f <- ffind(data,f_inc = 0.005)
      f <- ifelse(is.na(I_E_Edec_f[1,4])==TRUE, 1, I_E_Edec_f[nrow(I_E_Edec_f),4])
    }else{f <- 1}

    #add column with inc unflattened, plus different parameters
    unfl_data <- unflat_DI(data,f)
    data$inc_u <- unfl_data[,2]
    data$dec_r <- d2r(data$dec)
    data$inc_u_r <- d2r(data$inc_u)
    data$Slat_r <- rep(Slat_r)
    data$Slat_d <- rep(Slat_d)
    data$Slong_r <- rep(Slong_r)
    data$Slong_d <- rep(Slong_d)
    #calculate pole colatitude
    data$p_colat_r <- atan(2/tan(data$inc_u_r))
    data$p_colat_d <- r2d(data$p_colat_r)
    #calculate pole latitude
    data$PLat_r <- asin((sin(data$Slat_r)*cos(data$p_colat_r))+
                          (cos(data$Slat_r)*sin(data$p_colat_r)*cos(data$dec_r)))
    data$Plat_d<- r2d(data$PLat)
    #Longitudinal difference between site and pole
    data$LDist_r <- asin((sin(data$p_colat_r)*sin(data$dec_r))/cos(data$PLat_r))
    data$LDist_d <- r2d(data$LDist_r)
    #calculate longitude
    data$PLong_d <- ifelse(cos(data$p_colat_r)<(sin(data$Slat_r)*sin(data$PLat_r)),
                           data$Slong_d+180-data$LDist_d,
                           data$Slong_d+data$LDist_d)
    #long_lat in cartesian coordinates
    data$x <- cos(d2r(data$PLong_d))*cos(data$PLat)
    data$y <- sin(d2r(data$PLong_d))*cos(data$PLat)
    data$z <- sin(data$PLat)
    #average cartesian coordinates
    X_aver <- mean(data$x)
    Y_aver <- mean(data$y)
    Z_aver <- mean(data$z)
    #sum of all values along the axes
    X_sum <- sum(data$x)
    Y_sum <- sum(data$y)
    Z_sum <- sum(data$z)
    #magnitude of average
    B <- sqrt((X_aver^2)+(Y_aver^2)+(Z_aver^2))

    N <- as.numeric(length(data$dec))
    #calculate paleomagnetic pole
    Long_aver <- r2d(atan2(Y_aver,X_aver))
    #corrects for negative declination
    Long_aver <- ifelse(Long_aver<0,Long_aver+360,Long_aver)
    Lat_aver <- r2d(asin(Z_aver/B))
    #PPole long
    data$Pole_long <- rep(Long_aver)
    #PPole lat
    data$Pole_lat <- rep(Lat_aver)
    data$delta <- abs(data$PLong_d-data$Pole_long)
    #calculate angle between Pole and VGP
    data$diff <- r2d(acos((sin(d2r(data$Plat_d))*sin(d2r(data$Pole_lat)))+
                            (cos(d2r(data$Plat_d))*cos(d2r(data$Pole_lat))*cos(d2r(data$delta)))))
    if(VD==TRUE){
      #vandamme filtering calculation
      ASD <- sqrt(sum(((data$diff)^2)/(N-1)))
      A <- (1.8*ASD)+5
    }else{
      A <- cutoff
    }
    VGPcut <- as.numeric(which(data$diff>A), arr.ind = TRUE)
    if (length(VGPcut)!=0) {
      data <- data[-VGPcut,]
      #cut also lines from DIAP file if coordinates are geographic
      if(geo==TRUE){
        DIAP <- DIAP[-VGPcut,]
      }else{DI <- DI[-VGPcut,]}
    }
    data  <- data[,1:2]
    if (length(VGPcut)==0) break
  }
  #export cut file
  if(export==TRUE){
    if(geo==TRUE){write.csv(round(DIAP,digits = 2),paste(name,".csv"),row.names = FALSE)}else
    {write.csv(round(DI,digits = 2),paste(name,".csv"),row.names = FALSE)}
  }
  if(Shiny==F){
    cat(paste("Number of reiteration: ", n,"
"))
  }
  #return file
  ifelse(geo==TRUE,
         return(DIAP),
         return(DI))
}

#convert VGPs and site latitude and longitude in directions
DI_from_VGP <- function(VGPs, lat, long, export=FALSE,name="Directions") {
  #degree to radians and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  VGPs <- na.omit(VGPs)
  directions <- as.data.frame(matrix(nrow = 0,ncol = 2))
  for (i in 1:nrow(VGPs)){
    plong <- VGPs[i,1] %% 360
    plat <- VGPs[i,2]
    long <- long %% 360
    signdec <- 1
    Dphi <- abs(plong - long)
    if (Dphi != 0) signdec <- (plong - long) / Dphi
    if (lat == 90) lat <- 89.99
    thetaS <- d2r(90 - lat)
    thetaP <- d2r(90 - plat)
    Dphi <- d2r(Dphi)
    cosp <- cos(thetaS) * cos(thetaP) + sin(thetaS) * sin(thetaP) * cos(Dphi)
    thetaM <- acos(cosp)
    cosd <- (cos(thetaP) - cos(thetaM) * cos(thetaS)) / (sin(thetaM) * sin(thetaS))
    C <- abs(1 - cosd^2)
    dec <- ifelse(C != 0,-atan(cosd / sqrt(abs(C))) + (pi / 2),dec <- acos(cosd))
    if (-pi < signdec * Dphi && signdec < 0) dec <- 2 * pi - dec
    if (signdec * Dphi > pi) dec <- 2 * pi - dec
    dec <- r2d(dec) %% 360
    inc <- r2d(atan2(2 * cos(thetaM), sin(thetaM)))
    DecInc <- as.data.frame(t(c(dec,inc)))
    directions <- rbind(directions,DecInc)
  }
  colnames(directions) <- c("Dec","Inc")
  #export csv with directions if requested
  if(export==TRUE) {
    write.csv(round(directions,digits = 2),paste(name,".csv"),row.names = FALSE)
  }
  return(directions)
}

#function that calculate E_I couples of data and plot bootstrapped statistics
EI_boot <- function(DI,nb=1000,conf=95,export=TRUE, name="EI_boot_plot") {
  data <- DI
  data <- na.omit(data)
  data <- data[,1:2]
  colnames(data) <- c("dec", "inc")
  Inc_E_real <- inc_E_finder(data)
  Inc_E_real$V1inc <- abs(Inc_E_real$V1inc)
  Inc_E <- as.data.frame(matrix(ncol=3,nrow=0))
  cat(paste("Simulating",nb,"detasets.
"))
  for (i in 1:nb) {
    dataprov <- boots_DI(data)
    I_E_Ed <- inc_E_finder(dataprov)
    I_E_Ed$V1inc <- abs(I_E_Ed$V1inc)
    Inc_E <- rbind(Inc_E, I_E_Ed)
    if(((i%%50)==0)==TRUE) {
      cat(paste(i,"simulations out of",nb,"done
"))
    }
  }
  colnames(Inc_E) <- c("Inc","E","E_dec")
  Inc_E <- Inc_E[order(Inc_E$E),]
  Inc_E_bk <- Inc_E
  Inc_E <- Inc_E_bk

  confn <- conf/100
  num <- round((nb*(1-confn))/2,digits=0)
  Lconf <- num
  Uconf <- nb-num
  Inc_E <- Inc_E[Lconf:Uconf,]

  N <- as.character(length(data[,1]))
  Inc <- format(round(Inc_E_real$V1inc,1),nsmall=1)
  Ecut <- format(round(Inc_E_real$E,2),nsmall=2)
  V2 <- format(round(Inc_E_real$DV1V2,1),nsmall=1)

  #plot frame
  par(fig=c(0,0.7,0,1), new= FALSE)
  y_up <- ifelse(Inc_E_real$E>3.5, ceiling(Inc_E_real$E),3.5)

  plot(NA, xlim= c(0,90), ylim= c(1,y_up), xaxp= c(0,90, 9),
       xlab="Inclination (°)", ylab="Elongation")

  #Plot Bootstrapped data
  points(x=Inc_E$Inc,
         y=Inc_E$E,
         pch=16,
         col=rgb(0, 0, 1, 0.15),
         cex=0.6)

  #plot tk03.GAD model E-I
  x <- 0:90
  y <- tk03(x)
  points(x=x, y= y, type= "l", col="blue", lwd=3)

  #Plot real data E-I couple and values
  points(x=Inc_E_real$V1inc,y=Inc_E_real$E,pch=21,
         col="black", bg="white", cex=1.2)
  text <- paste("N:", N, "
Inc:", Inc, "
E:", Ecut,"
Edec:",V2)
  text(x=0, y=3.2,pos=4,text, cex= 0.8)

  #plot histogram of Dec top right
  par(fig=c(0.35,0.69,0.56,0.99), new=TRUE)

  hist(Inc_E$E_dec, xlim=c(-90,90), xaxp=c(-90,90,4),
       breaks = 90, xlab = "", ylab = "",
       main="", cex.axis=0.8, col= "blue", border="blue")
  #plot labels closer than standard to axes
  title(xlab = "Edec(°)", line=1.9, cex=0.3)
  title(ylab = "Frequency", line=2, cex=0.1)

  #plot confidence margin of declination if between -45 and 45
  if(Inc_E_real$DV1V2>-45 && Inc_E_real$DV1V2<45){
    Inc_E_Edec <- Inc_E_bk
    Inc_E_Edec <- Inc_E_Edec[order(Inc_E_Edec[,3]),]
    Inc_E_Edec <- Inc_E_Edec[Lconf:Uconf,]
    low_dec <- Inc_E_Edec[1,3]
    up_dec <- Inc_E_Edec[length(Inc_E_Edec[,3]),3]
    abline(v=low_dec,lwd=1,lty=2)
    abline(v=up_dec,lwd=1,lty=2)
  }
  #calculate low and high inclination error
  Inc_E <- Inc_E_bk
  Inc_E <- Inc_E[order(Inc_E[,1]),]
  Inc_E <- Inc_E[Lconf:Uconf,]
  low_inc <- Inc_E[1,1]
  up_inc <- Inc_E[length(Inc_E[,1]),1]

  if(bip_check(data)==TRUE){
    #plot equal area data two modes
    par(fig=c(0.55,1,0,0.6), new=TRUE)
    plot_DI(data,title="Original directions")
  }

  #plot equal area single mode
  par(fig=c(0.55,1,0.4,1), new=TRUE)
  plot_DI(data,single_mode = TRUE,title = "Single mode directions")

  if(Inc_E_real$DV1V2>-45 && Inc_E_real$DV1V2<45){
    results <- as.data.frame(matrix(ncol= 9, nrow=1))
    colnames(results) <- c("Inc","Low_inc","High_inc", "Elong","Low_E","High_E","E_dec","Low_E_dec","High_E_dec")
    results$Inc <- Inc
    results$Low_inc <- round(low_inc,digits = 1)
    results$High_inc <- round(up_inc,digits = 1)
    results$Elong <- Ecut
    results$Low_E <- round(min(Inc_E$E),digits= 2)
    results$High_E <- round(max(Inc_E$E),digits= 2)
    results$E_dec <- V2
    results$Low_E_dec <- round(low_dec,digits= 2)
    results$High_E_dec <- round(up_dec,digits=2)
  }else{
    results <- as.data.frame(matrix(ncol= 7, nrow=1))
    colnames(results) <- c("Inc","Low_inc","High_inc", "Elong","Low_E","High_E","E_dec")
    results$Inc <- Inc
    results$Low_inc <- round(low_inc,digits = 1)
    results$High_inc <- round(up_inc,digits = 1)
    results$Elong <- Ecut
    results$Low_E <- round(min(Inc_E$E),digits=2)
    results$High_E <- round(max(Inc_E$E),digits= 2)
    results$E_dec <- V2
  }
  print(results, row.names = FALSE)
  #reset screen
  par(fig=c(0,1,0,1))
  #export if requested
  if(export==TRUE){
    cat("
Results saved as ", paste(name,".csv"),"
Graph saved as", paste(name,".pdf"),"
")
    write.csv(results,file=paste(name,".csv"),row.names = FALSE)
    save_pdf(name=paste(name,".pdf"))
  }
}

#function that calculated DeltaDec and DeltaInc
ellips_DI <- function(DI,lat,long,export=FALSE){
  #degree to radians and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  #cut NA rows
  DI <- na.omit(DI)
  #calculate average dir and paleolat in radians
  aver_DI <- fisher(DI)
  aver_inc <- aver_DI[1,2]
  lat_r <- atan((tan(d2r(aver_inc)))/2)
  #calculate vgps
  poles <- VGP_DI(DI,in_file=FALSE,lat=lat,long=long,export=F,type="VGPs",Prnt=FALSE)
  #calculate A95
  PPole <- fisher(poles)
  A95 <- PPole[1,3]
  #calculated ∂Dec and ∂Inc
  D_dec <- r2d(asin((sin(d2r(A95)))/cos(lat_r)))
  D_inc <- r2d((2*d2r(A95))/(1+(3*(sin(lat_r))^2)))
  result <- as.data.frame(matrix(nrow = 1,ncol=4))
  colnames(result) <- c("dec","inc","delta_dec","delta_inc")
  result[1,1] <- aver_DI[1,1]
  result[1,2] <- aver_DI[1,2]
  result[1,3] <- D_dec
  result[1,4] <- D_inc
  result$N <- nrow(DI)
  #esport results if requested
  if(export==TRUE){
    write.csv(round(result,digits = 2),"confidence_ellipse.csv",row.names = FALSE)
    cat("Confidence ellipse:
")
    print(result)
  }
  return(result)
}

#plot bimodal elliptical confidence (calculated from A95 inversion) from dec_inc  and print results on console
ellips_plot <- function(DI,lat=0,long=0, plot=TRUE, on_plot=TRUE, col_d="red",col_u="white",col_l="black",symbol="c", text=FALSE,export=TRUE,save=FALSE,name="ellipse"){
  #degrees to radians and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  #functions converting inc(x) and dec(y) into equal area
  a2cx <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*sin(d2r(y))}
  a2cy <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*cos(d2r(y))}
  data <- DI
  #cut lines with empty cells
  data <- na.omit(data)
  data <- data[,1:2]
  colnames(data) <- c("dec", "inc")
  #directions in Cartesian coordinates
  data$x <- cos(d2r(data$dec))*cos(d2r(data$inc))
  data$y <- sin(d2r(data$dec))*cos(d2r(data$inc))
  data$z <- sin(d2r(data$inc))
  #averaged Cartesian coordinates
  x_av <- mean(data$x)
  y_av <- mean(data$y)
  z_av <- mean(data$z)
  #elements of the distribution matrix
  T_elements <- c(sum((data$x)*(data$x)),sum(data$x*data$y),sum(data$x*data$z),
                  sum(data$y*data$x),sum(data$y*data$y),sum(data$y*data$z),
                  sum(data$z*data$x),sum(data$z*data$y),sum(data$z*data$z))
  #distribution matrix
  T <- matrix(T_elements,nrow=3, byrow=TRUE)
  #calculate and copy eigenvalues and vectors
  T_e <- eigen(T,symmetric = TRUE)
  T_vec <- T_e$vectors
  T_val <- T_e$values
  #calculate dec inc of max variance
  V1inc <- r2d(asin(T_vec[3,1]/(sqrt((T_vec[1,1]^2)+(T_vec[2,1]^2)+(T_vec[3,1]^2)))))
  V1dec <- r2d(atan2(T_vec[2,1],T_vec[1,1]))
  V1dec <- V1dec%%360
  #next  calculates difference between dec_inc and average
  data$Dec_aver <- rep(V1dec)
  data$Inc_aver <- rep(V1inc)
  data$delta <- abs(data$dec-data$Dec_aver)
  data$diff <- r2d(acos((sin(d2r(data$inc))*sin(d2r(data$Inc_aver)))+
                          (cos(d2r(data$inc))*cos(d2r(data$Inc_aver))*cos(d2r(data$delta)))))
  #Isolate modes
  if(any(data$diff<=90)){
    mode1 <- as.data.frame(data$dec[data$diff<=90])
    mode1$inc <- data$inc[data$diff<=90]
    colnames(mode1) <- c("dec","inc")
  }
  if(any(data$diff>90)){
    mode2 <- as.data.frame(data$dec[data$diff>90])
    mode2$inc <- data$inc[data$diff>90]
    colnames(mode2) <- c("dec","inc")
  }
  #calculate ellipses
  if(exists("mode1")==TRUE) {ellips_M1 <- ellips_DI(mode1, lat=lat, long=long)}
  if(exists("mode2")==TRUE) {ellips_M2 <- ellips_DI(mode2, lat=lat, long=long)}

  if(plot==TRUE){
    if(on_plot==FALSE){
      par(fig=c(0,1,0,1), new= FALSE)
      plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1,
           xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
      equalarea()
    }
    if(exists("mode1")==TRUE){generate_ellips(ellips_M1[1,1],ellips_M1[1,2],ellips_M1[1,3],ellips_M1[1,4],
                                              on_plot = TRUE,symbol=symbol, col_d = col_d,
                                              col_u=col_u,col_l=col_l)}
    if(exists("mode2")==TRUE){generate_ellips(ellips_M2[1,1],ellips_M2[1,2],ellips_M2[1,3],ellips_M2[1,4],
                                              on_plot = TRUE,symbol=symbol, col_d = col_d,
                                              col_u=col_u,col_l=col_l)}
  }
  data_M12 <- common_DI(data)
  ellips_M12 <- ellips_DI(data_M12,lat=lat, long=long)
  #plot text with results
  N <- ellips_M12[1,5]
  Dec <- round(ellips_M12[1,1],digits=2)
  Inc <- round(ellips_M12[1,2],digits=2)
  Delta_dec <- round(ellips_M12[1,3],digits=2)
  Delta_inc <- round(ellips_M12[1,4],digits=2)

  if(any(data$diff<=90)) {
    cat("Ellipse Mode 1:
")
    print(round(ellips_M1, digits=2), row.names = FALSE)
    if(export==TRUE){write.csv(round(ellips_M1, digits=2),paste(name,"_mode_1.csv"), row.names = FALSE)}
  }
  if(any(data$diff>90)) {
    cat("Ellipse Mode 2:
")
    print(round(ellips_M2,digits=2), row.names = FALSE)
    if(export==TRUE){write.csv((round(ellips_M2,digits=2)),paste(name,"_mode_2.csv"), row.names = FALSE)}
  }
  if(any(data$diff>90)) {
    cat("Ellipse common mode:
")
    print(round(ellips_M12, digits=2), row.names = FALSE)
    if(export==TRUE){write.csv((round(ellips_M12, digits=2)),paste(name,"_mode_1&2.csv"), row.names = FALSE)}
  }
  if (text==TRUE){
    #plot text if true
    par(fig=c(0,1,0,1), new=TRUE)
    text <- paste("N: ",N,"
Dec: ", Dec,"
Inc: ", Inc,"
∆ dec: ", Delta_dec,"
∆ inc: ", Delta_inc)
    plot(NA, xlim=c(0,1), ylim=c(0,1),
         xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)

    text(x=0.79, y=0.02,pos=4,text, cex= 0.85)
    cat("\nDo not attempt to plot other directions or Fisher mean on the same diagram if text option is set TRUE.\n")
  }

  if(save==TRUE){save_pdf(name = paste(name,".pdf"),width = 8,height = 8)}
}

#function plotting equal area net
equalarea <- function(title="") {
  #functions converting degree and radians
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  #functions converting inc(x) and dec(y) into equal area
  a2cx <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*sin(d2r(y))}
  a2cy <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*cos(d2r(y))}
  #plot external circle
  frame_dec = 0:360
  frame_inc=rep(0,length(frame_dec))
  x = a2cx(frame_inc,frame_dec)
  y = a2cy(frame_inc,frame_dec)
  plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1,
       xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
  lines(x, y, col = "black")
  #plot 30deg and 60deg circles
  tdegD <- 0:360
  tdegI <- rep(30, length(tdegD))
  sdegI <- rep(60, length(tdegD))
  tdegx <- a2cx(tdegI,tdegD)
  tdegy <- a2cy(tdegI,tdegD)
  sdegx <- a2cx(sdegI,tdegD)
  sdegy <- a2cy(sdegI,tdegD)
  points(x=tdegx,y=tdegy,type="l", col="gray")
  points(x=sdegx,y=sdegy,type="l", col="gray")
  #plot lines every 45deg
  xd1_1 <- a2cx(0,225)
  xd1_2 <- a2cx(0,45)
  yd1_1 <- a2cy(0,225)
  yd1_2 <- a2cy(0,45)
  xd2_1 <- a2cx(0,-45)
  xd2_2 <- a2cx(0,135)
  yd2_1 <- a2cy(0,-45)
  yd2_2 <- a2cy(0,135)
  lines(c(-1, 1), c(0, 0), col = "gray")
  lines(c(0, 0), c(-1, 1), col = "gray")
  lines(c(xd1_1, xd1_2), c(yd1_1, yd1_2), col = "gray")
  lines(c(xd2_1, xd2_2), c(yd2_1, yd2_2), col = "gray")
  title(xlab = title, line=0.2, cex=0.1)
}

#Function that correct inclination shallowing after tk03.GAD model
ffind_boot <- function(DI,confidence=95,nb=1000, f_increment=0.01,export=TRUE,return=TRUE, name="Unflattened_dirs") {
  data <- DI[,1:2]
  data <- na.omit(data)
  N <- nrow(data)
  colnames(data) <- c("dec", "inc")
  #calculate E-I of real data
  Inc_E_R <- inc_E_finder(data)
  Inc_E_R$V1inc <- abs(Inc_E_R$V1inc)
  Inc <- round(Inc_E_R$V1inc, digits=1)
  Ecut <- round(Inc_E_R$E, digits=2)
  Edec <- round(Inc_E_R$DV1V2, digits=1)
  cat("Calculating precise inclination flattening of raw data.

")

  #calculate E-I correction sequence of real data
  Seq_I_E_R <- ffind(data,f_inc = 0.0005)
  colnames(Seq_I_E_R) <- c("V1inc","E","DV1V2","f")
  alert <- ifelse(length(Seq_I_E_R$V1inc)==1,"y","n")
  Ffinal <- round(Seq_I_E_R[length(Seq_I_E_R$f),4], digits=2)
  Inc_f <- round(Seq_I_E_R[length(Seq_I_E_R$V1inc),1], digits=1)
  Efinal <- round(Seq_I_E_R[length(Seq_I_E_R$E),2], digits=2)
  Edec_f <- round(Seq_I_E_R[length(Seq_I_E_R$DV1V2),3], digits=1)
  f <- min(Seq_I_E_R$f)
  if(alert=="y") f <- 1
  unf_data <- unflat_DI(data,f)
  colnames(unf_data) <- c("dec","inc")

  #plot frame
  par(fig=c(0,0.7,0,1), new= FALSE)
  plot(NA, xlim= c(0,90), ylim= c(1,3.5), xaxp= c(0,90, 9),
       xlab="Inclination (°)", ylab="Elongation")

  #plot tk03.GAD model E-I
  x <- 0:90
  y <- tk03(x)
  points(x=x, y= y, type= "l", col="blue", lwd=3)

  #Plot real data E-I and correction of real data
  points(x=Seq_I_E_R$V1inc, y= Seq_I_E_R$E, type= "l", col="red", lwd=3)
  points(x=Inc_E_R$V1inc,y=Inc_E_R$E,pch=21,
         col="black", bg="blue", cex=1.2)
  points(x=Inc_f,y=Efinal,
         pch=21, col="black", bg="red", cex=1.2)

  text <- paste("N:", N, "
Inc:", Inc,"
E:", Ecut,"
Edec:",Edec)

  text2 <- paste("f:", Ffinal, "
Inc_Unfl:", Inc_f, "
E_Unfl:", Efinal, "
Edec_Unfl:", Edec_f)

  text(x=0, y=3.2,pos=4,text, cex= 0.8)
  text(x=20, y=3.2, pos=4, text2, cex=0.8)
  if (alert=="y"){
    text3 <- "Distribution not flattened"
    text(x=0, y=3, pos=4,text3,cex=1)
  }

  #plot equal area data
  par(fig=c(0.55,1,0,0.6), new=TRUE)
  plot_DI(data,single_mode = TRUE,title = "Original directions")

  #plot equal area unflat data
  par(fig=c(0.55,1,0.4,1), new=TRUE)
  plot_DI(unf_data,single_mode = TRUE, col_d ="red",col_u = "pink",
          title="Unflattened directions")

  #create files for initial and final readings for the histograms
  init_E_I <- as.data.frame(matrix(ncol=3,nrow=0))
  final_E_I <- as.data.frame(matrix(ncol=3,nrow=0))
  colnames(init_E_I) <- c("Inc","E","E_dec")
  colnames(final_E_I)<- c("Inc","E","E_dec")
  if(alert=="y") par(fig=c(0,1,0,1))
  if(alert=="y") stop("
DISTRIBUTION NOT FLATTENED.")


  cat("Bootstrapping.
Simulation ends when", nb, "valid pseudosamples are saved.

")
  n <- 0
  par(fig=c(0,0.7,0,1), new=TRUE)
  plot(NA, xlim= c(0,90), ylim= c(1,3.5), xaxt="n",yaxt="n",
       xlab="", ylab="", axes=FALSE)
  repeat {
    n <- n+1
    Seq_I_E_B <- as.data.frame(matrix(ncol=3,nrow=0))
    dataprov <- boots_DI(data)
    Seq_I_E_B <- ffind(dataprov, f_inc = f_increment)
    #plot bootstrapped lines
    points(x=Seq_I_E_B$V1inc, y= Seq_I_E_B$E,
           type= "l", col=rgb(1, 0, 0, 0.15), lwd=1)
    i_E_I <- Seq_I_E_B[1,]
    f_E_I <- Seq_I_E_B[nrow(Seq_I_E_B),]
    colnames(i_E_I) <- c("Inc","E","E_dec")
    colnames(f_E_I)<- c("Inc","E","E_dec")

    #isolate initial and final readings for histograms
    init_E_I <- rbind(init_E_I,i_E_I)
    final_E_I <- rbind(final_E_I,f_E_I)
    init_E_I <- na.omit(init_E_I)
    final_E_I <- na.omit(final_E_I)
    if(((n%%50)==0)==TRUE) {
      cat(paste(n,"simulations done and",(nrow(final_E_I)),"pseudosamples saved
"))
    }

    if(nrow(final_E_I)==nb) {
      cat(paste("Saved",(length(final_E_I[,1])), "pseudosamples after", n,"simulations
"))
      break
    }
  }
  #replot real data with different color
  points(x=Seq_I_E_R$V1inc, y= Seq_I_E_R$E, type= "l", col="yellow", lwd=3)
  points(x=Inc_E_R$V1inc,y=Inc_E_R$E,pch=21,
         col="black", bg="blue", cex=1.2)
  points(x=x, y= y, type= "l", col="blue", lwd=3)
  points(x=Inc_f,y=Efinal,
         pch=21, col="black", bg="red", cex=1.2)


  #replot results in case covered by boostrapps
  text(x=0, y=3.2,pos=4,text, cex= 0.8)
  text(x=20, y=3.2, pos=4, text2, cex=0.8)

  colnames(init_E_I) <- c("Inc","E","E_dec","f")
  colnames(final_E_I)<- c("Inc","E","E_dec","f")
  final_E_I <- final_E_I[order(final_E_I$Inc),] #order final results by inclination
  final_E_Ibk <- final_E_I

  conf <- confidence/100
  num <- round((nb*(1-conf))/2,digits=0)
  Lconf <- num
  Uconf <- nb-num
  final_E_I <- final_E_I[Lconf:Uconf,]    #cut bootstrapped results for 95% confidence

  #draw two lines for 95% confidence margin
  arrows(x0=final_E_I[1,1],x1=final_E_I[1,1],
         y0=1,y1=final_E_I[1,2], length = 0,lty=2)

  arrows(x0=final_E_I[length(final_E_I$Inc),1],
         x1=final_E_I[length(final_E_I$Inc),1],
         y0= 1, y1= final_E_I[length(final_E_I$Inc),2],
         length = 0,lty=2)
  Inc_l95 <- round(final_E_I[1,1], digits= 1)
  Inc_u95 <- round(final_E_I[length(final_E_I$Inc),1], digits=1)

  text(x=final_E_I[1,1], y=1, pos=2, Inc_l95)
  text(x=final_E_I[length(final_E_I$Inc),1], y=1, pos=4, Inc_u95)


  #plot histogram of E_declination with respect V1 before and after correction
  par(fig=c(0.35,0.69,0.56,0.99), new=TRUE)
  hist(init_E_I$E_dec, xlim=c(-90,90), breaks= 90,
       axes=FALSE,xlab="",ylab="",col="blue", border="blue", main="")
  #plot lables closer than standard to axes
  title(xlab = "Edec(°)", line=1.9, cex=0.2)
  title(ylab = "Frequency", line=1.9,cex=0.2)
  #after
  par(fig=c(0.35,0.69,0.56,0.99), new=TRUE)
  hist(final_E_Ibk$E_dec, xlim=c(-90,90), xaxp=c(-90,90,4),
       breaks = 90, xlab = "", ylab = "",
       main="", cex.axis=0.8,col="red",border ="red")

  #recalculate boostrtapped confidence for Edec for plotting confidence margin
  final_E_I_Edec <- final_E_Ibk
  colnames(final_E_I_Edec) <- c("Inc","E","E_dec")
  final_E_I_Edec <- final_E_I_Edec[order(final_E_I_Edec$E_dec),]
  final_E_I_Edec <- final_E_I_Edec[Lconf:Uconf,]
  abline(v=final_E_I_Edec[1,3],lwd=1, lty=2)
  abline(v=final_E_I_Edec[length(final_E_I_Edec[,3]),3],lwd=1, lty=2)

  #plot histogram of inclination before and after correction
  par(fig=c(0.35,0.69,0.25,0.68), new=TRUE)
  hist(init_E_I$Inc,xlim=c(0,90), breaks=90,
       axes=FALSE,xlab="",ylab="",col="blue",border="blue", main="")
  par(fig=c(0.35,0.69,0.25,0.68), new=TRUE)
  hist(final_E_Ibk$Inc, xlim=c(0,90), xaxp=c(0,90,6),
       breaks=90,xlab = "", ylab = "",
       main="",cex.axis=0.8,col="red",border ="red")
  abline(v=final_E_I[1,1],lwd=1,lty=2)
  abline(v=final_E_I[length(final_E_I[,1]),1],lwd=1,lty=2)
  title(xlab = "Inc(°)", line=1.9, cex=0.2)
  title(ylab = "Frequency", line=1.9, cex=0.2)

  stat <- as.data.frame(matrix(ncol=1,nrow=1))
  colnames(stat) <- c("N")
  stat$N <- N
  stat$Inc <- Inc
  stat$E <- Ecut
  stat$Edec <- Edec
  stat$f <- round(f,digits=2)
  stat$Inc_unfl <- Inc_f
  stat$Low_inc <- Inc_l95
  stat$Hign_inc <- Inc_u95
  stat$E_unfl <- Efinal
  stat$Edec_unfl <- Edec_f
  stat$Edec_low <- round(final_E_I_Edec[1,3],digits = 1)
  stat$Edec_high <- round(final_E_I_Edec[length(final_E_I_Edec[,3]),3],digits = 1)

  par(fig=c(0,1,0,1))
  if(export==TRUE){
    cat("
Unflattened directions and statistics saved as .csv file
Graph saved as",paste(name,".pdf"),"

")
    write.csv(round(unf_data,digits = 2),paste(name,".csv"), row.names = FALSE)
    write.csv(stat,paste(name,"_statistic.csv"), row.names=FALSE)
    save_pdf(name=paste(name,".pdf"))
  }
  if(return==TRUE){return(unf_data)}
}

#flattening factor finder function from Dec Inc, results in Inc, E, and E declination
ffind <-function(DI, f_inc=0.005) {
  data <- DI[,1:2]
  data <- na.omit(data)
  colnames(data) <- c("dec", "inc")
  inc_E_seq <- as.data.frame(matrix(ncol=4,nrow=0))
  f <- 1.00
  #loop that stops only when tk03 line is crossed from below
  repeat {
    f <- f-f_inc
    data_unfl <- unflat_DI(data,f)
    inc_E_prov <- inc_E_finder(data_unfl)
    inc_E_prov$V1inc <- abs(inc_E_prov$V1inc)
    inc_E_prov$f <- f
    inc_E_seq <- rbind(inc_E_seq,inc_E_prov)
    if(nrow(inc_E_seq)>1){
      E_lim_low <- round(tk03(inc_E_seq[nrow(inc_E_seq)-1,1]), digits = 6)
      E_lim_high <- round(tk03(inc_E_seq[nrow(inc_E_seq),1]),digits = 6)
      if(round(max(inc_E_seq[nrow(inc_E_seq),1],digit=1)>89.5)) break
      if(inc_E_seq[nrow(inc_E_seq)-1,2]< E_lim_low && inc_E_seq[nrow(inc_E_seq),2]>=E_lim_high) break
    }
  }
  #next return data only when E goes below tk03 line
  Emin <- min(inc_E_seq$E)
  Imin <- inc_E_seq$V1inc[inc_E_seq$E==min(inc_E_seq$E)]
  Eminlim <- tk03(Imin)
  if(round(max(inc_E_seq$V1inc),digit=1)>89.5) {return(as.data.frame(t(c(NA, NA, NA, NA))))}
  if(Emin>Eminlim) {return(as.data.frame(t(c(NA, NA, NA, NA))))} else {return(inc_E_seq)}
}

#plot bimodal fisher from dec_inc and print results on console
fisher_plot <- function(DI, plot=TRUE, on_plot=TRUE,col_d="red",col_u="white",col_l="black",symbol="c",text=FALSE,export=TRUE,save=FALSE,name="Fisher_mean") {
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  data <- DI
  data <- na.omit(data)
  data <- data[,1:2]
  colnames(data) <- c("dec", "inc")
  #directions in Cartesian coordinates
  data$x <- cos(d2r(data$dec))*cos(d2r(data$inc))
  data$y <- sin(d2r(data$dec))*cos(d2r(data$inc))
  data$z <- sin(d2r(data$inc))
  #averaged Cartesian coordinates
  x_av <- mean(data$x)
  y_av <- mean(data$y)
  z_av <- mean(data$z)
  #elements of the distribution matrix
  T_elements <- c(sum((data$x)*(data$x)),sum(data$x*data$y),sum(data$x*data$z),
                  sum(data$y*data$x),sum(data$y*data$y),sum(data$y*data$z),
                  sum(data$z*data$x),sum(data$z*data$y),sum(data$z*data$z))
  #distribution matrix
  T <- matrix(T_elements,nrow=3, byrow=TRUE)
  #calculate and copy eigenvalues and vectors
  T_e <- eigen(T,symmetric = TRUE)
  T_vec <- T_e$vectors
  T_val <- T_e$values
  #calculate dec inc of max variance
  V1inc <- r2d(asin(T_vec[3,1]/(sqrt((T_vec[1,1]^2)+(T_vec[2,1]^2)+(T_vec[3,1]^2)))))
  V1dec <- r2d(atan2(T_vec[2,1],T_vec[1,1]))
  V1dec <- ifelse(V1dec<0,V1dec+360,V1dec)
  #next  calculates difference between dec_inc and average
  data$Dec_aver <- rep(V1dec)
  data$Inc_aver <- rep(V1inc)
  data$delta <- abs(data$dec-data$Dec_aver)
  data$diff <- r2d(acos((sin(d2r(data$inc))*sin(d2r(data$Inc_aver)))+
                          (cos(d2r(data$inc))*cos(d2r(data$Inc_aver))*cos(d2r(data$delta)))))
  #Isolate modes
  if(any(data$diff<=90)){
    mode1 <- as.data.frame(data$dec[data$diff<=90])
    mode1$inc <- data$inc[data$diff<=90]
    colnames(mode1) <- c("dec","inc")
  }
  if(any(data$diff>90)){
    mode2 <- as.data.frame(data$dec[data$diff>90])
    mode2$inc <- data$inc[data$diff>90]
    colnames(mode2) <- c("dec","inc")
  }
  if(exists("mode1")==TRUE) {fisher_M1 <- fisher(mode1)}
  if(exists("mode2")==TRUE) {fisher_M2 <- fisher(mode2)}
  if(plot==TRUE){
    if(on_plot==FALSE){
      par(fig=c(0,1,0,1), new= FALSE)
      plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1,
           xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
      equalarea()
      }
    if(exists("mode1")==TRUE){plot_a95(fisher_M1[1,1],fisher_M1[1,2],fisher_M1[1,3],
                                       on_plot = TRUE,symbol=symbol, col_d = col_d,
                                       col_u=col_u,col_l=col_l)}
    if(exists("mode2")==TRUE){plot_a95(fisher_M2[1,1],fisher_M2[1,2],fisher_M2[1,3],
                                       on_plot = TRUE,symbol=symbol, col_d = col_d,
                                       col_u=col_u,col_l=col_l)}
  }
  data_M12 <- common_DI(data)
  fisher_M12 <- fisher(data_M12)
  #plot text with results
  Dec <- round(fisher_M12[1,1],digits=2)
  Inc <- round(fisher_M12[1,2],digits=2)
  a <- round(fisher_M12[1,3],digits=2)
  N <- round(fisher_M12[1,4],digits=2)

  if(any(data$diff<=90)) {
    cat("fisher Mode 1:
")
print(round(fisher_M1, digits=2), row.names = FALSE)
if(export==TRUE){write.csv(round(fisher_M1, digits=2),paste(name,"_mode_1.csv"), row.names = FALSE)}
  }
  if(any(data$diff>90)) {
    cat("fisher Mode 2:
")
    print(round(fisher_M2,digits=2), row.names = FALSE)
    if(export==TRUE){write.csv((round(fisher_M2,digits=2)),paste(name,"_mode_2.csv"), row.names = FALSE)}
  }
  if(exists("fisher_M1")==TRUE | exists("fisher_M2")==TRUE) {
    cat("fisher common mode:
")
    print(round(fisher_M12, digits=2), row.names = FALSE)
    if(export==TRUE){write.csv((round(fisher_M12, digits=2)),paste(name,"_mode_1&2.csv"), row.names = FALSE)}
  }

  #plot text in figure if requested
  if (text==TRUE){
    #plot text if true
    par(fig=c(0,1,0,1), new=TRUE)
    text <- paste("N: ",N,"
Dec: ", Dec,"
Inc: ", Inc,"
a95%: ", a)
    plot(NA, xlim=c(0,1), ylim=c(0,1),
         xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)

    text(x=0.79, y=0.02,pos=4,text, cex= 0.85)
    cat("\nDo not attempt to plot other directions or Fisher mean on the same diagram if text option is set TRUE.\n")
  }

  if(save==TRUE){save_pdf(name = paste(name,".pdf"),width = 8,height = 8)}
}

#function that return fisher statistic from dec_inc
fisher <- function(DI, export=FALSE, name="fisher_mean"){
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  data <- DI[,1:2]
  data <- na.omit(data)
  colnames(data) <- c("dec", "inc")
  #add various columns to data, with conversions to rads, lat and long
  data$dec_r <- d2r(data$dec)
  data$inc_r <- d2r(data$inc)
  #long_lat in cartesian coordinates
  data$x <- cos(data$dec_r)*cos(data$inc_r)
  data$y <- sin(data$dec_r)*cos(data$inc_r)
  data$z <- sin(data$inc_r)
  #average cartesian coordinates
  X_aver <- mean(data$x)
  Y_aver <- mean(data$y)
  Z_aver <- mean(data$z)
  #sum of all values along the axes
  X_sum <- sum(data$x)
  Y_sum <- sum(data$y)
  Z_sum <- sum(data$z)
  #magnitude of average
  B <- sqrt((X_aver^2)+(Y_aver^2)+(Z_aver^2))
  #Fisher (1953) parameters R,K, a95
  N <- length(data$dec)
  R <- sqrt(X_sum^2+Y_sum^2+Z_sum^2)
  K <- (N-1)/(N-R)
  a95 <- r2d(acos(1-(((N-R)/R)*(((1/0.05)^(1/(N-1)))-1))))
  Dec_aver <- r2d(atan2(Y_aver,X_aver))
  #corrects for negative declination
  Dec_aver <- ifelse(Dec_aver<0,Dec_aver+360,Dec_aver)
  Inc_aver <- r2d(asin(Z_aver/B))
  result <- as.data.frame(matrix(ncol=6,nrow=1))
  colnames(result) <- c("dec","inc","a95","N","R","k")
  result[1,1:6] <- c(Dec_aver,Inc_aver,a95,N,R,K)
  if(export==TRUE){write.csv(round(result,digits = 2),paste(name,".csv"),row.names = FALSE)}
  return(result)
}

#flat directions with given f
flat_DI <- function(DI,f=1,export=FALSE,name="flattened_dirs") {
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  DI[,2]<- r2d(atan(tan(d2r(DI[,2]))*f))
  if(export==TRUE){write.csv(round(DI,digits = 2),paste(name,".csv"),row.names = FALSE)}
  return(DI)
}

#flip directions to antipodal
flip_DI <- function(DI,export=FALSE,name="flipped_dirs"){
  dat <- DI[,1:2]
  dat <- na.omit(dat)
  colnames(dat) <- c("dec","inc")
  dat_fl <- dat
  dat_fl$dec <- ifelse((dat$dec+180)>360,dat$dec-180,dat$dec+180)
  dat_fl$inc <- -dat$inc
  if(export==TRUE){write.csv(round(dat_fl,digits=2),paste(name, ".csv"),row.names = FALSE)}
  return(dat_fl)
}

#function that plots points on a KavrayskiyVII geographic map
geo_point <- function(S_file=FALSE,symbol="c",col="red",center=0,grid=30,A95=FALSE,fill_A=TRUE,export=TRUE,on_plot=FALSE){
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #functions converting long & lat to xy in KavrayskiyVII projection
  c2x <- function(lon,lat) {((3*d2r(lon))/2)*(sqrt((1/3)-((d2r(lat)/pi)^2)))}
  c2y <- function(lat) {d2r(lat)}

  #functions spherical (lon=x, lat=y) to Cartesian
  s2cx <- function(x,y) {cos(d2r(x))*cos(d2r(y))}
  s2cy <- function(x,y) {sin(d2r(x))*cos(d2r(y))}
  s2cz <- function(y) {sin(d2r(y))}

  #function that draw circle around point
  circle <- function(longitude,latitude,A,fill=FALSE){
    #empty new circle file
    circle <- as.data.frame(matrix(ncol=2,nrow=0))
    #loop that create A95 and rotate it around new coordinate (dec, inc)
    for (i in seq(0,360,0.5)){
      #temporary circle point
      circleP <- as.data.frame(matrix(ncol=2,nrow=1))
      S_lon <- longitude
      S_lat <- latitude
      A <- A

      #trigonometric parameters for rotation
      sbd <- -sin(d2r(S_lon))
      cbd <- cos(d2r(S_lon))
      sbi <- sin(d2r(90-S_lat))
      cbi <- cos(d2r(90-S_lat))

      # cartesian coordinates of the confidence circle
      x <- s2cx(i,(90-A))
      y <- s2cy(i,(90-A))
      z <- s2cz(90-A)

      #new rotated coordinate
      xn <- x*(sbd^2+cbd^2*cbi)+
        y*(cbd*sbd*(1-cbi))+
        z*sbi*cbd
      yn <- x*cbd*sbd*(1-cbi)+
        y*(cbd^2+sbd*sbd*cbi)-
        z*sbd*sbi
      zn <- -(x*cbd*sbi-
                y*sbi*sbd-
                z*cbi)
      #converted to spherical
      newlon <- r2d(atan2(yn,xn))
      newlon <- ifelse(newlon>180,newlon-360,newlon)
      newlon <- ifelse(newlon<(-180),newlon+360,newlon)
      newlat <- r2d(asin(zn))
      circleP[1,1:2] <- c(newlon,newlat)
      circle <- rbind(circle,circleP)
    }
    colnames(circle) <- c("lon","lat")

    #next is to divide circles into two polygons if crosses the end of map
    circle1 <- as.data.frame(matrix(nrow = 0,ncol = 2))
    colnames(circle1) <- c("lon","lat")
    #creates the breaking line
    breaker <- as.data.frame(t(c(NA,NA)))
    colnames(breaker) <- c("lon","lat")
    #when two points are on the different side of the map, based on 355° distance, it put break in between
    for(l in 2:nrow(circle)){
      if(abs(circle[l-1,1]-circle[l,1])>355){
        provv1 <- circle[l-1,]
        provv2 <- circle[l,]
        circle1 <- rbind(circle1,provv1,breaker,provv2)
      }else{circle1 <- rbind(circle1,circle[l-1,],circle[l,])}
    }

    circle1$x <- c2x(circle1$lon,circle1$lat)
    circle1$y <- c2y(circle1$lat)
    filling <- which(is.na(circle1))
    if(length(filling!=0)){
      lines(circle1$x,circle1$y, lwd=0.8)
    }
    else if(fill_A==FALSE) {lines(circle1$x,circle1$y, lwd=0.8)}
    else {polygon(circle1$x,circle1$y, lwd=0.8,col=rgb(1,0.9,0,0.25))}
  }

  #plot map
  if(on_plot==FALSE) {
    Map_KVII(grid=grid,center=center)
  }else{cat("
Double check the center meridian of the map!
")}


  #plot points from file
  if(S_file==TRUE){
    cat("Select file in .csv format")
    dat <- read.csv(file.choose())
    for(i in 1:nrow(dat)){
      S_lon <- dat[i,1]-center
      S_lon <- ifelse(S_lon<0,S_lon+360,S_lon)
      S_lon <- ifelse(S_lon>180,S_lon-360,S_lon)
      S_lat <- dat[i,2]
      if(A95==TRUE){
        circ <- dat[i,3]
        circle(S_lon,S_lat,circ,fill_A)
      }
      #select symbol
      if(symbol=="c") pch <- 21
      if(symbol=="s") pch <- 22
      if(symbol=="d") pch <- 23
      if(symbol=="t") pch <- 24

      site_x <- c2x(S_lon,S_lat)
      site_y <- c2y(S_lat)
      points(x=site_x,y=site_y,pch=pch, col="black",bg=col)
    }
  }
  if(S_file==FALSE){
    repeat{
      #plot point
      S_lon <- as.numeric(readline("Longitude -or press enter to exit-: "))
      if(is.na(S_lon)==TRUE) break
      S_lon <- S_lon-center
      S_lon <- ifelse(S_lon>180,
                      S_lon-360,S_lon)
      S_lon <- ifelse(S_lon<(-180),S_lon+360,S_lon)
      S_lat <- as.numeric(readline("Latitude: "))
      circ <- as.numeric(readline("Semi-angle of circle -press enter if no confidence angle is required-: "))
      if(is.na(circ)==FALSE){circle(S_lon,S_lat,circ,fill_A)}
      symbol <- readline("Symbol -press enter for circle-: ")
      #select symbol
      if(symbol=="") sym <- 21
      if(symbol=="c") sym <- 21
      if(symbol=="s") sym <- 22
      if(symbol=="d") sym <- 23
      if(symbol=="t") sym <- 24

      site_x <- c2x(S_lon,S_lat)
      site_y <- c2y(S_lat)
      color <- readline("Color -press enter for red-:")
      if(color=="") color <- "red"
      points(x=site_x,y=site_y,pch=sym, col="black",bg=color)
    }
  }
  if(export==TRUE){
    save_pdf(name="Map.pdf")
    cat("
Figure saved as Map.pdf
")
  }
}

#function that generate and plot confidence ellipses calculated from A95 back to directions, delta Inc > delta dec
generate_ellips <- function(D,I,delta_dec,delta_inc,col_d="red",col_u="white",col_l="black", symbol="c", on_plot=FALSE, save=FALSE, name="confidence_ellipse"){
  #degree to radians and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  #functions converting inc(x) and dec(y) into equal area
  a2cx <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*sin(d2r(y))}
  a2cy <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*cos(d2r(y))}
  #create circle with dummy bedding only for strain_DI to work
  circle <- as.data.frame(seq(0,360,2))
  circle$inc <- rep(90-delta_inc)
  colnames(circle) <- c("circ_dec","circ_inc")
  circle$dummy_az <- rep(0)
  circle$dummy_pl <- rep(0)
  #calculate parameter for adjusting delta_dec
  fol <- tan(d2r(delta_inc))/tan(d2r(delta_dec))
  #create matrix deforming circle
  M <- matrix_maker(Fol = fol,v1d = D,v1i = 0,v2d = 0,v2i = 90,v3d = ((D-90)%%360),v3i = 0, return_P=F)
  ell <- strain_DI(DIAP = circle,M = M)
  ell <- ell[,1:2]
  #uses bedding correction to place final ellipses in the right position
  ellipses <- bed_DI(ell,in_file = FALSE,bed_az = D,bed_plunge = 90-I,export = FALSE)
  colnames(ellipses) <- c("dec","inc")
  ellipses$x <- a2cx(abs(ellipses$inc),ellipses$dec)
  ellipses$y <- a2cy(abs(ellipses$inc),ellipses$dec)
  #restore screen
  par(fig=c(0,1,0,1))
  #standalone graph or on existing graph
  if (on_plot==FALSE) {
    plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1,
         xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
    equalarea()
  }
  #check for negative inclination
  inc <- I
  dec <- D
  UD <- ifelse(inc>0,"D","U")
  inc <- abs(inc)
  X <- a2cx(inc,dec)
  Y <- a2cy(inc,dec)
  if(symbol=="c") {pch <- 21}
  else if(symbol=="s") {pch <- 22}
  else if(symbol=="d") {pch <- 23}
  else if(symbol=="t") {pch <- 24}
  else{stop("Please select valid symbol. Check help for info.",call. = F)}
  if(UD=="D"){
    points(X,Y, pch=pch,cex=1.3, col="black",
           bg= col_d)
  }else{
    points(X,Y, pch=pch,cex=1.3, col="black",
           bg=col_u)
  }
  lines(ellipses$x,ellipses$y,lty=1, col=col_l, lwd=1.8)
  if(save==TRUE){save_pdf(name = paste(name,".pdf"),width = 8,height = 8)}
}

#function that take data Dec_Inc and return the average Inc, E, and E declination
inc_E_finder <- function(DI, export=FALSE, name="I_E_Edec") {
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  data <- DI[,1:2]
  data <- na.omit(data)
  colnames(data) <- c("dec", "inc")
  #directions in Cartesian coordinates
  data$x <- cos(d2r(data$dec))*cos(d2r(data$inc))
  data$y <- sin(d2r(data$dec))*cos(d2r(data$inc))
  data$z <- sin(d2r(data$inc))
  #averaged Cartesian coordinates
  x_av <- mean(data$x)
  y_av <- mean(data$y)
  z_av <- mean(data$z)
  #elements of the distribution matrix
  T_elements <- c(sum((data$x)*(data$x)),sum(data$x*data$y),sum(data$x*data$z),
                  sum(data$y*data$x),sum(data$y*data$y),sum(data$y*data$z),
                  sum(data$z*data$x),sum(data$z*data$y),sum(data$z*data$z))
  #distribution matrix
  T <- matrix(T_elements,nrow=3, byrow=TRUE)
  #calculate and copy eigenvalues and vectors
  T_e <- eigen(T)
  T_vec <- T_e$vectors
  T_val <- T_e$values
  #calculate dec inc of max variance
  V1inc <- r2d(asin(T_vec[3,1]/(sqrt((T_vec[1,1]^2)+(T_vec[2,1]^2)+(T_vec[3,1]^2)))))
  V1dec <- r2d(atan2(T_vec[2,1],T_vec[1,1]))
  V1dec <- ifelse(V1dec<0,V1dec+360,V1dec)

  #force V1 to positive
  V1inc <- abs(V1inc)
  V1dec <- ifelse(V1inc<0,(V1dec+180)%%360,V1dec)

  V2dec <- r2d(atan2(T_vec[2,2],T_vec[1,2]))
  V2dec <- ifelse(V2dec<0,V2dec+360,V2dec)
  V2inc <- r2d(asin(T_vec[3,2]/(sqrt((T_vec[1,2]^2)+(T_vec[2,2]^2)+(T_vec[3,2]^2)))))

  #Calculate difference between V1 and V2 to have the declination of V2 with respect to V1
  DV1V2 <- V1dec-V2dec
  DV1V2 <- ifelse(DV1V2<0,DV1V2+360,DV1V2)
  DV1V2 <- ifelse(DV1V2>90,ifelse(DV1V2<270,DV1V2-180,DV1V2),DV1V2)
  DV1V2 <- ifelse(DV1V2>270,DV1V2-360,DV1V2)

  E <- T_val[2]/T_val[3]

  inc_E <- as.data.frame(cbind(V1inc,E,DV1V2))
  #export result if requested
  if(export==TRUE){write.csv(round(inc_E,digits = 2),paste(name,".csv"),row.names = FALSE)}
  return(inc_E)
}

#Arason and Levi(2010) inclination only calculation
#adepted from the original fortran source code ARALEV available at http://hergilsey.is/arason/paleomag/aralev.txt
inc_only <- function(DI,dec=TRUE, print=TRUE,export=TRUE, name="Inclination_only",return=TRUE, arith_stat=FALSE) {
  #The Arason-Levi MLE Iteration Formula 1
  AL1 <- function(th, n, the, ak) {
    dr <- 0.0174532925199433  # Degrees to radians (pi/180)
    s <- 0
    c <- 0

    for (i in 1:n) {
      x <- ak * sin(the * dr) * sin(th[i] * dr)
      bessel_result <- bessel(x)
      bi1i0 <- bessel_result[3]

      s <- s + sin(th[i] * dr) * bi1i0
      c <- c + cos(th[i] * dr)
    }

    AL1 <- atan2(s, c) / dr
    AL1 <- ifelse(AL1 < 0.000001, 0.000001, AL1)
    AL1 <- ifelse(AL1 > 179.999999, 179.999999, AL1)

    return(AL1)
  }
  #The Arason-Levi MLE Iteration Formula 2
  AL2 <- function(th, n, the, ak) {
    dr <- 0.0174532925199433  # Degrees to radians (pi/180)
    dn <- n

    s <- 0
    c <- 0
    for (i in 1:n) {
      x <- ak * sin(the * dr) * sin(th[i] * dr)
      Btemp <- bessel(x)
      bi1i0 <- Btemp[3]

      s <- s + sin(th[i] * dr) * bi1i0
      c <- c + cos(th[i] * dr)
    }

    x <- dn * (1 / tanh(ak)) - cos(the * dr) * c - sin(the * dr) * s
    AL2 <- 1e10
    if (x / dn > 1e-10) {
      AL2 <- dn / x
    }
    if (AL2 < 1e-06) {
      AL2 <- 1e-06
    }

    return(AL2)
  }
  #Evaluation of the Hyperbolic Bessel functions I0(x), I1(x) and their ratio I1(x)/I0(x).
  bessel <- function(x) {
    p <- c(1.0, 3.5156229, 3.0899424, 1.2067492, 0.2659732, 0.360768e-1, 0.45813e-2)
    q <- c(0.39894228, 0.1328592e-1, 0.225319e-2, -0.157565e-2, 0.916281e-2, -0.2057706e-1, 0.2635537e-1, -0.1647633e-1, 0.392377e-2)
    u <- c(0.5, 0.87890594, 0.51498869, 0.15084934, 0.2658733e-1, 0.301532e-2, 0.32411e-3)
    v <- c(0.39894228, -0.3988024e-1, -0.362018e-2, 0.163801e-2, -0.1031555e-1, 0.2282967e-1, -0.2895312e-1, 0.1787654e-1, -0.420059e-2)

    if (abs(x) < 3.75) {
      t <- (x / 3.75)^2
      b0 <- p[1] + t * (p[2] + t * (p[3] + t * (p[4] + t * (p[5] + t * (p[6] + t * p[7])))))
      b1 <- x * (u[1] + t * (u[2] + t * (u[3] + t * (u[4] + t * (u[5] + t * (u[6] + t * u[7]))))))
      bi0e <- b0 / exp(abs(x))
      bi1e <- b1 / exp(abs(x))
      bi1i0 <- b1 / b0
    } else {
      t <- 3.75 / abs(x)
      b0 <- q[1] + t * (q[2] + t * (q[3] + t * (q[4] + t * (q[5] + t * (q[6] + t * (q[7] + t * (q[8] + t * q[9])))))))
      b1 <- v[1] + t * (v[2] + t * (v[3] + t * (v[4] + t * (v[5] + t * (v[6] + t * (v[7] + t * (v[8] + t * v[9])))))))
      if (x < 0) b1 <- -b1
      bi0e <- b0 / sqrt(abs(x))
      bi1e <- b1 / sqrt(abs(x))
      bi1i0 <- b1 / b0
    }

    return(c(bi0e, bi1e, bi1i0))
  }
  #Evaluation of the Hyperbolic Cotangens function coth(x)
  coth <- function(x) {
    if (x == 0) {
      return(0)
    }

    t <- abs(x)

    if (t < 0.001) {
      return(1 / t + t / 3 - t^3 / 45 + 2 * t^5 / 945)
    } else if (t <= 15) {
      ep <- exp(t)
      em <- exp(-t)
      return((ep + em) / (ep - em))
    } else {
      return(1)
    }

    if (x < 0) {
      return(-coth)
    }
  }
  #Evaluation of the Log-Likelihood function for inclination-only data.
  xlik <- function(th, n, the, ak) {
    dr <- 0.0174532925199433         # Degrees to radians (pi/180)
    pi <- 180 * dr
    dn <- n

    # Illegal use
    if (n < 1) {
      cat("ERROR: Data missing in xlik\n")
      return(-1e10)
    }

    if (n > 10000) {
      cat("ERROR: Too small dimension in xlik\n")
      return(-1e10)
    }

    # Uncomment the following lines if you want to check the range of ak
    # if (ak < 0) {
    #   cat("ERROR: Out of range in xlik\n")
    #   return(-1e10)
    # }

    # A1(k) = N ln(k) - N ln(sinh k) - N ln(2)
    a1 <- 0

    if (ak >= 0 && ak < 0.01) {
      q <- -ak * (1 - ak * (2 / 3 - ak * (1 / 3 - ak * (2 / 15 - ak * (8 / 45)))))
      a1 <- dn * (-log(2) - log(1 + q) - ak)
    } else if (ak >= 0.01 && ak <= 15) {
      a1 <- dn * (log(ak) - log(1 - exp(-2 * ak)) - ak)
    } else {
      a1 <- dn * (log(ak) - ak)
    }

    # A2(k,t,ti) = Sum(k cos t cos ti) + Sum(ln(BessIo(k sin t sin ti)))
    a2 <- 0

    for (i in 1:n) {
      x <- ak * sin(the * dr) * sin(th[i] * dr)
      bi0e <- bessel(x)[[1]]
      bi1e <- bessel(x)[[2]]
      bi1i0 <- bessel(x)[[3]]
      a2 <- a2 + ak * cos((th[i] - the) * dr) + log(bi0e)
    }

    # A3(ti) = Sum( ln(sin(ti)) ), Note: 0.000001 < ti < 179.999999
    a3 <- 0

    for (i in 1:n) {
      x <- th[i]
      if (x < 1e-6) x <- 1e-6
      if (x > 179.999999) x <- 179.999999
      a3 <- a3 + log(sin(x * dr))
    }

    # The log-likelihood function
    xlik <- a1 + a2 + a3

    return(xlik)
  }
  #Calculation of the arithmetic mean of inclination-only data
  armean <- function(xinc) {
    dr <- 0.01745329252         # Degrees to radians (pi/180)
    t63max <- 105.070062145     # 63 % of a sphere.
    a95max <- 154.158067237     # 95 % of a sphere.
    dn <- length(xinc)

    s <- sum(xinc)
    s2 <- sum(xinc^2)

    ainc <- s / dn

    sd <- 0
    ak <- -1

    if (dn > 1) {
      sd <- sqrt((s2 - s^2 / dn) / (dn - 1))
      ak <- (dn - 1) / ((s2 - s^2 / dn) * dr^2)
    }

    nf <- dn - 1

    tval_63 <- qt(0.63, df = nf)
    t63 <- tval_63 * sd

    tval_95 <- qt(0.95, df = nf)
    a95 <- tval_95 * sd / sqrt(dn)

    result <- as.data.frame(matrix(ncol=5, nrow=1))
    result[1] <- dn
    result[2] <- round(ainc, digits=2)
    result[3] <- round(ak, digits=2)
    result[4] <- round(t63,digits = 2)
    result[5] <- round(a95, digits = 2)
    colnames(result) <- c("N","Inc","Precision","Angular st.dev(63%)","a95")

    return(result)
  }

  #isolate inclination if file contains also declination
  DI <- na.omit(DI)
  if(dec==TRUE) {DI <- DI[,1:2]
  }else if(dec==F) {DI <- DI[,1]}
  if(dec==TRUE) {
    colnames(DI) <- c("dec","inc")
    xinc <- DI$inc
  }else{
    #convert the imported file in a list of numbers
    inc <- as.list(inc)
    xinc <- inc[[1]]
  }
  #generate only arithmetic statistic if requested
  if(arith_stat==TRUE){
    result <- armean(xinc)
    #print if request
    if(export==TRUE){write.csv(result,paste(name,".csv"), row.names = F)}

    if(print==TRUE){
      cat("Arithmetic average inclination result:

N:", result[1,1],"
Inclination:", result[1,2],"
Precision:",result[1,3],"
St.Dev_63:", result[1,4],"
alpha_95:", result[1,5],"

")
    }
    if(return==TRUE) {return(result)}
  }else{
    # main routine
    # Degrees to radians (pi/180)
    dr <- 0.0174532925199433
    # 63 % of a sphere
    t63max <- 105.070062145
    # 95 % of a sphere
    a95max <- 154.158067237
    n <- length(xinc)
    th <- numeric(n)
    dn <- n
    ierr <- 1
    #sdata checking and warnings
    if (n == 1) {
      stop("Only one observed inclination\n", call.=F)
    }
    if (n > 10000) {
      assign("inc_warn","Too many directions", envir = .GlobalEnv)
      stop("Too many directions, max=10000\n", call.=F)
    }
    if(length(unique(xinc))==1){
      assign("inc_warn","Directions are all identical", envir = .GlobalEnv)
      stop("Directions are all identical\n", call.=F)
    }
    if (any(xinc>90) | any(xinc<(-90))) {
      assign("inc_warn","Inclination must be between -90 and 90", envir = .GlobalEnv)
      stop("Inclination must be between -90 and 90\n", call.=F)
    }
    #file with co-incl
    for (i in 1:n) {
      th[i] <- 90 - xinc[i]
    }

    s <- sum(th)
    s2 <- sum(th^2)
    c <- sum(cos(th * dr)) / dn
    # initial theta guess (17)
    rt <- s / dn
    x <- (s2 - s^2 / dn) * dr^2
    rk <- ifelse(x / (dn - 1) > 1e-10, (dn - 1) / x, 1e10)
    rt1 <- rt
    rk1 <- rk
    ie1 <- 0

    the1 <- rt
    akap1 <- rk
    for (j in 1:10000) {
      rt <- AL1(th, n, rt, rk)
      rk <- AL2(th, n, rt, rk)
      dt <- abs(rt - the1)
      dk <- abs((rk - akap1) / rk)
      if (j > 10 && dt < 1e-6 && dk < 1e-6) break
      the1 <- rt
      akap1 <- rk
    }
    #likelihood for theta e k
    #ie1 <- 0
    the1 <- rt
    akap1 <- rk
    xl1 <- xlik(th, n, rt, rk)

    #likelihood for theta=0
    rt <- 0
    rk <- rk1
    akap2 <- rk
    ie2 <- 0
    for (j in 1:10000) {
      x <- coth(rk) - c
      if (x > 1e-10) {
        rk <- 1 / x
      } else {
        rk <- 1e10
      }
      dk <- abs((rk - akap2) / rk)
      if (j > 4 && dk < 1e-6) break
      if (rk < 1e-6) break
      akap2 <- rk
    }
    ie2 <- 1
    the2 <- 0
    akap2 <- rk
    xl2 <- xlik(th, n, rt, rk)

    #likelihood for theta=180
    rt <- 180
    rk <- rk1
    akap3 <- rk
    ie3 <- 0
    for (j in 1:10000) {
      x <- coth(rk) + c
      if (x > 1e-10) {
        rk <- 1 / x
      } else {
        rk <- 1e10
      }
      dk <- abs((rk - akap3) / rk)
      if (j > 4 && dk < 1e-6) break
      if (rk < 1e-6) break
      akap3 <- rk
    }
    ie3 <- 1
    the3 <- 180
    akap3 <- rk
    xl3 <- xlik(th, n, rt, rk)

    #likelihood for k=0
    rt <- 90
    rk <- 0
    the4 <- rt
    akap4 <- rk
    xl4 <- xlik(th, n, rt, rk)

    isol <- 1
    ierr <- ie1
    if (xl2 > xl1) {
      the1 <- the2
      akap1 <- akap2
      xl1 <- xl2
      isol <- 2
      ierr <- 1
    }

    #compares solutions
    if (xl3 > xl1) {
      the1 <- the3
      akap1 <- akap3
      xl1 <- xl3
      isol <- 3
      ierr <- 1
    }

    if (xl4 > xl1) {
      the1 <- the4
      akap1 <- akap4
      xl1 <- xl4
      isol <- 4
      ierr <- 0
    }

    ainc <- 90 - the1
    ak <- akap1
    if (ierr != 0) {
      assign("inc_warn","Convergence problems", envir = .GlobalEnv)
      cat("Convergence problems\n")
    }

    #Test of robustness with 16 surrounding points
    for (i in 1:16) {
      x <- i
      rt <- the1 + 0.01 * cos(22.5 * x * dr)
      if (rt < 0 || rt > 180) break
      rk <- akap1 * (1 + 0.001 * sin(22.5 * x * dr))
      xl <- xlik(th, n, rt, rk)
      if (xl > xl1) {
        ierr <- ierr + 2
        assign("inc_warn","Robustness failure", envir = .GlobalEnv)
        cat("Robustness failure\n")
      }
    }

    if (akap1 >= 20) {
      co <- 1 + log(1 - 0.63) / akap1
    } else if (akap1 > 0.1 && akap1 < 20) {
      co <- 1 + log(1 - 0.63 * (1 - exp(-2 * akap1))) / akap1
    } else if (akap1 <= 0.1) {
      co <- -0.26 + 0.4662 * akap1
    }

    #calculate confidences
    t63 <- 90 - (90*sign(co))
    if (abs(co) < 1) {
      t63 <- 90 - atan(co / sqrt(1 - co^2)) / dr
    }
    if (t63 > t63max) {
      t63 <- t63max
    }

    co <- 1 - (dn - 1) * (20^(1 / (dn - 1)) - 1) / (dn * (akap1 - 1) + 1)
    a95 <- 90 - (90*sign(co))
    if (abs(co) < 1) {
      a95 <- 90 - atan(co / sqrt(1 - co^2)) / dr
    }
    if (a95 > a95max) {
      a95 <- a95max
    }

    #calculates arith mean
    ari_mean <- armean(xinc)

    #compile result file
    result <- as.data.frame(matrix(ncol=5, nrow=1))
    result[1] <- n
    result[2] <- round(ainc, digits=2)
    result[3] <- round(ak, digits=2)
    result[4] <- round(t63,digits = 2)
    result[5] <- round(a95, digits = 2)
    result[6] <- round(ari_mean[1,2],digits = 2)
    colnames(result) <- c("N","Inc","Precision","Angular st.dev(63%)","a95","Aritm. mean")

    #print if request
    if(export==TRUE){write.csv(result,paste(name,".csv"), row.names = F)}

    if(print==TRUE){
      cat("Arason-Levi inclination only result:

N:", result[1,1],"
Inclination:", result[1,2],"
Precision:",result[1,3],"
St.Dev_63:", result[1,4],"
alpha_95:", result[1,5],"
Aritm. mean:",result[1,6],"

")
    }
  }
  if(return==TRUE) {return(result)}
}

#plot equal area of Arason and Levi(2010) inclination only calculation
inc_plot <- function(DI,dec=TRUE,plot=TRUE,bimodal=FALSE,on_plot=TRUE, col="black", print=TRUE,export=TRUE, save=TRUE,name="Inc_only", arith_stat=FALSE,Shiny=FALSE){
  #import dplyr for filter_ALL
  library(dplyr)
  #functions converting degree and radians
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  #functions converting inc(x) and dec(y) into equal area
  a2cx <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*sin(d2r(y))}
  a2cy <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*cos(d2r(y))}
  DI <- na.omit(DI)
  #splits modes
  if(Shiny==TRUE){
    if(arith_stat==FALSE){
      inconly_stat <- data.frame(matrix(nrow=0,ncol=6))
      colnames(inconly_stat) <- c("N","Inc","x","St.D.","a95","Ar. Mean")
    }else{
      inconly_stat <- data.frame(matrix(nrow=0,ncol=5))
      colnames(inconly_stat) <- c("N","Inc","x","St.D.","a95")
    }
  }
  if(bimodal==TRUE){
     if(dec==TRUE) {DI <- DI[,1:2]
    }else if(dec==F) {DI <- DI[,1]}
    dirs <- DI
    ifelse(dec==TRUE, colnames(dirs) <- c("dec","inc"), colnames(dirs) <- "inc")
    dirs_D <- filter_all(dirs, all_vars(inc>0))
    dirs_U <- filter_all(dirs, all_vars(inc<=0))
    #down_pointing
    if(print==TRUE){cat("Down-pointing\n")}
    inc_stat_D <- inc_only(DI = dirs_D,dec = dec, print = print,export=export, name=paste(name,"_down"), arith_stat=arith_stat)
    if(Shiny==TRUE){
      if(arith_stat==FALSE){
        colnames(inc_stat_D) <- c("N","Inc","x","St.D.","a95","Ar. Mean")
        rownames(inc_stat_D) <- "Mode 1"
        inconly_stat <- rbind(inconly_stat,inc_stat_D)
      }else{
        colnames(inc_stat_D) <- c("N","Inc","x","St.D.","a95")
        rownames(inc_stat_D) <- "Mode 1"
        inconly_stat <- rbind(inconly_stat,inc_stat_D)
      }
    }
    #up_pointing
    if(print==TRUE){cat("Up-pointing\n")}
    inc_stat_U <- inc_only(DI = dirs_U,dec = dec, print = print,export=export, name=paste(name,"_up"), arith_stat=arith_stat)
    if(Shiny==TRUE){
      if(arith_stat==FALSE){
        colnames(inc_stat_U) <- c("N","Inc","x","St.D.","a95","Ar. Mean")
        rownames(inc_stat_U) <- "Mode 2"
        inconly_stat <- rbind(inconly_stat,inc_stat_U)
      }else{
        colnames(inc_stat_U) <- c("N","Inc","x","St.D.","a95")
        rownames(inc_stat_U) <- "Mode 2"
        inconly_stat <- rbind(inconly_stat,inc_stat_U)
      }
    }
    #all_down_pointing
    if(print==TRUE){cat("All data\n")}
    dirs$inc <- abs(dirs$inc)
    inc_stat_ALL <- inc_only(DI = dirs,dec = dec, print = print,export=export, name=paste(name,"_all"), arith_stat=arith_stat)
    if(Shiny==TRUE){
      if(arith_stat==FALSE){
        colnames(inc_stat_ALL) <- c("N","Inc","x","St.D.","a95","Ar. Mean")
        rownames(inc_stat_ALL) <- "All"
        inconly_stat <- rbind(inconly_stat,inc_stat_ALL)
      }else{
        colnames(inc_stat_ALL) <- c("N","Inc","x","St.D.","a95")
        rownames(inc_stat_ALL) <- "All"
        inconly_stat <- rbind(inconly_stat,inc_stat_ALL)
      }
    }
  }else{
    dirs <- DI
    inc_stat_ALL <- inc_only(DI = dirs,dec = dec, print = print,export=export, name=paste(name,"_all"), arith_stat=arith_stat)
    if(Shiny==TRUE){
      if(arith_stat==FALSE){
        colnames(inc_stat_ALL) <- c("N","Inc","x","St.D.","a95","Ar. Mean")
        rownames(inc_stat_ALL) <- "All"
        inconly_stat <- rbind(inconly_stat,inc_stat_ALL)
      }else{
        colnames(inc_stat_ALL) <- c("N","Inc","x","St.D.","a95")
        rownames(inc_stat_ALL) <- "All"
        inconly_stat <- rbind(inconly_stat,inc_stat_ALL)
      }
    }
  }
  #plot if requested
  if(plot==TRUE){
    if(bimodal==TRUE){
      #create circles
      circle_D <- as.data.frame(seq(0,90,1))
      colnames(circle_D) <- "dec"
      circle_D$inc <- rep(inc_stat_D$Inc)
      #create a column with horizontal inclination for doubling the confidence area if crosses 0 inc
      if(inc_stat_D$Inc-inc_stat_D$a95 < 0){
        circle_D$x_hrz <- a2cx(0,circle_D$dec)
        circle_D$y_hrz <- a2cy(0,circle_D$dec)
      }
      circle_D$inc_l <- rep(abs(inc_stat_D$Inc-inc_stat_D$a95))
      circle_D$inc_h <- ifelse((inc_stat_D$Inc+inc_stat_D$a95)<90, rep(inc_stat_D$Inc+inc_stat_D$a95),rep(90))
      circle_D$x <- a2cx(circle_D$inc,circle_D$dec)
      circle_D$y <- a2cy(circle_D$inc,circle_D$dec)
      circle_D$x_l <- a2cx(circle_D$inc_l,circle_D$dec)
      circle_D$y_l <- a2cy(circle_D$inc_l,circle_D$dec)
      circle_D$x_h <- a2cx(circle_D$inc_h,circle_D$dec)
      circle_D$y_h <- a2cy(circle_D$inc_h,circle_D$dec)

      circle_U <- as.data.frame(seq(270,360,1))
      colnames(circle_U) <- "dec"
      circle_U$inc <- rep(abs(inc_stat_U$Inc))
      #create a column with horizontal inclination for doubling the confidence area if crosses 0 inc
      if(abs(inc_stat_U$Inc)-inc_stat_U$a95 < 0){
        circle_U$x_hrz <- a2cx(0,circle_U$dec)
        circle_U$y_hrz <- a2cy(0,circle_U$dec)
      }
      circle_U$inc_l <- rep(abs(abs(inc_stat_U$Inc)-inc_stat_U$a95))
      circle_U$inc_h <- ifelse((abs(inc_stat_U$Inc)+inc_stat_U$a95)<90, rep(abs(inc_stat_U$Inc)+inc_stat_U$a95), rep(90))
      circle_U$x <- a2cx(circle_U$inc,circle_U$dec)
      circle_U$y <- a2cy(circle_U$inc,circle_U$dec)
      circle_U$x_l <- a2cx(circle_U$inc_l,circle_U$dec)
      circle_U$y_l <- a2cy(circle_U$inc_l,circle_U$dec)
      circle_U$x_h <- a2cx(circle_U$inc_h,circle_U$dec)
      circle_U$y_h <- a2cy(circle_U$inc_h,circle_U$dec)
      circle_ALL <- as.data.frame(seq(90,270,1))
      colnames(circle_ALL) <- "dec"
    }else{
      #if biomdal is false draw a complete ALL circle
      circle_ALL <- as.data.frame(seq(0,360,1))
      colnames(circle_ALL) <- "dec"
    }
    #complete circle with all even if not bimodal
    circle_ALL$inc <- rep(abs(inc_stat_ALL$Inc))
    #create a column with horizontal inclination for doubling the confidence area if crosses 0 inc
    if(abs(inc_stat_ALL$Inc)-inc_stat_ALL$a95 < 0){
      circle_ALL$x_hrz <- a2cx(0,circle_ALL$dec)
      circle_ALL$y_hrz <- a2cy(0,circle_ALL$dec)
    }
    circle_ALL$inc_l <- rep(abs(abs(inc_stat_ALL$Inc)-inc_stat_ALL$a95))
    circle_ALL$inc_h <- ifelse((abs(inc_stat_ALL$Inc)+inc_stat_ALL$a95)<90,rep(abs(inc_stat_ALL$Inc)+inc_stat_ALL$a95),rep(90))
    circle_ALL$x <- a2cx(circle_ALL$inc,circle_ALL$dec)
    circle_ALL$y <- a2cy(circle_ALL$inc,circle_ALL$dec)
    circle_ALL$x_l <- a2cx(circle_ALL$inc_l,circle_ALL$dec)
    circle_ALL$y_l <- a2cy(circle_ALL$inc_l,circle_ALL$dec)
    circle_ALL$x_h <- a2cx(circle_ALL$inc_h,circle_ALL$dec)
    circle_ALL$y_h <- a2cy(circle_ALL$inc_h,circle_ALL$dec)

    #plot_circles
    if(on_plot==FALSE) equalarea()
    if(bimodal==TRUE){
      #Down
      lines(circle_D$x_l,circle_D$y_l,lty=2)
      lines(circle_D$x_h,circle_D$y_h,lty=2)
      #splits the confidence area in two if crosses zero inclination
      if(inc_stat_D$Inc-inc_stat_D$a95 < 0){
        conf1_D <- data.frame(cbind(circle_D$x_l,circle_D$y_l))
        confh_D <- data.frame(cbind(circle_D$x_hrz,circle_D$y_hrz))
        confh_D <- confh_D[nrow(confh_D):1,]
        first_D_area <- rbind(conf1_D,confh_D)
        polygon(first_D_area, col=rgb(0,0,1,0.30),border=NA)
        conf2_D <- data.frame(cbind(circle_D$x_h, circle_D$y_h))
        second_D_area <- rbind(conf2_D,confh_D)
        polygon(second_D_area, col=rgb(0,0,1,0.30),border=NA)

      }else{
        conf1_D <- data.frame(cbind(circle_D$x_l,circle_D$y_l))
        conf2_D <- data.frame(cbind(circle_D$x_h, circle_D$y_h))
        conf2_D <- conf2_D[nrow(conf2_D):1,]
        conf_D <- rbind(conf1_D,conf2_D)
        polygon(conf_D, col=rgb(0,0,1,0.30),border=NA)
      }
      lines(circle_D$x,circle_D$y,col=col, lwd=1.5)

      #Up
      lines(circle_U$x_l,circle_U$y_l,lty=2)
      lines(circle_U$x_h,circle_U$y_h,lty=2)
      #splits the confidence area in two if crosses zero inclination
      if(abs(inc_stat_U$Inc)-inc_stat_U$a95 < 0){
        conf1_U <- data.frame(cbind(circle_U$x_l,circle_U$y_l))
        confh_U <- data.frame(cbind(circle_U$x_hrz,circle_U$y_hrz))
        confh_U <- confh_U[nrow(confh_U):1,]
        first_U_area <- rbind(conf1_U,confh_U)
        polygon(first_U_area, col=rgb(0,0,1,0.30),border=NA)
        conf2_U <- data.frame(cbind(circle_U$x_h, circle_U$y_h))
        second_U_area <- rbind(conf2_U,confh_U)
        polygon(second_U_area, col=rgb(0,0,1,0.30),border=NA)
      }else{
        conf1_U <- data.frame(cbind(circle_U$x_l,circle_U$y_l))
        conf2_U <- data.frame(cbind(circle_U$x_h, circle_U$y_h))
        conf2_U <- conf2_U[nrow(conf2_U):1,]
        conf_U <- rbind(conf1_U,conf2_U)
        polygon(conf_U, col=rgb(0,1,1,0.30),border=NA)
      }
      lines(circle_U$x,circle_U$y,col=col, lwd=1.5 )
    }
    #ALL
    lines(circle_ALL$x_l,circle_ALL$y_l,lty=2)
    lines(circle_ALL$x_h,circle_ALL$y_h,lty=2)
    if(abs(inc_stat_ALL$Inc)-inc_stat_ALL$a95 < 0){
      conf1_ALL <- data.frame(cbind(circle_ALL$x_l,circle_ALL$y_l))
      confh_ALL <- data.frame(cbind(circle_ALL$x_hrz,circle_ALL$y_hrz))
      confh_ALL <- confh_ALL[nrow(confh_ALL):1,]
      first_ALL_area <- rbind(conf1_ALL,confh_ALL)
      polygon(first_ALL_area, col=rgb(0,0,1,0.30),border=NA)
      conf2_ALL <- data.frame(cbind(circle_ALL$x_h, circle_ALL$y_h))
      second_ALL_area <- rbind(conf2_ALL,confh_ALL)
      polygon(second_ALL_area, col=rgb(0,0,1,0.30),border=NA)

    }else{
      conf1_ALL <- data.frame(cbind(circle_ALL$x_l,circle_ALL$y_l))
      conf2_ALL <- data.frame(cbind(circle_ALL$x_h, circle_ALL$y_h))
      conf2_ALL <- conf2_ALL[nrow(conf2_ALL):1,]
      conf_ALL <- rbind(conf1_ALL,conf2_ALL)
      infill <- ifelse(inc_stat_ALL$Inc<0,rgb(0,1,1,0.30),rgb(0,0,1,0.30))
      if(bimodal==TRUE){infill <- rgb(1,0,0,0.30)}
      polygon(conf_ALL, col=infill,border=NA)
    }
    lines(circle_ALL$x,circle_ALL$y,col=col, lwd=1.5)
    if(save==TRUE){save_pdf(name = paste(name,".pdf"),width = 8,height = 8)}
  }
  if(Shiny==TRUE){return(inconly_stat[,1:5])} #cut last column with Arithmetic mean
}

#create matrix from fol,lin, and dec inc of vectors
matrix_maker <- function(Fol=1,Lin=1,v1d,v1i,v2d,v2i,v3d,v3i, export=FALSE, name="matrix",return_P=TRUE){
  library(matlib)
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  P <- Fol*Lin
  #eigenvalues
  v1 <- (3*Lin)/(Lin+(1/Fol)+1)
  v2 <- v1/Lin
  v3 <- v1/P
  A <- as.data.frame(matrix(c(v1d,v1i,
                              v2d,v2i,
                              v3d,v3i),
                            ncol=2,byrow=TRUE))
  colnames(A) <- c("dec","inc")
  A$x <- cos(d2r(A$dec))*cos(d2r(A$inc))
  A$y <- sin(d2r(A$dec))*cos(d2r(A$inc))
  A$z <- sin(d2r(A$inc))
  A_vec <- t(A[,3:5])
  A_val <- c(v1,0,0,0,v2,0,0,0,v3)
  A_val <- matrix(A_val,nrow=3,byrow=TRUE)
  M <- A_vec%*%A_val%*%inv(A_vec)
  if(return_P==TRUE){
    cat(paste("P:",P,"
"))
  }
  #export matrix if requested
  if(export==TRUE){write.csv(round(M,digits=5),paste(name,".csv"),row.names = TRUE)}
  return(round(M, digits=5))
}

#function that plots KavrayskiyVII geographic projection
Map_KVII <- function(grid=30, center=0, title="",seaCol="light cyan",landCol="light green",gridCol="gray") {
  library(rlist)
  if(center>180 | center<(-180)) stop("Please set center between -180° and 180°",call. = F)
  if(grid>90) stop("Please set the grid between 1° and 90°. If grid=0, grid is not plotted",call. = F)
  #functions converting degree and radians
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #functions converting long & lat to xy
  c2x <- function(lon,lat) {((3*d2r(lon))/2)*(sqrt((1/3)-((d2r(lat)/pi)^2)))}
  c2y <- function(lat) {d2r(lat)}

  #fix frame
  plot(NA, xlim=c(-2.6,2.6), ylim=c(-1.5,1.5), asp=1,
       xlab=title, xaxt="n",ylab="", yaxt="n", axes=FALSE)
  #import coastline from PmagDiR
  if (abs(center)==180){
    cl <- world_coastline_180
  }else{cl <- world_coastline}
  # set the coastline offset if longitude is not 0
  if(center!=0 | abs(center)!=180){
    #create the coastline-breaks file (breaks separating the single continental lines)
    sep <- as.data.frame(1)
    colnames(sep) <- "index"
    #finds all breaks in the coastline file and index them
    list <- as.data.frame(which(is.na(cl[,1]),arr.ind = TRUE))
    colnames(list) <- "index"
    sep <- rbind(sep, list)
    rm(list)
    #create empty list of all continents coastline
    conts <- list()
    #function that isolate coastlines and append to list
    for(i in 2:nrow(sep)){
      contour <- cl[sep[i-1,1]:sep[i,1],]
      contour <- na.omit(contour)
      conts <- list.append(conts,contour)
    }
    new_cl <- as.data.frame(matrix(nrow = 0,ncol = 2))
    colnames(new_cl) <- c("lon","lat")
    #for every coastline it changes the longitude depending on the new center and it fix it (-180<l<+180)
    for(i in 1:length(conts)){
      clc <- conts[[i]]
      clc$lon <- clc$lon-center
      clc$lon <- ifelse(clc$lon<(-180),clc$lon+360,clc$lon)
      clc$lon <- ifelse(clc$lon>180,clc$lon-360,clc$lon)
      #creates a new dataframe for single coastlines with breaks
      clc1 <- as.data.frame(matrix(nrow = 0,ncol = 2))
      colnames(clc1) <- c("lon","lat")
      #creates the breaking line
      breaker <- as.data.frame(t(c(NA,NA)))
      colnames(breaker) <- c("lon","lat")
      #when two points are on the different side of the map, based on 350° distance, it put break in between
      for(l in 2:nrow(clc)){
        if(abs(clc[l-1,1]-clc[l,1])>350){
          provv1 <- clc[l-1,]
          provv2 <- clc[l,]
          clc1 <- rbind(clc1,provv1,breaker,provv2)
        }else{clc1 <- rbind(clc1,clc[l-1,],clc[l,])}
      }
      #puts break after the new continent line
      clc1 <- rbind(clc1,breaker)
      #appends all new continent lines
      new_cl <- rbind(new_cl,clc1)
    }
  }
  #plot sea
  bord_left <- as.data.frame(matrix(ncol = 2, nrow = 181))
  bord_left[,2] <- -90:90
  bord_left[,1] <- rep(-180)
  bord_right <- as.data.frame(matrix(ncol = 2, nrow = 181))
  bord_right[,2] <- 90:-90
  bord_right[,1] <- rep(180)
  bord <- rbind(bord_left,bord_right)
  bord$x <- c2x(bord[,1], bord[,2])
  bord$y <- c2y(bord[,2])
  if(center==0 | abs(center)==180){
    polygon(bord$x,bord$y, col=seaCol,border = NA)
  }

  #set coastline if longitude is greenwich
  if(center==0 | abs(center)==180){
    new_cl <- cl
    polygon(x = c2x(new_cl$lon,new_cl$lat),
                        y = c2y(new_cl$lat), col=landCol, border=landCol)
    }else{
      lines(x = c2x(new_cl$lon,new_cl$lat),
        y = c2y(new_cl$lat), col="black")
    }

  #plot grid only if different from 0
  if(grid!=0){
    #plot_main_parallel
    #longitude circle
    lats <- seq(-(90-grid),(90-grid),grid)
    for(i in lats){
      lon_lat_p <-  as.data.frame(-180:180)
      lon_lat_p$lat <- rep(i)
      lon_lat_p$x <- c2x(lon_lat_p[,1],lon_lat_p[,2])
      lon_lat_p$y <- c2y(lon_lat_p[,2])
      lines(lon_lat_p$x,lon_lat_p$y,col=gridCol, pch=16, cex=0.3, lty=1)
    }
    #plot_main_meridians
    #fix meridians if center is not greenwich
    Gr <- (-center)
    LonLeft <- seq((Gr),-180,-grid)
    LonRight <- seq(Gr,180,grid)
    LonRight <- LonRight[-1]
    #plot meridians left of Greenwich
    for(i in LonLeft){
      lat_lon_m <- as.data.frame(seq(-89,89,1))
      lat_lon_m$lon <- rep(i)
      lat_lon_m$x <- c2x(lat_lon_m[,2],lat_lon_m[,1])
      lat_lon_m$y <- c2y(lat_lon_m[,1])
      lines(lat_lon_m$x,lat_lon_m$y,col=gridCol, pch=16, cex=0.3, lty=1)
    }
    #plot meridians right of Greenwich
    for(i in LonRight){
      lat_lon_m <- as.data.frame(seq(-89,89,1))
      lat_lon_m$lon <- rep(i)
      lat_lon_m$x <- c2x(lat_lon_m[,2],lat_lon_m[,1])
      lat_lon_m$y <- c2y(lat_lon_m[,1])
      lines(lat_lon_m$x,lat_lon_m$y,col=gridCol, pch=16, cex=0.3, lty=1)
    }
  }
  #plot contour
  polygon(bord$x,bord$y, col=NA)
}

#calculate PCA-derived direction and MAD from demagnetization steps
PCA_DI <- function(DII,anchor="f", export=FALSE,name="PCA") {
  #degree to radians and VV
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  data <- DII
  colnames(data) <- c("dec", "inc","int")
  #directions in Cartesian coordinates
  data$x <- data$int*cos(d2r(data$dec))*cos(d2r(data$inc))
  data$y <- data$int*sin(d2r(data$dec))*cos(d2r(data$inc))
  data$z <- data$int*sin(d2r(data$inc))
  #averaged Cartesian coordinates
  x_av <- mean(data$x)
  y_av <- mean(data$y)
  z_av <- mean(data$z)
  #copy coordinates for anchored directions
  if (anchor=="a"){
    data$xn <- data$x
    data$yn <- data$y
    data$zn <- data$z
  } else if(anchor=="i") {
    #includes origin and calculate new center of mass
    newrow <- c(0,0,0,0,0,0)
    data <- rbind(newrow,data)
    data$xn <- data$x-x_av
    data$yn <- data$y-y_av
    data$zn <- data$z-z_av
  } else {
    #calculate coordinates with new center of mass
    data$xn <- data$x-x_av
    data$yn <- data$y-y_av
    data$zn <- data$z-z_av
  }
  #elements of the distribution matrix
  T_elements <- c(sum((data$xn)*(data$xn)),sum(data$xn*data$yn),sum(data$xn*data$zn),
                  sum(data$yn*data$xn),sum(data$yn*data$yn),sum(data$yn*data$zn),
                  sum(data$zn*data$xn),sum(data$zn*data$yn),sum(data$zn*data$zn))

  T <- matrix(T_elements,nrow=3, byrow=TRUE)
  T_e <- eigen(T)
  T_vec <- T_e$vectors
  T_val <- T_e$value
  #calculate dec inc of max variance
  Vdec <- r2d(atan2(T_vec[2,1],T_vec[1,1]))
  Vdec <- ifelse(Vdec<0,Vdec+360,Vdec)
  Vinc <- r2d(asin(T_vec[3,1]/(sqrt((T_vec[1,1]^2)+(T_vec[2,1]^2)+(T_vec[3,1]^2)))))
  MAD <- r2d(atan(sqrt(((T_val[2])+(T_val[3]))/T_val[1])))
  N <- length(data[,1])

  dirs <- cbind(Vdec,Vinc,MAD,N)
  colnames(dirs) <- c("Dec", "Inc","MAD","N")
  if(export==TRUE){write.csv(round(dirs,digits = 2),paste(name,".csv"),row.names = FALSE)}
  return(dirs)
}

#function that plots a95 with dec,inc,a95
plot_a95 <- function(D,I,a, col_d="red",col_u="white",col_l="black", symbol="c", on_plot=FALSE, save=FALSE, name="F_a95"){
  library("dplyr", warn.conflicts = FALSE)
  #functions converting degree and radians
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  #functions converting inc(x) and dec(y) into equal area
  a2cx <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*sin(d2r(y))}
  a2cy <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*cos(d2r(y))}
  #functions spherical (Dec=x, Inc=y) to Cartesian
  s2cx <- function(x,y) {cos(d2r(x))*cos(d2r(y))}
  s2cy <- function(x,y) {sin(d2r(x))*cos(d2r(y))}
  s2cz <- function(y) {sin(d2r(y))}
  #save declination and inc and calculate new system for rotation
  dec <- D
  newSdec <- ifelse((D+180)>360,D-180,D+180)
  inc <- I
  newSinc <- 90-I
  newSdecr <- d2r(newSdec)
  newSincr <- d2r(newSinc)
  a95 <- a
  circle <- as.data.frame(matrix(ncol=2,nrow=0))
  #loop that create a95 and rotate it around new coordinate (dec, inc)
  for (i in seq(0,360,2)){
    circleP <- as.data.frame(matrix(ncol=2,nrow=1))
    x <- s2cx(i,(90-a95))
    y <- s2cy(i,(90-a95))
    z <- s2cz(90-a95)
    vec <- as.matrix(c(x,y,z))
    R_elements <- c(cos(newSincr)*cos(newSdecr), -sin(newSdecr), -sin(newSincr)*cos(newSdecr),
                    cos(newSincr)*sin(newSdecr), cos(newSdecr), -sin(newSincr)*sin(newSdecr),
                    sin(newSincr), 0, cos(newSincr))
    R <- matrix(R_elements,nrow=3, byrow=TRUE)
    newvec <- R%*%vec
    newdec <- r2d(atan2(newvec[2,1],newvec[1,1]))
    newdec <- ifelse(newdec<0,newdec+360,newdec)
    #absolute value avoid point outside the graph
    newinc <- abs(r2d(asin(newvec[3,1])))
    circleP[1,1:2] <- c(newdec,newinc)
    circle <- rbind(circle,circleP)
  }
  colnames(circle) <- c("dec","inc")
  circle$x <- a2cx(circle$inc,circle$dec)
  circle$y <- a2cy(circle$inc,circle$dec)
  #restore screen
  par(fig=c(0,1,0,1))
  #standalone graph or on existing graph
  if (on_plot==FALSE) {
    plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1,
         xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
    equalarea()
    }
  UD <- ifelse(inc>0,"D","U")
  inc <- abs(inc)
  X <- a2cx(inc,dec)
  Y <- a2cy(inc,dec)
  if(symbol=="c") {pch <- 21}
  else if(symbol=="s") {pch <- 22}
  else if(symbol=="d") {pch <- 23}
  else if(symbol=="t") {pch <- 24}
  else{stop("Please select valid symbol. Check help for info.",call. = F)}

  if(UD=="D"){
    points(X,Y, pch=pch,cex=1.3, col="black",
           bg= col_d)
  }else{
    points(X,Y, pch=pch,cex=1.3, col="black",
           bg=col_u)
  }
  lines(circle$x,circle$y,lty=1, col=col_l, lwd=1.8)
  if(save==TRUE){save_pdf(name = paste(name,".pdf"),width = 8,height = 8)}
}

#plot A95 on a spherical orthographic plot
plot_PA95 <- function(lon,lat,A,lon0=0,lat0=90,grid=30, col_f="red",col_b="white",col_l="black",col_A=rgb(1,0,0,0.30), symbol="c",size=1, coast=FALSE, on_plot=FALSE, save=FALSE, name="A95"){
  library("dplyr", warn.conflicts = FALSE)
  #functions converting degree and radians
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #functions converting long & lat to xy
  c2x <- function(lon,lat) {cos(d2r(lat))*sin(d2r(lon-lon0))}
  c2y <- function(lon,lat) {(cos(d2r(lat0))*sin(d2r(lat)))-(sin(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}
  #cut is cosin of c, when negative is behind projections, needs to be cut
  cut <- function(lon,lat) {(sin(d2r(lat0))*sin(d2r(lat)))+(cos(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}

  #functions spherical (lon=x, lat=y) to Cartesian
  s2cx <- function(x,y) {cos(d2r(x))*cos(d2r(y))}
  s2cy <- function(x,y) {sin(d2r(x))*cos(d2r(y))}
  s2cz <- function(y) {sin(d2r(y))}
  #save declination and inc and calculate new system for rotation
  newSlon <- ifelse((lon+180)>360,lon-180,lon+180)
  newSlat <- 90-lat
  newSlonr <- d2r(newSlon)
  newSlatr <- d2r(newSlat)
  a95 <- A
  circle <- as.data.frame(matrix(ncol=2,nrow=0))
  #loop that create a95 and rotate it around new coordinate (dec, inc)
  for (i in seq(0,360,2)){
    circleP <- as.data.frame(matrix(ncol=2,nrow=1))
    x <- s2cx(i,(90-a95))
    y <- s2cy(i,(90-a95))
    z <- s2cz(90-a95)
    vec <- as.matrix(c(x,y,z))
    R_elements <- c(cos(newSlatr)*cos(newSlonr), -sin(newSlonr), -sin(newSlatr)*cos(newSlonr),
                    cos(newSlatr)*sin(newSlonr), cos(newSlonr), -sin(newSlatr)*sin(newSlonr),
                    sin(newSlatr), 0, cos(newSlatr))
    R <- matrix(R_elements,nrow=3, byrow=TRUE)
    newvec <- R%*%vec
    newlon <- r2d(atan2(newvec[2,1],newvec[1,1]))
    newlon <- ifelse(newlon<0,newlon+360,newlon)
    #absolute value avoid point outside the graph
    newlat <- r2d(asin(newvec[3,1]))
    circleP[1,1:2] <- c(newlon,newlat)
    circle <- rbind(circle,circleP)
  }
  colnames(circle) <- c("lon","lat")
  circle$x <- c2x(circle$lon,circle$lat)
  circle$y <- c2y(circle$lon,circle$lat)
  circle$cut <- cut(circle$lon,circle$lat)
  #restore screen
  par(fig=c(0,1,0,1))
  #standalone graph or on existing graph
  if (on_plot==FALSE) {
    plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1,
         xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
    sph_ortho(lat = lat0,long = lon0,grid = grid, coast=coast)
    }

  X <- c2x(lon,lat)
  Y <- c2y(lon,lat)
  CUT <- cut(lon,lat)

  #plot alfa 95
  polygon(circle$x,circle$y, col=col_A, lwd=0.8,lty= ifelse(CUT>0,1,3))

  if(symbol=="c") {pch <- 21}
  else if(symbol=="s") {pch <- 22}
  else if(symbol=="d") {pch <- 23}
  else if(symbol=="t") {pch <- 24}
  else{stop("Please select valid symbol. Check help (type ?plot_PA95) for info.",call. = F)}

  if(CUT>0){
    points(X,Y, pch=pch,cex=size, col="black",
           bg= col_f)
  }else{
    points(X,Y, pch=pch,cex=size, col="black",
           bg=col_b)
  }
  if(save==TRUE){save_pdf(name = paste(name,".pdf"),width = 8,height = 8)}
}

#plot pole with A95 and Apparent polar wander path
plot_pole_APWP <- function(lon,lat,A,lon0=0,lat0=90,grid=30, col_f="red",col_b="white",col_l="black",col_A=rgb(1,0,0,0.30), symbol="c",coast=FALSE, on_plot=FALSE, save=FALSE, name="A95",APWP="V23", S_APWP=FALSE){
  #plot pole
  plot_PA95(lon=lon, lat = lat,A = A,lon0=lon0, lat0=lat0, grid=grid, col_f = col_f, col_b= col_b, col_l=col_l, col_A = col_A, symbol=symbol, coast=coast, on_plot = on_plot,save=FALSE)

  #plot APWP if requested during process
  pAPWP <- readline("Plot APWP? (y or n): ")
  if(pAPWP=="y"){
    plot_APWP(APWP=APWP, lon0 = lon0,lat0 = lat0,grid = grid,on_plot = T,S_APWP = S_APWP)
  }
  if(save==TRUE){save_pdf(name = paste(name,".pdf"),width = 8,height = 8)}
}

#plot only apparent polar wander path
plot_APWP <- function(APWP= "V23",lon0=0,lat0=90,grid=30,col="gray",symbol="c",size=0.6, coast=FALSE, on_plot=FALSE, save=FALSE, name="APWP",S_APWP=FALSE,Shiny=FALSE,Y=0,O=320,frame=1,Age_size=1){
  if (on_plot==FALSE) {
    plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1,
         xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
    sph_ortho(lat = lat0,long = lon0,grid = grid, coast=coast)
  }

  #plot APWP
  #functions converting degree and radians
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #functions converting long & lat to xy
  c2x <- function(lon,lat) {cos(d2r(lat))*sin(d2r(lon-lon0))}
  c2y <- function(lon,lat) {(cos(d2r(lat0))*sin(d2r(lat)))-(sin(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}
  #cut is cosin of c, when negative is behind projections, needs to be cut
  cut <- function(lon,lat) {(sin(d2r(lat0))*sin(d2r(lat)))+(cos(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}
  #questions on APW age and frame
  if(Shiny==FALSE){
    if(APWP=="V23"){
      cat("Frames:
(1) South Africa
(2) North America
(3) South America
(4) Europe
(5) India
(6) Australia
(7) Antarctica
(8) Pacific (0 to 80_Ma)
(9) Iberia (0 to 80 Ma)")
    } else if (APWP=="T12"){
      cat("Frames:
(1) South Africa
(2) North America
(3) Europe
(4) India
(5) Amazonia
(6) Australia
(7) East Antarctica")
    }
    frame <- as.numeric(readline("insert frame (number): "))
    cat("APWP range from 0 to 320 Ma every 10 Myr.
")
    Y <- round(as.numeric(readline("Insert younger age: ")),-1)
    O <- round(as.numeric(readline("Older age: ")),-1)
  }
  if(Shiny==TRUE){
    Y=Y
    O=O
    frame=frame
  }
  col1 <- (frame*2)+1
  col2 <- (frame*2)+2
  if(is.na(O)==TRUE){O <- 320}
  if(O>320){O <- 320}
  if(is.na(Y)==TRUE){Y <- 0}
  #if frame is 8 or 9 (only in V23) and age is too old it fixes it
  if(frame==8 && O>80) O <- 80
  if(frame==9 && O>80) O <- 80
  Y <- (Y/10)+1
  O <- (O/10)+1
  #select apwp file
  if(APWP=="V23") G <- V23_GAPWP
  if(APWP=="T12") G <- T12_GAPWP
  #flip if necessary
  if(S_APWP==FALSE) {G[,col1:col2] <- flip_DI(G[,col1:col2])}
  par(fig=c(0,1,0,1), new=TRUE)
  plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1,
       xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
  #line connecting APWP
  lin <- as.data.frame(c2x(G[Y:O,col1],G[Y:O,col2]))
  colnames(lin) <- "lx"
  lin$ly <- c2y(G[Y:O,col1],G[Y:O,col2])
  lin$cut <- cut(G[Y:O,col1],G[Y:O,col2])
  lines(lin$lx,lin$ly,cex=1)
  #plot poles APWP
  for (i in Y:O){
    plot_PA95(lon = G[i,col1],lat = G[i,col2],A = G[i,2],lon0 = lon0,lat0 = lat0,on_plot = T,col_f = col,symbol=symbol,size=size, col_l = "black",col_A=rgb(1,0.9,0,0.30))
  }
  text1 <- paste(G[Y,1],"Ma")
  text2 <- paste(G[O,1], "Ma")
  text(x=lin[1,1], y=lin[1,2],pos=4,substitute(paste(bold(text1))), cex= Age_size)
  text(x=lin[length(lin$lx),1], y=lin[length(lin$lx),2],pos=4,substitute(paste(bold(text2))), cex= Age_size)

  if(save==TRUE){save_pdf(name = paste(name,".pdf"),width = 8,height = 8)}
}

#function plotting directions
plot_DI <- function(DI,single_mode=FALSE, down=TRUE,symbol="c", col_d="blue",col_u="cyan",col_ext="black", on_plot=FALSE, title="",save=FALSE,name="Equal_area"){
  library("dplyr", warn.conflicts = FALSE)

  #functions converting degree and radians
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #functions converting inc(x) and dec(y) into equal area
  a2cx <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*sin(d2r(y))}
  a2cy <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*cos(d2r(y))}

  data <- DI[,1:2]
  data <- na.omit(data)
  data <- data[,1:2]
  colnames(data) <- c("dec","inc")
  if(single_mode==TRUE) {data <- common_DI(data,down=down)}
  data_U <- filter_all(data,all_vars(inc<0))
  data_D <- filter_all(data,all_vars(inc>=0))
  data_U$inc <- abs(data_U$inc)
  xU <- a2cx(data_U$inc,data_U$dec)
  yU <- a2cy(data_U$inc,data_U$dec)
  xD <- a2cx(data_D$inc,data_D$dec)
  yD <- a2cy(data_D$inc,data_D$dec)
  if(on_plot==FALSE){
    equalarea(title=title)
  }
  if(symbol=="c") {pch <- 21}
  else if(symbol=="s") {pch <- 22}
  else if(symbol=="d") {pch <- 23}
  else if(symbol=="t") {pch <- 24}
  else{stop("Please select valid symbol. Check help for info.",call. = F)}

  points(xD,yD, pch=pch,col=col_ext,
         bg= col_d)
  points(xU,yU, pch=pch,col=col_ext,
         bg=col_u)
  if(save==TRUE){save_pdf(name = paste(name,".pdf"),width = 8,height = 8)}
}

#plot plain intersection in equal area given declination and inclination of pole
plot_plane <- function(D,I, col_cD="black",col_cU="grey", pole=TRUE, col_d="red",col_u="white", symbol="s",on_plot=TRUE,save=FALSE,name="plane"){
  #functions degree to radians and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  #functions converting inc(x) and dec(y) into equal area
  a2cx <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*sin(d2r(y))}
  a2cy <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*cos(d2r(y))}
  #functions spherical (Dec=x, Inc=y) to Cartesian
  s2cx <- function(x,y) {cos(d2r(x))*cos(d2r(y))}
  s2cy <- function(x,y) {sin(d2r(x))*cos(d2r(y))}
  s2cz <- function(y) {sin(d2r(y))}
  #set intial declination of circle, to avoid lines when plotted
  #strike of pole dec
  data <- D-90
  data <- as.data.frame(seq(data,data+359,1))
  data$inc <- 0
  colnames(data) <- c("dec","inc")
  #define plane Azimuth and dip from pole
  Paz <- D+180
  Pdip <- 90-I
  #sines and cosines of plane coord
  sbd <- -sin(d2r(Paz))
  cbd <- cos(d2r(Paz))
  sbi <- sin(d2r(Pdip))
  cbi <- cos(d2r(Pdip))
  #create new rotated circle
  newDI <- as.data.frame(matrix(ncol=2,nrow=0))
  for(i in 1:length(data[,1])){
    newDI_p <- as.data.frame(matrix(ncol=2,nrow=1))
    x <- s2cx(data[i,1],data[i,2])
    y <- s2cy(data[i,1],data[i,2])
    z <- s2cz(data[i,2])
    xn <- x*(sbd^2+cbd^2*cbi)+
      y*(cbd*sbd*(1-cbi))+
      z*sbi*cbd
    yn <- x*cbd*sbd*(1-cbi)+
      y*(cbd^2+sbd*sbd*cbi)-
      z*sbd*sbi
    zn <- -(x*cbd*sbi-
              y*sbi*sbd-
              z*cbi)
    newdec <- r2d(atan2(yn,xn))
    newdec <- ifelse(newdec<0,newdec+360,newdec)
    newinc <- r2d(asin(zn))
    newDI_p[1,1:2] <- c(newdec,newinc)
    newDI <- rbind(newDI,newDI_p)
  }
  #convert inc to absolute
  newDI[,2] <- abs(newDI[,2])
  #positive inclination half circle
  circle_U <- newDI[1:180,]
  #negative inclination half circle
  circle_D <- newDI[181:360,]
  colnames(circle_U) <- c("dec","inc")
  colnames(circle_D) <- c("dec","inc")
  #convert to x y
  circle_U$x <- a2cx(circle_U$inc,circle_U$dec)
  circle_U$y <- a2cy(circle_U$inc,circle_U$dec)
  circle_D$x <- a2cx(circle_D$inc,circle_D$dec)
  circle_D$y <- a2cy(circle_D$inc,circle_D$dec)
  #standalone graph or on existing graph
  if (on_plot==FALSE) equalarea()
  UD <- ifelse(I>0,"D","U")
  I <- abs(I)
  X <- a2cx(I,D)
  Y <- a2cy(I,D)
  if(symbol=="c") {pch <- 21}
  else if(symbol=="s") {pch <- 22}
  else if(symbol=="d") {pch <- 23}
  else if(symbol=="t") {pch <- 24}
  else{stop("Please select valid symbol. Type ?plot_DI for info.",call. = F)}

  #plot pole only if pole==TRUE
  if(pole==TRUE){
    if(UD=="D"){
      points(X,Y, pch=pch,cex=1.3, col="black",
             bg= col_d)
    }else{
      points(X,Y, pch=pch,cex=1.3, col="black",
             bg=col_u)
    }
  }
  points(x=circle_U$x,y=circle_U$y,type="l", col=col_cU,lty=2)
  points(x=circle_D$x,y=circle_D$y,type="l", col=col_cD)
  if(save==TRUE){save_pdf(name = paste(name,".pdf"),width = 8,height = 8)}
}

#plot virtual geomagnetic poles
plot_VGP <- function(VGP,lat=90,long=0,grid=30, col="black", on_plot=FALSE,auto_cent=TRUE,exp=TRUE,coast=FALSE, title="",save=TRUE,A95=FALSE,name="VGP"){
  #functions converting degree and radians
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #functions converting long & lat to xy
  c2x <- function(lon,lat) {cos(d2r(lat))*sin(d2r(lon-lon0))}
  c2y <- function(lon,lat) {(cos(d2r(lat0))*sin(d2r(lat)))-(sin(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}
  #cut is cosin of c, when negative is behind projections, needs to be cut
  cut <- function(lon,lat) {(sin(d2r(lat0))*sin(d2r(lat)))+(cos(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}
  #manipulate data
  #VGP <- na.omit(VGP)
  colnames(VGP) <- c("lon","lat")
  vgpsN <- common_DI(VGP,down = ifelse(mean(VGP$lat)<0,FALSE,TRUE))
  PPole <- fisher(vgpsN)

  #fix point of view
  if(auto_cent==FALSE){
    #center of proj is Lon0 & Lat0
    lon0 <- long
    lat0 <- lat
  }else{
    lon0 <- PPole[1,1]
    lat0 <- PPole[1,2]
  }
  if(on_plot==FALSE){
    plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1,
         xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
    sph_ortho(lat=lat0,long=lon0,grid=grid,coast=coast, title=title)
    }

  coord <- as.data.frame(lon0)
  coord$lat0 <- lat0
  if(exp==TRUE){
    write.csv(round(coord, digits=2),file="center_coordinates.csv", row.names = F)
  }
  cat(paste("Center coordinates:

"))
  print(round(coord,digits=2),row.names=F)

  VGP$x <- c2x(VGP$lon,VGP$lat)
  VGP$y <- c2y(VGP$lon,VGP$lat)
  VGP$cut <- cut(VGP$lon,VGP$lat)

  points(VGP$x,VGP$y,pch=ifelse(VGP$cut>0,21,1),col="black",bg=col)
  if(A95==TRUE){
    plot_PA95(lon = PPole[1,1],lat = PPole[1,2],A = PPole[1,3],lon0 = lon0,lat0 = lat0,on_plot = TRUE,symbol = "d",col_l = "red")
    text <- paste("N: ",PPole[1,4],"
Long: ", round(PPole[1,1],digits=2),"
Lat: ", round(PPole[1,2], digits=2),"
A95: ", round(PPole[1,3], digits=2))
    text(x=0.75, y=-0.85,pos=4,text, cex= 0.85)
  }
  if(save==TRUE){save_pdf(name = paste(name,".pdf"),width = 8,height = 8)}
}

#reversal test boostrapped following Tauxe
revtest <- function(DI,nb=1000,export=TRUE, name="reversal_test"){
  #fucnctions deg to rads and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  data <- DI[,1:2]
  data <- na.omit(data)
  colnames(data) <- c("dec", "inc")

  #directions in Cartesian coordinates
  data$x <- cos(d2r(data$dec))*cos(d2r(data$inc))
  data$y <- sin(d2r(data$dec))*cos(d2r(data$inc))
  data$z <- sin(d2r(data$inc))

  #averaged Cartesian coordinates
  x_av <- mean(data$x)
  y_av <- mean(data$y)
  z_av <- mean(data$z)

  #elements of the distribution matrix
  T_elements <- c(sum((data$x)*(data$x)),sum(data$x*data$y),sum(data$x*data$z),
                  sum(data$y*data$x),sum(data$y*data$y),sum(data$y*data$z),
                  sum(data$z*data$x),sum(data$z*data$y),sum(data$z*data$z))

  #distribution matrix
  T <- matrix(T_elements,nrow=3, byrow=TRUE)

  #calculate and copy eigenvalues and vectors
  T_e <- eigen(T,symmetric = TRUE)
  T_vec <- T_e$vectors
  T_val <- T_e$values

  #calculate dec inc of max variance
  V1inc <- r2d(asin(T_vec[3,1]/(sqrt((T_vec[1,1]^2)+(T_vec[2,1]^2)+(T_vec[3,1]^2)))))
  V1dec <- r2d(atan2(T_vec[2,1],T_vec[1,1]))
  V1dec <- ifelse(V1dec<0,V1dec+360,V1dec)

  #flip V1 if negative
  V1dec <- ifelse(V1inc<0,ifelse((V1dec+180)>360,V1dec-180,V1dec+180),V1dec)
  V1inc <- ifelse(V1inc<0,-V1inc,V1inc)


  #next  calculates difference between dec_inc and average
  data$Dec_aver <- rep(V1dec)
  data$Inc_aver <- rep(V1inc)
  data$delta <- abs(data$dec-data$Dec_aver)
  data$diff <- r2d(acos((sin(d2r(data$inc))*sin(d2r(data$Inc_aver)))+
                          (cos(d2r(data$inc))*cos(d2r(data$Inc_aver))*cos(d2r(data$delta)))))
  #Isolate modes
  m1ind <- as.numeric(which(data$diff<=90), arr.ind = TRUE)
  m2ind <- as.numeric(which(data$diff>90), arr.ind = TRUE)

  #terminate if distribution is not bimodal
  if(length(m2ind)<1) stop("
DISTRIBUTION NOT BIMODAL")
  mode1 <- data[m1ind,1:2]
  mode2 <- data[m2ind,1:2]

  #flip mode 2 same as mode 1
  mode2$dec <- ifelse((mode2$dec+180)>360,mode2$dec-180,mode2$dec+180)
  mode2$inc <- -mode2$inc
  nb <- nb
  n <- 0
  mode1B <- as.data.frame(matrix(ncol=3, nrow=0))
  mode2B <- as.data.frame(matrix(ncol=3, nrow=0))

  #simulate pseudosamples of mode 1
  repeat{
    n <- n+1
    mode1B_p <- as.data.frame(matrix(ncol=3, nrow=1))
    Bdata <- boots_DI(mode1)
    Bdata$x <- cos(d2r(Bdata$dec))*cos(d2r(Bdata$inc))
    Bdata$y <- sin(d2r(Bdata$dec))*cos(d2r(Bdata$inc))
    Bdata$z <- sin(d2r(Bdata$inc))
    mode1B_p[1,1] <- mean(Bdata$x)
    mode1B_p[1,2] <- mean(Bdata$y)
    mode1B_p[1,3] <- mean(Bdata$z)
    mode1B <- rbind(mode1B,mode1B_p)
    if(((n%%50)==0)==TRUE){
      cat(paste(n,"simulations out of",nb,"of mode 1 done
"))
    }
    if(n==nb) break
  }
  colnames(mode1B) <- c("x","y","z")
  mode1B$dec <- r2d(atan2(mode1B$y,mode1B$x))
  mode1B$dec <- mode1B$dec%%360
  mode1B$inc <- r2d(asin(mode1B$z))
  n <- 0
  #simulate pseudosamples of mode 2
  repeat{
    n <- n+1
    mode2B_p <- as.data.frame(matrix(ncol=3, nrow=1))
    Bdata <- boots_DI(mode2)
    Bdata$x <- cos(d2r(Bdata$dec))*cos(d2r(Bdata$inc))
    Bdata$y <- sin(d2r(Bdata$dec))*cos(d2r(Bdata$inc))
    Bdata$z <- sin(d2r(Bdata$inc))
    mode2B_p[1,1] <- mean(Bdata$x)
    mode2B_p[1,2] <- mean(Bdata$y)
    mode2B_p[1,3] <- mean(Bdata$z)
    mode2B <- rbind(mode2B,mode2B_p)
    if(((n%%50)==0)==TRUE){
      cat(paste(n,"simulations out of",nb,"of mode 2 done
"))
    }
    if(n==nb) break
  }
  colnames(mode2B) <- c("x","y","z")
  mode2B$dec <- r2d(atan2(mode2B$y,mode2B$x))
  mode2B$dec <- ifelse(mode2B$dec<0,mode2B$dec+360,mode2B$dec)
  mode2B$inc <- r2d(asin(mode2B$z))

  #isolate components of models
  B1x <- sort(mode1B[,1])
  B1y <- sort(mode1B[,2])
  B1z <- sort(mode1B[,3])
  B2x <- sort(mode2B[,1])
  B2y <- sort(mode2B[,2])
  B2z <- sort(mode2B[,3])

  #define low and high boostrapped margins
  confn <- 0.95
  num <- round((nb*(1-confn))/2,digits=0)
  Lconf <- num
  Uconf <- nb-num
  B1x_l <- c(B1x[Lconf],B1x[Uconf])
  B2x_l <- c(B2x[Lconf],B2x[Uconf])
  B1y_l <- c(B1y[Lconf],B1y[Uconf])
  B2y_l <- c(B2y[Lconf],B2y[Uconf])
  B1z_l <- c(B1z[Lconf],B1z[Uconf])
  B2z_l <- c(B2z[Lconf],B2z[Uconf])

  #max and min values for graphs
  xmax <- round(max(c(B1x,B2x)), digits=1)+0.05
  xmin <- round(min(c(B1x,B2x)),digits=1)-0.05
  ymax <- round(max(c(B1y,B2y)), digits=1)+0.05
  ymin <- round(min(c(B1y,B2y)),digits=1)-0.05
  zmax <- round(max(c(B1z,B2z)), digits=1)+0.05
  zmin <- round(min(c(B1z,B2z)),digits=1)-0.05

  #function that extract intervals and counts from hist function and make cumulative curve
  cumulative_curve <- function(x){
    h <- hist(x, breaks=50,plot = FALSE)
    cnts <- h[["counts"]]
    t <- length(cnts)
    new_c <- as.data.frame(matrix(ncol=1,nrow = 1))
    for(i in 1:t){
      if(i==1){new_cp <- cnts[1]}
      if(i>1) {new_cp <- new_c[i,1]+cnts[i-1]}
      new_c <- rbind(new_c,new_cp)
    }
    new_c <- na.omit(new_c)
    breaks <- as.data.frame(h[["mids"]])
    cumul <- cbind(breaks,new_c)
    colnames(cumul) <- c("breaks","counts")
    cumul$counts <- cumul$counts/nb
    return(cumul)
  }
  cu1x <- cumulative_curve(B1x)
  cu2x <- cumulative_curve(B2x)
  cu1y <- cumulative_curve(B1y)
  cu2y <- cumulative_curve(B2y)
  cu1z <- cumulative_curve(B1z)
  cu2z <- cumulative_curve(B2z)
  text1 <- "Equal area projections"
  text2 <- "Normalized cumulative distributions"

  #clean screen to avoid figure over figure
  par(fig=c(0,1,0,1))
  plot(0, xlim=c(0,1), ylim=c(0,1),
       xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)

  #plot title for equal area
  par(fig=c(0,0.65,0.4,1))
  plot(NA, xlim=c(0,1), ylim=c(0,1),
       xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
  text(x=0.5,y=1,text1,cex=1.2)

  #plot title for cumulative distributions
  par(fig=c(0.55,1,0.5,1),new=TRUE)
  plot(NA, xlim=c(0,1), ylim=c(0,1),
       xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
  text(x=0.5,y=0.9,text2,cex=1.2)

  #plot equal area
  #real data
  par(fig=c(0,0.65,0.4,0.98),new=TRUE)
  plot_DI(mode1[,1:2],col_d = rgb(1,0,0,0.30),col_u=rgb(1,0.75,1,0.30),col_ext = NA, title = "Data common mode")
  plot_DI(mode2[,1:2],col_d = rgb(0,0,1,0.30),col_u=rgb(0,1,1,0.30),col_ext = NA,on_plot = TRUE)
  #pseudosamples mean
  par(fig=c(0,0.65,0,0.58), new=TRUE)
  plot_DI(mode1B[,4:5],col_d = rgb(1,0,0,0.20),col_u=rgb(1,0.75,1,0.30),col_ext = NA,title = "Pseudosample means")
  plot_DI(mode2B[,4:5],col_d =rgb(0,0,1,0.20),col_u=rgb(0,1,1,0.30),col_ext = NA,on_plot = TRUE)

  #plot cumulative distributions
  par(fig=c(0.55,1,0.58,0.91),new=TRUE)
  #plot x
  plot(0,type="n",xlim=c(xmin,xmax),ylim=c(0,1),yaxp= c(0,1,2), ylab="", xlab="",cex.axis=0.8)
  title(xlab = "x axis", line=2)
  rect(xleft = B1x_l[1],xright=B1x_l[2], ybottom=0,ytop=1, col=rgb(1,0,0,0.3), border = NA)
  rect(xleft = B2x_l[1],xright=B2x_l[2], ybottom=0,ytop=1, col=rgb(0,0,1,0.3), border = NA)
  points(cu1x, type="l", lwd=2, col="red")
  points(cu2x, type="l", lwd=2, col="blue")

  par(fig=c(0.55,1,0.34,0.67),new=TRUE)
  #plot y
  plot(0,type="n",xlim=c(ymin,ymax),ylim=c(0,1),yaxp= c(0,1,2), ylab="", xlab="",cex.axis=0.8)
  title(xlab = "y axis", line=2)
  rect(xleft = B1y_l[1],xright=B1y_l[2], ybottom=0,ytop=1, col=rgb(1,0,0,0.3), border = NA)
  rect(xleft = B2y_l[1],xright=B2y_l[2], ybottom=0,ytop=1, col=rgb(0,0,1,0.3), border = NA)
  points(cu1y, type="l", lwd=2, col="red")
  points(cu2y, type="l", lwd=2, col="blue")

  par(fig=c(0.55,1,0.10,0.43),new=TRUE)
  #plot z
  plot(0,type="n",xlim=c(zmin,zmax),ylim=c(0,1),yaxp= c(0,1,2), ylab="", xlab="",cex.axis=0.8)
  title(xlab = "z axis", line=2)
  rect(xleft = B1z_l[1],xright=B1z_l[2], ybottom=0,ytop=1, col=rgb(1,0,0,0.3), border = NA)
  rect(xleft = B2z_l[1],xright=B2z_l[2], ybottom=0,ytop=1, col=rgb(0,0,1,0.3), border = NA)
  points(cu1z, type="l", lwd=2, col="red")
  points(cu2z, type="l", lwd=2, col="blue")
  #reset screen
  par(fig=c(0,1,0,1))
  #export if requested
  if(export==TRUE){
    save_pdf(name=paste(name,".pdf"),height =8,width =11 )
    boot_stat <- matrix(c(B1x_l[1],B1x_l[2],B2x_l[1],B2x_l[2],
                          B1y_l[1],B1y_l[2],B2y_l[1],B2y_l[2],
                          B1z_l[1],B1z_l[2],B2z_l[1],B2z_l[2]),nrow = 3,byrow = TRUE)
    rownames(boot_stat) <- c("x","y","z")
    colnames(boot_stat) <- c("mode1_L","mode1_H","mode2_L","mode2_H")
    write.csv(round(boot_stat,digits = 3),paste(name,"bootstrap_stat.csv"),row.names = TRUE)
    cat("Figure saved as ",paste(name,".pdf"),"
")
  }
}

#function that rotate Longitude and Latitude around a Euler pole
rot_DI <- function(Lonlat,P_long=0,P_lat=90,rot=0){
  #degree to radians and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #functions spherical (long=x, lat=y) to Cartesian
  s2c1 <- function(x,y) {cos(d2r(x))*cos(d2r(y))}
  s2c2 <- function(x,y) {sin(d2r(x))*cos(d2r(y))}
  s2c3 <- function(y) {sin(d2r(y))}

  #functions converting cartesian to spherical
  c2sLon <- function(C1,C2) {r2d(atan2(C2,C1))}
  c2sLat <- function(C1,C2,C3) {r2d(asin(C3/(sqrt((C1^2)+(C2^2)+(C3^2)))))}

  #define cos and sin of rotation for simplicity
  cosrot <- cos(d2r(rot))
  sinrot <- sin(d2r(rot))
  #build rotation matrix
  R <- matrix(c((s2c1(P_long,P_lat)^2)*(1-cosrot)+cosrot,
                s2c1(P_long,P_lat)*s2c2(P_long,P_lat)*(1-cosrot)-s2c3(P_lat)*sinrot,
                s2c1(P_long,P_lat)*s2c3(P_lat)*(1-cosrot)+s2c2(P_long,P_lat)*sinrot,
                s2c2(P_long,P_lat)*s2c1(P_long,P_lat)*(1-cosrot)+s2c3(P_lat)*sinrot,
                (s2c2(P_long,P_lat)^2)*(1-cosrot)+cosrot,
                s2c2(P_long,P_lat)*s2c3(P_lat)*(1-cosrot)-s2c1(P_long,P_lat)*sinrot,
                s2c3(P_lat)*s2c1(P_long,P_lat)*(1-cosrot)-s2c2(P_long,P_lat)*sinrot,
                s2c3(P_lat)*s2c2(P_long,P_lat)*(1-cosrot)+s2c1(P_long,P_lat)*sinrot,
                (s2c3(P_lat)^2)*(1-cosrot)+cosrot),
              nrow = 3,ncol = 3,byrow = F)
  #creates result file
  Lonlat_R <- data.frame(matrix(ncol = 2,nrow = 0))
  #apply rotation to all data
  for(i in 1:nrow(Lonlat)){
    C <- matrix(c(s2c1(Lonlat[i,1],Lonlat[i,2]),
                  s2c2(Lonlat[i,1],Lonlat[i,2]),
                  s2c3(Lonlat[i,2])),
                nrow = 3,ncol = 1)
    C_rot <- R%*%C
    Lonlat_R_temp <- data.frame(t(c(c2sLon(C_rot[1,1],C_rot[2,1])%%360,
                                    c2sLat(C_rot[1,1],C_rot[2,1],C_rot[3,1]))))
    Lonlat_R <- rbind(Lonlat_R,Lonlat_R_temp)
  }
  colnames(Lonlat_R) <- c("Long_R","Lat_R")
  return(Lonlat_R)
}

#pdf printing standard size
save_pdf <- function(name="Figure.pdf",width=11,height=8){
  dev.print(pdf,name,width = width, height = height)
}

#plot spherical ortographic projection centered in specified coordinates
sph_ortho <- function(lat=90,long=0,grid=30,coast=FALSE, title="") {
  #functions converting degree and radians
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #functions converting long & lat to xy
  #center of proj is Lon0 & Lat0
  lon0 <- long
  lat0 <- lat
  c2x <- function(lon,lat) {cos(d2r(lat))*sin(d2r(lon-lon0))}
  c2y <- function(lon,lat) {(cos(d2r(lat0))*sin(d2r(lat)))-(sin(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}
  #cut is cosin of c, when negative is behind projections, needs to be cut
  cut <- function(lon,lat) {(sin(d2r(lat0))*sin(d2r(lat)))+(cos(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}

  #fix frame
  plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1,
       xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
  #plot coastline if true
  if(coast==TRUE){
    cst <- world_coastline
    colnames(cst) <- c("lon","lat")
    cst$x <- ifelse(cut(cst$lon,cst$lat)<0,NA,c2x(cst$lon,cst$lat))
    cst$y <- ifelse(cut(cst$lon,cst$lat)<0,NA,c2y(cst$lon,cst$lat))
    lines(cst$x,cst$y,col="black",lwd=0.75)
  }
  #if grid==0 does not plot parallels and meridians
  if(grid!=0){
    #plot_main_parallel
    #longitude circle
    lats <- seq(-(90-grid),(90-grid),grid)
    for(i in lats){
      lon_lat_p <-  as.data.frame(0:360)
      lon_lat_p$lat <- rep(i)
      lon_lat_p$x <- ifelse(cut(lon_lat_p[,1],lon_lat_p[,2])<0,NA,
                            c2x(lon_lat_p[,1],lon_lat_p[,2]))
      lon_lat_p$y <- ifelse(cut(lon_lat_p[,1],lon_lat_p[,2])<0,NA,
                            c2y(lon_lat_p[,1],lon_lat_p[,2]))
      lines(lon_lat_p$x,lon_lat_p$y,col="gray", pch=16, cex=0.3)
    }

    #plot_main_meridians
    lons <- seq(grid,360,grid)
    for(i in lons){
      lat_lon_m <- as.data.frame(seq(-89,89,1))
      lat_lon_m$lon <- rep(i)
      lat_lon_m$x <- ifelse(cut(lat_lon_m[,2],lat_lon_m[,1])<0,NA,
                            c2x(lat_lon_m[,2],lat_lon_m[,1]))
      lat_lon_m$y <- ifelse(cut(lat_lon_m[,2],lat_lon_m[,1])<0,NA,
                            c2y(lat_lon_m[,2],lat_lon_m[,1]))
      lines(lat_lon_m$x,lat_lon_m$y,col="gray", pch=16, cex=0.3)
    }
  }

  #plot black frame around globe
  a2cx <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*sin(d2r(y))}
  a2cy <- function(x,y) {sqrt(2)*sin((d2r(90-x))/2)*cos(d2r(y))}
  #plot external circle
  frame_dec = 0:360
  frame_inc=rep(0,length(frame_dec))
  x = a2cx(frame_inc,frame_dec)
  y = a2cy(frame_inc,frame_dec)
  lines(x, y, col = "black")
  title(xlab = title, line=0.2, cex=0.1)
}

#function that deform dec Inc from eigenvectors
strain_DI <- function(DIAP,M,export=FALSE,name="strained_dirs"){
  library(matlib)
  #degree to radians and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #load directions data
  data <- DIAP
  data <- na.omit(data)
  dirs <- data[,1:2]
  bed <- data[,3:4]
  colnames(dirs) <- c("dec", "inc")
  colnames(bed) <- c("b_az","b_plunge")

  #directions in Cartesian coordinates
  dirs$x <- cos(d2r(dirs$dec))*cos(d2r(dirs$inc))
  dirs$y <- sin(d2r(dirs$dec))*cos(d2r(dirs$inc))
  dirs$z <- sin(d2r(dirs$inc))

  #bedding in cartesian coordinates
  bed$x <- cos(d2r(bed$b_az))*cos(d2r(bed$b_plunge))
  bed$y <- sin(d2r(bed$b_az))*cos(d2r(bed$b_plunge))
  bed$z <- sin(d2r(bed$b_plunge))

  new_DI <- data.frame(matrix(ncol=2,nrow=0))
  new_bed <- data.frame(matrix(ncol=2,nrow=0))

  for (i in 1:nrow(data)){
    #strain dirs
    dircart <- t(as.matrix(dirs[i,3:5]))
    strain <- M%*%dircart
    NewInc <- r2d(asin(strain[3,1]/(sqrt((strain[1,1]^2)+(strain[2,1]^2)+(strain[3,1]^2)))))
    NewDec <- r2d(atan2(strain[2,1],strain[1,1]))
    NewDec <- NewDec%%360

    new_DI_p <- cbind(as.data.frame(NewDec),as.data.frame(NewInc))
    new_DI <- rbind(new_DI,new_DI_p)

    #strain bed
    bedcart <- t(as.matrix(bed[i,3:5]))
    b_strain <- M%*%bedcart
    New_pl <- r2d(asin(b_strain[3,1]/(sqrt((b_strain[1,1]^2)+(b_strain[2,1]^2)+(b_strain[3,1]^2)))))
    New_az <- r2d(atan2(b_strain[2,1],b_strain[1,1]))
    New_az <- New_az%%360

    new_bed_p <- cbind(as.data.frame(New_az),as.data.frame(New_pl))
    new_bed <- rbind(new_bed,new_bed_p)
  }
  colnames(new_DI) <- c("str_dec","str_inc")
  colnames(new_bed) <- c("str_B_az","str_Binc")
  str_data <- cbind(new_DI,new_bed)
  if(export==TRUE){write.csv(round(str_data,digits = 2),paste(name,".csv"),row.names = FALSE)}
  return(str_data)
}

#function generating E from TK03.GAD model from given inclination
tk03 <- function(I) {
  x <- sqrt(I^2)
  y <- 2.895-(0.01466*x)-(0.0003525*(x^2))+(0.00000316*(x^3)) #E-I equation
  return(y)
}

#unflat directions with given f
unflat_DI <- function(DI,f,export=FALSE,name="unflattened_dirs.csv") {
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  DI[,2]<- r2d(atan(tan(d2r(DI[,2]))/f))
  if(export==TRUE){write.csv(round(DI,digits = 2),paste(name,".csv"),row.names = FALSE)}
  return(DI)
}

#unstrain bootstrpped pseudosamples of directions
unstr_boot <- function(unstr_file,nb= 100,S_vec,Lin,Fol,ns=1,confidence=95,cross=FALSE,EdMAX=FALSE,EdMIN=FALSE,hist=TRUE,save=TRUE,name="Unstrain_bootstrap"){
  cat("
!!Unstrain of pseudosamples is SLOW!!

")
  #degree to radians and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #load directions data
  dat <- unstr_file[[1]]
  dat <- na.omit(dat)
  colnames(dat) <- c("dec", "inc","baz","binc")

  #load inc_e_dec result
  inc_e_dec <- unstr_file[[5]]

  #plot frame
  par(fig=c(0,0.7,0,1), new= FALSE)
  y_up <- ifelse(inc_e_dec[1,2]>3.5, ceiling(inc_e_dec[1,2]),3.5)

  plot(NA, xlim= c(0,90), ylim= c(1,y_up), xaxp= c(0,90, 9),
       xlab="Inclination (°)", ylab="Elongation")

  #plot tk03.GAD model E-I
  x <- 0:90
  y <- tk03(x)
  points(x=x, y= y, type= "l", col="blue", lwd=3)

  #plot unstrain curve
  points(x=inc_e_dec[,1], y= inc_e_dec[,2], type= "l", col="red", lwd=3)
  points(x=inc_e_dec[1,1],y=inc_e_dec[1,2],
         pch=21, col="black", bg="blue", cex=1.2)
  points(x=inc_e_dec[length(inc_e_dec[,1]),1],y=inc_e_dec[length(inc_e_dec[,1]),2],
         pch=21, col="black", bg="red", cex=1.2)

  #text for figure
  N <- length(dat[,1])
  Inc <- round(inc_e_dec[1,1],digits = 1)
  Ecut <- round(inc_e_dec[1,2],digits = 2)
  Edec <- round(inc_e_dec[1,3],digits = 1)
  text <- paste("N:", N, "
Inc:", Inc, "
E:", Ecut,"
Edec:",Edec)
  text(x=0, y=3.2,pos=4,text, cex= 0.8)

  #unstrained parameters for text2 of figure
  inc_nstr <- round(inc_e_dec[length(inc_e_dec[,1]),1],digits = 1)
  E_nstr <- round(inc_e_dec[length(inc_e_dec[,2]),2], digits=2)
  Edec_nstr <- round(inc_e_dec[length(inc_e_dec[,3]),3], digits=1)

  #put text unstrained in figure
  text2 <- paste("Unstrained","
Inc:", inc_nstr, "
E:", E_nstr, "
Edec:", Edec_nstr)
  text(x=20, y=3.2, pos=4, text2, cex=0.8)

  #create empty files
  init_E_I <- as.data.frame(matrix(ncol=3,nrow=0))
  final_E_I <- as.data.frame(matrix(ncol=3,nrow=0))
  colnames(init_E_I) <- c("Inc","E","E_dec")
  colnames(final_E_I)<- c("Inc","E","E_dec")

  #start bootstrapping
  n <- 0
  repeat{
    n <- n+1
    Seq_I_E_b <- as.data.frame(matrix(ncol=3,nrow=0))
    E_declin <- as.data.frame(matrix(ncol=1,nrow=0))
    data <- boots_DI(dat)
    dirs <- data[,1:2]
    bed <- data[,3:4]

    #directions in Cartesian coordinates
    dirs$x <- cos(d2r(dirs$dec))*cos(d2r(dirs$inc))
    dirs$y <- sin(d2r(dirs$dec))*cos(d2r(dirs$inc))
    dirs$z <- sin(d2r(dirs$inc))

    #bedding in Cartesian coordinates
    bed$x <- cos(d2r(bed$baz))*cos(d2r(bed$binc))
    bed$y <- sin(d2r(bed$baz))*cos(d2r(bed$binc))
    bed$z <- sin(d2r(bed$binc))

    #set parameters of deforming matrix
    Lincr <- (Lin-1)/ns
    Fincr <- (Fol-1)/ns
    L <- 1
    F <- 1

    #gradual unstrain of the n pseudosample
    repeat{
      #anisotropy degree
      P <- F*L
      #eigenvalues
      K1 <- (3*L)/(L+(1/F)+1)
      K2 <- K1/L
      K3 <- K1/P

      #matrix of new eigenvalue
      M <- c(K1,0,0,0,K2,0,0,0,K3)
      M <- matrix(M,nrow=3,byrow=T)

      #combines given eigenvalues with Strain directions
      S <- S_vec%*%M%*%inv(S_vec)
      new_DI <- data.frame(matrix(ncol=4,nrow=0))
      for (i in 1:length(data[,1])){
        #unstrain dirs
        dircart <- t(as.matrix(dirs[i,3:5]))
        unstr <- S%*%dircart
        NewInc <- r2d(asin(unstr[3,1]/(sqrt((unstr[1,1]^2)+(unstr[2,1]^2)+(unstr[3,1]^2)))))
        NewDec <- r2d(atan2(unstr[2,1],unstr[1,1]))
        NewDec <- ifelse(NewDec<0,NewDec+360,NewDec)

        #unstrain bedding
        bedcart <- t(as.matrix(bed[1,3:5]))
        bunstr <- S%*%bedcart
        Newbinc <- r2d(asin(bunstr[3,1]/(sqrt((bunstr[1,1]^2)+(bunstr[2,1]^2)+(bunstr[3,1]^2)))))
        Newbaz <- r2d(atan2(bunstr[2,1],bunstr[1,1]))
        Newbaz <- ifelse(Newbaz<0,Newbaz+360,Newbaz)

        new_DI_p <- cbind(as.data.frame(NewDec),as.data.frame(NewInc),
                          as.data.frame(Newbaz),as.data.frame(Newbinc))
        new_DI <- rbind(new_DI,new_DI_p)
      }
      colnames(new_DI) <- c("dec","inc","baz","binc")
      NBdecInc <- bed_DI(new_DI)
      inc_e_dec_p <- inc_E_finder(NBdecInc)
      Seq_I_E_b <- rbind(Seq_I_E_b,inc_e_dec_p)
      #take only absolute value of declination for the breaking conditions
      E_declin <- rbind(E_declin,abs(inc_e_dec_p[1,3]))

      #if cross== TRUE breaks loop when tk03.GAD is crossed
      if(cross==TRUE)  {if(any(Seq_I_E_b$E>tk03(Seq_I_E_b$V1inc))
                           && any(Seq_I_E_b$E<tk03(Seq_I_E_b$V1inc))) break}

      #if EdMAX==TRUE breaks loop when max of Edec is reached
      if(EdMAX==TRUE && length(E_declin[,1])>3){
        l <- length(E_declin[,1])
        if(E_declin[l,1]<E_declin[l-1,1] && E_declin[l-1,1]>E_declin[l-2,1]){    ##HERE
          E_declin <- E_declin[-l,]
          break
        }
      }
      #break loops when minimu Edec reached
      if(EdMIN==TRUE && length(E_declin[,1])>3){
        l <- length(E_declin[,1])
        if(E_declin[l,1]>E_declin[l-1,1] && E_declin[l-1,1]<E_declin[l-2,1]){
          E_declin <- E_declin[-l,]
          break
        }
      }

      #break loop if reaches defined F and L
      if(round(L,digits = 4) == round(Lin,digits = 4) &&
         round(F,digits = 4) == round(Fol,digits = 4)) break
      #lineation
      L <- L+Lincr
      #foliation
      F <- F+Fincr
    }

    #if cross==TRUE select only curves that cross tk03.GAD
    if(cross==TRUE){
      if(any(Seq_I_E_b$E>tk03(Seq_I_E_b$V1inc)) && any(Seq_I_E_b$E<tk03(Seq_I_E_b$V1inc))==TRUE){
        points(x=Seq_I_E_b$V1inc, y= Seq_I_E_b$E,
               type= "l", col=rgb(1, 0, 0, 0.15), lwd=1)
        i_E_I <- Seq_I_E_b[1,]
        f_E_I <- Seq_I_E_b[length(Seq_I_E_b[,1]),]
      }
      #if cross== FALSE, use any other condition to fill files
    }else{
      points(x=Seq_I_E_b$V1inc, y= Seq_I_E_b$E,
             type= "l", col=rgb(0, 0.5, 1, 0.15), lwd=1)
      i_E_I <- Seq_I_E_b[1,]
      f_E_I <- Seq_I_E_b[length(Seq_I_E_b[,1]),]
    }
    colnames(i_E_I) <- c("Inc","E","E_dec")
    colnames(f_E_I)<- c("Inc","E","E_dec")
    init_E_I <- rbind(init_E_I,i_E_I)
    final_E_I <- rbind(final_E_I,f_E_I)
    init_E_I <- na.omit(init_E_I)
    final_E_I <- na.omit(final_E_I)

    #message for bootstrapping update
    if(((n%%10)==0)==TRUE) {
      cat(paste(n,"simulations done and",(length(final_E_I[,1])),"pseudosamples saved
"))
    }
    if(length(final_E_I[,1])==nb) {
      cat(paste("Saved",(length(final_E_I[,1])), "pseudosamples after", n,"simulations
"))
      break
    }
  }
  #backup_files
  final_E_Ibk <- final_E_I
  init_E_Ibk <- init_E_I

  #Plot Bootstrapped data
  if(cross==FALSE){
    points(x=final_E_I$Inc,
           y=final_E_I$E,
           pch=16,
           col=rgb(1,0,0,0.7),
           cex=0.5)
  }

  #cut bootstrapped results for 95% (unless different) confidence
  conf <- confidence/100
  num <- round((nb*(1-conf))/2,digits=0)
  Lconf <- num
  Uconf <- nb-num
  final_E_I <- final_E_I[order(final_E_I$Inc),]
  final_E_I <- final_E_I[Lconf:Uconf,]

  points(x=x, y= y, type= "l", col="blue", lwd=3)

  #plot unstrain curve
  points(x=inc_e_dec[,1], y= inc_e_dec[,2], type= "l", col="yellow", lwd=2)
  points(x=inc_e_dec[1,1],y=inc_e_dec[1,2],
         pch=21, col="black", bg="blue", cex=1.2)
  #plot cross with EI line
  points(x=inc_e_dec[length(inc_e_dec[,1]),1],y=inc_e_dec[length(inc_e_dec[,1]),2],
         pch=21, col="black", bg=rgb(1,0.4,0.4,1), cex=1.2)

  #draw two lines for 95% confidence margin
  arrows(x0=final_E_I[1,1],x1=final_E_I[1,1],
         y0=1,y1=1.1, lwd=1.5, length = 0)

  arrows(x0=final_E_I[length(final_E_I$Inc),1],
         x1=final_E_I[length(final_E_I$Inc),1],
         y0= 1, y1= 1.1,lwd=1.5, length = 0)

  Inc_l <- round(final_E_I[1,1], digits= 1)
  Inc_u <- round(final_E_I[length(final_E_I$Inc),1], digits=1)

  text(x=final_E_I[1,1], y=1, pos=2, Inc_l)
  text(x=final_E_I[length(final_E_I$Inc),1], y=1, pos=4, Inc_u)

  #plot text results again if covered by red lines
  if(inc_e_dec[1,1]<40){
    text(x=0, y=3.2,pos=4,text, cex= 0.8)
    text(x=20, y=3.2, pos=4, text2, cex=0.8)
  }
  #recalculate boostrtapped confidence for Edec for plotting confidence margin
  final_E_I_Edec <- final_E_Ibk
  colnames(final_E_I_Edec) <- c("Inc","E","E_dec")
  final_E_I_Edec <- final_E_I_Edec[order(final_E_I_Edec$E_dec),]
  final_E_I_Edec <- final_E_I_Edec[Lconf:Uconf,]

  #plot histogram of E_declination with respect V1 before and after correction
  if(hist==TRUE){
    par(fig=c(0.35,0.69,0.56,0.99), new=TRUE)
    hist(init_E_Ibk$E_dec, xlim=c(-90,90), breaks= 90,
         axes=FALSE,xlab="",ylab="",col="blue", border="blue", main="")
    par(fig=c(0.35,0.69,0.56,0.99), new=TRUE)
    hist(final_E_Ibk$E_dec, xlim=c(-90,90), xaxp=c(-90,90,4),
         breaks = 90, xlab = "", ylab = "",
         main="", cex.axis=0.8,col="red",border ="red")
    abline(v=final_E_I_Edec[1,3],lwd=1, lty=2)
    abline(v=final_E_I_Edec[length(final_E_I_Edec[,3]),3],lwd=1, lty=2)

    #plot lables closer than standard to axes
    title(xlab = "Edec(°)", line=1.9)
    title(ylab = "Frequency", line=1.9)
  }

  #plot original directions single mode tilt corrected
  par(fig=c(0.55,1,0,0.6), new=TRUE)
  plot_DI(unstr_file[[2]],single_mode = TRUE, title="Original directions")

  #plot unstrained directions tilt corrected
  par(fig=c(0.55,1,0.4,1), new=TRUE)
  plot_DI(unstr_file[[3]], single_mode = TRUE, col_d = "red",col_u = "pink", title = "Unstrained directions")

  #build result file
  stat <- as.data.frame(matrix(ncol=1,nrow=1))
  colnames(stat) <- c("N")
  stat$N <- N
  stat$Inc <- Inc
  stat$E <- Ecut
  stat$Edec <- Edec
  stat$Inc_unstr <- inc_nstr
  stat$Low_inc <- round(final_E_I[1,1], digits = 1)
  stat$Hign_inc <- round(final_E_I[length(final_E_I[,1]),1], digits = 1)
  stat$E_unstr <- E_nstr
  stat$Edec_unstr <- Edec_nstr
  stat$Edec_low <- round(final_E_I_Edec[1,3],digits = 1)
  stat$Edec_high <- round(final_E_I_Edec[length(final_E_I_Edec[,3]),3],digits = 1)

  #save and export results if save=TRUE
  if(save==TRUE){
    cat("
Statistics saved as csv file.
Figure saved as pdf file.

      ")
  }
  #restore screen
  par(fig=c(0,1,0,1))
  if(save==TRUE){
    write.csv(stat,paste(name,"_statistic.csv"), row.names=FALSE)
    save_pdf(paste(name,".pdf"))
  }
}

#file DI require also 3 and 4 columns with Bed_az and bed_plunge, S_vec= strain matrix
unstr_DI <- function(DIAP,S_vec,Lin,Fol,n=1,cross=FALSE,EdMAX=FALSE,EdMIN=FALSE,save=TRUE,name="Unstrain"){
  #degree to radians and vice versa
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #load directions data
  data <- DIAP
  data <- na.omit(data)
  colnames(data) <- c("dec", "inc","baz","binc")
  dirs <- data[,1:2]
  bed <- data[,3:4]

  #directions in Cartesian coordinates
  dirs$x <- cos(d2r(dirs$dec))*cos(d2r(dirs$inc))
  dirs$y <- sin(d2r(dirs$dec))*cos(d2r(dirs$inc))
  dirs$z <- sin(d2r(dirs$inc))

  #bedding in Cartesian coordinates
  bed$x <- cos(d2r(bed$baz))*cos(d2r(bed$binc))
  bed$y <- sin(d2r(bed$baz))*cos(d2r(bed$binc))
  bed$z <- sin(d2r(bed$binc))

  inc_e_dec <- as.data.frame(matrix(ncol=3,nrow = 1))
  BdecInc <- bed_DI(data[,1:4])
  inc_e_dec[1,1:3] <- as.data.frame(inc_E_finder(BdecInc))
  inc_e_dec[1,3] <- abs(inc_e_dec[1,3])
  colnames(inc_e_dec) <- c("V1inc","E","DV1V2")

  #set parameters of deforming matrix
  Lincr <- (Lin-1)/n
  Fincr <- (Fol-1)/n
  L <- 1
  F <- 1
  repeat{
    #lineation
    L <- L+Lincr
    #foliation
    F <- F+Fincr
    #anisotropy degree
    P <- F*L
    #eigenvalues
    K1 <- (3*L)/(L+(1/F)+1)
    K2 <- K1/L
    K3 <- K1/P

    #matrix of new eigenvalue
    M <- c(K1,0,0,0,K2,0,0,0,K3)
    M <- matrix(M,nrow=3,byrow=T)

    #combines given eigenvalues with Strain directions
    S <- S_vec%*%M%*%inv(S_vec)
    S_e <- eigen(S,symmetric = T)
    S_vec <- S_e$vectors
    S_val <- S_e$values
    new_DI <- data.frame(matrix(ncol=4,nrow=0))
    for (i in 1:nrow(data)){
      #unstrain dirs
      dircart <- t(as.matrix(dirs[i,3:5]))
      unstr <- S%*%dircart
      NewInc <- r2d(asin(unstr[3,1]/(sqrt((unstr[1,1]^2)+(unstr[2,1]^2)+(unstr[3,1]^2)))))
      NewDec <- r2d(atan2(unstr[2,1],unstr[1,1]))
      NewDec <- NewDec%%360
      #unstrain bedding
      bedcart <- t(as.matrix(bed[1,3:5]))
      bunstr <- S%*%bedcart
      Newbinc <- r2d(asin(bunstr[3,1]/(sqrt((bunstr[1,1]^2)+(bunstr[2,1]^2)+(bunstr[3,1]^2)))))
      Newbaz <- r2d(atan2(bunstr[2,1],bunstr[1,1]))
      Newbaz <- Newbaz%%360

      new_DI_p <- cbind(as.data.frame(NewDec),as.data.frame(NewInc),
                        as.data.frame(Newbaz),as.data.frame(Newbinc))
      new_DI <- rbind(new_DI,new_DI_p)
    }
    colnames(new_DI) <- c("dec","inc","baz","binc")
    NBdecInc <- bed_DI(new_DI)
    inc_e_dec_p <- inc_E_finder(NBdecInc)
    inc_e_dec <- rbind(inc_e_dec,inc_e_dec_p)
    inc_e_dec[,3] <- abs(inc_e_dec[,3])

    #generate tk03 curve also for later plot
    tkx <- 0:90
    tky <- tk03(tkx)
    if (cross==TRUE){
      #break loop if crosses tk03 line
      if(any(inc_e_dec$E>tk03(inc_e_dec$V1inc)) && any(inc_e_dec$E<tk03(inc_e_dec$V1inc))){
        curve1 <- inc_e_dec[,1:2]
        curve1 <- curve1[order(curve1$E),]
        curve2 <- cbind(as.data.frame(tkx),as.data.frame(tky))
        crxy <- curve_cross(curve1,curve2)
        cat(paste(" Final L: ",round(L,digits = 5),"
", "Final F:",round(F,digits = 5),"
","Final P:",round(F*L,digits = 5)))
        break
      }
    }
    #break loops when maximum Edec reached
    if(EdMAX==TRUE && length(inc_e_dec[,3])>3){
      l <- length(inc_e_dec[,3])
      if(inc_e_dec[l,3]<inc_e_dec[l-1,3] && inc_e_dec[l-1,3]>inc_e_dec[l-2,3]){
        inc_e_dec <- inc_e_dec[-l,]
        break
      }
    }
    #break loops when minimu Edec reached
    if(EdMIN==TRUE && length(inc_e_dec[,3])>3){
      l <- length(inc_e_dec[,3])
      if(inc_e_dec[l,3]>inc_e_dec[l-1,3] && inc_e_dec[l-1,3]<inc_e_dec[l-2,3]){
        inc_e_dec <- inc_e_dec[-l,]
        break
      }
    }
    #break loop if reaches defined F and L
    if(round(L,digits = 4) == round(Lin,digits = 4) &&
       round(F,digits = 4) == round(Fol,digits = 4)) break
  }
  #set figure
  par(fig=c(0,0.7,0,1), new= FALSE)
  y_up <- ifelse(inc_e_dec[1,2]>3.5, ceiling(inc_e_dec[1,2]),3.5)

  plot(NA, xlim= c(0,90), ylim= c(1,y_up), xaxp= c(0,90, 9),
       xlab="Inclination (°)", ylab="Elongation")

  #text for figure
  N <- as.character(length(data[,1]))
  Inc <- round(inc_e_dec[1,1],digits = 1)
  Ecut <- round(inc_e_dec[1,2],digits = 2)
  V2 <- round(inc_e_dec[1,3],digits = 1)
  text <- paste("N:", N, "
Inc:", Inc, "
E:", Ecut,"
Edec:",V2)
  text(x=0, y=3.2,pos=4,text, cex= 0.8)

  #plot tk03.GAD model E-I
  points(x=tkx, y= tky, type= "l", col="blue", lwd=3)
  #plot unstrain curve
  points(x=inc_e_dec[,1], y= inc_e_dec[,2], type= "l", col="red", lwd=3)
  points(x=inc_e_dec[1,1],y=inc_e_dec[1,2],
         pch=21, col="black", bg="blue", cex=1.2)
  #plot cross with EI line if exists
  if(exists("crxy")==TRUE){
    points(x=crxy[1,1],y=crxy[1,2],
           pch=21, col="black", bg="red", cex=1.2)
  }else{
    points(x=inc_e_dec[length(inc_e_dec[,1]),1],
           y=inc_e_dec[length(inc_e_dec[,1]),2],
           pch=21, col="black", bg="red", cex=1.2)
  }
  #unstrained parameters for text2 of figure
  inc_nstr <- ifelse(exists("crxy")==TRUE,
                     round(crxy[1,1], digits=1),
                     round(inc_e_dec[length(inc_e_dec[,1]),1],digits = 1))

  E_nstr <- round(inc_e_dec[length(inc_e_dec[,2]),2], digits=2)
  Edec_nstr <- round(inc_e_dec[length(inc_e_dec[,3]),3], digits=1)

  #put text unstrained in figure
  text2 <- paste("Unstrained","
Inc:", inc_nstr, "
E:", E_nstr, "
Edec:", Edec_nstr)
  text(x=20, y=3.2, pos=4, text2, cex=0.8)

  #plot Edec during unstrain
  par(fig=c(0.35,0.69,0.56,0.99), new=TRUE)
  plot(x=inc_e_dec$V1inc, y=inc_e_dec$DV1V2, type="l", tck=-0.05,
       frame.plot=FALSE,cex.axis=0.8, mgp=c(1.5,0.5,0.1),
       xlab="",ylab="",yaxt="n",
       ylim=c(floor(min(inc_e_dec$DV1V2)),
                    ceiling(max(inc_e_dec$DV1V2))),
       xlim=c(floor(min(inc_e_dec$V1inc)),
              ceiling(max(inc_e_dec$V1inc))))
  axis(2, cex.axis=0.8, las=2)
  title(xlab = "Inc (°)", line=1.6)
  title(ylab = "Edec (°)", line=2.5)

  points(x=inc_e_dec[1,1],y=inc_e_dec[1,3],
         pch=21, col="black", bg="blue", cex=1.2)
  points(x=inc_e_dec[length(inc_e_dec[,1]),1],
         y=inc_e_dec[length(inc_e_dec[,1]),3],
         pch=21, col="black", bg="red", cex=1.2)

  #sets results file
  unstr_matrix <- S
  unstr_results <- list(data,BdecInc,NBdecInc, new_DI,inc_e_dec,unstr_matrix)
  names(unstr_results) <- c("Original dataset","original TC directions","unstrained TC directions",
                            "unstrained directions and bedding",
                            "inc, E, declination triplets","Unstrain matrix")
  #plot original direction single mode tilt corrected
  par(fig=c(0.55,1,0,0.6), new=TRUE)
  plot_DI(BdecInc,single_mode = TRUE, title="Original directions")

  #plot unstrained directions tilt corrected
  par(fig=c(0.55,1,0.4,1), new=TRUE)
  plot_DI(NBdecInc, single_mode = TRUE, col_d = "red",col_u = "pink", title = "Unstrained directions")

  #save fig
  if(save==TRUE){
    save_pdf(name = paste(name,".pdf"))
    write.csv(unstr_results[[3]],paste(name,"_unstrained TC directions.csv"),
              row.names = FALSE)
    cat("
Unstrained directions saved as csv file.
Figure saved as pdf file.

")
  }
  #restore screen
  par(fig=c(0,1,0,1))
  return(unstr_results)
}

#makes bootstrap stats of EI of initial and final distribution
unstr_stat <- function(unstr_file, nb=1000,confidence=95,hist=TRUE, export=TRUE,name="bootstrap_stat"){
  DI <- unstr_file[[2]]
  NDI <- unstr_file[[3]]
  inc_e_dec <- unstr_file[[5]]
  Inc_E <- as.data.frame(matrix(ncol=3,nrow=0))
  I_E_D0 <- inc_E_finder(DI)
  I_E_Df <- inc_E_finder(NDI)
  NInc_E <- as.data.frame(matrix(ncol=3,nrow=0))

  #bootstrap of original directions
  cat(paste("
",nb, "simulations of original directions:
"))
  for (i in 1:nb) {
    dataprov <- boots_DI(DI)
    I_E_Ed <- inc_E_finder(dataprov)
    I_E_Ed$V1inc <- abs(I_E_Ed$V1inc)
    Inc_E <- rbind(Inc_E, I_E_Ed)
    if(((i%%50)==0)==TRUE) {
      cat(paste(i,"simulations out of",nb,"done
"))
    }
  }
  #bootstrap of unstrained directions
  cat(paste("
",nb, "simulations of unstrained directions:
"))
  for (i in 1:nb) {
    dataprov <- boots_DI(NDI)
    I_E_Ed <- inc_E_finder(dataprov)
    I_E_Ed$V1inc <- abs(I_E_Ed$V1inc)
    NInc_E <- rbind(NInc_E, I_E_Ed)
    if(((i%%50)==0)==TRUE) {
      cat(paste(i,"simulations out of",nb,"done
"))
    }
  }
  #names column
  colnames(Inc_E) <- c("Inc","E","E_dec")
  colnames(NInc_E) <- c("Inc","E","E_dec")
  #makes backup
  Inc_E_bk <- Inc_E
  NInc_E_bk <- NInc_E
  #order for E
  Inc_E <- Inc_E[order(Inc_E$E),]
  NInc_E <- NInc_E[order(NInc_E$E),]
  #trims files for confidence
  conf <- confidence/100
  num <- round((nb*(1-conf))/2,digits=0)
  Lconf <- num
  Uconf <- nb-num
  Inc_E <- Inc_E[Lconf:Uconf,]
  NInc_E <- NInc_E[Lconf:Uconf,]

  #plot frame
  par(fig=c(0,0.7,0,1), new= FALSE)
  y_up <- ifelse(I_E_D0[1,2]>3.5, ceiling(I_E_D0[1,2]),3.5)

  plot(NA, xlim= c(0,90), ylim= c(1,y_up), xaxp= c(0,90, 9),
       xlab="Inclination (°)", ylab="Elongation")

  #text for figure
  N <- length(DI[,1])
  Inc <- round(inc_e_dec[1,1],digits = 1)
  Ecut <- round(inc_e_dec[1,2],digits = 2)
  V2 <- round(inc_e_dec[1,3],digits = 1)
  text <- paste("N:", N, "
Inc:", Inc, "
E:", Ecut,"
Edec:",V2)
  text(x=0, y=3.2,pos=4,text, cex= 0.8)

  #unstrained parameters for text2 of figure
  inc_nstr <- round(inc_e_dec[length(inc_e_dec[,1]),1],digits = 1)
  E_nstr <- round(inc_e_dec[length(inc_e_dec[,2]),2], digits=2)
  Edec_nstr <- round(inc_e_dec[length(inc_e_dec[,3]),3], digits=1)

  #put text unstrained in figure
  text2 <- paste("Unstrained","
Inc:", inc_nstr, "
E:", E_nstr, "
Edec:", Edec_nstr)
  text(x=20, y=3.2, pos=4, text2, cex=0.8)

  #plot tk03.GAD model E-I
  x <- 0:90
  y <- tk03(x)
  points(x=x, y= y, type= "l", col="blue", lwd=3)

  #Plot Bootstrapped data
  points(x=Inc_E$Inc,
         y=Inc_E$E,
         pch=16,
         col=rgb(0.5, 0, 1, 0.15),
         cex=0.6)

  points(x=NInc_E$Inc,
         y=NInc_E$E,
         pch=16,
         col=rgb(1, 0, 0, 0.15),
         cex=0.6)

  #plot unstrain curve
  points(x=inc_e_dec[,1], y= inc_e_dec[,2], type= "l", col="yellow", lwd=2)
  points(x=inc_e_dec[1,1],y=inc_e_dec[1,2],
         pch=21, col="black", bg="blue", cex=1.2)
  inc_e_final <- inc_e_dec[length(inc_e_dec[,1]),]
  unx <- inc_e_final[1,1]
  uny <- inc_e_final[1,2]
  points(x=unx,y=uny,
         pch=21, col="black", bg="red", cex=1.2)

  #plot original direction single mode tilt corrected
  par(fig=c(0.55,1,0.4,1), new=TRUE)
  plot_DI(NDI, single_mode = TRUE, col_d = "red",col_u = "pink", title = "Unstrained directions")

  #plot unstrained directions tilt corrected
  par(fig=c(0.55,1,0,0.6), new=TRUE)
  plot_DI(DI,single_mode = TRUE, title="Original directions")

  #plot histogram of Dec if requested
  if(hist==TRUE){
    par(fig=c(0.35,0.69,0.56,0.99), new=TRUE)
    hist(Inc_E$E_dec, xlim=c(-90,90), xaxp=c(-90,90,4),
         breaks = 90, xlab = "", ylab = "",
         main="", cex.axis=0.8, col= "blue", border="blue")
  }

  #unstrained declination
  if(hist==TRUE){
    par(fig=c(0.35,0.69,0.56,0.99), new=TRUE)
    hist(NInc_E$E_dec, xlim=c(-90,90), breaks= 90,
         axes=FALSE,xlab="",ylab="",col="red", border="red", main="")
    #plot labels closer than standard to axes
    title(xlab = "Edec(°)", line=1.9, cex=0.3)
    title(ylab = "Frequency", line=2, cex=0.1)
  }
  #plot confidence margin of declination if between -45 and 45
  if(inc_e_final[1,3]>-45 && inc_e_final[1,3]<45){
    Inc_E_Edec <- NInc_E_bk
    Inc_E_Edec <- Inc_E_Edec[order(Inc_E_Edec[,3]),]
    Inc_E_Edec <- Inc_E_Edec[Lconf:Uconf,]
    low_dec <- Inc_E_Edec[1,3]
    up_dec <- Inc_E_Edec[length(Inc_E_Edec[,3]),3]
    abline(v=low_dec,lwd=1,lty=2)
    abline(v=up_dec,lwd=1,lty=2)
  }
  #define low and high inc margins
  inc_conf <- NInc_E_bk
  inc_conf <- inc_conf[order(inc_conf$Inc),]
  inc_conf <- inc_conf[Lconf:Uconf,]
  low_inc <- inc_conf[1,1]
  up_inc <- inc_conf[length(inc_conf[,1]),1]

  if(inc_e_final[1,3]>-45 && inc_e_final[1,3]<45){
    results <- as.data.frame(matrix(ncol= 9, nrow=1))
    colnames(results) <- c("Inc","Low_inc","High_inc", "Elong","Low_E","High_E","E_dec","Low_E_dec","High_E_dec")
    results$Inc <- inc_nstr
    results$Low_inc <- round(low_inc,digits = 1)
    results$High_inc <- round(up_inc,digits = 1)
    results$Elong <- E_nstr
    results$Low_E <- round(min(NInc_E[,2]),digits= 2)
    results$High_E <- round(max(NInc_E[,2]),digits= 2)
    results$E_dec <- Edec_nstr
    results$Low_E_dec <- round(low_dec,digits= 2)
    results$High_E_dec <- round(up_dec,digits=2)
  }else{
    results <- as.data.frame(matrix(ncol= 7, nrow=1))
    colnames(results) <- c("Inc","Low_inc","High_inc", "Elong","Low_E","High_E","E_dec")
    results$Inc <- inc_nstr
    results$Low_inc <- round(low_inc,digits = 1)
    results$High_inc <- round(up_inc,digits = 1)
    results$Elong <- E_nstr
    results$Low_E <- round(min(NInc_E[,2]),digits= 2)
    results$High_E <- round(max(NInc_E[,2]),digits= 2)
    results$E_dec <- Edec_nstr
  }
  print(results, row.names = FALSE)
  par(fig=c(0,1,0,1))
  if(export==TRUE){
    cat("
Results saved as Unstrained_directions_statistic.csv
Graph saved as Unstrained_directions_plot.pdf

")

    #write results file
    write.csv(results,file=paste(name,".csv"),row.names = FALSE)

    #save figure as pdf
    save_pdf(paste(name,".pdf"))
  }
}

# #plot A95 from VGP data and compare with GAPWP
VGP_A95 <- function(VGP,lat=90,long=0,grid=30, auto_cent=TRUE, symbol="c",color="blue",col_A=rgb(1,0,0,0.3), coast=FALSE, on_plot=FALSE, save=FALSE, name="A95",APWP="V23", S_APWP=FALSE){
  library("dplyr", warn.conflicts = FALSE)

  #warning for on-plot, to avoid wrong coordinates
  if(on_plot==TRUE && auto_cent==TRUE) {
    stop("Please SPECIFY center coordinates when on_plot==TRUE",call. = F)
  }

  #functions converting degree and radians
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #functions converting long & lat to xy
  c2x <- function(lon,lat) {cos(d2r(lat))*sin(d2r(lon-lon0))}
  c2y <- function(lon,lat) {(cos(d2r(lat0))*sin(d2r(lat)))-(sin(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}
  #cut is cosin of c, when negative is behind projections, needs to be cut
  cut <- function(lon,lat) {(sin(d2r(lat0))*sin(d2r(lat)))+(cos(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}

  #functions spherical (lon=x, lat=y) to Cartesian
  s2cx <- function(x,y) {cos(d2r(x))*cos(d2r(y))}
  s2cy <- function(x,y) {sin(d2r(x))*cos(d2r(y))}
  s2cz <- function(y) {sin(d2r(y))}

  colnames(VGP) <- c("lon","lat")
  vgpsN <- common_DI(VGP,down = ifelse(mean(VGP$lat)<0,FALSE,TRUE))
  PPole <- fisher(vgpsN)
  Plon <- PPole[1,1]
  Plat <- PPole[1,2]
  A <- PPole[1,3]

  #fix point of view
  if(auto_cent==FALSE){
    #center of proj is Lon0 & Lat0
    lon0 <- long
    lat0 <- lat
  }else{
    lon0 <- Plon
    lat0 <- Plat
  }
  #save declination and inc and calculate new system for rotation
  newSlon <- ifelse((Plon+180)>360,Plon-180,Plon+180)
  newSlat <- 90-Plat
  newSlonr <- d2r(newSlon)
  newSlatr <- d2r(newSlat)
  a95 <- A
  circle <- as.data.frame(matrix(ncol=2,nrow=0))
  #loop that create a95 and rotate it around new coordinate (dec, inc)
  for (i in seq(0,360,2)){
    circleP <- as.data.frame(matrix(ncol=2,nrow=1))
    x <- s2cx(i,(90-a95))
    y <- s2cy(i,(90-a95))
    z <- s2cz(90-a95)
    vec <- as.matrix(c(x,y,z))
    R_elements <- c(cos(newSlatr)*cos(newSlonr), -sin(newSlonr), -sin(newSlatr)*cos(newSlonr),
                    cos(newSlatr)*sin(newSlonr), cos(newSlonr), -sin(newSlatr)*sin(newSlonr),
                    sin(newSlatr), 0, cos(newSlatr))
    R <- matrix(R_elements,nrow=3, byrow=TRUE)
    newvec <- R%*%vec
    newlon <- r2d(atan2(newvec[2,1],newvec[1,1]))
    newlon <- ifelse(newlon<0,newlon+360,newlon)
    #absolute value avoid point outside the graph
    newlat <- r2d(asin(newvec[3,1]))
    circleP[1,1:2] <- c(newlon,newlat)
    circle <- rbind(circle,circleP)
  }
  colnames(circle) <- c("lon","lat")
  circle$x <- c2x(circle$lon,circle$lat)
  circle$y <- c2y(circle$lon,circle$lat)
  circle$cut <- cut(circle$lon,circle$lat)
  #standalone graph or on existing graph
  if (on_plot==FALSE) {sph_ortho(lat = lat0,long = lon0,grid = grid,coast=coast)}

  if(symbol=="c") {sym <- 21}
  else if(symbol=="s") {sym <- 22}
  else if(symbol=="d") {sym <- 23}
  else if(symbol=="t") {sym <- 24}
  else{stop("Please select valid symbol. Check help for info.",call. = F)}

  #plot pole
  Px <- c2x(Plon,Plat)
  Py <- c2y(Plon,Plat)
  Pcut <- cut(Plon,Plat)

  #plot symbol open if behind
  if(Pcut>0) {points(Px,Py,pch=sym,col="black",bg=color)
  }else{points(Px,Py,pch=sym,col="black",bg="white")}

  polygon(circle$x,circle$y, col=col_A, lwd=1.2,lty= ifelse(Pcut>0,1,3))

  if(on_plot==FALSE){
    text <- paste("N: ",round(PPole[1,4],digits=2),"
Long: ", round(PPole[1,1],digits=2),"
Lat: ", round(PPole[1,2],digits=2),"
A95: ", round(PPole[1,3],digits=2))
    text(x=0.75, y=-0.85,pos=4,text, cex= 0.85)
  }

  #plot APWP if requested during process
  pAPWP <- readline("Plot APWP? (y or n): ")
  if(pAPWP=="y"){
    plot_APWP(APWP=APWP, lon0 = lon0,lat0 = lat0,grid = grid,on_plot = T,S_APWP = S_APWP)
  }
  if(save==TRUE){
    save_pdf(name = paste(name,".pdf"),width = 8,height = 8)
    cat("Figure saved as",name, ".pdf")
  }
}

#bootstrap of VGPs
VGP_boot <- function(VGP,nb=1000,lat=90,long=0,grid=30,auto_cent=TRUE,on_plot=FALSE,coast=FALSE,symbol="c",color= "blue",hist=TRUE,text=TRUE,save=FALSE, name="VGP_boot",APWP="V23", S_APWP=FALSE){

  #warning for on-plot, to avoid wrong coordinates
  if(on_plot==TRUE && auto_cent==TRUE) {
    stop("Please SPECIFY center coordinates when on_plot==TRUE",call. = F)
  }

  #functions converting degree and radians
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}

  #functions converting long & lat to xy
  c2x <- function(lon,lat) {cos(d2r(lat))*sin(d2r(lon-lon0))}
  c2y <- function(lon,lat) {(cos(d2r(lat0))*sin(d2r(lat)))-(sin(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}
  #cut is cosin of c, when negative is behind projections, needs to be cut
  cut <- function(lon,lat) {(sin(d2r(lat0))*sin(d2r(lat)))+(cos(d2r(lat0))*cos(d2r(lat))*cos(d2r(lon-lon0)))}
  #manipulate data
  colnames(VGP) <- c("lon","lat")
  vgpsN <- common_DI(VGP,down = ifelse(mean(VGP$lat)<0,FALSE,TRUE))
  PPole <- fisher(vgpsN)
  Plon <- PPole[1,1]
  Plat <- PPole[1,2]
  #fix point of view
  if(auto_cent==FALSE){
    #center of proj is Lon0 & Lat0
    lon0 <- long
    lat0 <- lat
  }else{
    lon0 <- Plon
    lat0 <- Plat
  }
  if(on_plot==FALSE){sph_ortho(lat=lat0,long=lon0,grid=grid, coast=coast)}
  cat("Bootstrapping.
Simulation ends when", nb, " pseudosamples are saved.

")
  bootlonlat <- as.data.frame(matrix(ncol = 2,nrow = 0))
  n <- 0
  repeat{
    n <- n+1
    VGPb <- boots_DI(VGP)
    VGPb_av <- fisher(VGPb)
    blon <- VGPb_av[1,1]
    blat <- VGPb_av[1,2]
    blonlat <- as.data.frame(t(c(blon,blat)))
    bootlonlat <- rbind(bootlonlat,blonlat)
    x <- c2x(blon,blat)
    y <- c2y(blon,blat)
    cutt <- cut(blon,blat)
    points(x,y,pch=ifelse(cutt>0,16,1),col=rgb(1,0,0,0.15))
    if(((n%%50)==0)==TRUE){

      cat(paste(n,"simulations out of",nb,"done
"))

      if(n==nb) break
    }
  }
  colnames(bootlonlat) <- c("vgp_lon","vgp_lat")
  bootlonlat$Plon <- rep(Plon)
  bootlonlat$Plat <- rep(Plat)
  bootlonlat$delta <- abs(bootlonlat$vgp_lon-bootlonlat$Plon)
  bootlonlat$diff <- r2d(acos((sin(d2r(bootlonlat$vgp_lat))*sin(d2r(bootlonlat$Plat)))+
                                (cos(d2r(bootlonlat$vgp_lat))*cos(d2r(bootlonlat$Plat))*cos(d2r(bootlonlat$delta)))))
  ang_dis <- as.data.frame(bootlonlat$diff)
  ang_dis <- (ang_dis[order(ang_dis[,1]),])
  conf <- 0.95
  Uconf <- round(nb*conf,digits=0)
  angular_conf <- ang_dis[Uconf]

  if(symbol=="c") {sym <- 21}
  else if(symbol=="s") {sym <- 22}
  else if(symbol=="d") {sym <- 23}
  else if(symbol=="t") {sym <- 24}
  else{stop("Please select valid symbol. Check help for info.",call. = F)}

  #plot pole
  Px <- c2x(Plon,Plat)
  Py <- c2y(Plon,Plat)
  Pcut <- cut(Plon,Plat)

  #plot symbol open if behind
  if(Pcut>0) {points(Px,Py,pch=sym,col="black",bg=color)
  }else{points(Px,Py,pch=sym,col="black",bg=NA)}

  #plot angular error estimation
  if(on_plot==FALSE && hist==TRUE){
    par(fig=c(0,0.5,0,0.5), new=TRUE)
    plot(NA, xlim=c(0,1), ylim=c(0,1),
         xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
    rect(xleft = 0,ybottom = 0,xright = 1,ytop = 1,col = "white",border = "black",cex=0.5)

    par(fig=c(0,0.5,0,0.5), new=TRUE)
    hist(x = bootlonlat$diff,xlab=NA,main="",ylab=NA,
         xlim=c(0,10),breaks=40,cex.axis=0.9,
         col="red",border ="red")
    abline(v=angular_conf,lwd=1,lty=2)
    title(xlab = "Angular distance (°)", line=2, cex=0.2)
    title(ylab = "Frequency", line=2, cex=0.2)
  }
  par(fig=c(0,1,0,1), new=TRUE)
  #plot text with results
  results <- as.data.frame(Plon)
  results$Plat <- Plat
  results$N <- length(VGP$lon)
  results$ang_conf <- angular_conf
  results <- round(results, digits=2)

  if (on_plot==FALSE && text==TRUE){
    text <- paste("N: ",results$N,"
Long: ", results$Plon,"
Lat: ", results$Plat,"
B95: ", results$ang_conf)
    plot(NA, xlim=c(0,1), ylim=c(0,1),
         xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)

    text(x=0.76, y=0,pos=4,text, cex= 0.85)
  }
  #plot APWP if requested during process
  #plot APWP if requested during process
  pAPWP <- readline("Plot APWP? (y or n): ")
  if(pAPWP=="y"){
    plot_APWP(APWP=APWP, lon0 = lon0,lat0 = lat0,grid = grid,on_plot = T,S_APWP = S_APWP)
  }
  #replot pole data on apwp
  par(fig=c(0,1,0,1), new=TRUE)
  plot(NA, xlim=c(-1,1), ylim=c(-1,1), asp=1,
       xlab="", xaxt="n",ylab="", yaxt="n", axes=FALSE)
  if(Pcut>0) {points(Px,Py,pch=sym,col="black",bg=color)
  }else{points(Px,Py,pch=sym,col="black",bg=NA)}

  if (save==TRUE){
    write.csv(results, file = paste(name,".csv"), row.names = FALSE)
    save_pdf(name = paste(name,".pdf"),width = 8,height = 8)
    cat("Figure saved as",name, ".pdf
","Result file saved as", name,".csv")
  }
}

#calculate virtual geomagnetic pole(s)
VGP_DI <- function(DI,in_file=FALSE,lat,long,export=TRUE,type="VGPsN",name="VGPs",Prnt=TRUE){
  #conversion functions
  d2r <- function(x) {x*(pi/180)}
  r2d <- function(x) {x*(180/pi)}
  #start data table
  data <- DI[,1:2]
  data <- na.omit(data)
  #add lat and long columns if not present in file
  if(in_file==FALSE){
    data$slat <- rep(lat)
    data$slong <- rep(long)
  }
  #fix col names
  colnames(data) <- c("dec","inc","slat","slong")
  #populate data table with used data
  data$dec_r <- d2r(data$dec)
  data$inc_r <- d2r(data$inc)
  data$slat_r <- d2r(data$slat)
  data$slong_r <- d2r(data$slong)
  #calculate pole colatitude
  data$p_colat_r <- atan(2/tan(data$inc_r))
  data$p_colat_d <- r2d(data$p_colat_r)
  #calculate pole latitude
  data$PLat_r <- asin((sin(data$slat_r)*cos(data$p_colat_r))+
                        (cos(data$slat_r)*sin(data$p_colat_r)*cos(data$dec_r)))
  data$Plat_d<- r2d(data$PLat)
  #Longitudinal difference between site and pole
  data$LDist_r <- asin((sin(data$p_colat_r)*sin(data$dec_r))/cos(data$PLat_r))
  data$LDist_d <- r2d(data$LDist_r)
  #calculate longitude
  data$PLong_d <- ifelse(cos(data$p_colat_r)<(sin(data$slat_r)*sin(data$PLat_r)),
                         data$slong+180-data$LDist_d,
                         data$slong+data$LDist_d)
  #adjust vgp for normal and reversed
  data$Pole_longitude <- ifelse(data$p_colat_d<0,ifelse((data$PLong_d+180)>360,
                                                        data$PLong_d-180,data$PLong_d+180),
                                data$PLong_d)
  data$Pole_latitude <- ifelse(data$p_colat_d<0,-data$Plat_d,data$Plat_d)
  #isolate VGPs with reversals
  VGPs <- data[,16:17]

  #isolate VGPs all normal
  #VGPsN <- data[,c(15,12)]
  VGPsN <- common_DI(VGPs)     #fixed 29.12.2023
  colnames(VGPsN) <- c("Plong_N","Plat_N")

  #calculate average VGP
  PmagPole <- fisher(VGPsN)

  #rename columns of PmagPole
  colnames(PmagPole)[1:2] <- c("long","lat")

  #set coordinates for rotation of pole
  ifelse(PmagPole[1,2]>0,
         RotAZ <- (PmagPole[1,1]+180)%%360,
         RotAZ <-  PmagPole[1,1])
  RotPL <- 90-abs(PmagPole[1,2])

  #rotate VGPs to North pole
  VGPsR <- bed_DI(DI = VGPs,in_file = F,bed_az = RotAZ,
                  bed_plunge = RotPL)
  #if(data[1,3]<0 && PmagPole$lat<0) VGPsR <- flip_DI(VGPsR)
  colnames(VGPsR) <- c("Plong_R","Plat_R")


  if(export==TRUE){
    write.csv(round(PmagPole,digits=2),file=paste(name,"_average_pole.csv"),row.names = F)
    write.csv(round(VGPsN,digits = 2),file=paste(name,"_single_mode.csv"),row.names = F)
    write.csv(round(VGPs, digits = 2),file=paste(name,"_bimodal.csv"),row.names = F)
    write.csv(round(VGPsR,digits = 2),file=paste(name,"_rotated.csv"),row.names = F)
    cat("File exported as csv within the working directory

")
  }
  if(Prnt==TRUE){
    cat("Paleomagnetic pole:

")
    print(round(PmagPole,digits=2),row.names=FALSE)
  }

  if(type=="VGPs"){return(VGPs)}
  if(type=="VGPsN"){return(VGPsN)}
  if(type=="VGPsR"){return(VGPsR)}
}

#plot decl, inc, VGP lat, polarity in stratigraphic depth, and directions and VGP plots if requested
#still under development!!
magstrat_DI <- function(DIP,lat=0,long=0,offset=0,col="red",name="polarity_plot",save=FALSE,plot_ext=TRUE,POLE=TRUE, E.A.=TRUE,cex.main=1,cex.lab=1,cex.axis=1,lwd.grid=1,h_grid=10, Shiny=FALSE){
  library(plyr, warn.conflicts = F)
  library(PmagDiR)
  dat <- na.omit(DIP)
  colnames(dat) <- c("dec","inc","posit")
  #calculate VGPs rotated
  dat[,4:5] <- VGP_DI(dat[,1:2],in_file = F,lat = lat,long = long,export = F,type = "VGPsR",Prnt = F)
  #associated polarity to VGPs, where normal=1, reversed=0
  dat$pol <- ifelse(dat$Plat_R>0,1,0)
  #create reversals empty data frame
  normals <- data.frame(matrix(ncol = 2,nrow = 0))
  colnames(normals) <- c("bottom","top")
  #populate top and bottom of normals table
  for(i in 2:nrow(dat)){
    if((dat[i,6]+dat[i-1,6])==1){
      pos <- (dat[i,3]+dat[i-1,3])/2
      newline <- data.frame(matrix(ncol = 2,nrow = 0))
      ifelse(dat[i,6]==1, newline[1,1] <- pos, newline[1,2] <- pos)
      colnames(newline) <- c("bottom","top")
      normals <- rbind(normals,newline)
    }
  }
  #fill first or last box of normals when empty
  if(is.na(normals[1,1])==TRUE) normals[1,1] <- min(dat$posit)
  if(is.na(normals[nrow(normals),2])==TRUE) normals[nrow(normals),2] <- max(dat$posit)
  #reduce table to lines with top and bottom
  if(nrow(normals)>1){
    for(l in 2:nrow(normals)){
      if(is.na(normals[l-1,2])==TRUE) {normals[l-1,2] <- normals[l,2]}
    }
  }
  #eliminate duplicates
  normals <- na.omit(normals)
  ymin <- plyr::round_any(min(dat$posit), 0.5, f= floor)
  ymax <- plyr::round_any(max(dat$posit), 0.5, f=ceiling)

  #fix declination if offset is not zero
  if(offset!=0){
    dat$dec <- (dat$dec+(abs(offset)))%%360
    dat$dec <- dat$dec-offset
  }

  ############## PLOT ##############
  #screen splitter matrix
  if(plot_ext==TRUE) {dev.new(width = 10,height = 7,noRStudioGD = T)}
  screen <- matrix(c(1,1,2,2,3,3,4),ncol=7,byrow = T)
  layout(screen)
  inclim <- round(max(abs(dat$inc)), digits = 0)
  #declination
  plot(dat$dec,dat$posit, type="o",
       pch=21,bg=col,ylab="Position (m)",
       xlim=c(-offset,(360-offset)),
       xaxp= c(-offset,(360-offset),4),
       ylim=c(ymin,ymax),
       xlab=NA,
       main="Declination (°)",
       cex.main=cex.main,
       cex.lab=cex.lab,
       cex.axis=cex.axis,
       panel.first= abline(v=c(seq(0-offset,360-offset,90)),
                           h=c(seq(round(min(dat$posit), digits = -1),
                                   round(max(dat$posit), digits = 0),h_grid)),
                           col="gray", lty="dotted",lwd=lwd.grid))
  plot(dat$inc,dat$posit,type="o",
       pch=21,bg=col,ylab=NA,xaxp= c(-90,90,6),
       xlim=c(-90,90),
       ylim=c(ymin,ymax),
       xlab=NA,
       main="Inclination (°)",
       cex.main=cex.main,
       cex.axis=cex.axis,
       panel.first= abline(v=c(seq(-90,90,30)),
                           h=c(seq(round(min(dat$posit), digits = -1),
                                   round(max(dat$posit), digits = 0),h_grid)),
                           col="gray", lty="dotted",lwd=lwd.grid))
  plot(dat$Plat_R,dat$posit,type="o",
       pch=21,bg=col,ylab=NA,xaxp= c(-90,90,4),
       xlim=c(-90,90),
       ylim=c(ymin,ymax),
       xlab=NA,
       main="VGP Lat. (°)",
       cex.main=cex.main,
       cex.axis=cex.axis,
       panel.first= abline(v=c(seq(-90,90,45)),
                           h=c(seq(round(min(dat$posit), digits = -1),
                                   round(max(dat$posit), digits= 0),h_grid)),
                           col="gray", lty="dotted",lwd=lwd.grid))

  #create frame for polarity
  plot(NA,
       xlim=c(0,1), xaxt="n",
       type="n", ylab=NA,
       xlab="",ylim=c(ymin,ymax),
       main="Polarity",
       cex.main=cex.main,
       cex.axis=cex.axis)
  rect(xleft=0,
       ybottom=normals$bottom,
       xright=1,
       ytop=normals$top,
       col=ifelse(nrow(normals==1) && any(dat[,6]==1), "black","white"),
       border=NA)
  if(save==TRUE){
    save_pdf(name =paste(name,".pdf"),width = 10,height = 8)
  }
  if(POLE==TRUE){
    dev.new(width = 7,height = 7,noRStudioGD = T)
    VGPsN <- VGP_DI(dat[,1:2],in_file = F,lat = lat,long = long,export = T,Prnt = T)
    plot_VGP(VGPsN, coast = T, A95 = T,save = save)
  }
  if(E.A.==TRUE){
    dev.new(width = 7,height = 7,noRStudioGD = T)
    plot_DI(dat[,1:2])
    fisher_plot(dat[,1:2],save = save, text=T)
  }
  if(Shiny==TRUE){
    Table_of_normal_polarity_zones <- round(normals,digits=2)
    return(Table_of_normal_polarity_zones)
  }
}




