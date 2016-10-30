
####### get threshold values (range over response / driver)
################from Scott

gamoramma <- function(driver=NA,response=NA,the.time=NA,plot.diag=FALSE,nboot=1000,deriv.choose=1)
{
  require(mgcv)
  #  require(time)
  
  ## Variables
  mytime <- the.time
  mydriver <- driver
  myresponse <- response
  
  time.order<-order(mytime)
  mydriver<-mydriver[time.order]
  myresponse<-myresponse[time.order]
  #n.pts<-length(myresponse)
  
  # Reorder according to driver
  driver.order<-order(mydriver)
  mydriver<-mydriver[driver.order]
  myresponse<-myresponse[driver.order]
  
  #####################
  # Fit the GAM model #
  #####################
  sp.len<-200
  #gam1<-gam(myresponse~s(mydriver, bs="cr"), se=T)
  
  gam1<-gam(myresponse~s(mydriver, bs="ts"), method="GCV.Cp", se=T)
  gam2<-gam(myresponse~mydriver, method="GCV.Cp", se=T)
  
  x1 <- as.data.frame(cbind(summary(gam1)$dev.expl,summary(gam1)$edf, summary(gam1)$sp.criterion, summary(gam1)$s.pv))   
  x1 <- rbind(x1,as.data.frame(cbind(summary(gam2)$dev.expl,NA, summary(gam2)$sp.criterion, summary(gam2)$p.pv[2])))
  names(x1)<-c("dev.expl","edf", "GCV", "p_value")
  if (x1[1,4]<0.05 && x1[1,3]<x1[2,3] && x1[,1,2]>1e-04) {
    
    
    pred <- predict.gam(gam1, se.fit=T)
    fds <- cbind(pred$fit,mydriver,pred$se.fit)
    fds2 <- order(fds[,2])
    response <- spline(pred$fit[fds2]~mydriver[fds2],n=sp.len)
    gam.data <- as.data.frame(cbind(response$y,response$x))
    colnames(gam.data)<-c("response","driver")
    
    set.seed(1122) # to get reproducibable results from bootstrap
    nb <- nboot #bootstrap_size  
    sp.len<-200 # Length of the spline for bootstrapping
    thresh <- matrix(nrow=sp.len, ncol=nb) ## create a matrix to be filled with bootstrapped splines (200 x 1000)
    #prev <- progressBar() ## just to display how long it remains for the bootstrap
    for (i in 1:nb) {
      #prev <- progressBar(i/nb, prev)
      bootsample <- sample(1:length(myresponse), replace=T)
      myresponsei <- myresponse[bootsample]
      mydriveri <- mydriver[bootsample]
      #######Computing threshold on new sample #######
      gam.object <- gam(myresponsei~s(mydriveri, bs="cr"), na.action='na.omit')
      prediction <- predict.gam(gam.object, se.fit=T)
      fds <- cbind(prediction$fit,mydriveri,prediction$se.fit,0,0)
      fds2 <- order(fds[,2])
      titu <- spline(prediction$fit[fds2]~mydriveri[fds2],n=sp.len)
      data <- cbind(titu$y,titu$x)
      thresh[,i]<- titu$y
    }
    
    # Matrix plot of bootstrap GAM replicates for visualization
    driv<-titu$x
    
    #############################
    # 1st Derivative Estimation #
    #############################
    dif1.line<-diff(gam.data$response, difference=1) # Actual 1st deriv estimate from original smoother
    
    dif1<- 1 ## Adjust which difference to take
    deriv.matrix.1<-matrix(nrow=sp.len-dif1, ncol=nb) ## create a matrix to for the 1st deriv estimates
    #prev <- progressBar() ## just to display how long it remains for the bootstrap
    for (i in 1:nb) {
      #prev <- progressBar(i/nb, prev)
      derivi<-thresh[,i]
      deriv.response<-diff(derivi, difference=dif1)
      driver.len<-length(driv)-dif1
      deriv.object <- cbind(deriv.response,driver.len)
      deriv.matrix.1[,i]<-deriv.object[,1]
    }
    
    
    #############################
    # 2nd Derivative Estimation #
    #############################
    dif2.line<-diff(gam.data$response, difference=6) # Actual 2nd deriv estimate from original smoother, 6th difference
    # may be better estimate for 2nd derivative... must investigate
    
    dif2<- 6 ## Adjust which difference to take
    deriv.matrix.2<-matrix(nrow=sp.len-dif2, ncol=nb) ## create a matrix to for the 2nd deriv estimates
    #prev <- progressBar() ## just to display how long it remains for the bootstrap
    for (i in 1:nb) {
      #prev <- progressBar(i/nb, prev)
      derivi<-thresh[,i]
      deriv.response<-diff(derivi, difference=dif2)
      driver.len<-length(driv)-dif2
      deriv.object <- cbind(deriv.response,driver.len)
      deriv.matrix.2[,i]<-deriv.object[,1]
    }
    
    # CI of GAM bootstrap
    ci<-matrix(nrow= 2, ncol= sp.len) ## create a matrix to be filled with bootstrapped CI
    rownames(ci)<-c("lower","upper")
    #prev <- progressBar() ## just to display how long it remains for the bootstrap
    for (i in 1:sp.len) {
      #prev <- progressBar(i/sp.len, prev)
      IC<-quantile(thresh[i,], c(0.025, 0.975))
      ci[,i]<-rbind(IC[1], IC[2])
    }
    
    # CI of 1st derivative 
    dif1.len<-sp.len-dif1
    ci1<-matrix(nrow= 2, ncol= dif1.len) ## create a matrix to be filled with bootstrapped CI
    rownames(ci1)<-c("lower","upper")
    #prev <- progressBar() ## just to display how long it remains for the bootstrap
    for (i in 1:dif1.len) {
      #prev <- progressBar(i/dif1.len, prev)
      IC<-quantile(deriv.matrix.1[i,], c(0.025, 0.975))
      ci1[,i]<-rbind(IC[1], IC[2])
    }
    
    # CI of 2nd derivative 
    dif2.len<-sp.len-dif2
    ci2<-matrix(nrow= 2, ncol= dif2.len) ## create a matrix to be filled with bootstrapped CI
    rownames(ci2)<-c("lower","upper")
    #prev <- progressBar() ## just to display how long it remains for the bootstrap
    for (i in 1:dif2.len) {
      #prev <- progressBar(i/dif2.len, prev)
      IC<-quantile(deriv.matrix.2[i,], c(0.025, 0.975))
      ci2[,i]<-rbind(IC[1], IC[2])
    }
    
    
    # Number of years and sequence of years
    Nmydriver <- length(gam.data$driver)
    mydrivervec <-seq(min(gam.data$driver),max(gam.data$driver))
    driver1<-gam.data$driver[-1]
    driver2<-gam.data$driver[-c(1:6)]
    
    # CI of response
    lower <- min(ci[1,  ])
    upper <- max(ci[2,  ])
    lower1 <- min(ci1[1,  ])
    upper1 <- max(ci1[2,  ])
    lower2 <- min(ci2[1,  ])
    upper2 <- max(ci2[2,  ])
    
    # Significant 1st derivative indicators
    changepts.lower1 <- seq(1, dif1.len)[ci1["lower",  ] > 0]
    changepts.upper1 <- seq(1, dif1.len)[ci1["upper",  ] < 0]
    
    # Significant second derivative indicators
    changepts.lower2 <- seq(1, dif2.len)[ci2["lower",  ] > 0]
    changepts.upper2 <- seq(1, dif2.len)[ci2["upper",  ] < 0]
    
    # Matrix of 0 values to make polygon for derivative graphs, update with length of significant indicators
    ze<-as.matrix(rep(0, length(changepts.upper1))) 
    
    results <- NULL
    if(deriv.choose==1)
    {
      results$driverchange <- c(gam.data$driver[changepts.lower1],gam.data$driver[changepts.upper1])
      results$responsechange <- c(gam.data$response[changepts.lower1],gam.data$response[changepts.upper1])
    }
    if(deriv.choose==2)
    {
      results$driverchange <- c(gam.data$driver[changepts.lower2],gam.data$driver[changepts.upper2])
      results$responsechange <- c(gam.data$response[changepts.lower2],gam.data$response[changepts.upper2])
    }
    return(results)
    
  }
  if(x1[1,4]>=0.05 || x1[1,3]>=x1[2,3] || x1[,1,2]<=1e-04)
  {
    results <- NULL
  }
  
  #end function gamoramma
}


