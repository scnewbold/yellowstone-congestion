#===============================================================================
# yellowstone-congestion.R
#===============================================================================

#-------------------------------------------------------------------------------
# PRELIMINARIES:
#-------------------------------------------------------------------------------
if(TRUE){ 
  
  # Clear environment to start fresh
  rm(list=ls())
  
  # Grab script name for file handling:
  script.name <- basename(rstudioapi::getSourceEditorContext()$path) 
  script.name <- gsub(".R","",script.name)
  
  # Packages
  list.of.packages <- c('tikzDevice','tinytex')
  
  new.packages <- list.of.packages[!(list.of.packages %in% 
                                       installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  lapply(list.of.packages,function(x){library(x,character.only=TRUE)})
  
  # Clear console:
  cat('\014');
  
  # Clear all plots:
  try(dev.off(dev.list()["RStudioGD"]),silent=TRUE)
  try(dev.off(),silent=TRUE)
  
  # Streamlined functions for formatting:
  s2  <- function(x){return(sprintf('%-.2f',x))}
  s3  <- function(x){return(sprintf('%-.3f',x))}
  s4  <- function(x){return(sprintf('%-.4f',x))}
  s5  <- function(x){return(sprintf('%-.5f',x))}
  s6  <- function(x){return(sprintf('%-.6f',x))}
  s10 <- function(x){return(sprintf('%-.10f',x))}
  
  # Define paths for file handling:
  this.dir <- dirname(parent.frame(2)$ofile) # source file dir
  setwd(this.dir)                            # set wd to source file dir
  code.path <- getwd()                       # define code path
  # setwd('..')                                # set wd up one level
  # setwd(paste(getwd(),'/output',sep=''))     # set wd to output folder
  # output.path <- getwd()                     # define output path
  # setwd('..')                                # set wd up one level
  # setwd(paste(getwd(),'/input',sep=''))      # set wd to input folder
  # input.path <- getwd()                      # define input path
  # setwd('..')                                # set wd up one level
  # setwd(code.path)                           # set wd to code folder
  
  # Single path all local:
  input.path  <- code.path
  output.path <- code.path
  
  # Create output file in working directory:
  date.time     <- gsub(" ","_",Sys.time())
  date.time     <- gsub("-","_",date.time)
  date.time     <- gsub(":","_",date.time)
  # out.file.name <- paste(output.path,'/',script.name,'-',date.time,'.out',sep='')
  out.file.name <- paste(output.path,'/',script.name,'.out',sep='')
  outfile       <- file.create(out.file.name) 
  
  # WRITE SOURCE FILE TO OUTPUT FILE:
  if(FALSE){
    source.file.name <- paste(code.path,'/',script.name,'.R',sep='')
    
    # Read lines of source file:
    Rscript <- readLines(source.file.name)
    
    # Write lines of source file to output file:
    for(i in 1:length(Rscript)){cat('\n',Rscript[i],file=out.file.name,append=TRUE)}
    cat('\n\n|---------------------------------------------------------------------------|',file=out.file.name,append=TRUE)
    cat('\n| R SCRIPT ABOVE                                                            |',file=out.file.name,append=TRUE)
    cat('\n|---------------------------------------------------------------------------|',file=out.file.name,append=TRUE)
    cat('\n| R OUTPUT BELOW                                                            |',file=out.file.name,append=TRUE)
    cat('\n|---------------------------------------------------------------------------|\n',file=out.file.name,append=TRUE)
    
  }
}

#===============================================================================
# FUNCTIONS:
#===============================================================================
{
  pause   <- function(sec){
    if(missing(sec)){
      readline(prompt='Paused. Press [enter] to continue.\n')
    }else{
      Sys.sleep(sec)
    }
  }
  
  tic     <- function(){
    start.time<<-proc.time()[3]
  }
  
  toc     <- function(){
    elapsed<-(proc.time()[3]-start.time)
    return(sprintf("elapsed time = %0.4f sec",elapsed))
  }
  
  data.fn <- function(){
    
    # IMPORT DATA:
    
    # visits.csv
    visits <- read.csv(paste(input.path,"/visits.csv",sep=""))
    
    Ytm <- c()
    for(i in 1:length(visits[,1])){
      if(sum(years==visits[i,1])>0){
        Ytm <- c(Ytm,visits$rec_visits[i])
      }
    }
    
    
    # states.csv
    states <- read.csv(paste(input.path,"/states.csv",sep=""))
    yj     <- states$visits_percap_2018
    N      <- states$pop_2018
    dist   <- states$dist
    med.inc <- states$median_income
    
    # untransformed travel cost (excluding entry fee):
    d <- 2*(dist * cost.per.mi / adults.per.vehicle + (dist/avg.speed) * 
              fraction.of.wage * med.inc / adults.per.hh / working.hrs) 
    
    return(list(N=N,d=d,Ytm=Ytm,yj=yj))
    
  }
  
  Y.fn    <- function(Y.inputs){
    
    N        <- Y.inputs$N
    beta     <- Y.inputs$beta
    phi.t    <- Y.inputs$phi.t
    alpha.m  <- Y.inputs$alpha.m
    theta    <- Y.inputs$theta
    lambda   <- Y.inputs$lambda
    sigma    <- Y.inputs$sigma
    c        <- Y.inputs$c
    
    J        <- length(N)
    
    adjuster <- .1
    
    epsilon  <- rnorm(J,-.5*sigma^2,sigma)
    
    Y <- matrix(0,J,1)
    done <- 0
    while(done==0){
      YY <- adjuster*N*exp(beta+phi.t+alpha.m+theta*(sum(Y)/10^6)+lambda*c+epsilon) + (1-adjuster)*Y
      if(max(abs(YY-Y)/Y)<0.00001){done <- 1}else{Y <- YY}
    }
    
    Ytmj <- Y
    Ytm  <- sum(Ytmj)
    
    return(list(Ytmj=Ytmj,Ytm=Ytm))
    
  }
  
  sim.fn  <- function(sim.inputs){
    
    N      <- sim.inputs$N
    beta   <- sim.inputs$beta
    phi    <- sim.inputs$phi
    alpha  <- sim.inputs$alpha
    theta  <- sim.inputs$theta
    lambda <- sim.inputs$lambda
    c      <- sim.inputs$c
    J      <- length(c)
    
    Y.inputs <- sim.inputs
    
    Y <- matrix(0,T*12,J)
    tm <- 0
    for(t in 1:T){
      for(m in 1:12){tm <- tm + 1
      Y.inputs$phi.t   <- phi[t]
      Y.inputs$alpha.m <- alpha[m]
      if(open[tm]==1){Y.inputs$theta <- theta}
      if(open[tm]==0){Y.inputs$theta <- 0    }
      Y[tm,] <- Y.fn(Y.inputs)$Ytmj 
      }
    }
    
    Yj <- Y
    Y  <- rowSums(Yj)
    
    return(list(Y=Y,Yj=Yj))
    
  }
  
  R2.fn   <- function(theta){
    depv  <- log(Y)-theta*(Y/10^6)*open
    outs  <- summary(lm(depv ~ D.yrs + D.mth + prcp + tmin + tmax))
    ehat  <- outs$residuals
    R2    <- outs$adj.r.squared
    return(R2)
  }
}

#===============================================================================
# MAIN PROGRAM:
#===============================================================================

#-------------------------------------------------------------------------------
# SIMULATION TEST [full year]:
#-------------------------------------------------------------------------------
if(TRUE){
  
  set.seed(123)
  
  Z <- 1000 # number of Monte Carlo iterations
  
  # Parameters:
  {
    sigma  <- .45
    
    beta   <- -8.5
    #                J   F   M   A   M   J   J   A   S   O   N   D
    alpha  <- log(c(.01,.01,.02,.02,.04,.16,.34,.27,.09,.02,.01,.01))
    alpha  <- alpha - alpha[1]       # Normalizes Jan alpha to 0
    
    T      <- 40                     # 40 years of data
    phi    <- seq(0,log(2),length=T) # year fixed effects increase linearly
    phi    <- phi - phi[1]           # Normalizes first phi to 0
    
    J      <- 49                     # Continental US states + DC
    
    lambda <- -0.01                  # travel cost coefficient 
    theta  <- -1.0                   # congestion coefficient
    
  }
  
  # Year and month and may-oct dummy variables:
  {
    D.yrs <- matrix(0,T*12,T)
    tm <- 0
    for(t in 1:T){ for(m in 1:12){ tm <- tm + 1; D.yrs[tm,t] <- 1 } }
    D.yrs <- D.yrs[,-1]
    
    D.mth <- matrix(0,T*12,12)
    D.opn <- matrix(0,T*12,1)
    tm <- 0
    for(t in 1:T){ for(m in 1:12){ tm <- tm + 1; D.mth[tm,m] <- 1 ; if(m>=5 & m<=10){D.opn[tm]<-1}} }
    D.mth <- D.mth[,-1]
    
    open <- matrix(1,T*12,1)
    # open <- D.opn
    
  }
  
  tau.star      <- matrix(0,Z,1)
  Y.star        <- matrix(0,Z,1)
  lambda.hat    <- matrix(0,Z,1)
  theta.hat     <- matrix(0,Z,1)
  theta.hat.v2  <- matrix(0,Z,1)
  theta.hat.v3  <- matrix(0,Z,1)
  tau.hat       <- matrix(0,Z,1)
  Y.hat         <- matrix(0,Z,1)
  
  start.time <- proc.time()[3]
  for(z in 1:Z){
    
    cat('\014',sprintf('Working on %-.0f of %-.0f Monte Carlo iterations [full year]\n',z,Z))
    cat('[Estimated time to completion =',sprintf('%-.1f min.]\n',(proc.time()[3]-start.time)/z*(Z-z)/60))
    
    N0    <- floor(rnorm(J,2e8/J,2e7/J)) # population in each state in year 1
    gN    <- .01 # mean population growth rate [per year]
    N     <- matrix(0,J,T)
    N[,1] <- N0
    for(t in 2:T){ N[,t] <- N[,t-1]*(1+gN*exp(rnorm(J)*.1-.5*.1^2)) }
    
    c0 <- runif(J,min=100,max=500)    # travel cost from each state in year 1
    gc <- 0 # mean travel cost growth rate [per year]
    c  <- matrix(0,J,T)
    c[,1] <- c0
    for(t in 2:T){ c[,t] <- c[,t-1]*(1+gc*exp(rnorm(J)*.1-.5*.1^2)) }
    
    # Optimal fee computation (using true parameters and exog variables 
    # N[,T] & c[,T]):
    {
      Y.inputs         <- list()
      Y.inputs$N       <- N[,T]
      Y.inputs$theta   <- theta
      Y.inputs$lambda  <- lambda
      
      Y.0       <- 4*10^6 / 12 # optimal fee for a specific month with given baseline vistation
      gamma     <- log(Y.0/(exp(theta*(Y.0/10^6))*sum(N[,T]*exp(lambda*c[,T]))))
      tau.range <- seq(0,2000,1)
      Y.tau     <- matrix(0,length(tau.range),1)
      V.tau     <- matrix(0,length(tau.range),1)
      k <- 0
      for(tau in tau.range){k <- k + 1
      Y.inputs$beta    <- gamma
      Y.inputs$phi.t   <- 0
      Y.inputs$alpha.m <- 0
      Y.inputs$sigma   <- sigma
      Y.inputs$c       <- c[,T]+tau
      Y.tau[k]         <- Y.fn(Y.inputs)$Ytm
      V.tau[k]         <- Y.tau[k]/(-lambda) + tau*Y.tau[k]
      }
      Y.star[z]   <- Y.tau[which(V.tau==max(V.tau))]
      tau.star[z] <- tau.range[which(V.tau==max(V.tau))]
    }
    
    sim.inputs         <- list()
    sim.inputs$N       <- N[,T]
    sim.inputs$beta    <- beta
    sim.inputs$phi     <- phi
    sim.inputs$alpha   <- alpha
    sim.inputs$theta   <- theta
    sim.inputs$lambda  <- lambda
    sim.inputs$sigma   <- sigma
    sim.inputs$c       <- c[,T]
    
    outs <- sim.fn(sim.inputs)
    Y    <- outs$Y  # [T*12 x 1] total trips in 12 months of T years
    Yj   <- outs$Yj # [T*12 x J] trips from each j in 12 months of T years
    
    # Step 1: travel cost coefficient estimation
    {
      # Regress ln of total visits from each state in final year on travel costs:
      Yj   <- colSums(Yj[(dim(Yj)[1]-11):dim(Yj)[1],]) # Now a [J x 1] vector
      outs <- lm(log(Yj) ~ c[,T])
      lambda.hat[z] <- outs$coefficients[2]
    }
    
    # Step 2: congestion coefficient estimation
    if(TRUE){
      theta.range <- seq(0,theta*3,length=100)
      R2 <- matrix(0,length(theta.range),1)
      k  <- 0
      for(theta.k in theta.range){k <- k + 1
      depv  <- log(Y) - theta.k*(Y/10^6)*open
      outs  <- summary(lm(depv ~ D.yrs + D.mth))
      ehat  <- outs$residuals
      R2[k] <- outs$adj.r.squared
      }
      theta.hat[z] <- theta.range[which.max(R2)] 
      
      depv  <- log(Y)-theta.hat[z]*(Y/10^6)*open
      outs  <- summary(lm(depv ~ D.yrs + D.mth))
      ehat  <- outs$residuals
    }
    
    # Step 2b: congestion coefficient estimation (naive--endogenous Y on rhs)
    if(TRUE){
      depv  <- log(Y)
      YY    <- Y/10^6 * open
      outs  <- summary(lm(depv ~ D.yrs + D.mth + YY))
      theta.hat.v2[z] <- tail(outs$coeff,1)[1]
    }
    
    # Step 2c: congestion coefficient estimation (lagged Y on rhs)
    if(TRUE){
      depv  <- log(Y[13:length(Y)])
      YY    <- Y[1:(length(Y)-12)]/10^6 * open[1:(length(Y)-12)]
      outs  <- summary(lm(depv ~ D.yrs[13:length(Y)] + D.mth[13:length(Y)] + YY))
      theta.hat.v3[z] <- tail(outs$coeff,1)[1]
    }
    
    # Optimal fee estimation (using estimated parameters and exog variables
    # N[,T] & c[,T]) using theta.hat:
    if(TRUE){
      Y.inputs         <- list()
      Y.inputs$N       <- N[,T]
      Y.inputs$theta   <- theta.hat[z]
      Y.inputs$lambda  <- lambda.hat[z]
      
      Y.0   <- 4*10^6/12
      gamma <- log(Y.0/(exp(theta.hat[z]*(Y.0/10^6))*sum(N[,T]*exp(lambda.hat[z]*c[,T]))))
      tau.range <- seq(0,2000,1)
      Y.tau     <- matrix(0,length(tau.range),1)
      V.tau     <- matrix(0,length(tau.range),1)
      k <- 0
      for(tau in tau.range){k <- k + 1
      Y.inputs$beta    <- gamma
      Y.inputs$phi.t   <- 0
      Y.inputs$alpha.m <- 0
      Y.inputs$sigma   <- sigma
      Y.inputs$c       <- c[,T]+tau
      Y.tau[k]         <- Y.fn(Y.inputs)$Ytm
      V.tau[k]         <- Y.tau[k]/(-lambda.hat[z]) + tau*Y.tau[k]
      }
      Y.hat[z]   <- Y.tau[which(V.tau==max(V.tau))]
      tau.hat[z] <- tau.range[which(V.tau==max(V.tau))]
    }
    
    # Histograms:
    if(z>1){
      
      par(mfrow=c(2,3))
      
      hist(lambda.hat[1:z],xlim=c(2*lambda,0),main="travel cost coefficient",xlab="",ylab="")
      lines(x=c(lambda,lambda),y=c(0,Z),col='red')
      lines(x=c(mean(lambda.hat[1:z]),mean(lambda.hat[1:z])),y=c(0,500),col='blue',lty=1)
      lines(x=c(mean(lambda.hat[1:z])+1.96*sd(lambda.hat[1:z])/sqrt(z),
                mean(lambda.hat[1:z])+1.96*sd(lambda.hat[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
      lines(x=c(mean(lambda.hat[1:z])-1.96*sd(lambda.hat[1:z])/sqrt(z),
                mean(lambda.hat[1:z])-1.96*sd(lambda.hat[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
      
      hist(theta.hat[1:z],30,xlim=c(3*theta,-3*theta),main="congestion coefficient",xlab="",ylab="")
      lines(x=c(theta,theta),y=c(0,Z),col='red')
      lines(x=c(mean(theta.hat[1:z]),mean(theta.hat[1:z])),y=c(0,Z),col='blue',lty=1)
      lines(x=c(mean(theta.hat[1:z])+1.96*sd(theta.hat[1:z])/sqrt(z),
                mean(theta.hat[1:z])+1.96*sd(theta.hat[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
      lines(x=c(mean(theta.hat[1:z])-1.96*sd(theta.hat[1:z])/sqrt(z),
                mean(theta.hat[1:z])-1.96*sd(theta.hat[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
      
      hist(tau.hat[1:z],30,xlim=c(-mean(tau.star[1:z]),3*mean(tau.star[1:z])),main="optimal entry fee",xlab="",ylab="")
      lines(x=c(mean(tau.star[1:z]),mean(tau.star[1:z])),y=c(0,Z),col='red')
      lines(x=c(mean(tau.hat[1:z]),mean(tau.hat[1:z])),y=c(0,500),col='blue',lty=1)
      lines(x=c(mean(tau.hat[1:z])+1.96*sd(tau.hat[1:z])/sqrt(z),
                mean(tau.hat[1:z])+1.96*sd(tau.hat[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
      lines(x=c(mean(tau.hat[1:z])-1.96*sd(tau.hat[1:z])/sqrt(z),
                mean(tau.hat[1:z])-1.96*sd(tau.hat[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
      
      hist(theta.hat.v2[1:z],30,xlim=c(3*theta,-3*theta),main="congestion coefficient (uncorrected)",xlab="",ylab="")
      lines(x=c(theta,theta),y=c(0,Z),col='red')
      lines(x=c(mean(theta.hat.v2[1:z]),mean(theta.hat.v2[1:z])),y=c(0,500),col='blue',lty=1)
      lines(x=c(mean(theta.hat.v2[1:z])+1.96*sd(theta.hat.v2[1:z])/sqrt(z),
                mean(theta.hat.v2[1:z])+1.96*sd(theta.hat.v2[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
      lines(x=c(mean(theta.hat.v2[1:z])-1.96*sd(theta.hat.v2[1:z])/sqrt(z),
                mean(theta.hat.v2[1:z])-1.96*sd(theta.hat.v2[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
      
      hist(theta.hat.v3[1:z],30,main="congestion coefficient (lagged Y)",xlab="",ylab="")
      lines(x=c(theta,theta),y=c(0,Z),col='red')
      lines(x=c(mean(theta.hat.v3[1:z]),mean(theta.hat.v3[1:z])),y=c(0,500),col='blue',lty=1)
      lines(x=c(mean(theta.hat.v3[1:z])+1.96*sd(theta.hat.v3[1:z])/sqrt(z),
                mean(theta.hat.v3[1:z])+1.96*sd(theta.hat.v3[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
      lines(x=c(mean(theta.hat.v3[1:z])-1.96*sd(theta.hat.v3[1:z])/sqrt(z),
                mean(theta.hat.v3[1:z])-1.96*sd(theta.hat.v3[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
      
      plot(0,0)
      
    }
    
  }
  
  # Table of summary results:
  {
    cat('\nSimulation test summary results [full year]:\n',file=out.file.name,append=TRUE)
    cat('=========================================\n',file=out.file.name,append=TRUE)
    cat('Parameter      True    E[Bias]    E[s.e.]\n',file=out.file.name,append=TRUE)
    cat('-----------------------------------------\n',file=out.file.name,append=TRUE)
    cat(sprintf('\\hs Travel cost coefficient, $\\lambda$ &  %10.3f & %10.5f & %10.5f\\\\\n',lambda,mean(lambda.hat)-lambda,sd(lambda.hat)),file=out.file.name,append=TRUE)
    cat(sprintf('\\hs Congestion coefficient, $\\theta$   &  %10.2f & %10.5f & %10.5f\\\\\n',theta,mean(theta.hat)-theta,sd(theta.hat)),file=out.file.name,append=TRUE)
    cat(sprintf('\\hs Optimal entry fee, $\\tau^\\star$ [\\$/visitor]   &  %10.2f & %10.1f & %10.1f\\\\\n',mean(tau.star),mean(tau.hat-tau.star),sd(tau.hat)),file=out.file.name,append=TRUE)
    cat(sprintf("\\hs Uncorrected congestion coefficient, $\\theta'$ & %10.2f & %10.5f & %10.5f\\\\\n",theta,mean(theta.hat.v2)-theta,sd(theta.hat.v2)),file=out.file.name,append=TRUE)
    cat('=========================================\n',file=out.file.name,append=TRUE)
  }
  
  # pdf figure: simulation histograms
  if(TRUE){
    fig.name <- paste(output.path,'/sim-histograms-fullyear',sep='')
    tikz(paste(fig.name,'.tex',sep=''),
         width      = 6.5,
         height     = 6.5,
         pointsize  = 12,
         standAlone = TRUE)
    
    par(mfrow = c(2,2),        # 2x2 layout
        oma   = c(0, 0, 1, 0), # rows of text at the outer bottom, left, top, right
        mar   = c(4, 4, 2, 1), # rows of text at ticks and to separate plots
        mgp   = c(2, 1, 0) )   # axis label at 2 rows distance, tick labels at 1 row
    
    
    hist(lambda.hat[1:z],xlim=c(2*lambda,0),main="travel cost coefficient",xlab="",ylab="")
    lines(x=c(lambda,lambda),y=c(0,Z),col='red')
    lines(x=c(mean(lambda.hat[1:z]),mean(lambda.hat[1:z])),y=c(0,500),col='blue',lty=1)
    lines(x=c(mean(lambda.hat[1:z])+1.96*sd(lambda.hat[1:z])/sqrt(z),
              mean(lambda.hat[1:z])+1.96*sd(lambda.hat[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
    lines(x=c(mean(lambda.hat[1:z])-1.96*sd(lambda.hat[1:z])/sqrt(z),
              mean(lambda.hat[1:z])-1.96*sd(lambda.hat[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
    
    hist(theta.hat[1:z],30,xlim=c(2*theta,-2*theta),main="congestion coefficient",xlab="",ylab="")
    lines(x=c(theta,theta),y=c(0,Z),col='red')
    lines(x=c(mean(theta.hat[1:z]),mean(theta.hat[1:z])),y=c(0,Z),col='blue',lty=1)
    lines(x=c(mean(theta.hat[1:z])+1.96*sd(theta.hat[1:z])/sqrt(z),
              mean(theta.hat[1:z])+1.96*sd(theta.hat[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
    lines(x=c(mean(theta.hat[1:z])-1.96*sd(theta.hat[1:z])/sqrt(z),
              mean(theta.hat[1:z])-1.96*sd(theta.hat[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
    
    hist(tau.hat[1:z],30,xlim=c(-mean(tau.star[1:z]),3*mean(tau.star[1:z])),main="optimal entry fee",xlab="",ylab="")
    lines(x=c(mean(tau.star[1:z]),mean(tau.star[1:z])),y=c(0,Z),col='red')
    lines(x=c(mean(tau.hat[1:z]),mean(tau.hat[1:z])),y=c(0,500),col='blue',lty=1)
    lines(x=c(mean(tau.hat[1:z])+1.96*sd(tau.hat[1:z])/sqrt(z),
              mean(tau.hat[1:z])+1.96*sd(tau.hat[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
    lines(x=c(mean(tau.hat[1:z])-1.96*sd(tau.hat[1:z])/sqrt(z),
              mean(tau.hat[1:z])-1.96*sd(tau.hat[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
    
    hist(theta.hat.v2[1:z],30,xlim=c(2*theta,-2*theta),main="congestion coefficient (uncorrected)",xlab="",ylab="")
    lines(x=c(theta,theta),y=c(0,Z),col='red')
    lines(x=c(mean(theta.hat.v2[1:z]),mean(theta.hat.v2[1:z])),y=c(0,500),col='blue',lty=1)
    lines(x=c(mean(theta.hat.v2[1:z])+1.96*sd(theta.hat.v2[1:z])/sqrt(z),
              mean(theta.hat.v2[1:z])+1.96*sd(theta.hat.v2[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
    lines(x=c(mean(theta.hat.v2[1:z])-1.96*sd(theta.hat.v2[1:z])/sqrt(z),
              mean(theta.hat.v2[1:z])-1.96*sd(theta.hat.v2[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
    
    dev.off()
    tinytex::latexmk(paste(fig.name,".tex",sep=""))
    file.remove(paste(fig.name,".tex",sep=""))
  }
  
}

#-------------------------------------------------------------------------------
# SIMULATION TEST [may-oct]:
#-------------------------------------------------------------------------------
if(TRUE){
  
  set.seed(123)
  
  Z <- 1000 # number of Monte Carlo iterations
  
  # Parameters:
  if(FALSE){
    sigma  <- .45
    
    beta   <- -8.5
    #                J   F   M   A   M   J   J   A   S   O   N   D
    alpha  <- log(c(.01,.01,.02,.02,.04,.16,.34,.27,.09,.02,.01,.01))
    alpha  <- alpha - alpha[1]       # Normalizes Jan alpha to 0
    
    T      <- 40                     # 40 years of data
    phi    <- seq(0,log(2),length=T) # year fixed effects increase linearly
    phi    <- phi - phi[1]           # Normalizes first phi to 0
    
    J      <- 49                     # Continental US states + DC
    
    lambda <- -0.01                  # travel cost coefficient 
    theta  <- -1.0                   # congestion coefficient
    
  }
  
  # Year and month and may-oct dummy variables:
  {
    D.yrs <- matrix(0,T*12,T)
    tm <- 0
    for(t in 1:T){ for(m in 1:12){ tm <- tm + 1; D.yrs[tm,t] <- 1 } }
    D.yrs <- D.yrs[,-1]
    
    D.mth <- matrix(0,T*12,12)
    D.opn <- matrix(0,T*12,1)
    tm <- 0
    for(t in 1:T){ for(m in 1:12){ tm <- tm + 1; D.mth[tm,m] <- 1 ; if(m>=5 & m<=10){D.opn[tm]<-1}} }
    D.mth <- D.mth[,-1]
    
    # open <- matrix(1,T*12,1)
    open <- D.opn
    
  }
  
  tau.star      <- matrix(0,Z,1)
  Y.star        <- matrix(0,Z,1)
  lambda.hat    <- matrix(0,Z,1)
  theta.hat     <- matrix(0,Z,1)
  theta.hat.v2  <- matrix(0,Z,1)
  tau.hat       <- matrix(0,Z,1)
  Y.hat         <- matrix(0,Z,1)
  
  start.time <- proc.time()[3]
  for(z in 1:Z){
    
    cat('\014',sprintf('Working on %-.0f of %-.0f Monte Carlo iterations [may-oct]\n',z,Z))
    cat('[Estimated time to completion =',sprintf('%-.1f min.]\n',(proc.time()[3]-start.time)/z*(Z-z)/60))
    
    N0    <- floor(rnorm(J,2e8/J,2e7/J)) # population in each state in year 1
    gN    <- .01 # mean population growth rate [per year]
    N     <- matrix(0,J,T)
    N[,1] <- N0
    for(t in 2:T){ N[,t] <- N[,t-1]*(1+gN*exp(rnorm(J)*.1-.5*.1^2)) }
    
    c0 <- runif(J,min=100,max=500)    # travel cost from each state in year 1
    gc <- 0 # mean travel cost growth rate [per year]
    c  <- matrix(0,J,T)
    c[,1] <- c0
    for(t in 2:T){ c[,t] <- c[,t-1]*(1+gc*exp(rnorm(J)*.1-.5*.1^2)) }
    
    # Optimal fee computation (using true parameters and exog variables 
    # N[,T] & c[,T]):
    {
      Y.inputs         <- list()
      Y.inputs$N       <- N[,T]
      Y.inputs$theta   <- theta
      Y.inputs$lambda  <- lambda
      
      Y.0       <- 4*10^6 / 12 # optimal fee for a specific month with given baseline vistation
      gamma     <- log(Y.0/(exp(theta*(Y.0/10^6))*sum(N[,T]*exp(lambda*c[,T]))))
      tau.range <- seq(0,2000,1)
      Y.tau     <- matrix(0,length(tau.range),1)
      V.tau     <- matrix(0,length(tau.range),1)
      k <- 0
      for(tau in tau.range){k <- k + 1
      Y.inputs$beta    <- gamma
      Y.inputs$phi.t   <- 0
      Y.inputs$alpha.m <- 0
      Y.inputs$sigma   <- sigma
      Y.inputs$c       <- c[,T]+tau
      Y.tau[k]         <- Y.fn(Y.inputs)$Ytm
      V.tau[k]         <- Y.tau[k]/(-lambda) + tau*Y.tau[k]
      }
      Y.star[z]   <- Y.tau[which(V.tau==max(V.tau))]
      tau.star[z] <- tau.range[which(V.tau==max(V.tau))]
    }
    
    sim.inputs         <- list()
    sim.inputs$N       <- N[,T]
    sim.inputs$beta    <- beta
    sim.inputs$phi     <- phi
    sim.inputs$alpha   <- alpha
    sim.inputs$theta   <- theta
    sim.inputs$lambda  <- lambda
    sim.inputs$sigma   <- sigma
    sim.inputs$c       <- c[,T]
    
    outs <- sim.fn(sim.inputs)
    Y    <- outs$Y  # [T*12 x 1] total trips in 12 months of T years
    Yj   <- outs$Yj # [T*12 x J] trips from each j in 12 months of T years
    
    # plot(Y,type='l')
    # pause()
    
    # Step 1: travel cost coefficient estimation
    {
      # Regress ln of total visits from each state in final year on travel costs:
      Yj   <- colSums(Yj[(dim(Yj)[1]-11):dim(Yj)[1],]) # Now a [J x 1] vector
      outs <- lm(log(Yj) ~ c[,T])
      lambda.hat[z] <- outs$coefficients[2]
    }
    
    # Step 2: congestion coefficient estimation
    if(TRUE){
      theta.range <- seq(0,theta*3,length=100)
      R2 <- matrix(0,length(theta.range),1)
      k  <- 0
      for(theta.k in theta.range){k <- k + 1
      depv  <- log(Y) - theta.k*(Y/10^6)*open
      outs  <- summary(lm(depv ~ D.yrs + D.mth))
      ehat  <- outs$residuals
      R2[k] <- outs$adj.r.squared
      }
      theta.hat[z] <- theta.range[which.max(R2)] 
      
      depv  <- log(Y)-theta.hat[z]*(Y/10^6)*open
      outs  <- summary(lm(depv ~ D.yrs + D.mth))
      ehat  <- outs$residuals
    }
    
    # Step 2b: congestion coefficient estimation (naive--endogenous Y on rhs)
    if(TRUE){
      depv  <- log(Y)
      YY    <- Y/10^6 * open
      outs  <- summary(lm(depv ~ D.yrs + D.mth + YY))
      theta.hat.v2[z] <- tail(outs$coeff,1)[1]
    }
    
    # Optimal fee estimation (using estimated parameters and exog variables
    # N[,T] & c[,T]) using theta.hat:
    if(TRUE){
      Y.inputs         <- list()
      Y.inputs$N       <- N[,T]
      Y.inputs$theta   <- theta.hat[z]
      Y.inputs$lambda  <- lambda.hat[z]
      
      Y.0   <- 4*10^6/12
      gamma <- log(Y.0/(exp(theta.hat[z]*(Y.0/10^6))*sum(N[,T]*exp(lambda.hat[z]*c[,T]))))
      tau.range <- seq(0,2000,1)
      Y.tau     <- matrix(0,length(tau.range),1)
      V.tau     <- matrix(0,length(tau.range),1)
      k <- 0
      for(tau in tau.range){k <- k + 1
      Y.inputs$beta    <- gamma
      Y.inputs$phi.t   <- 0
      Y.inputs$alpha.m <- 0
      Y.inputs$sigma   <- sigma
      Y.inputs$c       <- c[,T]+tau
      Y.tau[k]         <- Y.fn(Y.inputs)$Ytm
      V.tau[k]         <- Y.tau[k]/(-lambda.hat[z]) + tau*Y.tau[k]
      }
      Y.hat[z]   <- Y.tau[which(V.tau==max(V.tau))]
      tau.hat[z] <- tau.range[which(V.tau==max(V.tau))]
    }
    
    # Histograms:
    if(z>1){
      
      par(mfrow=c(2,2))
      
      hist(lambda.hat[1:z],xlim=c(2*lambda,0),main="travel cost coefficient",xlab="",ylab="")
      lines(x=c(lambda,lambda),y=c(0,Z),col='red')
      lines(x=c(mean(lambda.hat[1:z]),mean(lambda.hat[1:z])),y=c(0,500),col='blue',lty=1)
      lines(x=c(mean(lambda.hat[1:z])+1.96*sd(lambda.hat[1:z])/sqrt(z),
                mean(lambda.hat[1:z])+1.96*sd(lambda.hat[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
      lines(x=c(mean(lambda.hat[1:z])-1.96*sd(lambda.hat[1:z])/sqrt(z),
                mean(lambda.hat[1:z])-1.96*sd(lambda.hat[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
      
      hist(theta.hat[1:z],30,xlim=c(3*theta,-3*theta),main="congestion coefficient",xlab="",ylab="")
      lines(x=c(theta,theta),y=c(0,Z),col='red')
      lines(x=c(mean(theta.hat[1:z]),mean(theta.hat[1:z])),y=c(0,Z),col='blue',lty=1)
      lines(x=c(mean(theta.hat[1:z])+1.96*sd(theta.hat[1:z])/sqrt(z),
                mean(theta.hat[1:z])+1.96*sd(theta.hat[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
      lines(x=c(mean(theta.hat[1:z])-1.96*sd(theta.hat[1:z])/sqrt(z),
                mean(theta.hat[1:z])-1.96*sd(theta.hat[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
      
      hist(tau.hat[1:z],30,xlim=c(-mean(tau.star[1:z]),3*mean(tau.star[1:z])),main="optimal entry fee",xlab="",ylab="")
      lines(x=c(mean(tau.star[1:z]),mean(tau.star[1:z])),y=c(0,Z),col='red')
      lines(x=c(mean(tau.hat[1:z]),mean(tau.hat[1:z])),y=c(0,500),col='blue',lty=1)
      lines(x=c(mean(tau.hat[1:z])+1.96*sd(tau.hat[1:z])/sqrt(z),
                mean(tau.hat[1:z])+1.96*sd(tau.hat[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
      lines(x=c(mean(tau.hat[1:z])-1.96*sd(tau.hat[1:z])/sqrt(z),
                mean(tau.hat[1:z])-1.96*sd(tau.hat[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
      
      hist(theta.hat.v2[1:z],30,xlim=c(3*theta,-3*theta),main="congestion coefficient (uncorrected)",xlab="",ylab="")
      lines(x=c(theta,theta),y=c(0,Z),col='red')
      lines(x=c(mean(theta.hat.v2[1:z]),mean(theta.hat.v2[1:z])),y=c(0,500),col='blue',lty=1)
      lines(x=c(mean(theta.hat.v2[1:z])+1.96*sd(theta.hat.v2[1:z])/sqrt(z),
                mean(theta.hat.v2[1:z])+1.96*sd(theta.hat.v2[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
      lines(x=c(mean(theta.hat.v2[1:z])-1.96*sd(theta.hat.v2[1:z])/sqrt(z),
                mean(theta.hat.v2[1:z])-1.96*sd(theta.hat.v2[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
      
    }
    
  }
  
  # Table of summary results:
  {
    cat('\nSimulation test summary results [May-Oct]:\n',file=out.file.name,append=TRUE)
    cat('=========================================\n',file=out.file.name,append=TRUE)
    cat('Parameter      True    E[Bias]    E[s.e.]\n',file=out.file.name,append=TRUE)
    cat('-----------------------------------------\n',file=out.file.name,append=TRUE)
    cat(sprintf('\\hs Travel cost coefficient, $\\lambda$ &  %10.3f & %10.5f & %10.5f\\\\\n',lambda,mean(lambda.hat)-lambda,sd(lambda.hat)),file=out.file.name,append=TRUE)
    cat(sprintf('\\hs Congestion coefficient, $\\theta$   &  %10.2f & %10.5f & %10.5f\\\\\n',theta,mean(theta.hat)-theta,sd(theta.hat)),file=out.file.name,append=TRUE)
    cat(sprintf('\\hs Optimal entry fee, $\\tau^\\star$ [\\$/visitor]   &  %10.2f & %10.1f & %10.1f\\\\\n',mean(tau.star),mean(tau.hat-tau.star),sd(tau.hat)),file=out.file.name,append=TRUE)
    cat(sprintf("\\hs Uncorrected congestion coefficient, $\\theta'$ & %10.2f & %10.5f & %10.5f\\\\\n",theta,mean(theta.hat.v2)-theta,sd(theta.hat.v2)),file=out.file.name,append=TRUE)
    cat('=========================================\n',file=out.file.name,append=TRUE)
  }
  
  # pdf figure: simulation histograms
  if(TRUE){
    fig.name <- paste(output.path,'/sim-histograms-mayoct',sep='')
    tikz(paste(fig.name,'.tex',sep=''),
         width      = 6.5,
         height     = 6.5,
         pointsize  = 12,
         standAlone = TRUE)
    
    par(mfrow = c(2,2),        # 2x2 layout
        oma   = c(0, 0, 1, 0), # rows of text at the outer bottom, left, top, right
        mar   = c(4, 4, 2, 1), # rows of text at ticks and to separate plots
        mgp   = c(2, 1, 0) )   # axis label at 2 rows distance, tick labels at 1 row
    
    
    hist(lambda.hat[1:z],xlim=c(2*lambda,0),main="travel cost coefficient",xlab="",ylab="")
    lines(x=c(lambda,lambda),y=c(0,Z),col='red')
    lines(x=c(mean(lambda.hat[1:z]),mean(lambda.hat[1:z])),y=c(0,500),col='blue',lty=1)
    lines(x=c(mean(lambda.hat[1:z])+1.96*sd(lambda.hat[1:z])/sqrt(z),
              mean(lambda.hat[1:z])+1.96*sd(lambda.hat[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
    lines(x=c(mean(lambda.hat[1:z])-1.96*sd(lambda.hat[1:z])/sqrt(z),
              mean(lambda.hat[1:z])-1.96*sd(lambda.hat[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
    
    hist(theta.hat[1:z],30,xlim=c(2*theta,-2*theta),main="congestion coefficient",xlab="",ylab="")
    lines(x=c(theta,theta),y=c(0,Z),col='red')
    lines(x=c(mean(theta.hat[1:z]),mean(theta.hat[1:z])),y=c(0,Z),col='blue',lty=1)
    lines(x=c(mean(theta.hat[1:z])+1.96*sd(theta.hat[1:z])/sqrt(z),
              mean(theta.hat[1:z])+1.96*sd(theta.hat[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
    lines(x=c(mean(theta.hat[1:z])-1.96*sd(theta.hat[1:z])/sqrt(z),
              mean(theta.hat[1:z])-1.96*sd(theta.hat[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
    
    hist(tau.hat[1:z],30,xlim=c(-mean(tau.star[1:z]),3*mean(tau.star[1:z])),main="optimal entry fee",xlab="",ylab="")
    lines(x=c(mean(tau.star[1:z]),mean(tau.star[1:z])),y=c(0,Z),col='red')
    lines(x=c(mean(tau.hat[1:z]),mean(tau.hat[1:z])),y=c(0,500),col='blue',lty=1)
    lines(x=c(mean(tau.hat[1:z])+1.96*sd(tau.hat[1:z])/sqrt(z),
              mean(tau.hat[1:z])+1.96*sd(tau.hat[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
    lines(x=c(mean(tau.hat[1:z])-1.96*sd(tau.hat[1:z])/sqrt(z),
              mean(tau.hat[1:z])-1.96*sd(tau.hat[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
    
    hist(theta.hat.v2[1:z],30,xlim=c(2*theta,-2*theta),main="congestion coefficient (uncorrected)",xlab="",ylab="")
    lines(x=c(theta,theta),y=c(0,Z),col='red')
    lines(x=c(mean(theta.hat.v2[1:z]),mean(theta.hat.v2[1:z])),y=c(0,500),col='blue',lty=1)
    lines(x=c(mean(theta.hat.v2[1:z])+1.96*sd(theta.hat.v2[1:z])/sqrt(z),
              mean(theta.hat.v2[1:z])+1.96*sd(theta.hat.v2[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
    lines(x=c(mean(theta.hat.v2[1:z])-1.96*sd(theta.hat.v2[1:z])/sqrt(z),
              mean(theta.hat.v2[1:z])-1.96*sd(theta.hat.v2[1:z])/sqrt(z)),y=c(0,Z),col='blue',lty=2)
    
    
    dev.off()
    tinytex::latexmk(paste(fig.name,".tex",sep=""))
    file.remove(paste(fig.name,".tex",sep=""))
  }
  
}

#-------------------------------------------------------------------------------
# ESTIMATION:
#-------------------------------------------------------------------------------
if(TRUE){
  
  set.seed(123)
  
  BS                  <- 1e3 # number of iterations for bootstrap
  
  adults.per.vehicle  <- 2.7   # 2.7 adults
  persons.per.vehicle <- 3.3   # plus 0.6 kids
  fraction.of.wage    <- (1/3 + 1/2) / 2
  avg.speed           <- 55   # [mi/hr]
  adults.per.hh       <- 2.09 # in 2018 (253815197 / 121520180)
  working.hrs         <- 2080 # [hr/yr]
  cost.per.mi         <- .25  # [$/mi]
  
  # For CV vs Y fig:
  years   <- 1979:2023
  data    <- data.fn()
  Ytm.all <- data$Ytm
  
  # For analysis:
  years <- c(1979:2019,2021,2023)
  # years <- 1979:2023
  T <- length(years)
  
  # IMPORT DATA:
  {
    data <- data.fn()
    N    <- data$N
    d    <- data$d
    Ytm  <- data$Ytm
    yj   <- data$yj
  }
  
  # YEAR AND MONTH AND MAY-OCT DUMMY VARIABLES:
  {
    D.yrs <- matrix(0,T*12,T)
    tm <- 0
    for(t in 1:T){ 
      for(m in 1:12){ 
        tm <- tm + 1
        D.yrs[tm,t] <- 1 
      } 
    }
    D.yrs <- D.yrs[,-40] # Exclude 2018
    
    D.mth <- matrix(0,T*12,12)
    D.opn <- matrix(0,T*12,1) # open season (may-oct) dummy
    tm <- 0
    for(t in 1:T){
      for(m in 1:12){
        tm <- tm + 1
        D.mth[tm,m] <- 1
        if(m>=5 & m<=10){D.opn[tm] <- 1}
      }
    }
    D.mth <- D.mth[,-7]            # Exclude July
    
    open <- matrix(1,T,1)  # include all months
    open <- D.opn                  # include only open months (may-oct)
  }
  
  # pdf figure: CV vs Y
  if(TRUE){
    
    CV    <- vector(length=length(Ytm.all)/12)
    Ytm.sum <- vector(length=length(Ytm.all)/12)
    for(i in 1:(length(Ytm.all)/12)){
      zzz <- ((i-1)*12+1):(i*12)       # Pulls out correct months for a year
      zzz <- zzz[-c(1,2,3,4,11,12)]    # Drop jan-apr and nov-dec
      CV[i] <- sd(Ytm.all[zzz])/mean(Ytm.all[zzz]) # Calculates CV for each year across months
      Ytm.sum[i] <- sum(Ytm.all[zzz])
    }
    
    outs <- lm(CV ~ Ytm.sum)
    
    fig.name <- paste(output.path,'/cv-vs-y',sep='')
    tikz(paste(fig.name,'.tex',sep=''),
         width      = 4.5,
         height     = 4.5,
         pointsize  = 12,
         standAlone = TRUE)
    
    #par(mfrow=c(1,2))
    
    par(mar=c(3.5,3.5,1,1)) # figure border whitespace [bottom,left,top,right]
    
    bottom <- 0.0
    top    <- 1.0
    left   <- 1500000
    right  <- 5000000
    
    plot(Ytm.sum,CV,
         type='p',
         lwd=1,
         main='',
         xlab='',
         ylab='',
         axes=FALSE,
         cex=1,
         xlim=c(left,right),
         ylim=c(bottom,top),
         cex.axis=1,
         cex.lab =1)
    
    abline(lm(CV ~ Ytm.sum),lty=3,lwd=1,col='black')
    
    axis(1,las=1,at=seq(1.5,5,.5)*1e6,labels=seq(1.5,5,.5),cex.axis=1) # Draw x axis
    axis(2,las=1,labels=TRUE,cex.axis=1) # Draw y axis
    
    mtext('$CV$',side=2, line=2.5, cex=1, las=0) # y-axis label
    mtext('$Y$ [millions]' ,side=1,line=2.5, cex=1, las=1) # x-axis label
    
    text(left+(right-left)*.95,bottom+(top-bottom)*.95,sprintf('$b$ = %-.3g',summary(outs)$coefficients[2,1]),pos=2)
    text(left+(right-left)*.95,bottom+(top-bottom)*.90,sprintf('$s.e.$ = %-.3g',summary(outs)$coefficients[2,2]),pos=2)
    text(left+(right-left)*.95,bottom+(top-bottom)*.85,sprintf('$R^2$ = %-.4f',summary(outs)$r.squared),pos=2)
    
    years.all <- 1979:2023
    zz <- which(years.all>=2020 & years.all<=2022)
    text(Ytm.sum[zz], CV[zz], labels=years.all[zz], cex= 1, pos=4)
    zz <- which(years.all==2023)
    text(Ytm.sum[zz], CV[zz], labels=years.all[zz], cex= 1, pos=1)
    
    dev.off()
    tinytex::latexmk(paste(fig.name,".tex",sep=""))
    file.remove(paste(fig.name,".tex",sep=""))
  } 
  
  # CONGESTION ESTIMATION:
  {
    # Weather data:
    {
      # Weather data downloaded from: https://www.ncei.noaa.gov/cdo-web/search
      weather <- read.csv(paste(input.path,"/weather.csv",sep=""))
      prcp    <- c()
      tmin    <- c()
      tmax    <- c()
      yr      <- weather$YEAR
      mo      <- weather$MONTH
      yrs     <- sort(unique(yr))
      mos     <- sort(unique(mo))
      t <- 0
      for(yr in years){
        for(mo in mos){
          t <- t + 1
          prcp[t] <- mean(weather$PRCP[which(weather$YEAR==yr & weather$MONTH==mo)],na.rm=TRUE)
          tmin[t] <- mean(weather$TMIN[which(weather$YEAR==yr & weather$MONTH==mo)],na.rm=TRUE)
          tmax[t] <- mean(weather$TMAX[which(weather$YEAR==yr & weather$MONTH==mo)],na.rm=TRUE)
        }
      }
      
    }
    
    theta.range <- -seq(0,4,length=500)
    R2 <- matrix(0,length(theta.range),1)
    k <- 0
    for(theta.k in theta.range){k <- k + 1
    depv  <- log(Ytm)-theta.k*(Ytm/10^6) * open
    outs  <- summary(lm(depv ~ D.yrs + D.mth + prcp + tmin + tmax))
    ehat  <- outs$residuals
    R2[k] <- outs$adj.r.squared
    }
    theta.hat <- theta.range[which.max(R2)] 
    
    depv      <- log(Ytm)-theta.hat*(Ytm/10^6)*open
    outs      <- lm(depv ~ D.yrs + D.mth + prcp + tmin + tmax)
    depv.hat  <- predict(outs)
    Ytm.hat   <- Ytm * 0
    for(tm in 1:T){
      YYtm <- seq(1,1e6,10)
      e2   <- (depv.hat[tm]-log(YYtm)+theta.hat*(YYtm/10^6)*open[t])^2
      Ytm.hat[tm] <- YYtm[which.min(e2)]
    }
    # Counterfactual (w/o congestion) predicted visits:
    Ytmc <- Ytm*exp(-theta.hat*(Ytm/10^6)*open)
    
    # Step 2 estimation [version 2]:
    if(TRUE){
      depv <- log(Ytm)
      YYtm <- Ytm/10^6
      outs <- summary(lm(depv ~ D.yrs + D.mth + prcp + tmin + tmax + YYtm*open))
      theta.hat.v2 <- tail(outs$coeff,1)[1]
    }
    
    # pdf figure:
    {
      fig.name <- paste(output.path,'/congestion-peak',sep='')
      tikz(paste(fig.name,'.tex',sep=''),
           width      = 5.5,
           height     = 5.5,
           pointsize  = 12,
           standAlone = TRUE)
      
      par(mfrow = c(1,1),        # 1x1 layout
          oma   = c(0, 0, 0, 0), # rows of text at the outer bottom, left, top, right
          mar   = c(4, 4.5, 1, 1), # rows of text at ticks and to separate plots
          mgp   = c(2, 1, 0) )   # axis label at 2 rows distance, tick labels at 1 row
      
      left   <- min(theta.range)
      right  <- max(theta.range)
      bottom <- min(R2)* 0 + 0.977
      top    <- max(R2)* 0 + 0.981
      
      plot(theta.range,R2,
           type='l',
           main='',
           cex.main = 1,
           xlab='',
           ylab='',
           axes=FALSE,
           ylim = c(bottom,top),
           xlim = c(left,right),
           cex.axis=1,
           cex.lab =1)
      
      lines(c(theta.hat,theta.hat),c(0,max(R2)),lty=3)
      
      axis(1,las=1,labels=TRUE,cex.axis=1) # Draw x axis
      axis(2,las=1,labels=TRUE,cex.axis=1) # Draw y axis
      
      mtext('$\\tilde{R}^2$'    , side=2, line=3.5, cex=1, las=0) # y-axis label
      mtext('$\\theta$ [$\\times 10^6$]', side=1, line=2.5, cex=1, las=1) # x-axis label
      
      #text(left+(right-left)*.95,bottom+(top-bottom)*.95,sprintf('$\\tilde{R}^2$  = %-.3f',adjR2.stage1),pos=2)
      #text(left+(right-left)*.95,bottom+(top-bottom)*.90,sprintf('$\\hat{\\lambda}$ = %-.5f',lambda.hat),pos=2)
      
      dev.off()
      tinytex::latexmk(paste(fig.name,".tex",sep=""))
      file.remove(paste(fig.name,".tex",sep=""))
    }
    
  }
  
  # RESIDUALS BOOTSTRAP AND BIAS CORRECTION:
  if(TRUE){
    
    depv0  <- log(Ytm)-theta.hat*(Ytm/10^6)*open
    outs0  <- lm(depv0 ~ D.yrs + D.mth + prcp + tmin + tmax)
    ehat0  <- outs0$residuals
    hats0  <- predict(outs0)
    
    beta.hat  <- outs0$coefficients[1]
    phi.hat   <- outs0$coefficients[2:(length(years))]
    alpha.hat <- outs0$coefficients[(1+length(years)):(length(years)+11)]
    omega.hat <- tail(outs0$coefficients,3)
    
    Ytmhat0 <- c()
    for(t in 1:length(depv0)){
      YYtm <- seq(100,4e6,100)
      e2   <- (hats0[t]-log(YYtm)+theta.hat*(YYtm/10^6)*open[t])^2
      Ytmhat0[t] <- YYtm[which.min(e2)]
    }
    
    plot(Ytm,Ytmhat0);lines(sort(Ytm),sort(Ytm),lty=3)
    ehatYtm0 <- Ytmhat0-Ytm
    
    beta.BS  <- matrix(0,BS,1)
    phi.BS   <- matrix(0,BS,length(years)-1)
    alpha.BS <- matrix(0,BS,11)
    omega.BS <- matrix(0,BS,3)
    theta.BS <- matrix(0,BS,1)
    start.time <- proc.time()[3]
    for(bs in 1:BS){
      
      # Re-sample residuals clustered by month:
      {
        depv <- matrix(0,length(depv0),1)
        for(month in 1:6){
          tt <- seq(month,length(depv0),12)
          depv[tt] <- hats0[tt] + sample(ehat0[tt],length(tt),replace=TRUE)
        }
        
        depv <- hats0 + sample(ehat0,length(depv),replace=TRUE)
        
        Ytmbs <- c()
        for(t in 1:length(depv)){
          YYtm <- seq(100,4e6,100)
          e2 <- (depv[t]-log(YYtm)+theta.hat*(YYtm/10^6)*open[t])^2
          Ytmbs[t] <- YYtm[which.min(e2)]
        }
        
      }
      
      # Re-estimate congestion coefficient:
      {
        R2 <- matrix(0,length(theta.range),1)
        k <- 0
        for(theta.k in theta.range){k <- k + 1
        depv.k <- log(Ytmbs)-theta.k*(Ytmbs/10^6)*open
        outs   <- summary(lm(depv.k ~ D.yrs + D.mth + prcp + tmin + tmax))
        R2[k]  <- outs$adj.r.squared
        }
        theta.BS[bs] <- theta.range[which.max(R2)] 
        depv.k       <- log(Ytmbs)-theta.BS[bs]*(Ytmbs/10^6)*open
        outs         <- lm(depv.k ~ D.yrs + D.mth + prcp + tmin + tmax)
        beta.BS[bs]  <- outs$coefficients[1]
        phi.BS[bs,]  <- outs$coefficients[2:(length(years))]
        alpha.BS[bs,]<- outs$coefficients[(1+length(years)):(length(years)+11)]
        omega.BS[bs,]<- tail(outs$coefficients,3)
      }
      
      par(mfrow=c(2,2))
      plot(hats0,depv); lines(sort(depv0),sort(depv0),lty=3)
      plot(Ytmhat0,Ytmbs);  lines(sort(Ytmhat0),sort(Ytmhat0),lty=3)
      hist(theta.BS[1:bs])
      lines(c(theta.hat,theta.hat),c(0,bs),col='red')
      lines(c(mean(theta.BS),mean(theta.BS)),c(0,bs),col='blue')
      lines(c(mean(theta.BS)-1.96*sd(theta.BS)/sqrt(bs),mean(theta.BS)-1.96*sd(theta.BS)/sqrt(bs)),c(0,bs),col='blue',lty=3)
      lines(c(mean(theta.BS)+1.96*sd(theta.BS)/sqrt(bs),mean(theta.BS)+1.96*sd(theta.BS)/sqrt(bs)),c(0,bs),col='blue',lty=3)
      plot(0,0)
      
      pause(.01)
      
      # Display progress in console:
      if(bs>1){
        cat('\014Working on bootstrap rep',sprintf('%-.0f',bs),'of',sprintf('%-.0f\n',BS))
        cat(sprintf('So far: lambda.hat.BC = %-.4f, se = %-.4f\n',theta.hat-(mean(theta.BS[1:bs])-theta.hat),sd(theta.BS[1:bs])))
        cat('[Estimated time to completion =',sprintf('%-.1f min.]\n',(proc.time()[3]-start.time)/bs*(BS-bs)/60))
      }
      
    }
    
    # Table of bootstrap results:
    {
      
      beta.se  <- sd(beta.BS)
      phi.se   <- apply(phi.BS,2,sd)
      alpha.se <- apply(alpha.BS,2,sd)
      omega.se <- apply(omega.BS,2,sd)
      theta.se <- sd(theta.BS)
      
      cat('\nResiduals bootstrap results:\n',file=out.file.name,append=TRUE)
      cat('============================================================\n',file=out.file.name,append=TRUE)
      cat('Parameter    Estimate      E[Bias]      E[s.e.]   t-stat[BC]\n',file=out.file.name,append=TRUE)
      cat('------------------------------------------------------------\n',file=out.file.name,append=TRUE)
      cat(sprintf("$\\beta$ &  %10.2f & %10.5f & %10.5f & %10.5f\\\\\n",
                  beta.hat ,mean(beta.BS)-beta.hat,beta.se,(beta.hat-(mean(beta.BS-beta.hat)))/beta.se),
          file=out.file.name,append=TRUE)
      for(k in 1:length(phi.hat)){
        cat(sprintf("$\\phi$%-.0f & %10.2f & %10.5f & %10.5f & %10.5f\\\\\n",
                    k,phi.hat[k],mean(phi.BS[,k])-phi.hat[k],sd(phi.BS[,k]),(phi.hat[k]-(mean(phi.BS[,k]-phi.hat[k])))/phi.se[k]),
            file=out.file.name,append=TRUE)
      }
      for(k in 1:length(alpha.hat)){
        cat(sprintf("$\\alpha$%-.0f & %10.2f & %10.5f & %10.5f & %10.5f\\\\\n",
                    k,alpha.hat[k],mean(alpha.BS[,k])-alpha.hat[k],sd(alpha.BS[,k]),(alpha.hat[k]-(mean(alpha.BS[,k]-alpha.hat[k])))/alpha.se[k]),
            file=out.file.name,append=TRUE)
      }
      for(k in 1:length(omega.hat)){
        cat(sprintf("$\\omega$%-.0f & %10.5f & %10.5f & %10.5f & %10.5f\\\\\n",
                    k,omega.hat[k],mean(omega.BS[,k])-omega.hat[k],sd(omega.BS[,k]),(omega.hat[k]-(mean(omega.BS[,k]-omega.hat[k])))/omega.se[k]),
            file=out.file.name,append=TRUE)
      }
      cat(sprintf("$\\theta$ &  %10.5f & %10.5f & %10.5f & %10.5f\\\\\n",
                  theta.hat,mean(theta.BS)-theta.hat,sd(theta.BS),(theta.hat-(mean(theta.BS-theta.hat)))/theta.se),
          file=out.file.name,append=TRUE)
      cat('============================================================\n',file=out.file.name,append=TRUE)
    }
    
    # pdf figure: bootstrap histogram
    {
      fig.name <- paste(output.path,'/bootstrap-histogram',sep='')
      tikz(paste(fig.name,'.tex',sep=''),
           width      = 6.5,
           height     = 6.5,
           pointsize  = 12,
           standAlone = TRUE)
      
      par(mfrow = c(1,1),        # 1x1 layout
          oma   = c(0, 0, 1, 0), # rows of text at the outer bottom, left, top, right
          mar   = c(4, 4, 2, 1), # rows of text at ticks and to separate plots
          mgp   = c(2, 1, 0) )   # axis label at 2 rows distance, tick labels at 1 row
      
      hist(theta.BS[1:bs],xlim=c(3*theta.hat,0),main='Residuals bootstrap distribution of $\\hat{\\theta}$',xlab='',ylab='')
      lines(c(theta.hat,theta.hat),c(0,bs),col='red')
      lines(c(mean(theta.BS),mean(theta.BS)),c(0,bs),col='blue')
      lines(c(mean(theta.BS)-1.96*sd(theta.BS)/sqrt(bs),mean(theta.BS)-1.96*sd(theta.BS)/sqrt(bs)),c(0,bs),col='blue',lty=2)
      lines(c(mean(theta.BS)+1.96*sd(theta.BS)/sqrt(bs),mean(theta.BS)+1.96*sd(theta.BS)/sqrt(bs)),c(0,bs),col='blue',lty=2)
      
      dev.off()
      tinytex::latexmk(paste(fig.name,".tex",sep=""))
      file.remove(paste(fig.name,".tex",sep=""))
    }
    
    # pdf figure: plot of estimated year fixed effects
    if(TRUE){
      fig.name <- paste(output.path,'/year-fixed-effects',sep='')
      tikz(paste(fig.name,'.tex',sep=''),
           width      = 6.5,
           height     = 3.5,
           pointsize  = 12,
           standAlone = TRUE)
      
      par(mfrow = c(1,1),        # 1x1 layout
          oma   = c(0, 0, 1, 0), # rows of text at the outer bottom, left, top, right
          mar   = c(4, 4, 2, 1), # rows of text at ticks and to separate plots
          mgp   = c(2, 1, 0) )   # axis label at 2 rows distance, tick labels at 1 row
      
      top <- max(phi.hat+1.96*phi.se)
      bot <- min(phi.hat-1.96*phi.se)
      
      phi.years <- years[-which(years==2018)]
      plot(years,c(phi.hat[which(phi.years<2018)],0,phi.hat[which(phi.years>2018)]),
           pch=1,
           xaxt="n",
           xlab="",
           ylab="$\\hat{\\phi}$",
           ylim=c(bot,top))
      points(2018,0,pch=1)
      points(c(years[which(years<2018)],years[which(years>2018)]),
             c(phi.hat[which(phi.years<2018)],phi.hat[which(phi.years>2018)]),pch=16)
      for(k in 1:length(phi.hat)){
        lines(c(phi.years[k],phi.years[k]),
              c(phi.hat[k]-1.96*phi.se[k],phi.hat[k]+1.96*phi.se[k]))
      }
      axis(1,at=years,labels=years,las=2,cex.axis=0.8)
      
      dev.off()
      tinytex::latexmk(paste(fig.name,".tex",sep=""))
      file.remove(paste(fig.name,".tex",sep=""))
    }
    
    # pdf figure: plot of estimated month fixed effects
    if(TRUE){
      fig.name <- paste(output.path,'/month-fixed-effects',sep='')
      tikz(paste(fig.name,'.tex',sep=''),
           width      = 6.5,
           height     = 3.5,
           pointsize  = 12,
           standAlone = TRUE)
      
      par(mfrow = c(1,1),        # 1x1 layout
          oma   = c(0, 0, 1, 0), # rows of text at the outer bottom, left, top, right
          mar   = c(4, 4, 2, 1), # rows of text at ticks and to separate plots
          mgp   = c(2, 1, 0) )   # axis label at 2 rows distance, tick labels at 1 row
      
      Months <-  c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
      alpha.Months <- c('Jan','Feb','Mar','Apr','May','Jun','Aug','Sep','Oct','Nov','Dec')
      
      top <- max(1.96*alpha.se)
      bot <- min(alpha.hat-1.96*alpha.se)
      
      plot(1:12,c(alpha.hat[1:6],0,alpha.hat[7:11]),
           pch=1,
           xaxt="n",
           xlab="",
           ylab="$\\hat{\\alpha}$",
           ylim=c(bot,top))
      points(7,0,pch=1)
      points(c(1:6,8:12),
             c(alpha.hat[1:6],alpha.hat[7:11]),pch=16)
      for(k in 1:6){
        lines(c(k,k),c(alpha.hat[k]-1.96*alpha.se[k],alpha.hat[k]+1.96*alpha.se[k]))
      }
      for(k in 7:11){
        lines(c(k+1,k+1),c(alpha.hat[k]-1.96*alpha.se[k],alpha.hat[k]+1.96*alpha.se[k]))
      }
      axis(1,at=(1:12),labels=Months,las=2,cex.axis=0.8)
      
      dev.off()
      tinytex::latexmk(paste(fig.name,".tex",sep=""))
      file.remove(paste(fig.name,".tex",sep=""))
    }
    
    # Bias-corrected theta estimate:
    theta.hat.BC <- theta.hat-(mean(theta.BS)-theta.hat)
    
    # Counterfactual (w/o congestion) predicted visits:
    Ytmc <- Ytm*exp(-theta.hat.BC*(Ytm/10^6)*open)
    
  }
  
  # TRAVEL COST ESTIMATION:
  for(transform in c('none','linear','exponential')){
    
    # IMPORT DATA:
    {
      data <- data.fn()
      N    <- data$N
      d    <- data$d
      Ytm  <- data$Ytm
      yj   <- data$yj
    }
    
    tau0 <- 35/persons.per.vehicle
    
    # YEAR AND MONTH DUMMY VARIABLES:
    {
      D.yrs <- matrix(0,T*12,T)
      tm <- 0
      for(t in 1:T){ for(m in 1:12){ tm <- tm + 1; D.yrs[tm,t] <- 1 } }
      D.yrs <- D.yrs[,-1]
      
      D.mth <- matrix(0,T*12,12)
      tm <- 0
      for(t in 1:T){ for(m in 1:12){ tm <- tm + 1; D.mth[tm,m] <- 1 } }
      D.mth <- D.mth[,-1]
      
    }
    
    if(transform=='none'){
      
      c <- d + tau0
      
    }
    
    if(transform=='linear'){
      
      # kappa optimized by hand 
      # (needs to be redone whenever travel cost assumptions are changed)
      # kappa <- 0.00070  # R2 = 0.8171
      kappa <- 0.00071  # R2 = 0.8173
      # kappa <- 0.00072  # R2 = 0.8158
      
      c     <- d * (1-kappa*(d-min(d))) + tau0
      
    }
    
    if(transform=='exponential'){
      
      # delta optimized by hand
      # (needs to be redone whenever travel cost assumptions are changed)
      # delta <- 0.00135 # R2 = 0.8346
      delta <- 0.00140 # R2 = 0.8355
      # delta <- 0.00145 # R2 = 0.8320
      
      c     <- d * exp(-delta*(d-min(d))) + tau0
      
    }
    
    # -- regression:
    { 
      outs         <- lm(log(yj) ~ c) # c includes tau0
      lambda.hat   <-  outs$coefficients[2]
      lambda.se    <-  summary(outs)$coefficients[2,2]
      logyj.hat    <-  fitted(outs)
      adjR2.stage1 <-  summary(outs)$adj.r.squared
      
      # pdf figure:
      {
        fig.name <- paste(output.path,'/lny-vs-c-',transform,'',sep='')
        tikz(paste(fig.name,'.tex',sep=''),
             width      = 4.5,
             height     = 2.5,
             pointsize  = 12,
             standAlone = TRUE)
        
        par(mfrow = c(1,1),        # 1x1 layout
            oma   = c(0, 0, 0, 0), # rows of text at the outer bottom, left, top, right
            mar   = c(4, 4, 1, 1), # rows of text at ticks and to separate plots
            mgp   = c(2, 1, 0) )   # axis label at 2 rows distance, tick labels at 1 row
        
        left   <- 0
        right  <- 1100
        bottom <- -6
        top    <- -2
        
        plot(c,log(yj),
             main='',
             cex.main = 1,
             xlab='',
             ylab='',
             axes=FALSE,
             ylim = c(bottom,top),
             xlim = c(left,right),
             cex.axis=1,
             cex.lab =1)
        
        axis(1,las=1,at=seq(0,1100,100),labels=seq(0,1100,100),cex.axis=1) # Draw x axis
        axis(2,las=1,labels=TRUE,cex.axis=1) # Draw y axis
        
        mtext('$\\ln{(y)}$', side=2, line=2.5, cex=1, las=0) # y-axis label
        
        if(transform=='none')       {mtext('$c=d$'                              , side=1, line=2.5, cex=1, las=1)} # x-axis label
        if(transform=='exponential'){mtext('$c=de^{-\\delta(d-\\underline{d})}$', side=1, line=2.5, cex=1, las=1)} # x-axis label
        if(transform=='linear')     {mtext('$c=d(1-\\kappa(d-\\underline{d}))$' , side=1, line=2.5, cex=1, las=1)} # x-axis label
        
        order <- sort(c,index.return=TRUE)$ix
        
        lines(c[order],logyj.hat[order],lty=2)
        
        text(left+(right-left)*.95,bottom+(top-bottom)*.95,sprintf('$\\tilde{R}^2$  = %-.4f',adjR2.stage1),pos=2)
        text(left+(right-left)*.95,bottom+(top-bottom)*.85,sprintf('$\\hat{\\lambda}$ = %-.5f',lambda.hat),pos=2)
        text(left+(right-left)*.95,bottom+(top-bottom)*.75,sprintf('$s.e.$ = %-.5f',lambda.se),pos=2)
        if(transform=='linear'){
          text(left+(right-left)*.95,bottom+(top-bottom)*.65,sprintf('$\\kappa$ = %-.6f',kappa),pos=2)
        }
        if(transform=='exponential'){
          text(left+(right-left)*.95,bottom+(top-bottom)*.65,sprintf('$\\delta$ = %-.5f',delta),pos=2)
        }
        
        dev.off()
        tinytex::latexmk(paste(fig.name,".tex",sep=""))
        file.remove(paste(fig.name,".tex",sep=""))
      }
      
      # remove tau0 from c
      c <- c - tau0
      
    }
    
    # OPTIMAL FEE ESTIMATION:
    if(TRUE){
      
      months  <- seq(1,12,1)
      Y.hat   <- matrix(0,length(months),1)
      tau.hat <- matrix(0,length(months),1)
      for(month in months){
        
        Y.inputs         <- list()
        Y.inputs$N       <- N
        Y.inputs$theta   <- theta.hat.BC * (month %in% 5:10)
        Y.inputs$lambda  <- lambda.hat
        Y.inputs$sigma   <- 0
        
        Y.0       <- Ytm[468+month] # total visits in month of 2018
        gamma     <- log(Y.0/(exp(theta.hat.BC*(Y.0/10^6))*sum(N*exp(lambda.hat*(c+tau0)))))
        tau.range <- seq(0,500,1)
        Y.tau     <- matrix(0,length(tau.range),1)
        V.tau     <- matrix(0,length(tau.range),1)
        k <- 0
        for(tau in tau.range){k <- k + 1
          Y.inputs$beta    <- gamma
          Y.inputs$phi.t   <- 0
          Y.inputs$alpha.m <- 0
          Y.inputs$c       <- c + tau
          Y.tau[k]         <- Y.fn(Y.inputs)$Ytm
          V.tau[k]         <- -Y.tau[k]/lambda.hat + tau*Y.tau[k]
        }
        Y.hat[month]   <- Y.tau[which.max(V.tau)]
        tau.hat[month] <- tau.range[which.max(V.tau)]
        
        Y.inputs$c <- c + tau.hat[month]
        outs       <- Y.fn(Y.inputs)
        Yj         <- outs$Ytmj
        avg.c      <- sum((c+tau0)*Yj)/sum(Yj)
        elas       <- lambda.hat * avg.c
        
      }
      
      # Table of observed and optimal visits and surplus:
      {
        cat("\nBaseline and counterfactual visits and surplus (transform =",transform,")\n",file=out.file.name,append=TRUE)
        cat("=============================================================================================\n",file=out.file.name,append=TRUE)
        cat("        Observed                                    Optimal                                  \n",file=out.file.name,append=TRUE)
        cat("Month   tau        Y          tau Y      TS         tau        Y          tau Y      TS      \n",file=out.file.name,append=TRUE)
        cat('-----   --------   --------   --------   --------   --------   --------   --------   --------\n',file=out.file.name,append=TRUE)
        m <- 0
        for(month in c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')){
          m <- m + 1
          cat(month,sprintf('  & %8.1f & %8.3f & %8.3f & %8.3f & %8.1f & %8.3f & %8.3f & %8.3f \\\\\n',
                            tau0,
                            Ytm[468+m] / 1e6,
                            tau0 * Ytm[468+m] / 1e6,
                            Ytm[468+m]/(-lambda.hat) / 1e6 + tau0 * Ytm[468+m] / 1e6,
                            tau.hat[m],
                            Y.hat[m] / 1e6,
                            tau.hat[m] * Y.hat[m] / 1e6,
                            Y.hat[m]/(-lambda.hat) / 1e6 + tau.hat[m] * Y.hat[m] / 1e6),file=out.file.name,append=TRUE)
        }
        cat('\\hline\n',file=out.file.name,append=TRUE)
        cat('Total',sprintf('&          & %8.3f & %8.3f & %8.3f &          & %8.3f & %8.3f & %8.3f \\\\\n',
                            sum(Ytm[468+1:12]) / 1e6,
                            sum(Ytm[468+1:12] * tau0) / 1e6,
                            sum(Ytm[468+1:12]/(-lambda.hat)) / 1e6 + sum(Ytm[468+1:12] * tau0) / 1e6,
                            sum(Y.hat) / 1e6,
                            sum(Y.hat * tau.hat) / 1e6,
                            sum(Y.hat/(-lambda.hat)) / 1e6 + sum(Y.hat * tau.hat) / 1e6),file=out.file.name,append=TRUE)
        cat("=============================================================================================\n",file=out.file.name,append=TRUE)
        cat('[Price elasticity of demand =',sprintf('%-.2f.]\n\n',lambda.hat*avg.c),file=out.file.name,append=TRUE)
      }
      
      # pdf figure: optimal fees by month
      {
        fig.name <- paste(output.path,'/monthly-fees-',transform,'',sep='')
        tikz(paste(fig.name,'.tex',sep=''),
             width      = 3.5,
             height     = 3.5,
             pointsize  = 12,
             standAlone = TRUE)

        par(mar=c(3.5,3.5,1,1)) # figure border whitespace [bottom,left,top,right]
        
        bottom <- 0
        top    <- max(tau.hat)
        left   <- 0
        right  <- 12
        
        plot( 1:12,tau.hat,
              type='l',
              lty=2,
              lwd=1,
              main='',
              xlab='',
              ylab='',
              axes=FALSE,
              cex=.5,
              xlim=c(left,right),
              ylim=c(bottom,top),
              cex.axis=.8,
              cex.lab =.5)
        
        points(1:12,tau.hat)
        
        axis(1,las=1,labels=TRUE,cex.axis=.65) # Draw x axis
        axis(2,las=1,labels=TRUE,cex.axis=.65) # Draw y axis
        
        mtext('$\\hat{\\tau}^{\\star}$',side=2, line=2.7, cex=.65, las=0, srt=90) # y-axis label
        mtext('month of year'          ,side=1, line=2  , cex=.65, las=1, srt=0 ) # x-axis label
        
        dev.off()
        tinytex::latexmk(paste(fig.name,".tex",sep=""))
        file.remove(paste(fig.name,".tex",sep=""))
      }
      
    }
    
  }
  
  # TAU VS THETA:
  if(TRUE){
    theta.range <- seq(-2,0,length=100)
    tau.theta   <- matrix(0,length(theta.range),1)
    j <- 0
    for(theta in theta.range){j <- j + 1
    
    Y.inputs         <- list()
    Y.inputs$N       <- N
    Y.inputs$theta   <- theta
    Y.inputs$lambda  <- lambda.hat
    Y.inputs$sigma   <- 0
    
    Y.0       <- Ytm[469:480][7] # visits in July 2018
    gamma     <- log(Y.0/(exp(theta*(Y.0/10^6))*sum(N*exp(lambda.hat*(c+tau0)))))
    tau.range <- seq(0,250,1)
    Y.tau     <- matrix(0,length(tau.range),1)
    V.tau     <- matrix(0,length(tau.range),1)
    k <- 0
    for(tau in tau.range){k <- k + 1
      Y.inputs$beta    <- gamma
      Y.inputs$phi.t   <- 0
      Y.inputs$alpha.m <- 0
      Y.inputs$c       <- c + tau
      Y.tau[k]         <- Y.fn(Y.inputs)$Ytm
      V.tau[k]         <- Y.tau[k]/(-lambda.hat) + tau*Y.tau[k]
    }
    Y.hat        <- Y.tau[which(V.tau==max(V.tau))]
    tau.theta[j] <- tau.range[which(V.tau==max(V.tau))]
    
    }
    
    plot(theta.range,tau.theta,type='l')
    
    # pdf figure: tau.star vs theta
    {
      fig.name <- paste(output.path,'/tau-vs-theta',sep='')
      tikz(paste(fig.name,'.tex',sep=''),
           width      = 5.5,
           height     = 5.5,
           pointsize  = 12,
           standAlone = TRUE)
      
      #par(mfrow=c(1,2))
      
      par(mar=c(3.5,3.5,1,1)) # figure border whitespace [bottom,left,top,right]
      
      bottom <- 0
      top    <- 100
      left   <- min(theta.range)
      right  <- max(theta.range)
      
      plot(theta.range,tau.theta,
           type='l',
           lty=1,
           lwd=1,
           main='',
           xlab='',
           ylab='',
           axes=FALSE,
           cex=.5,
           xlim=c(left,right),
           ylim=c(bottom,top),
           cex.axis=1,
           cex.lab =1)
      
      lines(c(theta.hat.BC,theta.hat.BC),c(0,tau.theta[which.min((theta.range-theta.hat.BC)^2)]),lty=3)
      lines(c(left,theta.range[which.min((theta.range-theta.hat.BC)^2)]),c(tau.theta[which.min((theta.range-theta.hat.BC)^2)],tau.theta[which.min((theta.range-theta.hat.BC)^2)]),lty=3)
      
      axis(1,las=1,labels=TRUE,cex.axis=1) # Draw x axis
      axis(2,las=1,labels=TRUE,cex.axis=1) # Draw y axis
      
      mtext('$\\tau^{\\star}$ [\\$/visitor]',side=2, line=2.5, cex=1, las=0, srt=-90) # y-axis label
      mtext('$\\theta$'                     ,side=1, line=2  , cex=1, las=1) # x-axis label
      
      dev.off()
      tinytex::latexmk(paste(fig.name,".tex",sep=""))
      file.remove(paste(fig.name,".tex",sep=""))
    }
    
  }
  
  # CS, R, TS vs TAU:
  if(TRUE){
    
    tau.range <- seq(0,200,1)
    CS <- c()
    TS <- c()
    R  <- c()
    
    Y.inputs         <- list()
    Y.inputs$N       <- N
    Y.inputs$theta   <- theta.hat.BC
    Y.inputs$lambda  <- lambda.hat
    Y.inputs$sigma   <- 0
    
    Y.0   <- Ytm[469:480][7] # visits in July 2018
    gamma <- log(Y.0/(exp(theta.hat.BC*(Y.0/10^6))*sum(N*exp(lambda.hat*(c+tau0)))))
    Y.tau     <- matrix(0,length(tau.range),1)
    V.tau     <- matrix(0,length(tau.range),1)
    k <- 0
    for(tau in tau.range){k <- k + 1
    Y.inputs$beta    <- gamma
    Y.inputs$phi.t   <- 0
    Y.inputs$alpha.m <- 0
    Y.inputs$c       <- c+tau
    Y.tau[k]         <- Y.fn(Y.inputs)$Ytm
    TS[k]            <- -Y.tau[k]/lambda.hat + tau*Y.tau[k]
    CS[k]            <- -Y.tau[k]/lambda.hat
    R[k]             <- tau*Y.tau[k]
    }
    
    tau.TS <- tau.range[which.max(TS)]
    tau.R  <- tau.range[which.max(R)]
    
    fig.name <- paste(output.path,'/surplus-vs-tau',sep='')
    tikz(paste(fig.name,'.tex',sep=''),
         width      = 5.5,
         height     = 5.5,
         pointsize  = 12,
         standAlone = TRUE)
    
    par(mar=c(3.5,5,1,1)) # figure border whitespace [bottom,left,top,right]
    
    bottom <- min(c(CS,R,TS))*0.9
    top    <- max(c(CS,R,TS))*1.1 * 0 + 100*1e6
    left   <- min(tau.range)
    right  <- max(tau.range)
    
    plot( tau.range,CS,
          type='l',
          lty=2,
          lwd=1,
          col='black',
          main='',
          xlab='',
          ylab='',
          axes=FALSE,
          cex=1,
          xlim=c(left,right),
          ylim=c(bottom,top),
          cex.axis=1,
          cex.lab =1)
    
    lines(tau.range,R                  ,lty=6,lwd=1,col='black')
    lines(tau.range,TS                 ,lty=1,lwd=1,col='black')
    lines(c(tau.TS,tau.TS),c(0,max(TS)),lty=3,      col='black')
    lines(c(tau0,tau0),c(0,TS[which.min((tau.range-tau0)^2)]),lty=3,col='black')
    lines(c(tau.R,tau.R)  ,c(0,max(R)) ,lty=3,      col='blue')
    
    axis(1,las=1,labels=TRUE,cex.axis=1) # Draw x axis
    axis(2,las=1,labels=seq(0,100,10),at=seq(0,100e6,10e6),cex.axis=1) # Draw y axis
    
    mtext('[$10^6$ \\$]'                    ,side=2, line=2.5, cex=1, las=0) # y-axis label
    mtext('Entry fee, $\\tau$ [\\$/visitor]',side=1, line=2.5, cex=1, las=1) # x-axis label
    
    legend(left+(right-left)*.62,
           bottom+(top-bottom)*1.0,
           legend=c('Consumer surplus','Revenues','Total surplus'),
           col=c('black','black','black'),
           lty=c(2,6,1),
           cex=.85,
           lwd=1.15,
           seg.len=3.5)
    
    dev.off()
    tinytex::latexmk(paste(fig.name,".tex",sep=""))
    file.remove(paste(fig.name,".tex",sep=""))
  } 
  
  # ACTUAL, PREDICTED, AND COUNTERFACTUAL VISITS:
  if(TRUE){
    
    # pdf figure: Yseries
    {
      Yt    <- colSums(matrix(Ytm,12,length(Ytm)/12))
      Ythat <- colSums(matrix(Ytm.hat,12,length(Ytm)/12))
      Yct   <- colSums(matrix(Ytmc,12,length(Ytm)/12))
      
      fig.name <- paste(output.path,'/Yseries',sep='')
      tikz(paste(fig.name,'.tex',sep=''),
           width      = 5.5,
           height     = 3.5,
           pointsize  = 12,
           standAlone = TRUE)
      
      par(mar=c(3.5,3.5,1,1)) # figure border whitespace [bottom,left,top,right]
      
      bottom <- 0
      top    <- max(c(Yt,Yct))*1.1
      left   <- min(years)
      right  <- max(years)
      
      plot( years,Yt,
            type='p',
            lwd=1,
            main='',
            xlab='',
            ylab='',
            axes=FALSE,
            cex=.5,
            xlim=c(left,right),
            ylim=c(bottom,top),
            cex.axis=.8,
            cex.lab =.5)
      
      lines(years,Yct  ,lty=1   ,lwd=1,col='darkgray')
      lines(years,Ythat,lty=1   ,lwd=1,col='black')
      
      axis(1,at=years,las=1,labels=TRUE,cex.axis=.65) # Draw x axis
      axis(2,         las=1,labels=TRUE,cex.axis=.65) # Draw y axis
      
      mtext('Visitors',side=2, line=2.7, cex=.65, las=0) # y-axis label
      mtext('Year'    ,side=1, line=2  , cex=.65, las=1) # x-axis label
      
      legend(left+(right-left)*0,
             bottom+(top-bottom)*.95,
             legend=c('$Y$','$Y_0$','$\\hat{Y}$'),
             col=c('black','darkgray','black'),
             lty=c(NA,1,1),
             pch=c(1,NA,NA),
             cex=.65)
      
      dev.off()
      tinytex::latexmk(paste(fig.name,".tex",sep=""))
      file.remove(paste(fig.name,".tex",sep=""))
    }
    
    # pdf figure: Yratio
    {
      fig.name <- paste(output.path,'/Yratio',sep='')
      tikz(paste(fig.name,'.tex',sep=''),
           width      = 5.5,
           height     = 3.5,
           pointsize  = 12,
           standAlone = TRUE)
      
      par(mar=c(3.5,3.5,1,1)) # figure border whitespace [bottom,left,top,right]
      
      bottom <- 1
      top    <- 2.5
      left   <- min(years)
      right  <- max(years)
      
      plot( years,Yct/Yt,
            type='l',
            lwd=1,
            main='',
            xlab='',
            ylab='',
            axes=FALSE,
            cex=.5,
            xlim=c(left,right),
            ylim=c(bottom,top),
            cex.axis=1,
            cex.lab =1)
      
      axis(1,at=years,las=1,labels=TRUE,cex.axis=1) # Draw x axis
      axis(2,         las=1,labels=TRUE,cex.axis=1) # Draw y axis
      
      mtext('Counterfactual visitor ratio',side=2, line=2.5, cex=1, las=0, srt=90) # y-axis label
      mtext('Year'                        ,side=1, line=2.5, cex=1, las=1, srt=0 ) # x-axis label
      
      dev.off()
      tinytex::latexmk(paste(fig.name,".tex",sep=""))
      file.remove(paste(fig.name,".tex",sep=""))
    }
    
  }
  
}

#-------------------------------------------------------------------------------
# POLICY EXPERIMENTS:
#-------------------------------------------------------------------------------

cat("\nWorking on policy experiments... ")

# Max peak season (jun-aug) entry fee w $0 off season entry fee holding total 
# visits constant in 2018:
if(TRUE){
  
  # theta.hat.BC <- -1.0
  
  # Fix entry fees at $0 for Jan-May and Oct-Dec. Then find fixed fee for 
  # Jun-Sep that maxes surplus while keeping total visitation equal to total 
  # annual visitation in 2018.
  
  # IMPORT DATA:
  {
    data <- data.fn()
    N    <- data$N
    d    <- data$d
    Ytm  <- data$Ytm
    yj   <- data$yj
  }
  
  # Confirm tau = 0 gives 2018 observed values:
  {
    Ytm2018 <- Ytm[469:480]
    
    Y.inputs         <- list()
    Y.inputs$N       <- N
    Y.inputs$lambda  <- lambda.hat
    Y.inputs$sigma   <- 0
    
    Y.tau0 <- 0
    V.tau0 <- 0
    R.tau0 <- 0
    
    for(m in 1:12){
      Y.0   <- Ytm2018[m] 
      gamma <- log(Y.0/(exp((m>=5 & m<=10)*theta.hat.BC*(Y.0/10^6))*sum(N*exp(lambda.hat*(c+tau0)))))
      
      Y.inputs$beta    <- gamma
      Y.inputs$phi.t   <- 0
      Y.inputs$alpha.m <- 0
      Y.inputs$theta   <- (m>=5 & m<=10)*theta.hat.BC 
      Y.inputs$c       <- c+tau0
      outs             <- Y.fn(Y.inputs)
      Y.tau0           <- Y.tau0 + outs$Ytm
      R.tau0           <- R.tau0 + tau0*outs$Ytm
      V.tau0           <- V.tau0 - outs$Ytm/lambda.hat + tau0*outs$Ytm
    }
    
  }
  
  # find optimal peak fee for jun-sep:
  tau.range <- seq(-100,100,.05)
  Y.tau     <- matrix(0,length(tau.range),1)
  V.tau     <- matrix(0,length(tau.range),1)
  R.tau     <- matrix(0,length(tau.range),1)
  
  for(m in 1:12){
    
    Y.inputs         <- list()
    Y.inputs$N       <- N
    Y.inputs$lambda  <- lambda.hat
    Y.inputs$sigma   <- 0
    
    Y.0   <- Ytm[468+m] # total visits in month m of 2018
    gamma <- log(Y.0/(exp((m %in% 5:10)*theta.hat.BC*(Y.0/10^6))*sum(N*exp(lambda.hat*(c+tau0)))))
    
    k <- 0
    for(tau in tau.range){k <- k + 1
      Y.inputs$beta    <- gamma
      Y.inputs$phi.t   <- 0
      Y.inputs$alpha.m <- 0
      Y.inputs$theta   <- (m %in% 5:10)*theta.hat.BC 
      Y.inputs$c       <- c + tau * (m %in% 6:9)
      outs             <- Y.fn(Y.inputs)
      Y.tau[k]         <- Y.tau[k] + outs$Ytm
      V.tau[k]         <- V.tau[k] - outs$Ytm/lambda.hat + tau*outs$Ytm * (m %in% 6:9)
      R.tau[k]         <- R.tau[k] + tau*outs$Ytm * (m %in% 6:9)
    }
    
  }
  
  k        <- which.min((Y.tau-sum(Ytm2018))^2)
  tau.peak <- tau.range[k]
  Y.peak   <- Y.tau[k]
  V.peak   <- V.tau[k]
  R.peak   <- R.tau[k]
}

# Lower entry fees for local states (Wyoming, Montana, Idaho):
# (Set local entry to a small amount, e.g., $5---then find non-local entry fee
# to max total surplus)
if(TRUE){
  
  # Wyoming: j = 48
  # Montana: j = 24
  # Idaho:   j = 10
  
  local    <- c(10,24,48)
  nonlocal <- c(1:9,11:23,25:47,49)
  
  # Local entry fee [$/visitor]:
  tau.L <- 5
  
  # IMPORT DATA:
  {
    data <- data.fn()
    N    <- data$N
    d    <- data$d
    Ytm  <- data$Ytm
    yj   <- data$yj
  }
  
  # find optimal peak fee for non-local visitors:
  tauN.range <- seq(0,200,.1)
  Y.tauN     <- matrix(0,length(tauN.range),1)
  V.tauN     <- matrix(0,length(tauN.range),1)
  R.tauN     <- matrix(0,length(tauN.range),1)
  
  for(m in 1:12){
    
    Y.inputs         <- list()
    Y.inputs$N       <- N
    Y.inputs$lambda  <- lambda.hat
    Y.inputs$sigma   <- 0
    
    Y.0   <- Ytm[468+m] # total visits in month m of 2018
    gamma <- log(Y.0/(exp((m %in% 5:10)*theta.hat.BC*(Y.0/10^6))*sum(N*exp(lambda.hat*(c+tau0)))))
    
    k <- 0
    for(tauN in tauN.range){k <- k + 1
      Y.inputs$beta    <- gamma
      Y.inputs$phi.t   <- 0
      Y.inputs$alpha.m <- 0
      Y.inputs$theta   <- (m %in% 5:10)*theta.hat.BC
      c.m <- c; c.m[nonlocal] <- c.m[nonlocal] + tauN
      Y.inputs$c       <- c.m
      outs             <- Y.fn(Y.inputs)
      Y.tauN[k]        <- Y.tauN[k] + outs$Ytm
      V.tauN[k]        <- V.tauN[k] - outs$Ytm/lambda.hat + tauN*sum(outs$Ytmj[nonlocal])
      R.tauN[k]        <- R.tauN[k] + tauN*sum(outs$Ytmj[nonlocal])
    }
    
  }
  
  tau.localV <- tauN.range[which.max(V.tauN)]
  Y.localV   <- Y.tauN[which.max(V.tauN)]
  V.localV   <- V.tauN[which.max(V.tauN)]
  R.localV   <- R.tauN[which.max(V.tauN)]
  tau.localR <- tauN.range[which.max(R.tauN)]
  Y.localR   <- Y.tauN[which.max(R.tauN)]
  V.localR   <- V.tauN[which.max(R.tauN)]
  R.localR   <- R.tauN[which.max(R.tauN)]
  
}

# Table of results:
{
  cat('\nPolicy experiments summary statistics:\n',file=out.file.name,append=TRUE)
  cat('==============================================================================================\n',file=out.file.name,append=TRUE)
  cat('Policy                           Fee ($/visitor)  Visitors (M)    Surplus (M$)   Revenues (M$)\n',file=out.file.name,append=TRUE)
  cat('----------------------------------------------------------------------------------------------\n',file=out.file.name,append=TRUE)
  cat(sprintf('Status quo                     & %13.2f & %13.3f & %13.1f & %13.1f\\\\\n',tau0,sum(Ytm2018)/1e6,V.tau0/1e6  ,R.tau0/1e6  ),file=out.file.name,append=TRUE)
  cat(sprintf('Peak pricing (constant visits) & %13.2f & %13.3f & %13.1f & %13.1f\\\\\n',tau.peak,Y.peak/1e6    ,V.peak/1e6  ,R.peak/1e6  ),file=out.file.name,append=TRUE)
  cat(sprintf('Local discount (max surplus)   & %13.2f & %13.3f & %13.1f & %13.1f\\\\\n',tau.localV,Y.localV/1e6  ,V.localV/1e6,R.localV/1e6),file=out.file.name,append=TRUE)
  cat(sprintf('Local discount (max revenue)   & %13.2f & %13.3f & %13.1f & %13.1f\\\\\n',tau.localR,Y.localR/1e6  ,V.localR/1e6,R.localR/1e6),file=out.file.name,append=TRUE)
  cat('==============================================================================================\n',file=out.file.name,append=TRUE)
}

cat("Done.\n")

# WRITE SOURCE FILE TO OUTPUT FILE:
if(TRUE){
  source.file.name <- paste(code.path,'/',script.name,'.R',sep='')
  
  # Read lines of source file:
  Rscript <- readLines(source.file.name)
  
  # Write lines of source file to output file:  
  cat('\n\n|---------------------------------------------------------------------------|',file=out.file.name,append=TRUE)
  cat('\n| R OUTPUT ABOVE                                                            |',file=out.file.name,append=TRUE)
  cat('\n|---------------------------------------------------------------------------|',file=out.file.name,append=TRUE)
  cat('\n| R SCRIPT BELOW                                                            |',file=out.file.name,append=TRUE)
  cat('\n|---------------------------------------------------------------------------|\n',file=out.file.name,append=TRUE)
  for(i in 1:length(Rscript)){cat('\n',Rscript[i],file=out.file.name,append=TRUE)}
  
}

