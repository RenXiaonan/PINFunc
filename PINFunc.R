#############################################################
#####                                                   #####
#####                                                   #####
#####                    Data Generation                #####
#####                                                   #####
#####                                                   #####
#############################################################

################## Generate Data (single data) ##############
### There are four data distributions for four Cases.     ###

####### Case I - Linear under Homogeneity model  #######
# Named as DataHoLinear.

DataHoLinear <- function(n1, n2, n3, b1, b2, b3, p, corval){
  # n1: sample size of dataset 1
  # n2: sample size of dataset 2
  # n3: sample size of dataset 3
  # b1: model coefficients of dataset 1
  # b2: model coefficients of dataset 1
  # b3: model coefficients of dataset 1
  # p : dimensions
  # corval: correlation coefficient of power structure
  # 
  # Output: 
  # This function will result three datasets with size n1, n2,and n3 respectively.
  mean=matrix(c(rep(0,p)),nrow=p)
  sigma=diag(p)
  for(sid in 1:p){
    for(sidd in 1:p){
      sigma[sid,sidd]=corval^abs(sid-sidd) # Power structure
    }
  }
  
  x1<-mvrnorm(n1,mean,sigma)
  x2<-mvrnorm(n2,mean,sigma)
  x3<-mvrnorm(n3,mean,sigma)
  
  SNR=3
  
  s1=x1%*%b1
  s2=x2%*%b2
  s3=x3%*%b3
  
  r1<-rnorm(n1,0,1)
  r2<-rnorm(n2,0,1)
  r3<-rnorm(n3,0,1)
  
  k1 = sqrt(var(s1)/(SNR*var(r1)))
  k2 = sqrt(var(s2)/(SNR*var(r2)))
  k3 = sqrt(var(s3)/(SNR*var(r3)))
  
  y1<- s1 +r1*c(k1)
  y2<- s2 +r2*c(k2)
  y3<- s3 +r3*c(k3)
  
  data <- list(x1=x1, x2=x2, x3=x3, y1=y1, y2=y2, y3=y3)
  return(data)		
}
###############################################################


####### Case II - Non-linear under Homogeneity model  #######
# Named as DataHoNonLinear.

DataHoNonLinear <- function(n1, n2, n3, b1, b2, b3, p, corval){
  # The parameter details are the same as Case I
  mean=matrix(c(rep(0,p)),nrow=p)
  sigma=diag(p)
  for(sid in 1:p){
    for(sidd in 1:p){
      sigma[sid,sidd]=corval^abs(sid-sidd)
    }
  }
  
  x1<-mvrnorm(n1,mean,sigma)
  x2<-mvrnorm(n2,mean,sigma)
  x3<-mvrnorm(n3,mean,sigma)
  
  x1[,3]=pnorm(x1[,3])
  x1[,4]=pnorm(x1[,4])
  x1[,5]=pnorm(x1[,5])
  
  x2[,3]=pnorm(x2[,3])
  x2[,4]=pnorm(x2[,4])
  x2[,5]=pnorm(x2[,5])
  
  x3[,3]=pnorm(x3[,3])
  x3[,4]=pnorm(x3[,4])
  x3[,5]=pnorm(x3[,5])
  
  SNR=3
  
  s1=x1[,-c(3,4,5)]%*%b1[-c(3,4,5),]+sin(4*pi*x1[,4])*b1[4,]+sin(4*pi*x1[,3])*b1[3,]+exp(2.5*x1[,5])*b1[5,]
  s2=x2[,-c(3,4,5)]%*%b2[-c(3,4,5),]+sin(4*pi*x2[,4])*b2[4,]+sin(4*pi*x2[,3])*b2[3,]+exp(2.5*x2[,5])*b2[5,]
  s3=x3[,-c(3,4,5)]%*%b3[-c(3,4,5),]+sin(4*pi*x3[,4])*b3[4,]+sin(4*pi*x3[,3])*b3[3,]+exp(2.5*x3[,5])*b3[5,]
  
  r1<-rnorm(n1,0,1)
  r2<-rnorm(n2,0,1)
  r3<-rnorm(n3,0,1)
  
  k1 = sqrt(var(s1)/(SNR*var(r1)))
  k2 = sqrt(var(s2)/(SNR*var(r2)))
  k3 = sqrt(var(s3)/(SNR*var(r3)))
  
  y1<- s1 +r1*c(k1)
  y2<- s2 +r2*c(k2)
  y3<- s3 +r3*c(k3)
  
  data <- list(x1=x1, x2=x2, x3=x3, y1=y1, y2=y2, y3=y3)
  return(data)		
}
###############################################################


####### Case III - linear under heterogeneity models  #######
# Named as DataHeLinear.

DataHeLinear <- function(n1, n2, n3, b1, b2, b3, p, corval){
  # The parameter details are the same as Case I
  mean=matrix(c(rep(0,p)),nrow=p)
  sigma=diag(p)
  
  for(sid in 1:p){
    for(sidd in 1:p){
      sigma[sid,sidd]=corval^abs(sid-sidd)
    }
  }
  
  x1<-mvrnorm(n1,mean,sigma)
  x2<-mvrnorm(n2,mean,sigma)
  x3<-mvrnorm(n3,mean,sigma)
  
  SNR=3
  
  s1=x1%*%b1
  s2=x2%*%b2
  s3=x3%*%b3
  
  r1<-rnorm(n1,0,1)
  r2<-rnorm(n2,0,1)
  r3<-rnorm(n3,0,1)
  
  k1 = sqrt(var(s1)/(SNR*var(r1)))
  k2 = sqrt(var(s2)/(SNR*var(r2)))
  k3 = sqrt(var(s3)/(SNR*var(r3)))
  
  y1<- s1 +r1*c(k1)
  y2<- s2 +r2*c(k2)
  y3<- s3 +r3*c(k3)
  
  data <- list(x1=x1, x2=x2, x3=x3, y1=y1, y2=y2, y3=y3)
  return(data)		
}
###############################################################


####### Case IV - Nonlinear under heterogeneity models  #######
# Named as DataHeNonLinear.

DataHeNonLinear <- function(n1, n2, n3, b1, b2, b3, p, corval){
  # The parameter details are the same as Case I
  mean=matrix(c(rep(0,p)),nrow=p)
  sigma=diag(p)
  for(sid in 1:p){
    for(sidd in 1:p){
      sigma[sid,sidd]=corval^abs(sid-sidd)
    }
  }
  
  x1<-mvrnorm(n1,mean,sigma)
  x2<-mvrnorm(n2,mean,sigma)
  x3<-mvrnorm(n3,mean,sigma)
  
  SNR=3
  
  x1[,3]=pnorm(x1[,3])
  x1[,4]=pnorm(x1[,4])
  x1[,5]=pnorm(x1[,5])
  x1[,6]=pnorm(x1[,6])
  
  x2[,3]=pnorm(x2[,3])
  x2[,4]=pnorm(x2[,4])
  x2[,5]=pnorm(x2[,5])
  x2[,7]=pnorm(x2[,7])
  
  x3[,3]=pnorm(x3[,3])
  x3[,4]=pnorm(x3[,4])
  x3[,5]=pnorm(x3[,5])
  x3[,8]=pnorm(x3[,8])
  
  s1=x1[,-c(3,4,5,6)]%*%b1[-c(3,4,5,6),]+sin(4*pi*x1[,3])*b1[3,]+sin(4*pi*x1[,4])*b1[4,]+exp(2.5*x1[,5])*b1[5,]+exp(2.5*x1[,6])*b1[6,]
  s2=x2[,-c(3,4,5,7)]%*%b2[-c(3,4,5,7),]+sin(4*pi*x2[,3])*b2[3,]+sin(4*pi*x2[,4])*b2[4,]+exp(2.5*x2[,5])*b2[5,]+exp(2.5*x2[,7])*b2[7,]
  s3=x3[,-c(3,4,5,8)]%*%b3[-c(3,4,5,8),]+sin(4*pi*x3[,3])*b3[3,]+sin(4*pi*x3[,4])*b3[4,]+exp(2.5*x3[,5])*b3[5,]+exp(2.5*x3[,8])*b3[8,]
  
  r1<-rnorm(n1,0,1)
  r2<-rnorm(n2,0,1)
  r3<-rnorm(n3,0,1)
  
  k1 = sqrt(var(s1)/(SNR*var(r1)))
  k2 = sqrt(var(s2)/(SNR*var(r2)))
  k3 = sqrt(var(s3)/(SNR*var(r3)))
  
  y1<- s1 +r1*c(k1)
  y2<- s2 +r2*c(k2)
  y3<- s3 +r3*c(k3)
  
  data <- list(x1=x1, x2=x2, x3=x3, y1=y1, y2=y2, y3=y3)
  return(data)		
}
###############################################################




#############################################################
#####                                                   #####
#####                                                   #####
#####                    Data Estimation                #####
#####                                                   #####
#####                                                   #####
#############################################################


##################### Data division #####################
DataDiv <- function(data, divrate){
  # Divide each dataset into training set and test set. The training set is used to train a PIN network, and the test set is used to select important variables for new data.
  x1=data$x1
  x2=data$x2
  x3=data$x3
  
  y1=data$y1
  y2=data$y2
  y3=data$y3
  
  nn1=trainrate*n1
  nn2=trainrate*n2
  nn3=trainrate*n3
  
  train_id_1=sample(1:n1,nn1)
  test_id_1=setdiff(1:n1,train_id_1)
  
  train_id_2=sample(1:n2,nn2)
  test_id_2=setdiff(1:n2,train_id_2)
  
  train_id_3=sample(1:n3,nn3)
  test_id_3=setdiff(1:n3,train_id_3)
  
  x_train_yu=list(x1[train_id_1,],x2[train_id_2,],x3[train_id_3,])
  x_test_yu=list(x1[test_id_1,],x2[test_id_2,],x3[test_id_3,])
  
  y_train_yu=list(y1[train_id_1,],y2[train_id_2,],y3[train_id_3,])
  y_test_yu=list(y1[test_id_1,],y2[test_id_2,],y3[test_id_3,])
  
  x=list(x1,x2,x3)
  y=list(y1,y2,y3)
  
  dataDiv <- list(x_train_yu, x_test_yu, y_train_yu, y_test_yu, x, y)
  return(dataDiv)		
}
###############################################################


################## Data Preprocessing #################
DataC <- function(dataDiv){
  # Preprocess the X and Y in the training set and test set of each data set separately.
  
  x_train_yu=dataDiv[[1]]
  x_test_yu=dataDiv[[2]]
  y_train_yu=dataDiv[[3]]
  y_test_yu=dataDiv[[4]]
  x=dataDiv[[5]]
  y=dataDiv[[6]]
  
  # For Y
  y_train=list()
  for(i in 1:length(y_train_yu)){
    y_train[[i]]<-(y_train_yu[[i]]-min(y_train_yu[[i]]))/(max(y_train_yu[[i]])-min(y_train_yu[[i]]))
  }
  y_test=list()
  for(i in 1:length(y_test_yu)){
    y_test[[i]]<-(y_test_yu[[i]]-min(y_test_yu[[i]]))/(max(y_test_yu[[i]])-min(y_test_yu[[i]]))
  }
  
  # For X
  x_train=list()
  for(i in 1:length(x_train_yu)){
    x_train[[i]]<-scale(x_train_yu[[i]],center = T,scale=T)
  }
  x_test=list()
  for(i in 1:length(x_test_yu)){
    x_test[[i]]<-scale(x_test_yu[[i]],center = T,scale=T)
  }
  
  # Return
  dataC <- list(x_train, y_train, x_test, y_test)
  return(dataC)
}
###############################################################


############ Derivation of activation function 'tanh' #########
dtanh=function(x){
  y=1-tanh(x)^2
  return(y)
}
###############################################################




##################### Algorithm - HoPIN #######################
####### This function will give the selected feature   ########
hopin=function(x,y,m,varepsilon,maxiter,alpha,eta,lambda,Lambda,a,threshold,node){
  ##### Input: #####
  # x: a list contain all datasets' X
  #	y: a list contain all datasets' response Y
  #	m: number of layers
  # node: number of hidden nodes
  # varepsilon: condition of convergence
  #	maxiter: maximum number of iterations
  # alpha: momentum parameter
  # eta: learning rate
  # lambda: tuning parameter for variable selection in the objective function
  # Lambda: tuning parameter for preventing overfitting in the objective function
  #	a: smoothing parameter for the smoothing function in the objective function
  #	threshold: variable selection threshold
  #
  ##### Output: #####
  # w, W, b: the values of w, W, and b in the fitted network
  #	index: selected important variables for each dataset (a vector)
  #	SSE: total SSE between fitted values and true values
  #	tt: number of final iterations
    
    p=dim(x[[1]])[2] #p is the number of variables
    q=node #q is the number of hidden nodes
    d=length(x) #d is the number of datasets
    tt=0 #t is the record of iteration
    
    w=list()
    for(i in 1:d){
      w[[i]]=diag(p)
    } # w is initialization weight matrix list for input
    
    Vw=list()
    for(i in 1:d){
      Vw[[i]]=matrix(0,nrow=p,ncol=p)
    } # Vw is initialization 1st moment matrix list for w
    
    dw=list()
    for(i in 1:d){
      dw[[i]]=matrix(0,nrow=p,ncol=p)
    } # dw initialization
    
    diffw=list()
    for(i in 1:d){
      diffw[[i]]=matrix(0,nrow=p,ncol=p)
    } # w difference initialization diffw  
    
    if(m==1){
      
      W=list()
      for(i in 1:d){
        W[[i]]=list()
        W[[i]][[1]]=matrix(runif(q*p,min = - sqrt(6/(q+p)), max = sqrt(6/(q+p)) ),nrow=q)
      }
      
      for(i in 1:d){
        W[[i]][[m+1]]=matrix(runif(q,min = - sqrt(6/(q+1)), max = sqrt(6/(q+1)) ),nrow=1)
      } # W is initialization weight matrix list
      
      b=list()
      for(i in 1:d){
        b[[i]]=list()
        for(l in 1:m){
          b[[i]][[l]]=matrix(0,nrow=q)
        }
      }
      
      for(i in 1:d){
        b[[i]][[m+1]]=matrix(0,nrow=1)
      } # b is initialization bias matrix list
      
      O=list()
      for(i in 1:d){
        O[[i]]=matrix(rep(0,dim(x[[i]])[1]),ncol = 1)
      }
      
      for(i in 1:d){
        for(j in 1:(dim(x[[i]])[1])){
          l=1
          h=t( tanh( W[[i]][[l]] %*% w[[i]] %*% t(t(x[[i]][j,])) + b[[i]][[l]] ) )
          O[[i]][j]=tanh( W[[i]][[m+1]] %*% t(h) + b[[i]][[m+1]] )
          
        }
      } # compute Output O
      
      dO=list()
      for(i in 1:d){
        dO[[i]]= matrix(rep(0,dim(x[[i]])[1]),ncol = 1)
      }
      for(i in 1:d){
        for(j in 1:(dim(x[[i]])[1])){
          dO[[i]][j]=dtanh( W[[i]][[m+1]] %*% t(h) + b[[i]][[m+1]] )
        }
      } # dO initialization
      
      VW=list()
      for(i in 1:d){
        VW[[i]]=list()
        VW[[i]][[1]]=matrix(rep(0,q*p),nrow=q)
      }
      
      for(i in 1:d){
        VW[[i]][[m+1]]=matrix(rep(0,q),nrow=1)
      } # VW is initialization 1st moment matrix list for W
      
      Vb=list()
      for(i in 1:d){
        Vb[[i]]=list()
        for(l in 1:m){
          Vb[[i]][[l]]=matrix(0,nrow=q)
        }
      }
      
      for(i in 1:d){
        Vb[[i]][[m+1]]=matrix(0,nrow=1)
      } # Vb is initialization 1st moment matrix list for b
      
      H=list()
      for(i in 1:d){
        H[[i]]=list()
      }
      for(i in 1:d){
        for(l in 1:m){
          H[[i]][[l]]=matrix( rep( 0 , dim(x[[i]])[1] * q ),ncol=q)
        }
      } #H initialization
      
      dH=list()
      for(i in 1:d){
        dH[[i]]=list()
      }
      for(i in 1:d){
        for(l in 1:m){
          dH[[i]][[l]]=matrix( rep( 0 , dim(x[[i]])[1] * q ),ncol=q)
        }
      } # dH initialization
      
      C=diag(p) # C initialization
      
      dW=list()
      for(i in 1:d){
        dW[[i]]=list()
        dW[[i]][[1]]=matrix( rep(0,q*p) ,nrow=q)
      }
      
      for(i in 1:d){
        dW[[i]][[m+1]]=matrix( rep(0,q) ,nrow=1)
      } # dW initialization
      
      db=list()
      for(i in 1:d){
        db[[i]]=list()
        for(l in 1:m){
          db[[i]][[l]]=matrix(0,nrow=q)
        }
      }
      
      for(i in 1:d){
        db[[i]][[m+1]]=matrix(0,nrow=1)
      } # db initialization
      
      diffW=list()
      for(i in 1:d){
        diffW[[i]]=list()
        diffW[[i]][[1]]=matrix( rep(0,q*p) ,nrow=q)
      }
      
      for(i in 1:d){
        diffW[[i]][[m+1]]=matrix( rep(0,q) ,nrow=1)
      } #W difference initialization diffW
    }else{
      W=list()
      for(i in 1:d){
        W[[i]]=list()
        W[[i]][[1]]=matrix(runif(q*p,min = - sqrt(6/(q+p)), max = sqrt(6/(q+p)) ),nrow=q)
      }
      
      for(i in 1:d){
        for(l in 2:m){
          W[[i]][[l]]=matrix(runif(q*q,min = - sqrt(6/(q+q)), max = sqrt(6/(q+q)) ),nrow=q)
        }
      }
      
      for(i in 1:d){
        W[[i]][[m+1]]=matrix(runif(q,min = - sqrt(6/(q+1)), max = sqrt(6/(q+1)) ),nrow=1)
      } #W is initialization weight matrix list
      
      b=list()
      for(i in 1:d){
        b[[i]]=list()
        for(l in 1:m){
          b[[i]][[l]]=matrix(0,nrow=q)
        }
      }
      
      for(i in 1:d){
        b[[i]][[m+1]]=matrix(0,nrow=1)
      } #b is initialization bias matrix list
      
      O=list()
      for(i in 1:d){
        O[[i]]=matrix(rep(0,dim(x[[i]])[1]),ncol = 1)
      }
      
      for(i in 1:d){
        for(j in 1:(dim(x[[i]])[1])){
          
          l=1
          h=t( tanh( W[[i]][[l]] %*% w[[i]] %*% t(t(x[[i]][j,])) + b[[i]][[l]] ) )
          repeat{
            l=l+1
            h=t (tanh( W[[i]][[l]] %*% t(h) + b[[i]][[l]] ) )
            if(l == m)
              break
          }
          O[[i]][j]=tanh( W[[i]][[m+1]] %*% t(h) + b[[i]][[m+1]] )
          
        }
      } #compute Output O
      
      dO=list()
      for(i in 1:d){
        dO[[i]]= matrix(rep(0,dim(x[[i]])[1]),ncol = 1)
      }
      for(i in 1:d){
        for(j in 1:(dim(x[[i]])[1])){
          l=1
          h=t( tanh( W[[i]][[l]] %*% w[[i]] %*% t(t(x[[i]][j,])) + b[[i]][[l]] ) )
          repeat{
            l=l+1
            h=t (tanh( W[[i]][[l]] %*% t(h) + b[[i]][[l]] ) )
            if(l == m)
              break
          }
          dO[[i]][j]=dtanh( W[[i]][[m+1]] %*% t(h) + b[[i]][[m+1]] )        
        }
      } #dO initialization
      
      VW=list()
      for(i in 1:d){
        VW[[i]]=list()
        VW[[i]][[1]]=matrix(rep(0,q*p),nrow=q)
      }
      
      for(i in 1:d){
        for(l in 2:m){
          VW[[i]][[l]]=matrix(rep(0,q*q),nrow=q)
        }
      }
      
      for(i in 1:d){
        VW[[i]][[m+1]]=matrix(rep(0,q),nrow=1)
      } #VW is initialization 1st moment matrix list for W
      
      Vb=list()
      for(i in 1:d){
        Vb[[i]]=list()
        for(l in 1:m){
          Vb[[i]][[l]]=matrix(0,nrow=q)
        }
      }
      
      for(i in 1:d){
        Vb[[i]][[m+1]]=matrix(0,nrow=1)
      } #Vb is initialization 1st moment matrix list for b
    
      H=list()
      for(i in 1:d){
        H[[i]]=list()
      }
      for(i in 1:d){
        for(l in 1:m){
          H[[i]][[l]]=matrix( rep( 0 , dim(x[[i]])[1] * q ),ncol=q)
        }
      } #H initialization
      
      dH=list()
      for(i in 1:d){
        dH[[i]]=list()
      }
      for(i in 1:d){
        for(l in 1:m){
          dH[[i]][[l]]=matrix( rep( 0 , dim(x[[i]])[1] * q ),ncol=q)
        }
      } #dH initialization
      
      C=diag(p) #C initialization
      
      dW=list()
      for(i in 1:d){
        dW[[i]]=list()
        dW[[i]][[1]]=matrix( rep(0,q*p) ,nrow=q)
      }
      
      for(i in 1:d){
        for(l in 2:m){
          dW[[i]][[l]]=matrix( rep(0,q*q) ,nrow=q)
        }
      }
      
      for(i in 1:d){
        dW[[i]][[m+1]]=matrix( rep(0,q) ,nrow=1)
      } #dW initialization
      
      db=list()
      for(i in 1:d){
        db[[i]]=list()
        for(l in 1:m){
          db[[i]][[l]]=matrix(0,nrow=q)
        }
      }
      
      for(i in 1:d){
        db[[i]][[m+1]]=matrix(0,nrow=1)
      } #db initialization
      
      diffW=list()
      for(i in 1:d){
        diffW[[i]]=list()
        diffW[[i]][[1]]=matrix( rep(0,q*p) ,nrow=q)
      }
      
      for(i in 1:d){
        for(l in 2:m){
          diffW[[i]][[l]]=matrix( rep(0,q*q) ,nrow=q)
        }
      }
      
      for(i in 1:d){
        diffW[[i]][[m+1]]=matrix( rep(0,q) ,nrow=1)
      } #W difference initialization diffW
    }
    
    
    plot(1:maxiter, seq(0,0.01,length.out = maxiter), xlim = c(1,maxiter), ylim = c (0,0.01), type = "n", axes = F, xlab = "iteration times", ylab="maximum difference")
    box()
    axis(side = 1, at = c(seq(1,maxiter+1,length.out=51)))
    axis(side = 2, at = c(seq(0,0.01,length.out = 21)))
    abline(h=0.001,col="red")
    #generate a null plot
    
    repeat{
      w_C<-w
      for(k in 1:p){
        sum_w=0
        for(i in 1:d){
          sum_w = sum_w + (w_C[[i]][k,k])^2
        }
        if(sum_w <= d * 10^(-20)){
          for(i in 1:d){
            w_C[[i]][k,k]=10^(-10)
          }
        }
      }#prevent small value 
      
      if(m==1){
        for(i in 1:d){
          for(j in 1:(dim(x[[i]])[1])){
            H[[i]][[1]][j,]= t( tanh( W[[i]][[1]] %*% w[[i]] %*% t(t(x[[i]][j,])) + b[[i]][[1]] ) ) 
          }
        }
      }else{
        for(i in 1:d){
          for(j in 1:(dim(x[[i]])[1])){
            H[[i]][[1]][j,]= t( tanh( W[[i]][[1]] %*% w[[i]] %*% t(t(x[[i]][j,])) + b[[i]][[1]] ) ) 
            
          }
        }
        for(i in 1:d){
          for(l in 2:m){
            for(j in 1:(dim(x[[i]])[1])){
              H[[i]][[l]][j,] = t( tanh( W[[i]][[l]] %*% t(t(H[[i]][[l-1]][j,])) + b[[i]][[l]] ) )
            }
          }
        }#compute H
      }
      
      if(m==1){
        for(i in 1:d){
          for(j in 1:(dim(x[[i]])[1])){
            dH[[i]][[1]][j,]= t( dtanh( W[[i]][[1]] %*% w[[i]] %*% t(t(x[[i]][j,])) + b[[i]][[1]] ) ) 
          }
        }
      }else{
        for(i in 1:d){
          for(j in 1:(dim(x[[i]])[1])){
            dH[[i]][[1]][j,]= t( dtanh( W[[i]][[1]] %*% w[[i]] %*% t(t(x[[i]][j,])) + b[[i]][[1]] ) ) 
          }
        }
      
        for(i in 1:d){
          for(l in 2:m){
            for(j in 1:(dim(x[[i]])[1])){
              dH[[i]][[l]][j,] = t( dtanh( W[[i]][[l]] %*% t(t(H[[i]][[l-1]][j,])) + b[[i]][[l]] ) )
            }
          }
        }#compute dH
      }
      
      
      for(k in 1:p){
        hat_w=0
        for(i in 1:d){
          hat_w = hat_w + (w_C[[i]][k,k])^2 
        }
        hat_w=sqrt(hat_w)
        if(hat_w >= a){
          C[k,k]=1/sqrt(hat_w)
        }else{
          C[k,k]=(hat_w)*(-1/(2*a^2)) + 3/(2*a)
        }
      }#compute C for t iteration by smoothing function
      
      tt=tt+1
      #iteration count
      
      if(m==1){
        
        for(i in 1:d){
          sum_dw = 0
          for(j in 1:(dim(x[[i]])[1])){
            sum_dw_j = c(dO[[i]][j]) * W[[i]][[2]]
            sum_dw_j = sum_dw_j %*% diag( dH[[i]][[1]][j,] ) %*% W[[i]][[1]]
            sum_dw_j = t(t( x[[i]][j,] )) %*% sum_dw_j  
            sum_dw_j = t( sum_dw_j ) * (O[[i]][j]-y[[i]][j])
            sum_dw = sum_dw + sum_dw_j
          }
          dw[[i]] = sum_dw / (dim(x[[i]])[1]) + lambda * w[[i]] * C 
        }#w
        
        for(i in 1:d){
          sum_dW = 0
          for(j in 1:(dim(x[[i]])[1])){
            sum_dW_j = c(dO[[i]][j]) * W[[i]][[2]]
            sum_dW_j = sum_dW_j %*% diag( dH[[i]][[1]][j,] )
            sum_dW_j = w[[i]] %*% t(t( x[[i]][j,] )) %*% sum_dW_j  
            sum_dW_j = t( sum_dW_j ) * (O[[i]][j]-y[[i]][j])
            sum_dW = sum_dW + sum_dW_j
          }
          dW[[i]][[1]] = sum_dW / (dim(x[[i]])[1]) + Lambda * W[[i]][[1]]   
        }#1
        
        for(i in 1:d){
          sum_dW = 0
          for(j in 1:(dim(x[[i]])[1])){
            sum_dW_j = t(H[[i]][[m]][j,]) * c(dO[[i]][j]) * (O[[i]][j]-y[[i]][j])  
            sum_dW = sum_dW + sum_dW_j
          }
          dW[[i]][[m+1]] = sum_dW / (dim(x[[i]])[1]) + Lambda * W[[i]][[m+1]] 
        }#m+1 #update dW 1:m+1
        
      }else{
       
        for(i in 1:d){
          sum_dw = 0
          for(j in 1:(dim(x[[i]])[1])){
            sum_dw_j = c(dO[[i]][j]) * W[[i]][[m+1]]
            for(l in m:2){
              sum_dw_j = sum_dw_j %*% diag( dH[[i]][[l]][j,] ) %*% W[[i]][[l]]
            }
            sum_dw_j = sum_dw_j %*% diag( dH[[i]][[1]][j,] ) %*% W[[i]][[1]]
            sum_dw_j = t(t( x[[i]][j,] )) %*% sum_dw_j  
            sum_dw_j = t( sum_dw_j ) * (O[[i]][j]-y[[i]][j])
            sum_dw = sum_dw + sum_dw_j
          }
          dw[[i]] = sum_dw / (dim(x[[i]])[1]) + lambda * w[[i]] * C 
        } #w
        
        for(i in 1:d){
          sum_dW = 0
          for(j in 1:(dim(x[[i]])[1])){
            sum_dW_j = c(dO[[i]][j]) * W[[i]][[m+1]]
            for(l in m:2){
              sum_dW_j = sum_dW_j %*% diag( dH[[i]][[l]][j,] ) %*% W[[i]][[l]]
            }
            sum_dW_j = sum_dW_j %*% diag( dH[[i]][[1]][j,] )
            sum_dW_j = w[[i]] %*% t(t( x[[i]][j,] )) %*% sum_dW_j  
            sum_dW_j = t( sum_dW_j ) * (O[[i]][j]-y[[i]][j])
            sum_dW = sum_dW + sum_dW_j
          }
          dW[[i]][[1]] = sum_dW / (dim(x[[i]])[1]) + Lambda * W[[i]][[1]] 
        } #1
        
        for(i in 1:d){
          
          if(m==2){
            sum_dW=0
            for(j in 1:(dim(x[[i]])[1])){
              sum_dW_j = c(dO[[i]][j]) * W[[i]][[m+1]]
              sum_dW_j = sum_dW_j %*% diag( dH[[i]][[m]][j,] )
              sum_dW_j = t(t( H[[i]][[m-1]][j,] )) %*% sum_dW_j  
              sum_dW_j = t( sum_dW_j ) * (O[[i]][j]-y[[i]][j])
              sum_dW = sum_dW + sum_dW_j
            }
            dW[[i]][[m]] = sum_dW / (dim(x[[i]])[1]) + Lambda * W[[i]][[m]]  #m
          }else{
            for(l in 2:(m-1)){
              sum_dW = 0
              for(j in 1:(dim(x[[i]])[1])){
                sum_dW_j = c(dO[[i]][j]) * W[[i]][[m+1]]
                for(ll in m:(l+1)){
                  sum_dW_j = sum_dW_j %*% diag( dH[[i]][[ll]][j,] ) %*% W[[i]][[ll]]
                }
                sum_dW_j = sum_dW_j %*% diag( dH[[i]][[l]][j,] )
                sum_dW_j = t(t( H[[i]][[l-1]][j,] )) %*% sum_dW_j  
                sum_dW_j = t( sum_dW_j ) * (O[[i]][j]-y[[i]][j])
                sum_dW = sum_dW + sum_dW_j
              }
              dW[[i]][[l]] = sum_dW / (dim(x[[i]])[1]) + Lambda * W[[i]][[l]]  
            }
          }#2:m-1 
          sum_dW=0
          for(j in 1:(dim(x[[i]])[1])){
            sum_dW_j = c(dO[[i]][j]) * W[[i]][[m+1]]
            sum_dW_j = sum_dW_j %*% diag( dH[[i]][[m]][j,] )
            sum_dW_j = t(t( H[[i]][[m-1]][j,] )) %*% sum_dW_j  
            sum_dW_j = t( sum_dW_j ) * (O[[i]][j]-y[[i]][j])
            sum_dW = sum_dW + sum_dW_j
          }
          dW[[i]][[m]] = sum_dW / (dim(x[[i]])[1]) + Lambda * W[[i]][[m]]  
          #m
        }#2:m
        
        for(i in 1:d){
          sum_dW = 0
          for(j in 1:(dim(x[[i]])[1])){
            sum_dW_j = t(H[[i]][[m]][j,]) * c(dO[[i]][j]) * (O[[i]][j]-y[[i]][j])  
            sum_dW = sum_dW + sum_dW_j
          }
          dW[[i]][[m+1]] = sum_dW / (dim(x[[i]])[1]) + Lambda * W[[i]][[m+1]] 
        }#m+1
        #update dW 1:m+1
      }
      
      if(m==1){
        for(i in 1:d){
          sum_db = 0
          for(j in 1:(dim(x[[i]])[1])){
            sum_db_j = c(dO[[i]][j]) * W[[i]][[2]]
            sum_db_j = sum_db_j %*% diag( dH[[i]][[1]][j,] )
            sum_db_j = t( sum_db_j ) * (O[[i]][j]-y[[i]][j])
            sum_db = sum_db + sum_db_j
          }
          db[[i]][[1]] = sum_db / (dim(x[[i]])[1])  
        } #1
        
        for(i in 1:d){
          sum_db = 0
          for(j in 1:(dim(x[[i]])[1])){
            sum_db_j = c(dO[[i]][j]) * (O[[i]][j]-y[[i]][j])  
            sum_db = sum_db + sum_db_j
          }
          
          db[[i]][[m+1]] = sum_db / (dim(x[[i]])[1])  
        }#m+1
        #update db 1:m+1
        
      }else{
        
        for(i in 1:d){
          sum_db = 0
          for(j in 1:(dim(x[[i]])[1])){
            sum_db_j = c(dO[[i]][j]) * W[[i]][[m+1]]
            for(l in m:2){
              sum_db_j = sum_db_j %*% diag( dH[[i]][[l]][j,] ) %*% W[[i]][[l]]
            }
            sum_db_j = sum_db_j %*% diag( dH[[i]][[1]][j,] )
            sum_db_j = t( sum_db_j ) * (O[[i]][j]-y[[i]][j])
            sum_db = sum_db + sum_db_j
          }
          db[[i]][[1]] = sum_db / (dim(x[[i]])[1])  
        }#1
        
        for(i in 1:d){
          if(m==2){
            sum_db=0
            for(j in 1:(dim(x[[i]])[1])){
              sum_db_j = c(dO[[i]][j]) * W[[i]][[m+1]]
              sum_db_j = sum_db_j %*% diag( dH[[i]][[m]][j,] )
              sum_db_j = t( sum_db_j ) * (O[[i]][j]-y[[i]][j])
              sum_db = sum_db + sum_db_j
            }
            db[[i]][[m]] = sum_db / (dim(x[[i]])[1]) 
            #m
          }else{
            for(l in 2:(m-1)){
              sum_db = 0
              for(j in 1:(dim(x[[i]])[1])){
                sum_db_j = c(dO[[i]][j]) * W[[i]][[m+1]]
                for(ll in m:(l+1)){
                  sum_db_j = sum_db_j %*% diag( dH[[i]][[ll]][j,] ) %*% W[[i]][[ll]]
                }
                sum_db_j = sum_db_j %*% diag( dH[[i]][[l]][j,] )
                sum_db_j = t( sum_db_j ) * (O[[i]][j]-y[[i]][j])
                sum_db = sum_db + sum_db_j
              }
              db[[i]][[l]] = sum_db / (dim(x[[i]])[1]) 
            } #2:m-1 
            
            sum_db=0
            for(j in 1:(dim(x[[i]])[1])){
              sum_db_j = c(dO[[i]][j]) * W[[i]][[m+1]]
              sum_db_j = sum_db_j %*% diag( dH[[i]][[m]][j,] )
              sum_db_j = t( sum_db_j ) * (O[[i]][j]-y[[i]][j])
              sum_db = sum_db + sum_db_j
            }
            db[[i]][[m]] = sum_db / (dim(x[[i]])[1]) 
            #m
          }
        }
        #2:m
        
        for(i in 1:d){
          sum_db = 0
          for(j in 1:(dim(x[[i]])[1])){
            sum_db_j = c(dO[[i]][j]) * (O[[i]][j]-y[[i]][j])  
            sum_db = sum_db + sum_db_j
          }
          db[[i]][[m+1]] = sum_db / (dim(x[[i]])[1])
        }
        #m+1
        #update db 1:m+1
        
      }
      
      w_old = w
      #backup w
      for(i in 1:d){
        Vw[[i]] = alpha * Vw[[i]] - eta * dw[[i]]
        
        w[[i]] = w[[i]] + Vw[[i]] 
      }
      #backpropagation by Momentum for w
      
      for(i in 1:d){
        for(j in 1:p){
          for(k in 1:p){
            if(j!=k){
              w[[i]][j,k]=0
            }
            
          }
        }
      }
      #adjust non-diag
      
      W_old = W
      #backup W
      
      for(i in 1:d){
        for(l in 1:(m+1)){
          VW[[i]][[l]] = alpha * VW[[i]][[l]] - eta * dW[[i]][[l]]
          Vb[[i]][[l]] = alpha * Vb[[i]][[l]] - eta * db[[i]][[l]]
          
          W[[i]][[l]] = W[[i]][[l]] + VW[[i]][[l]] 
          b[[i]][[l]] = b[[i]][[l]] + Vb[[i]][[l]] 
        }
      }
      #backpropagation by Momentum for W b
      
      for(i in 1:d){
        for(l in 1:(m+1)){
          diffW[[i]][[l]] = W[[i]][[l]] - W_old[[i]][[l]]
        }
      }
      
      diffW_max = matrix(0,d,m+1)
      for(i in 1:d){
        for(l in 1:(m+1)){
          diffW_max[i,l]=max(abs(diffW[[i]][[l]]))
        }
      }
      
      ###
      for(i in 1:d){
        diffw[[i]] = w[[i]] - w_old[[i]]
      }
      
      diffw_max = rep(0,d)
      for(i in 1:d){
        diffw_max[i]=max(abs(diffw[[i]]))
      }
      
      maximum1 = max(diffW_max)
      maximum2 = max(diffw_max)
      maximum = max(maximum1,maximum2)
      #compute maximum difference value 
      
      points(tt,maximum)
      #plot points
      
      if(maximum < varepsilon | tt == maxiter){
        break
      }
      #stop criterion
      
      if(m==1){
        
        for(i in 1:d){
          for(j in 1:(dim(x[[i]])[1])){
            
            l=1
            h=t( tanh( W[[i]][[l]] %*% w[[i]] %*% t(t(x[[i]][j,])) + b[[i]][[l]] ) )
            O[[i]][j]= tanh( W[[i]][[m+1]] %*% t(h) + b[[i]][[m+1]] )
            
          }
        }
        #compute Output O
        #update O 
        
      }else{
        
        for(i in 1:d){
          for(j in 1:(dim(x[[i]])[1])){
            
            l=1
            h=t( tanh( W[[i]][[l]] %*% w[[i]] %*% t(t(x[[i]][j,])) + b[[i]][[l]] ) )
            repeat{
              l=l+1
              h=t (tanh( W[[i]][[l]] %*% t(h) + b[[i]][[l]] ) )
              if(l == m)
                break
            }
            O[[i]][j]= tanh( W[[i]][[m+1]] %*% t(h) + b[[i]][[m+1]] )
            
          }
        }
        #compute Output O
        #update O 
        
      }
      
    }
    #HoPIN training
  
  index=NULL
  #generate variable index
  
  sum_w=diag(p)
  for(i in 1:d){
    sum_w = sum_w * w[[i]]
  }
  for(k in 1:p){
    sum_w[k,k]<-abs(sum_w[k,k])^(1/d)*sign(sum_w[k,k])
  }
  
  max_w=max(abs(sum_w))
  for(k in 1:p){
    if( abs(sum_w[k,k]) >= max_w * threshold ){
      index=cbind(index,k)
    }
  }
  
  SSE=0
  for(i in 1:d){
    for(j in 1:(dim(x[[i]])[1])){
      SSE=SSE+(y[[i]][j]-O[[i]][j])^2
    }
  } #compute SSE among d datasets
  
  return(list(m=m,w=w,W=W,b=b,index=index,SSE=SSE,iteration=tt,yhat=O))
}
###############################################################




##################### Algorithm - HePIN #######################
####### This function will give the selected feature   ########
hepin=function(x,y,m,varepsilon,maxiter,alpha,eta,lambda,beta,Lambda,a,sig,threshold,node){
  ##### Input: #####
  # x: a list contain all datasets' X
  #	y: a list contain all datasets' response Y
  #	m: number of layers
  # node: number of hidden nodes
  # varepsilon: condition of convergence
  #	maxiter: maximum number of iterations
  # alpha: momentum parameter
  # eta: learning rate
  # lambda, beta: tunning parameters for variable selection in the objective function
  # Lambda: tunning parameter for preventing overfitting in the objective function
  #	a, sig: smoothing parameters for the smoothing function in the objective function
  #	threshold: threshold: variable selection threshold
  #
  ##### Output: #####
  # w, W, b: the values of w, W, and b in the fitted network
  #	index: selected important variables for each dataset (a list)
  #	SSE: total SSE between fitted values and true values
  #	tt: number of final iterations

  p=dim(x[[1]])[2]#p is the number of variables
  q=node#q is the number of hidden nodes
  d=length(x)#d is the number of datasets
  tt=0#t is the record of iteration
  
  w=list()
  for(i in 1:d){
    #w[[i]]=diag(runif(p))
    w[[i]]=diag(p)
  }
  #w is initialization weight matrix list for input
  
  Vw=list()
  for(i in 1:d){
    Vw[[i]]=matrix(0,nrow=p,ncol=p)
  }
  #Vw is initialization 1st moment matrix list for w
  
  
  dw=list()
  for(i in 1:d){
    dw[[i]]=matrix(0,nrow=p,ncol=p)
  }
  
  #dw initialization
  
  diffw=list()
  for(i in 1:d){
    diffw[[i]]=matrix(0,nrow=p,ncol=p)
  }
  #w difference initialization diffw  
  
  if(m==1){
    
    W=list()
    for(i in 1:d){
      W[[i]]=list()
      W[[i]][[1]]=matrix(runif(q*p,min = - sqrt(6/(q+p)), max = sqrt(6/(q+p)) ),nrow=q)
    }
    
    for(i in 1:d){
      W[[i]][[m+1]]=matrix(runif(q,min = - sqrt(6/(q+1)), max = sqrt(6/(q+1)) ),nrow=1)
    }
    #W is initialization weight matrix list
    
    b=list()
    for(i in 1:d){
      b[[i]]=list()
      for(l in 1:m){
        b[[i]][[l]]=matrix(0,nrow=q)
      }
    }
    
    for(i in 1:d){
      b[[i]][[m+1]]=matrix(0,nrow=1)
    }
    #b is initialization bias matrix list
    
    O=list()
    for(i in 1:d){
      O[[i]]=matrix(rep(0,dim(x[[i]])[1]),ncol = 1)
    }
    
    for(i in 1:d){
      for(j in 1:(dim(x[[i]])[1])){
        l=1
        h=t( tanh( W[[i]][[l]] %*% w[[i]] %*% t(t(x[[i]][j,])) + b[[i]][[l]] ) )
        O[[i]][j]=tanh( W[[i]][[m+1]] %*% t(h) + b[[i]][[m+1]] )
        
      }
    }
    #compute Output O
    
    dO=list()
    for(i in 1:d){
      dO[[i]]= matrix(rep(0,dim(x[[i]])[1]),ncol = 1)
    }
    for(i in 1:d){
      for(j in 1:(dim(x[[i]])[1])){
        dO[[i]][j]=dtanh( W[[i]][[m+1]] %*% t(h) + b[[i]][[m+1]] )
      }
    }
    #dO initialization
    
    VW=list()
    for(i in 1:d){
      VW[[i]]=list()
      VW[[i]][[1]]=matrix(rep(0,q*p),nrow=q)
    }
    
    for(i in 1:d){
      VW[[i]][[m+1]]=matrix(rep(0,q),nrow=1)
    }
    #VW is initialization 1st moment matrix list for W
    
    Vb=list()
    for(i in 1:d){
      Vb[[i]]=list()
      for(l in 1:m){
        Vb[[i]][[l]]=matrix(0,nrow=q)
      }
    }
    
    for(i in 1:d){
      Vb[[i]][[m+1]]=matrix(0,nrow=1)
    }
    #Vb is initialization 1st moment matrix list for b
    
    
    H=list()
    for(i in 1:d){
      H[[i]]=list()
    }
    for(i in 1:d){
      for(l in 1:m){
        H[[i]][[l]]=matrix( rep( 0 , dim(x[[i]])[1] * q ),ncol=q)
      }
    }
    #H initialization
    
    dH=list()
    for(i in 1:d){
      dH[[i]]=list()
    }
    for(i in 1:d){
      for(l in 1:m){
        dH[[i]][[l]]=matrix( rep( 0 , dim(x[[i]])[1] * q ),ncol=q)
      }
    }
    #dH initialization
    
    C=diag(p)
    #C initialization
    CC=list()
    for (i in 1:d){
      CC[[i]]=diag(p)
    }
    #CC initialization
    
    dW=list()
    for(i in 1:d){
      dW[[i]]=list()
      dW[[i]][[1]]=matrix( rep(0,q*p) ,nrow=q)
    }
    
    for(i in 1:d){
      dW[[i]][[m+1]]=matrix( rep(0,q) ,nrow=1)
    }
    #dW initialization
    
    db=list()
    for(i in 1:d){
      db[[i]]=list()
      for(l in 1:m){
        db[[i]][[l]]=matrix(0,nrow=q)
      }
    }
    
    for(i in 1:d){
      db[[i]][[m+1]]=matrix(0,nrow=1)
    }
    #db initialization
    
    diffW=list()
    for(i in 1:d){
      diffW[[i]]=list()
      diffW[[i]][[1]]=matrix( rep(0,q*p) ,nrow=q)
    }
    
    for(i in 1:d){
      diffW[[i]][[m+1]]=matrix( rep(0,q) ,nrow=1)
    }
    #W difference initialization diffW
  }else{
    
    W=list()
    for(i in 1:d){
      W[[i]]=list()
      W[[i]][[1]]=matrix(runif(q*p,min = - sqrt(6/(q+p)), max = sqrt(6/(q+p)) ),nrow=q)
    }
    
    for(i in 1:d){
      for(l in 2:m){
        W[[i]][[l]]=matrix(runif(q*q,min = - sqrt(6/(q+q)), max = sqrt(6/(q+q)) ),nrow=q)
      }
    }
    
    for(i in 1:d){
      W[[i]][[m+1]]=matrix(runif(q,min = - sqrt(6/(q+1)), max = sqrt(6/(q+1)) ),nrow=1)
    }
    #W is initialization weight matrix list
    
    b=list()
    for(i in 1:d){
      b[[i]]=list()
      for(l in 1:m){
        b[[i]][[l]]=matrix(0,nrow=q)
      }
    }
    
    for(i in 1:d){
      b[[i]][[m+1]]=matrix(0,nrow=1)
    }
    #b is initialization bias matrix list
    
    O=list()
    for(i in 1:d){
      O[[i]]=matrix(rep(0,dim(x[[i]])[1]),ncol = 1)
    }
    
    for(i in 1:d){
      for(j in 1:(dim(x[[i]])[1])){
        
        l=1
        h=t( tanh( W[[i]][[l]] %*% w[[i]] %*% t(t(x[[i]][j,])) + b[[i]][[l]] ) )
        repeat{
          l=l+1
          h=t (tanh( W[[i]][[l]] %*% t(h) + b[[i]][[l]] ) )
          if(l == m)
            break
        }
        O[[i]][j]=tanh( W[[i]][[m+1]] %*% t(h) + b[[i]][[m+1]] )
        
      }
    }
    #compute Output O
    
    dO=list()
    for(i in 1:d){
      dO[[i]]= matrix(rep(0,dim(x[[i]])[1]),ncol = 1)
    }
    for(i in 1:d){
      for(j in 1:(dim(x[[i]])[1])){
        l=1
        h=t( tanh( W[[i]][[l]] %*% w[[i]] %*% t(t(x[[i]][j,])) + b[[i]][[l]] ) )
        repeat{
          l=l+1
          h=t (tanh( W[[i]][[l]] %*% t(h) + b[[i]][[l]] ) )
          if(l == m)
            break
        }
        dO[[i]][j]=dtanh( W[[i]][[m+1]] %*% t(h) + b[[i]][[m+1]] )        
      }
    }
    #dO initialization
    
    VW=list()
    for(i in 1:d){
      VW[[i]]=list()
      VW[[i]][[1]]=matrix(rep(0,q*p),nrow=q)
    }
    
    for(i in 1:d){
      for(l in 2:m){
        VW[[i]][[l]]=matrix(rep(0,q*q),nrow=q)
      }
    }
    
    for(i in 1:d){
      VW[[i]][[m+1]]=matrix(rep(0,q),nrow=1)
    }
    #VW is initialization 1st moment matrix list for W
    
    Vb=list()
    for(i in 1:d){
      Vb[[i]]=list()
      for(l in 1:m){
        Vb[[i]][[l]]=matrix(0,nrow=q)
      }
    }
    
    for(i in 1:d){
      Vb[[i]][[m+1]]=matrix(0,nrow=1)
    }
    #Vb is initialization 1st moment matrix list for b
    
    
    H=list()
    for(i in 1:d){
      H[[i]]=list()
    }
    for(i in 1:d){
      for(l in 1:m){
        H[[i]][[l]]=matrix( rep( 0 , dim(x[[i]])[1] * q ),ncol=q)
      }
    }
    #H initialization
    
    dH=list()
    for(i in 1:d){
      dH[[i]]=list()
    }
    for(i in 1:d){
      for(l in 1:m){
        dH[[i]][[l]]=matrix( rep( 0 , dim(x[[i]])[1] * q ),ncol=q)
      }
    }
    #dH initialization
    
    C=diag(p)
    #C initialization
    CC=list()
    for (i in 1:d){
      CC[[i]]=diag(p)
    }
    #CC initialization
    
    dW=list()
    for(i in 1:d){
      dW[[i]]=list()
      dW[[i]][[1]]=matrix( rep(0,q*p) ,nrow=q)
    }
    
    for(i in 1:d){
      for(l in 2:m){
        dW[[i]][[l]]=matrix( rep(0,q*q) ,nrow=q)
      }
    }
    
    for(i in 1:d){
      dW[[i]][[m+1]]=matrix( rep(0,q) ,nrow=1)
    }
    #dW initialization
    
    db=list()
    for(i in 1:d){
      db[[i]]=list()
      for(l in 1:m){
        db[[i]][[l]]=matrix(0,nrow=q)
      }
    }
    
    for(i in 1:d){
      db[[i]][[m+1]]=matrix(0,nrow=1)
    }
    #db initialization
    
    diffW=list()
    for(i in 1:d){
      diffW[[i]]=list()
      diffW[[i]][[1]]=matrix( rep(0,q*p) ,nrow=q)
    }
    
    for(i in 1:d){
      for(l in 2:m){
        diffW[[i]][[l]]=matrix( rep(0,q*q) ,nrow=q)
      }
    }
    
    for(i in 1:d){
      diffW[[i]][[m+1]]=matrix( rep(0,q) ,nrow=1)
    }
    #W difference initialization diffW
    
  }
  
  
  plot(1:maxiter, seq(0,0.01,length.out = maxiter), xlim = c(1,maxiter), ylim = c (0,0.01), type = "n", axes = F, xlab = "iteration times", ylab="maximum difference")
  box()
  axis(side = 1, at = c(seq(1,maxiter+1,length.out=51)))
  axis(side = 2, at = c(seq(0,0.01,length.out = 21)))
  abline(h=0.001,col="red")
  #generate a null plot
  
  repeat{
    
    w_C<-w
    
    for(k in 1:p){
      
      sum_w=0
      for(i in 1:d){
        sum_w = sum_w + (w_C[[i]][k,k])^2
      }
      
      if(sum_w <= d * 10^(-20)){
        for(i in 1:d){
          w_C[[i]][k,k]=10^(-10)
        }
      }
    }
    #prevent small value 
    
    if(m==1){
      
      for(i in 1:d){
        for(j in 1:(dim(x[[i]])[1])){
          H[[i]][[1]][j,]= t( tanh( W[[i]][[1]] %*% w[[i]] %*% t(t(x[[i]][j,])) + b[[i]][[1]] ) ) 
        }
      }
      
    }else{
      
      for(i in 1:d){
        for(j in 1:(dim(x[[i]])[1])){
          H[[i]][[1]][j,]= t( tanh( W[[i]][[1]] %*% w[[i]] %*% t(t(x[[i]][j,])) + b[[i]][[1]] ) ) 
        }
      }
      
      for(i in 1:d){
        for(l in 2:m){
          for(j in 1:(dim(x[[i]])[1])){
            H[[i]][[l]][j,] = t( tanh( W[[i]][[l]] %*% t(t(H[[i]][[l-1]][j,])) + b[[i]][[l]] ) )
          }
        }
      }
      #compute H
    }
    
    
    if(m==1){
      
      for(i in 1:d){
        for(j in 1:(dim(x[[i]])[1])){
          dH[[i]][[1]][j,]= t( dtanh( W[[i]][[1]] %*% w[[i]] %*% t(t(x[[i]][j,])) + b[[i]][[1]] ) ) 
        }
      }
      
    }else{
      
      for(i in 1:d){
        for(j in 1:(dim(x[[i]])[1])){
          dH[[i]][[1]][j,]= t( dtanh( W[[i]][[1]] %*% w[[i]] %*% t(t(x[[i]][j,])) + b[[i]][[1]] ) ) 
        }
      }
      
      for(i in 1:d){
        for(l in 2:m){
          for(j in 1:(dim(x[[i]])[1])){
            dH[[i]][[l]][j,] = t( dtanh( W[[i]][[l]] %*% t(t(H[[i]][[l-1]][j,])) + b[[i]][[l]] ) )
          }
        }
      }
      #compute dH
    }
    
    for(k in 1:p){
      hat_w=0
      for(i in 1:d){
        hat_w = hat_w + (w_C[[i]][k,k])^2 
      }
      hat_w=sqrt(hat_w)
      if(hat_w >= a){
        C[k,k]=1/sqrt(hat_w)
      }else{
        C[k,k]=(hat_w)*(-1/(2*a^2)) + 3/(2*a)
      }
      
    }
    #compute C for t iteration by smoothing function
    
    for(i in 1:d){
      for(k in 1:p){
        if(abs(w_C[[i]][k,k]) >= sig){
          CC[[i]][k,k]=1/abs(w_C[[i]][k,k])
        }else{
          CC[[i]][k,k]=1/sqrt((w_C[[i]][k,k])^2 + sig)
        }
      }
    }
    #compute CC for t iteration by smoothing function
    
    tt=tt+1
    #iteration count
    
    if(m==1){
      
      for(i in 1:d){
        sum_dw = 0
        for(j in 1:(dim(x[[i]])[1])){
          sum_dw_j = c(dO[[i]][j]) * W[[i]][[2]]
          sum_dw_j = sum_dw_j %*% diag( dH[[i]][[1]][j,] ) %*% W[[i]][[1]]
          sum_dw_j = t(t( x[[i]][j,] )) %*% sum_dw_j  
          sum_dw_j = t( sum_dw_j ) * (O[[i]][j]-y[[i]][j])
          sum_dw = sum_dw + sum_dw_j
        }
        dw[[i]] = sum_dw / (dim(x[[i]])[1]) + (1-beta) * lambda * w[[i]] * C + beta * lambda * w[[i]] * CC[[i]] 
      }
      #w
      
      for(i in 1:d){
        sum_dW = 0
        for(j in 1:(dim(x[[i]])[1])){
          sum_dW_j = c(dO[[i]][j]) * W[[i]][[2]]
          sum_dW_j = sum_dW_j %*% diag( dH[[i]][[1]][j,] )
          sum_dW_j = w[[i]] %*% t(t( x[[i]][j,] )) %*% sum_dW_j  
          sum_dW_j = t( sum_dW_j ) * (O[[i]][j]-y[[i]][j])
          sum_dW = sum_dW + sum_dW_j
        }
        dW[[i]][[1]] = sum_dW / (dim(x[[i]])[1]) + Lambda * W[[i]][[1]]   
      }
      #1
      
      for(i in 1:d){
        sum_dW = 0
        for(j in 1:(dim(x[[i]])[1])){
          sum_dW_j = t(H[[i]][[m]][j,]) * c(dO[[i]][j]) * (O[[i]][j]-y[[i]][j])  
          sum_dW = sum_dW + sum_dW_j
        }
        dW[[i]][[m+1]] = sum_dW / (dim(x[[i]])[1]) + Lambda * W[[i]][[m+1]] 
      }
      #m+1
      #update dW 1:m+1
      
    }else{
      
      for(i in 1:d){
        sum_dw = 0
        for(j in 1:(dim(x[[i]])[1])){
          sum_dw_j = c(dO[[i]][j]) * W[[i]][[m+1]]
          for(l in m:2){
            sum_dw_j = sum_dw_j %*% diag( dH[[i]][[l]][j,] ) %*% W[[i]][[l]]
          }
          sum_dw_j = sum_dw_j %*% diag( dH[[i]][[1]][j,] ) %*% W[[i]][[1]]
          sum_dw_j = t(t( x[[i]][j,] )) %*% sum_dw_j  
          sum_dw_j = t( sum_dw_j ) * (O[[i]][j]-y[[i]][j])
          sum_dw = sum_dw + sum_dw_j
        }
        dw[[i]] = sum_dw / (dim(x[[i]])[1]) + (1-beta) * lambda * w[[i]] * C + beta * lambda * w[[i]] * CC[[i]]  
      }
      #w
      
      for(i in 1:d){
        sum_dW = 0
        for(j in 1:(dim(x[[i]])[1])){
          sum_dW_j = c(dO[[i]][j]) * W[[i]][[m+1]]
          for(l in m:2){
            sum_dW_j = sum_dW_j %*% diag( dH[[i]][[l]][j,] ) %*% W[[i]][[l]]
          }
          sum_dW_j = sum_dW_j %*% diag( dH[[i]][[1]][j,] )
          sum_dW_j = w[[i]] %*% t(t( x[[i]][j,] )) %*% sum_dW_j  
          sum_dW_j = t( sum_dW_j ) * (O[[i]][j]-y[[i]][j])
          sum_dW = sum_dW + sum_dW_j
        }
        dW[[i]][[1]] = sum_dW / (dim(x[[i]])[1]) + Lambda * W[[i]][[1]]  
      }
      #1
      
      for(i in 1:d){
        if(m==2){
          sum_dW=0
          for(j in 1:(dim(x[[i]])[1])){
            sum_dW_j = c(dO[[i]][j]) * W[[i]][[m+1]]
            sum_dW_j = sum_dW_j %*% diag( dH[[i]][[m]][j,] )
            sum_dW_j = t(t( H[[i]][[m-1]][j,] )) %*% sum_dW_j  
            sum_dW_j = t( sum_dW_j ) * (O[[i]][j]-y[[i]][j])
            sum_dW = sum_dW + sum_dW_j
          }
          dW[[i]][[m]] = sum_dW / (dim(x[[i]])[1]) + Lambda * W[[i]][[m]]  
          #m
        }else{
          for(l in 2:(m-1)){
            sum_dW = 0
            for(j in 1:(dim(x[[i]])[1])){
              sum_dW_j = c(dO[[i]][j]) * W[[i]][[m+1]]
              for(ll in m:(l+1)){
                sum_dW_j = sum_dW_j %*% diag( dH[[i]][[ll]][j,] ) %*% W[[i]][[ll]]
              }
              sum_dW_j = sum_dW_j %*% diag( dH[[i]][[l]][j,] )
              sum_dW_j = t(t( H[[i]][[l-1]][j,] )) %*% sum_dW_j  
              sum_dW_j = t( sum_dW_j ) * (O[[i]][j]-y[[i]][j])
              sum_dW = sum_dW + sum_dW_j
            }
            dW[[i]][[l]] = sum_dW / (dim(x[[i]])[1]) + Lambda * W[[i]][[l]]  

          }
        }
        #2:m-1 
        sum_dW=0
        for(j in 1:(dim(x[[i]])[1])){
          sum_dW_j = c(dO[[i]][j]) * W[[i]][[m+1]]
          sum_dW_j = sum_dW_j %*% diag( dH[[i]][[m]][j,] )
          sum_dW_j = t(t( H[[i]][[m-1]][j,] )) %*% sum_dW_j  
          sum_dW_j = t( sum_dW_j ) * (O[[i]][j]-y[[i]][j])
          sum_dW = sum_dW + sum_dW_j
        }
        dW[[i]][[m]] = sum_dW / (dim(x[[i]])[1]) + Lambda * W[[i]][[m]]  
        #m
      }
      #2:m
      
      for(i in 1:d){
        sum_dW = 0
        for(j in 1:(dim(x[[i]])[1])){
          sum_dW_j = t(H[[i]][[m]][j,]) * c(dO[[i]][j]) * (O[[i]][j]-y[[i]][j])  
          sum_dW = sum_dW + sum_dW_j
        }
        dW[[i]][[m+1]] = sum_dW / (dim(x[[i]])[1]) + Lambda * W[[i]][[m+1]] 
      }
      #m+1
      #update dW 1:m+1
    }
    
    if(m==1){
      for(i in 1:d){
        sum_db = 0
        for(j in 1:(dim(x[[i]])[1])){
          sum_db_j = c(dO[[i]][j]) * W[[i]][[2]]
          sum_db_j = sum_db_j %*% diag( dH[[i]][[1]][j,] )
          sum_db_j = t( sum_db_j ) * (O[[i]][j]-y[[i]][j])
          sum_db = sum_db + sum_db_j
        }
        db[[i]][[1]] = sum_db / (dim(x[[i]])[1])  
      }
      #1
      for(i in 1:d){
        sum_db = 0
        for(j in 1:(dim(x[[i]])[1])){
          sum_db_j = c(dO[[i]][j]) * (O[[i]][j]-y[[i]][j])  
          sum_db = sum_db + sum_db_j
        }
        db[[i]][[m+1]] = sum_db / (dim(x[[i]])[1])  
      }
      #m+1
      #update db 1:m+1
    }else{
      for(i in 1:d){
        sum_db = 0
        for(j in 1:(dim(x[[i]])[1])){
          sum_db_j = c(dO[[i]][j]) * W[[i]][[m+1]]
          for(l in m:2){
            sum_db_j = sum_db_j %*% diag( dH[[i]][[l]][j,] ) %*% W[[i]][[l]]
          }
          sum_db_j = sum_db_j %*% diag( dH[[i]][[1]][j,] )
          sum_db_j = t( sum_db_j ) * (O[[i]][j]-y[[i]][j])
          sum_db = sum_db + sum_db_j
        }
        db[[i]][[1]] = sum_db / (dim(x[[i]])[1])  
      }
      #1
      
      for(i in 1:d){
        if(m==2){
            sum_db=0
          for(j in 1:(dim(x[[i]])[1])){
            sum_db_j = c(dO[[i]][j]) * W[[i]][[m+1]]
            sum_db_j = sum_db_j %*% diag( dH[[i]][[m]][j,] )
            sum_db_j = t( sum_db_j ) * (O[[i]][j]-y[[i]][j])
            sum_db = sum_db + sum_db_j
          }
          db[[i]][[m]] = sum_db / (dim(x[[i]])[1]) 
          #m
          
        }else{
          for(l in 2:(m-1)){
            sum_db = 0
            for(j in 1:(dim(x[[i]])[1])){
              sum_db_j = c(dO[[i]][j]) * W[[i]][[m+1]]
              for(ll in m:(l+1)){
                sum_db_j = sum_db_j %*% diag( dH[[i]][[ll]][j,] ) %*% W[[i]][[ll]]
              }
              sum_db_j = sum_db_j %*% diag( dH[[i]][[l]][j,] )
              sum_db_j = t( sum_db_j ) * (O[[i]][j]-y[[i]][j])
              sum_db = sum_db + sum_db_j
            }
            db[[i]][[l]] = sum_db / (dim(x[[i]])[1]) 
          }
          #2:m-1 
          
          sum_db=0
          for(j in 1:(dim(x[[i]])[1])){
            sum_db_j = c(dO[[i]][j]) * W[[i]][[m+1]]
            sum_db_j = sum_db_j %*% diag( dH[[i]][[m]][j,] )
            sum_db_j = t( sum_db_j ) * (O[[i]][j]-y[[i]][j])
            sum_db = sum_db + sum_db_j
          }
          db[[i]][[m]] = sum_db / (dim(x[[i]])[1]) 
          #m
        }
      }
      #2:m
      
      for(i in 1:d){
        sum_db = 0
        for(j in 1:(dim(x[[i]])[1])){
          sum_db_j = c(dO[[i]][j]) * (O[[i]][j]-y[[i]][j])  
          sum_db = sum_db + sum_db_j
        }
        db[[i]][[m+1]] = sum_db / (dim(x[[i]])[1])  
      }
      #m+1
      #update db 1:m+1
    }
    
    w_old = w
    #backup w
    for(i in 1:d){
      Vw[[i]] = alpha * Vw[[i]] - eta * dw[[i]]
      w[[i]] = w[[i]] + Vw[[i]] 
    }
    #backpropagation by Momentum for w
    
    for(i in 1:d){
      for(j in 1:p){
        for(k in 1:p){
          if(j!=k){
            w[[i]][j,k]=0
          }
        }
      }
    }
    #adjust non-diag
    
    W_old = W
    #backup W
    
    for(i in 1:d){
      for(l in 1:(m+1)){
        VW[[i]][[l]] = alpha * VW[[i]][[l]] - eta * dW[[i]][[l]]
        Vb[[i]][[l]] = alpha * Vb[[i]][[l]] - eta * db[[i]][[l]]
        W[[i]][[l]] = W[[i]][[l]] + VW[[i]][[l]] 
        b[[i]][[l]] = b[[i]][[l]] + Vb[[i]][[l]] 
      }
    }
    #backpropagation by Momentum for W b
    
    for(i in 1:d){
      for(l in 1:(m+1)){
        diffW[[i]][[l]] = W[[i]][[l]] - W_old[[i]][[l]]
      }
    }
    
    diffW_max = matrix(0,d,m+1)
    for(i in 1:d){
      for(l in 1:(m+1)){
        diffW_max[i,l]=max(abs(diffW[[i]][[l]]))
      }
    }
    
    ###
    for(i in 1:d){
      diffw[[i]] = w[[i]] - w_old[[i]]
    }
    
    diffw_max = rep(0,d)
    for(i in 1:d){
      diffw_max[i]=max(abs(diffw[[i]]))
    }
    
    maximum1 = max(diffW_max)
    maximum2 = max(diffw_max)
    maximum = max(maximum1,maximum2)
    #compute maximum difference value 
    
    points(tt,maximum)
    #plot points
    
    if(maximum < varepsilon | tt == maxiter){
      break
    }
    #stop criterion
    
    if(m==1){
      
      for(i in 1:d){
        for(j in 1:(dim(x[[i]])[1])){
          l=1
          h=t( tanh( W[[i]][[l]] %*% w[[i]] %*% t(t(x[[i]][j,])) + b[[i]][[l]] ) )
          O[[i]][j]= tanh( W[[i]][[m+1]] %*% t(h) + b[[i]][[m+1]] )
          
        }
      }
      #compute Output O
      #update O 
      
    }else{
      
      for(i in 1:d){
        for(j in 1:(dim(x[[i]])[1])){
          l=1
          h=t( tanh( W[[i]][[l]] %*% w[[i]] %*% t(t(x[[i]][j,])) + b[[i]][[l]] ) )
          repeat{
            l=l+1
            h=t (tanh( W[[i]][[l]] %*% t(h) + b[[i]][[l]] ) )
            if(l == m)
              break
          }
          O[[i]][j]= tanh( W[[i]][[m+1]] %*% t(h) + b[[i]][[m+1]] )
        }
      }
      #compute Output O
      #update O 
    }
  }
  #HePIN training
  
  lastw=list()
  choose_w = rep(0,d)
  for(i in 1:d){
    choose_w[i]<-max(abs(diag(w[[i]])))
    lastw[[i]]<-diag(w[[i]])
  }
  
  ############
  index=list()
  for(i in 1:d){
    index[[i]]=rep(0,p)
  }
  
  for(i in 1:d){
    for(k in 1:p){
      if( abs(w[[i]][k,k]) >= choose_w[i] * threshold ){
        index[[i]][k]=k
      }
    }
  }
  
  for(i in 1:d){
    index[[i]]=index[[i]][index[[i]]>0]
  }
  
  ############SSE
  SSE=0
  for(i in 1:d){
    for(j in 1:(dim(x[[i]])[1])){
      SSE=SSE+(y[[i]][j]-O[[i]][j])^2
    }
  }#compute SSE among d datasets
  
  return(list(m=m,w=w,W=W,b=b,index=index,SSE=SSE,iteration=tt,yhat=O))
}
###################################################################


#############################################################
#####                                                   #####
#####                                                   #####
#####                    Data prediction                #####
#####                                                   #####
#####                                                   #####
#############################################################

pre_pin=function(newx,pin){
  ##### Input: #####
  # newx: List of new values for X at which predictions are to be made.
  #	pin: Fitted "PIN" object.
  ##### Output: #####
  # A list containing the predicted values of all datasets given newx.
  
  m=pin$m
  b=pin$b
  w=pin$w
  W=pin$W
  d=length(newx)
  O_pre=list()
  for(i in 1:d){
    O_pre[[i]]=matrix(rep(0,dim(newx[[i]])[1]),ncol = 1)
  }
  
  if(m==1){
    
    for(i in 1:d){
      for(j in 1:(dim(newx[[i]])[1])){
        
        l=1
        h=t(tanh( W[[i]][[l]] %*% w[[i]] %*% t(t(newx[[i]][j,])) + b[[i]][[l]] ) )
        O_pre[[i]][j]=tanh( W[[i]][[m+1]] %*% t(h) + b[[i]][[m+1]] )
        
      }
    }
    
  }else{
    
    for(i in 1:d){
      for(j in 1:(dim(newx[[i]])[1])){
        
        l=1
        h=t(tanh( W[[i]][[l]] %*% w[[i]] %*% t(t(newx[[i]][j,])) + b[[i]][[l]] ) )
        repeat{
          l=l+1
          h=t (tanh( W[[i]][[l]] %*% t(h) + b[[i]][[l]] ) )
          if(l == m)
            break
        }
        O_pre[[i]][j]=tanh( W[[i]][[m+1]] %*% t(h) + b[[i]][[m+1]] )
        
      }
    }
  }
  return(O_pre)
}
###################################################################

#############################################################
#####                                                   #####
#####                                                   #####
#####       Evaluation Index for the i-th dataset       #####
#####                                                   #####
#####                                                   #####
#############################################################
### Include: SEN, SPE, GM, CCR.                         #####
#############################################################


###################### Evaluation Index #####################
Eva <- function(index, b){
  ##### Input: #####
  # index: important variables selected by the network for the i-th dataset
  #	b: model coefficients of the i-th dataset
  ##### Output: #####
  # Evaluation Index (SEN, SPE, GM, CCR) for the i-th dataset.
  
  trueindex=which(b!=0)
  TP=0
  FP=0
  FN=0
  TN=0
  TP=length(intersect(index,trueindex))
  FP=length(setdiff(index,trueindex))
  FN=length(trueindex)-TP
  TN=p-length(trueindex)-FP
  SEN=TP/(TP+FN)
  SPE=TN/(TN+FP)
  GM=sqrt(SEN*SPE)
  MR=(FP+FN)/(TP+FN+TN+FP)
  CCR=1-MR
  
  eva_index=return(c(SEN,SPE,GM,CCR))
}
###################################################################
