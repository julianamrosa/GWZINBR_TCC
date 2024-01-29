cv <- function(h){
  #N, wt, x, y, G, yhat, yhat2, pihat, hv, coord, _dist_, seq, offset, alphai, S, Si, parg, pargg, ujg, bg, lambdag, pos0, pos1, pos02, nvar  
  if(method=="adaptiven"){
    #hv <- matrix(0,1,1)
    hv <- 0
    #yhat <- matrix(0,1,1)
    yhat <- 0
    # create &output from hv[colname='h'];
  }
  for (i in 1:N){
    for (i in 1:N){
      #seqi <- matrix(i, N, 1)
      seqi <- rep(i, N)
      distan <- cbind(seqi, t(sequ), distance[,i])
      if (distancekm){
        distan[,3] <- distan[,3]*111
      }
    }
    u <- nrow(distan)
    #w <- matrix(0, u, 1)
    w <- rep(0, u)
    for (jj in 1:u){
      w[jj] <- exp(-0.5*(distan[jj,3]/h)^2)
      if (method=="fixed_bsq" | method=="adaptiven"){
        w[jj] <- (1-(distan[jj,3]/h)^2)^2
      }
      if (bandwidth=="cv"){
        w[i] <- 0
      }
    }
    if (method=="fixed_bsq" | method=="adaptiven"){
      position <- which(distan[,3]<=h)
      w[position] <- 0
    }
    if (method=="adaptive_bsq"){
      distan <- distan[order(distan[, 3]), ]
      distan <- cbind(distan, 1:nrow(distan))
      w <- matrix(0, N, 2)
      hn <- distan[h,3]
      for (jj in 1:N){
        if (distan[jj,4]<=h){
          w[jj,1] <- (1-(distan[jj,3]/hn)^2)^2
        }
        else{
          w[jj,1] <- 0
        }
        w[jj,2] <- distan[jj,2]
      }
      if (bandwidth=="cv"){
        w[which(w[,2]==i)] <- 0
      }
      w <- w[order(w[, 2]), ]
      w <- w[,1]
    }
    b <- bg
    nj <- x%*%b+Offset #checar multiplicador
    uj <- exp(nj)
    par <- parg
    lambda <- lambdag
    njl <- G*lambda
    njl <- ifelse(njl>maxg, maxg, njl)
    njl <- ifelse(njl<-maxg, -maxg, njl)
    if (model!="zip" & model!="zinb"){
      zk <- 0
    }
    else{
      lambda0 <- (ncol(pos0)-sum((parg/(uj+parg))^parg))/N
      if (lambda0>0){
        lambda0 <- log(lambda0/(1-lambda0))
        #lambda <- rbind(lambda0, matrix(0, ncol(G)-1, 1))
        lambda <- rbind(lambda0, rep(0, ncol(G)-1))
        njl <- G*lambda
      }
      zk <- 1/(1+exp(-njl)*(par/(par+uj))^par)
      zk <- ifelse(y>0, 0, zk)
    }
    dllike <- 1
    llike <- 0
    j <- 1
    while (abs(dllike)>0.00001 & j<=600){
      ddpar <- 1
      #while (abs(ddpar)>0.000001)
      dpar <- 1
      parold <- par
      aux1 <- 1
      aux2 <- 1
      int <- 1
      if (model=="zip" | model=="poisson"){
        alpha <- E^-6
        par <- 1/alpha
      }
      else{
        if (par<=E^-5){
          if (i>1){
            par <- 1/alphai[i-1,2]
          }
        }
        if (par>=E^6){
          par <- E^6
          dpar <- 0
          alpha <- 1/par
          b <- bg
          uj <- exp(x%*%b+Offset)
          lambda <- lambdag
          njl <- G*lambda
          njl <- ifelse(njl>maxg, maxg, njl)
          njl <- ifelse(njl<-maxg, -maxg, njl)
          zk <- 1/(1+exp(-njl)*(parg/(parg+uj))^parg)
          zk <- ifelse(y>0, 0, zk)
          if (any(lambda)==0){
            zk <- 0
          }
        }
        while (abs(dpar)>0.000001 & aux2<200){
          par <- ifelse(par<E^-10, E^-10, par)
          gf <- sum(w*wt*(1-zk)*(digamma(par+y)-digamma(par)+log(par)+1-log(par+uj)-(par+y)/(par+uj)))
          hess <- sum(w*wt*(1-zk)*(trigamma(par+y)-trigamma(par)+1/par-2/(par+uj)+(y+par)/(par+uj)^2))
          hess <- ifelse(hess==0, E^-23, hess)
          par0 <- par
          par <- par0-solve(hess)%*%gf #multiplicador
          dpar <- par-par0
          if (par>=E^6){
            par <- E^6
            dpar <- 0
            alpha <- 1/par
            b <- bg
            uj <- exp(x%*%b+Offset)
            lambda <- lambdag
            njl <- G*lambda
            njl <- ifelse(njl>maxg, maxg, njl)
            njl <- ifelse(njl<-maxg,-maxg,njl)
            zk <- 1/(1+exp(-njl)*(parg/(parg+uj))^parg)
            zk <- ifelse(y>0,0,zk)
            if (any(lambda)==0){
              zk <- 0
            }
          }
          aux2 <- aux2+1
        }
        if (par<=E^-5){
          par <- E^6
          b <- bg
          uj <- exp(x%*%b+Offset)
          lambda <- lambdag
          njl <- G*lambda
          njl <- ifelse(njl>maxg, maxg, njl)
          njl <- ifelse(njl<-maxg, -maxg, njl)
          zk <- 1/(1+exp(-njl)*(parg/(parg+uj))^parg)
          zk <- ifelse(y>0, 0, zk)
          if (any(lambda)==0){
            zk <- 0
          }
        }
        alpha <- 1/par
        #print(c("i", "j", "b", "lambda", "par", "alpha", "aux2"))
        #print(c(i, j, b, lambda, par, alpha, aux2))
      }
      dev <- 0
      ddev <- 1
      nj <- x%*%b+Offset
      nj <- ifelse(nj>700, 700, nj)
      nj <- ifelse(nj<-700, -700, nj)
      uj <- exp(nj)
      while (abs(ddev)>0.000001 & aux1<100){
        uj <- ifelse(uj>E^100,E^100,uj)
        Ai <- (1-zk)*((uj/(1+alpha*uj)+(y-uj)*(alpha*uj/(1+2*alpha*uj+alpha^2*uj^2))))
        Ai <- ifelse(Ai<=0,E^-5,Ai)
        uj <- ifelse(uj<E^-150,E^-150,uj)
        denz <- (((uj/(1+alpha*uj)+(y-uj)*(alpha*uj/(1+2*alpha*uj+alpha^2*uj^2))))*(1+alpha*uj))
        denz <- ifelse(denz==0,E^-5,denz)
        zj <- (nj+(y-uj)/denz)-Offset
        if (det(t(x)%*%(w*Ai*x*wt))==0){
          #b <- matrix(0, nvar, 1)
          b <- rep(0, nvar)
        }
        else{
          b <- solve(t(x)%*%(w*Ai*x*wt))%*%t(x)%*%(w*Ai*wt*zj)
        }
        nj <- x%*%b+Offset
        nj <- ifelse(nj>700, 700, nj)
        nj <- ifelse(nj<-700, -700, nj)
        uj <- exp(nj)
        olddev <- dev
        uj <- ifelse(uj>E^10, E^10, uj)
        uj <- ifelse(uj==0, E^-10, uj)
        if (par==E^6){
          #gamma1=/*(gamma(par+y)/(gamma(y+1)#gamma(par)))#*/(uj/(uj+par))##y#exp(-uj)
          gamma1 <- (uj/(uj+par))^y*exp(-uj)
        }
        else{
          #gamma1=/*(gamma(par+y)/(gamma(y+1)#gamma(par)))#*/(uj/(uj+par))##y#(par/(uj+par))##par
          gamma1 <- (uj/(uj+par))^y*(par/(uj+par))^par
        }
        gamma1 <- ifelse(gamma1<=0, E^-10, gamma1)
        dev <- sum((1-zk)*(log(gamma1)))
        ddev <- dev-olddev
        #print(c("b", "par", "aux1", "dev", "olddev", "ddev"))
        #print(c(b, par, aux1, dev, olddev, ddev))
        aux1 <- aux1+1
      }
      ddpar <- par-parold
      #end
      ## PARAMOS AQUI ##
      if (model=="zip" | model=="zinb"){
        if (j==1){
          alphatemp <- alpha
          lambdatemp <- lambda[1]
        }
        else{
          alphatemp <- rbind(alphatemp, alpha)
          lambdatemp <- rbind(lambdatemp, lambda[1])
        }
        alphatemp <- round(alphatemp, 7)
        lambdatemp <- round(lambdatemp, 7)
        #print(c('i', 'j', 'alphatemp', 'nrow(alphatemp)', 'ncol(unique(alphatemp))'))
        #print(c(i, j, alphatemp, nrow(alphatemp), ncol(unique(alphatemp))))
        if (model=="zinb"){
          condition <- (j>300 & nrow(alphatemp)>ncol(unique(alphatemp)) & nrow(lambdatemp)>ncol(unique(lambdatemp)))
        }
        else if (model=="zip"){
          #print i j lambdatemp (nrow(lambdatemp)) (ncol(unique(lambdatemp)));
          condition <- (j>300 & nrow(lambdatemp)>ncol(unique(lambdatemp)))
        }
        if (condition){
          #lambda <- matrix(0, ncol(G), 1)
          lambda <- rep(0, ncol(G))
          njl <- G*lambda
          #zk <- matrix(0, n, 1)
          zk <- rep(0, n)
        }
        else{
          aux3 <- 1
          dev <- 0
          ddev <- 1
          njl <- G*lambda
          njl <- ifelse(njl>maxg, maxg, njl)
          njl <- ifelse(njl<-maxg, -maxg, njl)
          pi <- exp(njl)/(1+exp(njl))
          while (abs(ddev)>0.000001 & aux3<100){
            Aii <- pi*(1-pi)
            Aii <- ifelse(Aii<=0, E^-5, Aii)	
            zj <- njl+(zk-pi)/Aii
            if (det(t(G*Aii*w*wt)%*%G)=0){ #multiplicador
              lambda <- matrix(0, ncol(G), 1)
            }
            else{
              lambda <- solve(t(G*Aii*w*wt)%*%G)%*%t(G*Aii*w*wt)%*%zj
            }
            njl <- G*lambda
            njl <- ifelse(njl>maxg, maxg, njl)
            njl <- ifelse(njl<-maxg, -maxg, njl)
            pi <- exp(njl)/(1+exp(njl))
            olddev <- dev
            dev <- sum(zk*njl-log(1+exp(njl)))
            ddev <- dev-olddev
            #print(c("lambda", "aux3", "dev", "olddev", "ddev"))
            #print(c(lambda, aux3, dev, olddev, ddev))
            aux3 <- aux3+1
          }
        }
      }
      njl <- G*lambda
      njl <- ifelse(njl>maxg, maxg, njl)
      njl <- ifelse(njl<-maxg, -maxg, njl)
      zk <- 1/(1+exp(-njl)*(par/(par+uj))^par)
      zk <- ifelse(y>0, 0, zk)
      if (any(lambda)==0){
        #zk <- matrix(0, n, 1)
        zk <- rep(0, n)
      }
      if (model!="zip" & model!="zinb"){
        zk <- 0
      }
      oldllike <- llike
      llike <- sum(zk*(njl)-log(1+exp(njl))+(1-zk)*(log(gamma1)))
      dllike <- llike-oldllike
      #print(c("i", "j", "b", "alpha", "lambda", "llike", "dllike"))
      #print(c(i, j, b, alpha, lambda, llike, dllike))
      j <- j+1
    }
    if (method=="fixed_g" | method=="fixed_bsq" | method=="adaptive_bsq"){
    yhat[i]=uj[i];
    pihat[i]=njl[i];
    alphai[i]= alpha;
    if det(x`*(w#Ai#x#wt))=0 then S[i]=0;
    else S[i]= (x[i,]*inv(x`*(w#Ai#x#wt))*(x#w#Ai#wt)`)[i];
    %if %upcase(&model)=ZIP or %upcase(&model)=ZINB %then %do;
    yhat[i]=(uj#(1-exp(njl)/(1+exp(njl))))[i];
    yhat2[i]=uj[i];
    if det(G`*(w#Aii#G#wt))=0 then Si[i]=0;
    else Si[i]= (G[i,]*inv(G`*(w#Aii#G#wt))*(G#w#Aii#wt)`)[i];
    if any(lambda)=0 then Si[i]=0;
    %end;
  }
  CV=((y-yhat)#wt)`*(y-yhat);
  
  _par_=1/alphai;
  
  %IF %UPCASE(&model)=ZINB or %UPCASE(&model)=ZIP %THEN %DO;
  if any(lambda)=0 then do;
  ll=sum(-log(0+exp(pihat[pos0]))+log(0*exp(pihat[pos0])+(_par_[pos0]/(_par_[pos0]+yhat2[pos0]))##_par_[pos0]))+
  sum(-log(0+exp(pihat[pos1]))+lgamma(_par_[pos1]+y[pos1,])-lgamma(y[pos1,]+1)-lgamma(_par_[pos1])+
  y[pos1]#log(yhat2[pos1]/(_par_[pos1]+yhat2[pos1]))+_par_[pos1]#log(_par_[pos1]/(_par_[pos1]+yhat2[pos1])));
  llnull1=sum(-log(1+zk[pos0])+log(zk[pos0]+(_par_[pos0]/(_par_[pos0]+y[pos0,]))##_par_[pos0]))+
  sum(-log(1+zk[pos1])+lgamma(_par_[pos1]+y[pos1,])-lgamma(y[pos1,]+1)-lgamma(_par_[pos1])+
  y[pos1]#log(y[pos1,]/(_par_[pos1]+y[pos1,]))+_par_[pos1]#log(_par_[pos1]/(_par_[pos1]+y[pos1,])));
  llnull2=sum(-log(1+0)+log(0+(_par_/(_par_+y[:]))##_par_))+
  sum(-log(1+0)+lgamma(_par_+y)-lgamma(y+1)-lgamma(_par_)+
  y#log(y[:]/(_par_+y[:]))+_par_#log(_par_/(_par_+y[:])));
  end;
  else do;
  ll=sum(-log(1+exp(pihat[pos0]))+log(exp(pihat[pos0])+(_par_[pos0]/(_par_[pos0]+yhat2[pos0]))##_par_[pos0]))+
  sum(-log(1+exp(pihat[pos1]))+lgamma(_par_[pos1]+y[pos1])-lgamma(y[pos1]+1)-lgamma(_par_[pos1])+
  y[pos1]#log(yhat2[pos1]/(_par_[pos1]+yhat2[pos1]))+_par_[pos1]#log(_par_[pos1]/(_par_[pos1]+yhat2[pos1])));
  llnull1=sum(-log(1+zk[pos0])+log(zk[pos0]+(_par_[pos0]/(_par_[pos0]+y[pos0]))##_par_[pos0]))+
  sum(-log(1+zk[pos1])+lgamma(_par_[pos1]+y[pos1])-lgamma(y[pos1]+1)-lgamma(_par_[pos1])+
  y[pos1]#log(y[pos1]/(_par_[pos1]+y[pos1]))+_par_[pos1]#log(_par_[pos1]/(_par_[pos1]+y[pos1])));
  end;
  dev=2*(llnull1-ll);
  npar=sum(S)+sum(Si);
  AIC=2*npar-2*ll;
  AICc=AIC+2*(npar*(npar+1)/(n-npar-1));
  %IF %UPCASE(&model)=ZINB %THEN %DO;
  AIC=2*(npar+npar/(ncol(x)+ncol(G)))-2*ll;
  AICc=AIC+2*((npar+npar/(ncol(x)+ncol(G)))*((npar+npar/(ncol(x)+ncol(G)))+1)/(n-(npar+npar/(ncol(x)+ncol(G)))-1));
  %END;
  %END;
  %IF %UPCASE(&model)=POISSON or %UPCASE(&model)=NEGBIN %THEN %DO;
  if ncol(pos02)=0 then do;
  pos0=pos1;
  pos0x=1;
  pos0xl=1;
  end;
  else do;
  pos0x=(_par_[pos0]/(_par_[pos0]+yhat[pos0]))##_par_[pos0];
  pos0xl=(_par_[pos0]/(_par_[pos0]+y[pos0,]))##_par_[pos0];
  end;
  ll=sum(-log(0+exp(pihat[pos0]))+log(0*exp(pihat[pos0])+pos0x))+
  sum(-log(0+exp(pihat[pos1]))+lgamma(_par_[pos1]+y[pos1,])-lgamma(y[pos1,]+1)-lgamma(_par_[pos1])+
  y[pos1]#log(yhat[pos1]/(_par_[pos1]+yhat[pos1]))+
  _par_[pos1]#log(_par_[pos1]/(_par_[pos1]+yhat[pos1])));
  llnull1=sum(-log(1+zk)+log(zk+pos0xl))+
  sum(-log(1+zk)+lgamma(_par_[pos1]+y[pos1,])-lgamma(y[pos1,]+1)-lgamma(_par_[pos1])+
  y[pos1]#log(y[pos1,]/(_par_[pos1]+y[pos1,]))+_par_[pos1]#log(_par_[pos1]/(_par_[pos1]+y[pos1,])));
  dev=2*(llnull1-ll);
  npar=sum(S);
  AIC=2*npar-2*ll;
  AICc=AIC+2*(npar*(npar+1)/(n-npar-1));
  %IF %UPCASE(&model)=NEGBIN %THEN %DO;
  AIC=2*(npar+npar/ncol(x))-2*ll;
  AICc=AIC+2*(npar+npar/ncol(x))*(npar+npar/ncol(x)+1)/(n-(npar+npar/ncol(x))-1);
  %END;
  %END;
  
  %if %UPCASE(&bandwidth)=AIC %THEN %DO;CV=AICC;%END;
  %END;
  %ELSE %IF %UPCASE(&method)=ADAPTIVEN %THEN %DO;
  yhat[1]=x[i,]*b;
  CV=((y[i]-yhat)#wt)`*(y[i]-yhat);
  %END;
  
  free dist w;
  res=cv||npar;
  return (res);
}