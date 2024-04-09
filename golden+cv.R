Golden <- function(data, formula, xvarinf, weight,
                   lat, long, globalmin=TRUE,
                   method, model="zinb", bandwidth="cv", offset, 
                   force=FALSE, maxg=100, distancekm=FALSE){
  output <- list()
  E <- 10
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf)
  mt <- attr(mf, "terms")
  XVAR <- attr(mt, "term.labels")
  #y <<- model.extract(mf, "response")
  y <- model.extract(mf, "response")
  #N <<- length(y)
  N <- length(y)
  #x <<- model.matrix(mt, mf)
  x <- model.matrix(mt, mf)
  if (is.null(xvarinf)){
    #G <<- matrix(1, N, 1)
    G <- matrix(1, N, 1)
    # G <<- rep(1, N) #matriz coluna
    #lambdag <<- matrix(0, ncol(G), 1) #ncol(G) em vez de length(G)
    lambdag <- matrix(0, ncol(G), 1) #ncol(G) em vez de length(G)
    # lambdag <<- rep(0, length(G))
  }
  else{
    #G[colname=varnamezip]
    #G <<- unlist(data[, xvarinf])
    G <- unlist(data[, xvarinf])
    #G <<- cbind(rep(1, N), G)
    G <- cbind(rep(1, N), G)
  }
  #wt <<- matrix(1, N, 1)
  #wt <<- rep(1, N)
  wt <- rep(1, N)
  if (!is.null(weight)){
    #wt <<- unlist(data[, weight])
    wt <- unlist(data[, weight])
  }
  #Offset <- matrix(0, N, 1)
  #Offset <<- rep(0, N)
  Offset <- rep(0, N)
  if (!is.null(offset)){
    #Offset <<- unlist(data[, offset])
    Offset <- unlist(data[, offset])
  }
  # x <<- cbind(rep(1, N), x)
  #nvar <<- ncol(x)
  nvar <- ncol(x)
  #yhat <<- matrix(0, N, 1)
  #yhat <<- rep(0, N)
  yhat <- rep(0, N)
  #yhat2 <<- matrix(0, N, 1)
  #yhat2 <<- rep(0, N)
  yhat2 <- rep(0, N)
  #pihat <<- matrix(0, N, 1)
  #pihat <<- rep(0, N)
  pihat <- rep(0, N)
  #alphai <<- matrix(0, N, 1)
  #alphai <<- rep(0, N)
  alphai <- rep(0, N)
  #S <<- matrix(0, N, 1)
  #S <<- rep(0, N)
  S <- rep(0, N)
  #Si <<- matrix(0, N, 1)
  #Si <<- rep(0, N)
  Si <- rep(0, N)
  Iy <- ifelse(y>0, 1, y)
  Iy <- 1-Iy
  #pos0 <<- which(y==0)
  pos0 <- which(y==0)
  #pos02 <<- which(y==0)
  pos02 <- which(y==0)
  #pos1 <<- which(y>0)
  pos1 <- which(y>0)
  
  # /**** global estimates ****/
  uj <- (y+mean(y))/2 
  nj <- log(uj)
  #parg <<- sum((y-uj)^2/uj)/(N-nvar)
  parg <- sum((y-uj)^2/uj)/(N-nvar)
  ddpar <- 1
  cont <- 1
  cont3 <- 0
  while (abs(ddpar)>0.000001 & cont<100){
    dpar <- 1
    parold <- parg
    cont1 <- 1
    if (model == "zip" | model == "poisson"){
      #parg <<- 1/E^(-6)
      parg <- 1/E^(-6)
      alphag <- 1/parg
    }  
    if (model == "zinb" | model == "negbin"){
      if (cont>1){ 
        #parg <<- 1/(sum((y-uj)^2/uj)/(N-nvar))
        parg <- 1/(sum((y-uj)^2/uj)/(N-nvar))
      }  
      while (abs(dpar)>0.0001 & cont1<200){
        if (parg<0){
          #parg <<- 0.00001
          parg <- 0.00001
        }
        #parg <<- ifelse(parg<E^-10, E^-10, parg)
        parg <- ifelse(parg<E^-10, E^-10, parg)
        gf <- sum(digamma(parg+y)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+y)/(parg+uj))
        hess <- sum(trigamma(parg+y)-trigamma(parg)+1/parg-2/(parg+uj)+(y+parg)/(parg+uj)^2)
        hess <- ifelse(hess==0, E^-23, hess)
        par0 <- parg
        #parg <<- par0-solve(hess)%*%gf
        #parg <<- par0-as.vector(solve(hess))*gf
        parg <- par0-as.vector(solve(hess))*gf
        if (parg>E^5){
          dpar <- 0.0001
          cont3 <- cont3+1
          if (cont3==1){
            #parg <<- 2
            parg <- 2
          } 
          else if (cont3==2) {
            #parg <<- E^5
            parg <- E^5
          }
          else if (cont3==3){
            #parg <<- 0.0001
            parg <- 0.0001
          } 
        }
        else{
          dpar <- parg-par0
          cont1 <- cont1+1
        }
        if (parg>E^6){
          #parg <<- E^6
          parg <- E^6
          dpar <- 0
        }
      }
      #alphag <<- 1/parg
      alphag <- 1/parg
    }
    devg <- 0
    ddev <- 1
    cont2 <- 0
    while (abs(ddev)>0.000001 & cont2<100){
      Ai <- (uj/(1+alphag*uj))+(y-uj)*(alphag*uj/(1+2*alphag*uj+alphag^2*uj*uj))
      Ai <- ifelse(Ai<=0,E^-5,Ai)
      zj <- nj+(y-uj)/(Ai*(1+alphag*uj))-Offset
      if (det(t(x)%*%(Ai*x))==0) {
        # bg <<- matrix(0, ncol(x),1)
        #bg <<- rep(0,ncol(x))
        bg <- rep(0,ncol(x))
      } 
      else{
        #bg <<- solve(t(x)%*%(Ai*x))%*%t(x)%*%(Ai*zj)
        bg <- solve(t(x)%*%(Ai*x))%*%t(x)%*%(Ai*zj)
      }
      nj <- as.vector(x%*%bg+Offset)
      nj <- ifelse(nj>700,700,nj)
      uj <- exp(nj)
      olddev <- devg
      uj <- ifelse(uj<E^-150,E^-150,uj)
      uj <- ifelse(uj>100000,100000,uj)
      tt <- y/uj
      tt <- ifelse(tt==0,E^-10,tt)
      devg <- 2*sum(y*log(tt)-(y+1/alphag)*log((1+alphag*y)/(1+alphag*uj)))
      if (cont2>100){
        ddev <- 0.0000001
      } 
      else{
        ddev <- devg-olddev
      } 
      cont2 <- cont2+1
    }
    cont <-cont+1
    ddpar <- parg-parold
  }
  if (!is.null(xvarinf)){
    lambda0 <- (length(pos0)-sum((parg/(uj+parg))^parg))/N
    if (lambda0 <= 0) {
      #lambdag <<- matrix(0, ncol(G),1) 
      #lambdag <<- rep(0, ncol(G))
      lambdag <- rep(0, ncol(G))
    }
    else {
      #lambda0 <<- log(lambda0/(1-lambda0))
      lambda0 <- log(lambda0/(1-lambda0))
      #lambdag <<- rbind(lambda0, matrix(0,ncol(G)-1,1))
      #lambdag <<- rbind(lambda0, rep(0,ncol(G)-1))
      lambdag <- rbind(lambda0, rep(0,ncol(G)-1))
    }
    #pargg <<- parg
    pargg <- parg
    #ujg <<- uj
    ujg <- uj
    if (is.null(nrow(pos0)) | any(lambdag)==0){ 
      if (is.null(nrow(pos0))) {
        #pos0 <<- pos1
        pos0 <- pos1
        if (model=="zinb" | model=="zip"){
          model <- "negbin"
        }
      }
      if (!force){
        model <- "negbin"
      }
    }
  } #para o teste que estamos rodando, o SAS e o R nao entram aqui
  njl <- G%*%lambdag
  if (model!="zip" & model!="zinb"){
    zkg <- 0
  }
  else{
    zkg <- 1/(1+exp(-G%*%lambdag)*(parg/(parg+uj))^parg)
    zkg <- ifelse(y>0,0,zkg)
  }
  dllike <- 1
  llikeg <- 0
  j <- 0
  contador <- 0
  while (abs(dllike)>0.00001 & j<600){
    contador <- contador + 1
    ddpar <- 1
    cont <- 1
    contador2 <- 0
    while (abs(ddpar)>0.000001 & cont<100){
      contador2 <- contador2+1
      dpar <- 1
      parold <- parg
      aux1 <- 1
      aux2 <- 1
      aux3 <- 1
      cont3 <- 0
      int <- 1
      if (model=="zip" | model=="poisson"){
        alphag <- E^-6
        #parg <<- 1/alphag
        parg <- 1/alphag
      }
      else{
        if (j>0){
          #parg <<- 1/(sum((y-uj)^2/uj)/(N-nvar))
          parg <- 1/(sum((y-uj)^2/uj)/(N-nvar))
        }
        while (abs(dpar)>0.0001 & aux2<200){
          if (parg<0){
            #parg <<- 0.00001
            parg <- 0.00001
          }
          #parg <<- ifelse(parg<E^-10, E^-10, parg)
          parg <- ifelse(parg<E^-10, E^-10, parg)
          gf <- sum((1-zkg)*(digamma(parg+y)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+y)/(parg+uj)))
          hess <- sum((1-zkg)*(trigamma(parg+y)-trigamma(parg)+1/parg-2/(parg+uj)+(y+parg)/(parg+uj)^2))
          hess <- ifelse(hess==0, E^-23, hess)
          par0 <- parg
          #parg <<- as.vector(par0-solve(hess)%*%gf) #multiplicador
          parg <- as.vector(par0-solve(hess)%*%gf) #multiplicador
          if (aux2>50 & parg>E^5){
            dpar <- 0.0001
            cont3 <- cont3+1
            if (cont3==1){
              #parg <<- 2
              parg <- 2
            }
            else if (cont3==2){
              #parg <<- E^5
              parg <- E^5
            }
            else if (cont3==3){
              #parg <<- 0.0001
              parg <- 0.0001
            }
          }
          else{
            dpar <- parg-par0
          }
          if (parg>E^6){
            #parg <<- E^6
            parg <- E^6
            dpar <- 0
          }
          aux2 <- aux2+1
        }
        alphag <- 1/parg
      }
      devg <- 0
      ddev <- 1
      nj <- x%*%bg+Offset
      uj <- exp(nj)
      while (abs(ddev)>0.000001 & aux1<100){
        uj <- ifelse(uj>E^100, E^100, uj)
        Ai <- as.vector((1-zkg)*((uj/(1+alphag*uj)+(y-uj)*(alphag*uj/(1+2*alphag*uj+alphag^2*uj^2)))))
        Ai <- ifelse(Ai<=0, E^-5, Ai)
        uj <- ifelse(uj<E^-150, E^-150, uj)
        zj <- (nj+(y-uj)/(((uj/(1+alphag*uj)+(y-uj)*(alphag*uj/(1+2*alphag*uj+alphag^2*uj^2))))*(1+alphag*uj)))-Offset
        if (det(t(x)%*%(Ai*x))==0){
          #bg <<- matrix(0, nvar, 1)
          #bg <<- rep(0, nvar)
          bg <- rep(0, nvar)
        }
        else{
          #bg <<- solve(t(x)%*%(Ai*x))%*%t(x)%*%(Ai*zj)
          bg <- solve(t(x)%*%(Ai*x))%*%t(x)%*%(Ai*zj)
        }
        nj <- x%*%bg+Offset
        nj <- ifelse(nj>700, 700, nj)
        nj <- ifelse(nj<(-700), -700, nj)
        uj <- exp(nj)
        olddev <- devg
        gamma1 <- (uj/(uj+parg))^y*(parg/(uj+parg))^parg #(gamma(par+y)/(gamma(y+1)#gamma(par)))#
        gamma1 <- ifelse(gamma1<=0, E^-10, gamma1)
        devg <- sum((1-zkg)*(log(gamma1)))
        ddev <- devg-olddev
        #print(c('bg', 'aux1', 'devg', 'olddev', 'ddev'))
        #print(c(bg, aux1, devg, olddev, ddev))
        #prints comentados
        aux1 <- aux1+1
      }
      ddpar <- parg-parold
      cont <- cont+1
    }
    if (model == "zip" |model == 'zinb'){
      devg <- 0
      ddev <- 1
      njl <- G%*%lambdag
      njl <- ifelse(njl > maxg, maxg, njl)
      njl <- ifelse(njl < (-maxg),-maxg, njl)
      pig <- exp(njl)/(1+exp(njl))
      contador3 <- 0
      while (abs(ddev)>0.000001 & aux3<100){
        contador3 <- contador3 + 1
        Ai <- pig*(1-pig)
        Ai <- ifelse(Ai<=0, E^-5, Ai)
        zj <- njl+(zkg-pig)*1/Ai
        if (det(t(G*Ai)%*%G)==0){
          #lambdag <<- matrix(0, ncol(G), 1)
          lambdag <- matrix(0, ncol(G), 1)
          #lambdag <<- rep(0, ncol(G))
        }  
        else {
          #lambdag <<- solve(t(G*Ai)%*%G)%*%t(G*Ai)%*%zj
          lambdag <- solve(t(G*Ai)%*%G)%*%t(G*Ai)%*%zj
        }
        njl <- G%*%lambdag
        njl <- ifelse(njl > maxg, maxg, njl)
        njl <- ifelse(njl < (-maxg),-maxg, njl)
        pig <- exp(njl)/(1+exp(njl))
        olddev <- devg
        devg <- sum(zkg*njl-log(1+exp(njl)))
        ddev <- devg-olddev
        #print(c('lambdag', 'devg', 'olddev', 'ddev'))
        #print(c(lambdag, devg, olddev, ddev))
        #prints comentados
        aux3 <- aux3+1
      }
    }
    zkg <- 1/(1+exp(-njl)*(parg/(parg+uj))^parg)
    zkg <- ifelse(y>0, 0, zkg)
    if (model != 'zip' & model != 'zinb'){
      zkg <- 0
    } 
    oldllike <- llikeg
    llikeg <- sum(zkg*(njl)-log(1+exp(njl))+(1-zkg)*(log(gamma1)))
    dllike <- llikeg-oldllike
    #print(c('j', 'bg', 'alphag', 'lambdag', 'llikeg', 'dllike'))
    #print(c(j, bg, alphag, lambdag, llikeg, dllike))
    #prints comentados
    j <- j+1
  }
  long <- unlist(data[, long])
  lat <- unlist(data[, lat])
  #COORD <<- matrix(c(long, lat), ncol=2)
  COORD <- matrix(c(long, lat), ncol=2)
  # distance <<- dist(COORD,"euclidean") # _dist_ <- distance
  #sequ <<- 1:N # seq <- sequ
  sequ <- 1:N # seq <- sequ
  cv <- function(h){
    #N, wt, x, y, G, yhat, yhat2, pihat, hv, coord, _dist_, seq, offset, alphai, S, Si, parg, pargg, ujg, bg, lambdag, pos0, pos1, pos02, nvar  
    # if(method=="adaptiven"){
    #   #hv <- matrix(0,1,1)
    #   hv <- 0
    #   #yhat <- matrix(0,1,1)
    #   yhat <- 0
    #   # create &output from hv[colname='h'];
    # }
    for (i in 1:N){
      for (j in 1:N){
        #seqi <- matrix(i, N, 1)
        seqi <- rep(i, N)
        dx <- sp::spDistsN1(COORD,COORD[i,])
        distan <- cbind(seqi, sequ, dx)
        if (distancekm){
          distan[,3] <- distan[,3]*111
        }
      }
      u <- nrow(distan)
      #w <- matrix(0, u, 1)
      w <- rep(0, u)
      for (jj in 1:u){
        w[jj] <- exp(-0.5*(distan[jj,3]/h)^2)
        #if (i==172 & jj==1){print(-0.5*(distan[jj,3]/h)^2)}
        if (method=="fixed_bsq"){ # | method=="adaptiven"
          w[jj] <- (1-(distan[jj,3]/h)^2)^2
        }
        if (bandwidth=="cv"){
          w[i] <- 0
        }
      }
      #if (i==172){print(w)}
      if (method=="fixed_bsq"){ # | method=="adaptiven"
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
      njl <- G%*%lambda
      njl <- ifelse(njl>maxg, maxg, njl)
      njl <- ifelse(njl<(-maxg), -maxg, njl)
      if (model!="zip" & model!="zinb"){
        zk <- 0
      }
      else{
        lambda0 <- (length(pos0)-sum((parg/(uj+parg))^parg))/N
        if (lambda0>0){
          lambda0 <- log(lambda0/(1-lambda0))
          #lambda <- rbind(lambda0, matrix(0, ncol(G)-1, 1))
          lambda <- rbind(lambda0, rep(0, ncol(G)-1))
          njl <- G%*%lambda
        }
        zk <- 1/(1+exp(-njl)*(par/(par+uj))^par)
        zk <- ifelse(y>0, 0, zk)
      }
      dllike <- 1
      llike <- 0
      j <- 1
      contador4 <- 0
      while (abs(dllike)>0.00001 & j<=600){
        contador4 <- contador4+1
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
            njl <- G%*%lambda
            njl <- ifelse(njl>maxg, maxg, njl)
            njl <- ifelse(njl<(-maxg), -maxg, njl)
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
            #if (i==172 & contador4==1){print(w)}
            #if (i==172 & contador4==1){print(c(sum(w), sum(wt), sum(zk), par, sum(y), sum(uj)))}
            #w e par com problema --> origem em w
            #wt, zk, y e uj ok
            hess <- ifelse(hess==0, E^-23, hess)
            par0 <- par
            #print(c(i, contador4)) #aqui
            #if (i==172 & contador4==1){print(hess)} erro!
            par <- as.vector(par0-solve(hess)%*%gf) #multiplicador
            #par <- as.vector(par0-MASS::ginv(hess)%*%gf) #multiplicador
            dpar <- par-par0
            if (par>=E^6){
              par <- E^6
              dpar <- 0
              alpha <- 1/par
              b <- bg
              uj <- exp(x%*%b+Offset)
              lambda <- lambdag
              njl <- G%*%lambda
              njl <- ifelse(njl>maxg, maxg, njl)
              njl <- ifelse(njl<(-maxg),-maxg,njl)
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
            njl <- G%*%lambda
            njl <- ifelse(njl>maxg, maxg, njl)
            njl <- ifelse(njl<(-maxg), -maxg, njl)
            zk <- 1/(1+exp(-njl)*(parg/(parg+uj))^parg)
            zk <- ifelse(y>0, 0, zk)
            if (any(lambda)==0){
              zk <- 0
            }
          }
          alpha <- 1/par
          # if (i==244){
          #   print(c("i", "j", "b", "lambda", "par", "alpha", "aux2"))
          #   print(c(i, j, b, lambda, par, alpha, aux2))
          # }
          #prints comentados
        }
        dev <- 0
        ddev <- 1
        nj <- x%*%b+Offset
        nj <- ifelse(nj>700, 700, nj)
        nj <- ifelse(nj<(-700), -700, nj)
        uj <- exp(nj)
        contador5 <- 0
        while (abs(ddev)>0.000001 & aux1<100){
          contador5 <- contador5+1
          uj <- ifelse(uj>E^100,E^100,uj)
          Ai <- as.vector((1-zk)*((uj/(1+alpha*uj)+(y-uj)*(alpha*uj/(1+2*alpha*uj+alpha^2*uj^2)))))
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
            b <- solve(t(x)%*%(w*Ai*x*wt), tol=E^-60)%*%t(x)%*%(w*Ai*wt*zj)
            #b <- MASS::ginv(t(x)%*%(w*Ai*x*wt))%*%t(x)%*%(w*Ai*wt*zj)
          }
          nj <- x%*%b+Offset
          nj <- ifelse(nj>700, 700, nj)
          nj <- ifelse(nj<(-700), -700, nj)
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
          # if (i==244){
          #   print(c("b", "par", "aux1", "dev", "olddev", "ddev"))
          #   print(c(b, par, aux1, dev, olddev, ddev))
          # }
          #prints comentados
          aux1 <- aux1+1
        }
        ddpar <- par-parold
        if (model=="zip" | model=="zinb"){
          if (j==1){
            alphatemp <- alpha
            lambdatemp <- lambda[1] 
          }
          else{
            alphatemp <- c(alphatemp, alpha)
            lambdatemp <- c(lambdatemp, lambda[1])
          }
          alphatemp <- round(alphatemp, 7)
          ambdatemp <- round(lambdatemp, 7)
          # if (i==244){
          #   print(c('i', 'j', 'alphatemp', 'length(alphatemp)', 'length(unique(alphatemp))'))
          #   print(c(i, j, alphatemp, length(alphatemp), length(unique(alphatemp))))
          # }
          #prints comentados
          if (model=="zinb"){
            condition <- (j>300 & length(alphatemp)>length(unique(alphatemp)) & length(lambdatemp)>length(unique(lambdatemp)))
          }
          else if (model=="zip"){
            #print('i', 'j', 'lambdatemp', '(length(lambdatemp))', '(length(unique(lambdatemp)))')
            #print(i, j, lambdatemp, (length(lambdatemp)), (length(unique(lambdatemp))))
            #prints comentados
            condition <- (j>300 & length(lambdatemp)>length(unique(lambdatemp)))
          }
          if (condition){
            #lambda <- matrix(0, ncol(G), 1)
            lambda <- rep(0, ncol(G))
            njl <- G%*%lambda
            #zk <- matrix(0, n, 1)
            zk <- rep(0, N)
          }
          else{
            aux3 <- 1
            dev <- 0
            ddev <- 1
            njl <- G%*%lambda
            njl <- ifelse(njl>maxg, maxg, njl)
            njl <- ifelse(njl<(-maxg), -maxg, njl)
            pi <- exp(njl)/(1+exp(njl))
            contador6 <- 0
            while (abs(ddev)>0.000001 & aux3<100){
              contador6 <- contador6+1
              Aii <- pi*(1-pi)
              Aii <- ifelse(Aii<=0, E^-5, Aii)	
              zj <- njl+(zk-pi)/Aii
              if (det(t(G*Aii*w*wt)%*%G)==0){ #multiplicador
                lambda <- matrix(0, ncol(G), 1)
              }
              else{
                lambda <- solve(t(G*Aii*w*wt)%*%G)%*%t(G*Aii*w*wt)%*%zj
                #lambda <- MASS::ginv(t(G*Aii*w*wt)%*%G)%*%t(G*Aii*w*wt)%*%zj
              }
              njl <- G%*%lambda
              njl <- ifelse(njl>maxg, maxg, njl)
              njl <- ifelse(njl<(-maxg), -maxg, njl)
              pi <- exp(njl)/(1+exp(njl))
              olddev <- dev
              dev <- sum(zk*njl-log(1+exp(njl)))
              ddev <- dev-olddev
              # if (i==244){
              #   print(c("lambda", "aux3", "dev", "olddev", "ddev"))
              #   print(c(lambda, aux3, dev, olddev, ddev))
              # }
              #prints comentados
              aux3 <- aux3+1
            }
          }
        }
        njl <- G%*%lambda
        njl <- ifelse(njl>maxg, maxg, njl)
        njl <- ifelse(njl<(-maxg), -maxg, njl)
        zk <- 1/(1+exp(-njl)*(par/(par+uj))^par)
        zk <- ifelse(y>0, 0, zk)
        if (any(lambda)==0){
          #zk <- matrix(0, n, 1)
          zk <- rep(0, N)
        }
        if (model!="zip" & model!="zinb"){
          zk <- 0
        }
        oldllike <- llike
        llike <- sum(zk*(njl)-log(1+exp(njl))+(1-zk)*(log(gamma1)))
        if(i==244){
        }
        dllike <- llike-oldllike
        # if (i==244){
        #   print(c("i", "j", "b", "alpha", "lambda", "llike", "dllike"))
        #   print(c(i, j, b, alpha, lambda, llike, dllike))
        # }
        #prints comentados
        j <- j+1
      }
      yhat[i] <- uj[i]
      pihat[i] <- njl[i]
      #alphai[i] <<-  alpha
      alphai_ <- get("alphai")
      alphai_[i] <- alpha
      assign("alphai", alphai_, envir=parent.frame())
      alphai[i] <-  alpha
      if (det(t(x)%*%(w*Ai*x*wt))==0){
        S[i] <- 0
      }
      else {
        #S[i] <- (x[i,]*solve(t(x)%*%(w*Ai*x*wt))%*%t(x*w*Ai*wt))[i]
        S[i] <- (x[i,]%*%solve(t(x)%*%(w*Ai*x*wt), tol=E^-60)%*%t(x*w*Ai*wt))[i]
        #S[i] <- (x[i,]%*%MASS::ginv(t(x)%*%(w*Ai*x*wt))%*%t(x*w*Ai*wt))[i]
      }
      if(model=="zip" | model=="zinb"){
        yhat[i] <- (uj*(1-exp(njl)/(1+exp(njl))))[i]
        yhat2[i] <- uj[i]
        if (det(t(G)%*%(w*Aii*G*wt))==0){
          Si[i] <- 0
        }
        else{
          Si[i] <- (G[i,]%*%solve(t(G)%*%(w*Aii*G*wt))%*%t(G*w*Aii*wt))[i]
          if (any(lambda)==0){
            Si[i] <- 0
          }
        }
      }
      # if(i==1){
      #   max_dist <<- max(dx) 
      # }
      # max_dist <<- max(max_dist,max(dx))
    }
    CV <- t((y-yhat)*wt)%*%(y-yhat)
    par_ <- 1/alphai
    if(model=="zinb" | model == "zip"){
      if (any(lambda)==0){
        ll <- sum(-log(0+exp(pihat[pos0]))+log(0*exp(pihat[pos0])+(par_[pos0]/(par_[pos0]+yhat2[pos0]))^par_[pos0]))+
          sum(-log(0+exp(pihat[pos1]))+lgamma(par_[pos1]+y[pos1])-lgamma(y[pos1]+1)-lgamma(par_[pos1])+
                y[pos1]*log(yhat2[pos1]/(par_[pos1]+yhat2[pos1]))+par_[pos1]*log(par_[pos1]/(par_[pos1]+yhat2[pos1]))) #flag -> y vetor
        
        llnull1 <- sum(-log(1+zk[pos0])+log(zk[pos0]+(par_[pos0]/(par_[pos0]+y[pos0]))^par_[pos0]))+
          sum(-log(1+zk[pos1])+lgamma(par_[pos1]+y[pos1])-lgamma(y[pos1]+1)-lgamma(par_[pos1])+
                y[pos1]*log(y[pos1]/(par_[pos1]+y[pos1]))+par_[pos1]*log(par_[pos1]/(par_[pos1]+y[pos1]))) #flag -> y vetor
        
        llnull2 <- sum(-log(1+0)+log(0+(par_/(par_+mean(y)))^par_))+
          sum(-log(1+0)+lgamma(par_+y)-lgamma(y+1)-lgamma(par_)+y*log(mean(y)/(par_+mean(y)))+par_*log(par_/(par_+mean(y)))) 
      }
      else{
        ll <- sum(-log(1+exp(pihat[pos0]))+log(exp(pihat[pos0])+(par_[pos0]/(par_[pos0]+yhat2[pos0]))^par_[pos0]))+
          sum(-log(1+exp(pihat[pos1]))+lgamma(par_[pos1]+y[pos1])-lgamma(y[pos1]+1)-lgamma(par_[pos1])+
                y[pos1]*log(yhat2[pos1]/(par_[pos1]+yhat2[pos1]))+par_[pos1]*log(par_[pos1]/(par_[pos1]+yhat2[pos1])))
        
        llnull1 <- sum(-log(1+zk[pos0])+log(zk[pos0]+(par_[pos0]/(par_[pos0]+y[pos0]))^par_[pos0]))+
          sum(-log(1+zk[pos1])+lgamma(par_[pos1]+y[pos1])-lgamma(y[pos1]+1)-lgamma(par_[pos1])+
                y[pos1]*log(y[pos1]/(par_[pos1]+y[pos1]))+par_[pos1]*log(par_[pos1]/(par_[pos1]+y[pos1])))
      }
      dev <- 2*(llnull1-ll)
      npar <- sum(S)+sum(Si)
      AIC <- 2*npar-2*ll
      AICc <- AIC+2*(npar*(npar+1)/(N-npar-1))
      if(model=="zinb"){
        AIC <- 2*(npar+npar/(ncol(x)+ncol(G)))-2*ll
        AICc <- AIC+2*((npar+npar/(ncol(x)+ncol(G)))*((npar+npar/(ncol(x)+ncol(G)))+1)/(N-(npar+npar/(ncol(x)+ncol(G)))-1))
      }
    }
    else if(model=="poisson" | model=="negbin"){
      if (length(pos02)==0){ #flag ncol -> length
        pos0 <- pos1
        pos0x <- 1
        pos0xl <- 1
      } 
      else {
        pos0x <- (par_[pos0]/(par_[pos0]+yhat[pos0]))^par_[pos0]
        pos0xl <- (par_[pos0]/(par_[pos0]+y[pos0]))^par_[pos0] #flag -> y vetor
      }
      ll <- sum(-log(0+exp(pihat[pos0]))+log(0*exp(pihat[pos0])+pos0x))+
        sum(-log(0+exp(pihat[pos1]))+lgamma(par_[pos1]+y[pos1])-lgamma(y[pos1]+1)-lgamma(par_[pos1])+
              y[pos1]*log(yhat[pos1]/(par_[pos1]+yhat[pos1]))+
              par_[pos1]*log(par_[pos1]/(par_[pos1]+yhat[pos1]))) #flag -> y vetor
      
      llnull1 <- sum(-log(1+zk)+log(zk+pos0xl))+ sum(-log(1+zk)+lgamma(par_[pos1]+y[pos1])-lgamma(y[pos1]+1)-lgamma(par_[pos1])+
                                                       y[pos1]*log(y[pos1]/(par_[pos1]+y[pos1]))+par_[pos1]*log(par_[pos1]/(par_[pos1]+y[pos1]))) #flag -> y vetor
      dev <- 2*(llnull1-ll)
      npar <- sum(S)
      AIC <- 2*npar-2*ll
      AICc <- AIC+2*(npar*(npar+1)/(N-npar-1))
      if(model == "negbin"){
        AIC <- 2*(npar+npar/ncol(x))-2*ll
        AICc <-AIC+2*(npar+npar/ncol(x))*(npar+npar/ncol(x)+1)/(N-(npar+npar/ncol(x))-1)
      }
    }
    if(bandwidth == "aic"){
      CV <- AICc
    }
    res <- cbind(CV, npar)
    return (res)
  }
  
  # DEFINING GOLDEN SECTION SEARCH PARAMETERS 
  if(method=="fixed_g" | method=="fixed_bsq"){
    ax <- 0
    bx <- as.integer(max(dist(COORD))+1)
    #print(bx)
    if(distancekm){ #flag --> estava distacekm=="yes"
      bx <- bx*111
    }
  }  
  if(method=="adaptive_bsq"){
    ax <- 5
    bx <- N
  }
  r <- 0.61803399
  tol <- 0.1
  if (!globalmin){
    lower <- ax
    upper <- bx
    xmin <- matrix(0, 1, 2)
    #xmin <- rep(0, 2)
  }
  #do GMY=1 to 1;
  #ax1=lower[GMY];
  #bx1=upper[GMY];
  else{
    lower <- cbind(ax, (1-r)*bx, r*bx)
    upper <- cbind((1-r)*bx, r*bx, bx)
    xmin <- matrix(0, 3, 2)
  }
  #print(upper) #errado
  for (GMY in 1:3){
    ax1 <- lower[GMY]
    bx1 <- upper[GMY]
    h0 <- ax1
    h3 <- bx1
    h1 <- bx1-r*(bx1-ax1)
    #print(c(bx1, r, ax1)) #bx1 errado, o resto ok
    h2 <- ax1+r*(bx1-ax1)
    # print(c('h0', 'h1', 'h2', 'h3'))
    # print(c(h0, h1, h2, h3))
    if (GMY==1){
      h_values <- data.frame('h0'=h0, 'h1'=h1, 'h2'=h2, 'h3'=h3) #flag saída
    }
    else{
      h_values <- rbind(h_values, c(h0, h1, h2, h3))
    }
    ################################
    #print(h1) #errado
    res1 <- cv(h1)
    CV1 <- res1[1]
    res2 <- cv(h2)
    CV2 <- res2[1]
    if (GMY==1){
      band <- data.frame('GSS_count'=GMY, 'h1'=h1, 'cv1'=CV1, 'h2'=h2, 'cv2'=CV2) #flag saída
    }
    else{
      band <- rbind(band, c(GMY, h1, CV1, h2, CV2))
    }
    int <- 1
    while(abs(h3-h0) > tol*(abs(h1)+abs(h2)) & int<200){
      if (CV2<CV1){
        h0 <- h1
        h1 <- h3-r*(h3-h0)
        h2 <- h0+r*(h3-h0)
        CV1 <- CV2
        res2 <- cv(h2)
        CV2 <- res2[1]
      }
      else{
        h3 <- h2
        h1 <- h3-r*(h3-h0)
        h2 <- h0+r*(h3-h0)
        CV2 <- CV1
        res1 <- cv(h1)
        CV1 <- res1[1]
      }
      band <- rbind(band, c(GMY, h1, CV1, h2, CV2))
      int <- int+1
    }
    if (CV1<CV2){
      golden <- CV1
      xmin[GMY,1] <- golden
      xmin[GMY,2] <- h1
      npar <- res1[2]
      if (method=="adaptive_bsq"){
        xmin[GMY,2] <- floor(h1)
      }
    }
    else{
      golden <- CV2
      xmin[GMY,1] <- golden
      xmin[GMY,2] <- h2
      npar <- res2[2]
      if (method=="adaptive_bsq"){
        xmin[GMY,2] <- floor(h2)
      }
    }
    if (bandwidth=="aic"){
      # print(c('golden', 'xmin', 'npar'))
      # print(c(golden, xmin[GMY,2], npar))
      # if (GMY==1){
      #   gss_results <- data.frame('golden'=golden, 'xmin'=xmin[GMY,2], 'npar'=npar) #flag saída
      # }
      # else{
      #   gss_results <- rbind(gss_results, c(golden, xmin[GMY,2], npar))
      # }
    }
    else{
      # print(c('golden', 'xmin'))
      # print(c(golden, xmin[GMY,2]))
      # if (GMY==1){
      #   gss_results <- data.frame('golden'=golden, 'xmin'=xmin[GMY,2]) #flag saída
      # }
      # else{
      #   gss_results <- rbind(gss_results, c(golden, xmin[GMY,2]))
      # }
    }
    if (!globalmin){
      break
    }
  }
  min_bandwidth <- as.data.frame(xmin)
  names(min_bandwidth) <- c(bandwidth, 'bandwidth') #flag saída
  output <- append(output, list(h_values))
  names(output)[length(output)] <- "h_values"
  # output <- append(output, list(gss_results))
  # names(output)[length(output)] <- "gss_results"
  output <- append(output, list(band))
  names(output)[length(output)] <- "iterations"
  output <- append(output, list(min_bandwidth))
  names(output)[length(output)] <- "gss_results" #troca de min_bandwidth para gss_results
  if (globalmin){
    # print(c('golden', 'bandwidth'))
    # print(xmin)
    # xming <- xmin[which(xmin[,1]==min(xmin[,1])), 2]
    # print('(Da Silva and Mendes, 2018)')
    # print(c('Global Minimum', xming))
    message('Global Minimum (Da Silva and Mendes, 2018)') #flag saída
  }
  hh <- min_bandwidth[which(unlist(min_bandwidth[, bandwidth])==min(unlist(min_bandwidth[, bandwidth]))), 'bandwidth'] #flag saída
  output <- append(output, list(hh))
  names(output)[length(output)] <- "min_bandwidth"
  message('Bandwidth: ', hh) #flag saída
  # print('band')
  # print(band)
  # print('min_bandwidth')
  # print(min_bandwidth)
  return(output)
} # fecha golden
