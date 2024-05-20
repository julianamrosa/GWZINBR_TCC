gwzinbr <- function(data, formula, xvarinf, weight=NULL,
                    lat, long, grid=NULL, method, model = "zinb",
                    offset=NULL, distancekm=FALSE, force=TRUE, int_inf=TRUE,
                    maxg=100, h=NULL){
  
  E <- 10
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf)
  mt <- attr(mf, "terms")
  XVAR <- attr(mt, "term.labels")
  y <<- model.extract(mf, "response")
  N <<- length(y) 
  x <<- model.matrix(mt, mf)
  if (is.null(xvarinf)){
    G <<- matrix(1, N, 1) #certo
    # G <<- rep(1, N) #matriz coluna
    lambdag <<- matrix(0, ncol(G), 1) #ncol(G) em vez de length(G)
    # lambdag <<- rep(0, length(G))
  }
  else{
    #G <<- unlist(data[, xvarinf])
    G <<- as.matrix(data[, xvarinf])
  }
  if (int_inf){ #o que e int_inf? 
    G <<- cbind(rep(1, N), G)
  }
  # x <<- cbind(rep(1, N), x)
  yhat <<- rep(0, N)
  yhat2 <<- rep(0, N)
  pihat <<- rep(0, N)
  nvar <<- ncol(x)
  wt <<- rep(1, N)
  if (!is.null(weight)){
    wt <<- unlist(data[, weight])
  }
  Offset <<- rep(0, N)
  if (!is.null(offset)){
    Offset <<- unlist(data[, offset])
  }
  Iy <- ifelse(y>0, 1, y)
  Iy <- 1-Iy
  Iy2 <- Iy
  pos0 <<-  which(y==0)
  pos0 <<- t(as.matrix(pos0))
  pos02 <<- which(y==0)
  pos02 <<- t(as.matrix(pos02))
  pos1 <<- which(y>0)
  pos1 <<- t(as.matrix(pos1))
  
  # /**** global estimates ****/ 
  uj <- (y+mean(y))/2 
  nj <- log(uj)
  parg <<- sum((y-uj)^2/uj)/(N-nvar)
  ddpar <- 1
  cont <- 1
  cont3 <- 0
  while (abs(ddpar)>0.000001 & cont<100){
    dpar <- 1
    parold <- parg
    cont1 <- 1
    if (model == "zip" | model == "poisson"){
      parg <<- 1/E^(-6)
      alphag <- 1/parg
    }  
    if (model == "zinb" | model == "negbin"){
      if (cont>1){ 
        parg <<- 1/(sum((y-uj)^2/uj)/(N-nvar))
      }  
      while (abs(dpar)>0.0001 & cont1<200){
        if (parg<0){
          parg <<- 0.00001
        }
        parg <<- ifelse(parg<E^-10, E^-10, parg)
        gf <- sum(digamma(parg+y)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+y)/(parg+uj))
        hessg <- sum(trigamma(parg+y)-trigamma(parg)+1/parg-2/(parg+uj)+(y+parg)/(parg+uj)^2)
        hessg <- ifelse(hessg==0, E^-23, hessg)
        par0 <- parg
        #parg <<- par0-solve(hess)%*%gf
        parg <<- par0-as.vector(solve(hessg))*gf
        if (parg>E^5){
          dpar <- 0.0001
          cont3 <- cont3+1
          if (cont3==1){
            parg <<- 2
          } 
          else if (cont3==2) {
            parg <<- E^5
          }
          else if (cont3==3){
            parg <<- 0.0001  
          } 
        }
        else{
          dpar <- parg-par0
          cont1 <- cont1+1
        }
        if (parg>E^6){
          parg <<- E^6
          dpar <- 0
        }
      }
      alphag <<- 1/parg
    }
    devg <- 0
    ddev <- 1
    cont2 <- 0
    while (abs(ddev)>0.000001 & cont2<200){
      Ai <- (uj/(1+alphag*uj))+(y-uj)*(alphag*uj/(1+2*alphag*uj+alphag^2*uj*uj))
      Ai <- ifelse(Ai<=0,E^-5,Ai)
      zj <- nj+(y-uj)/(Ai*(1+alphag*uj))-Offset
      if (det(t(x)%*%(Ai*x))==0) {
        # bg <<- matrix(0, ncol(x),1)
        bg <<- rep(0,ncol(x))
      } 
      else{
        bg <<- solve(t(x)%*%(Ai*x))%*%t(x)%*%(Ai*zj)
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
    cont <- cont+1
    ddpar <- parg-parold
  }
  if (!is.null(xvarinf)){
    lambda0 <- (length(pos0)-sum((parg/(uj+parg))^parg))/N
    if (lambda0 <= 0) {
      #lambdag <<- matrix(0, ncol(G),1) 
      lambdag <<- rep(0, ncol(G))
      cat("NOTE: Expected number of zeros (", round(sum((parg/(uj+parg))^parg), 2), 
          ") >= number of zeros (", ncol(pos0), "). No Need of Zero Model.\n")
    }
    else{
      cat("NOTE: Expected number of zeros (", round(sum((parg/(uj+parg))^parg), 2), 
          ") < number of zeros (", ncol(pos0), "). Zero Model Used.\n")
      lambda0 <<- log(lambda0/(1-lambda0))
      lambdag <<- rbind(lambda0, rep(0, ncol(G)-1))
    }
    pargg <<- parg
    ujg <<- uj
    if (is.null(nrow(pos0)) | any(lambdag)==0){ 
      if (is.null(nrow(pos0))) {
        pos0 <<- pos1
        if (model=="zinb" | model=="zip"){
          model <- "negbin"
        }
      }
      if (!force){
        model <- "negbin"
      }
    }
  }
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
  while (abs(dllike)>0.00001 & j<600){
    ddpar <- 1
    cont <- 1
    while (abs(ddpar)>0.000001 & cont<100){
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
        while (abs(dpar)>0.0001 & aux2<200){
          if (parg<0){
            #parg <<- 0.00001
            parg <- 0.00001
          }
          parg <- ifelse(parg<E^-10, E^-10, parg)
          gf <- sum((1-zkg)*(digamma(parg+y)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+y)/(parg+uj)))
          hessg <- sum((1-zkg)*(trigamma(parg+y)-trigamma(parg)+1/parg-2/(parg+uj)+(y+parg)/(parg+uj)^2))
          hessg <- ifelse(hessg==0, E^-23, hessg)
          par0 <- parg
          parg <- as.vector(par0-solve(hessg)%*%gf)
          if ( aux2 > 50 |parg > E^5) {
            dpar <- 0.0001
            cont3 <- cont3+1
            if (cont3==1) {
              parg <- 2
            } 
            else if (cont3==2) {
              parg <- E^5
            } 
            else if (cont3==3) {
              parg <- 0.0001
            }
          } 
          else {
            dpar <- parg-par0
          }
          if (parg>E^6){
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
      nj <- ifelse(nj>700, 700, nj)
      nj <- ifelse(nj<(-700), -700, nj)
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
    if (model == "zip" | model == "zinb"){
      devg <- 0
      ddev <- 1
      njl <- G%*%lambdag
      njl <- ifelse(njl > maxg, maxg, njl)
      njl <- ifelse(njl < (-maxg),-maxg, njl)
      pig <- exp(njl)/(1+exp(njl))
      while (abs(ddev)>0.000001 & aux3<100){ 
        Ai <- as.vector(pig*(1-pig))
        Ai <- ifelse(Ai<=0, E^-5, Ai)
        zj <- njl+(zkg-pig)*1/Ai
        if(det(t(G)%*%(Ai*G))==0){ 
          lambdag <- matrix(0, ncol(G), 1)
        }  
        else {
          #lambdag <<- solve(t(G*Ai)%*%G)%*%t(G*Ai)%*%zj
          lambdag <- solve(t(G)%*%(Ai*G))%*%t(G)%*%(Ai*zj)
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
    j <- j+1
  }
  g1x <- parg/(parg+uj)
  g2x <- uj/(parg+uj)
  hgx <- exp(njl)+g1x^parg
  daa <- zkg*((g1x^parg*(log(g1x)+g2x))^2*(1-1/hgx)/hgx+g1x^parg*(g2x^2/parg)/hgx)+(1-zkg)*(trigamma(parg+y)-trigamma(parg)-2/(uj+parg)+1/parg+(y+parg)/(uj+parg))^2
  dab <- zkg*(g1x^(2*parg+1)*uj*(log(g1x)+g2x)/hgx^2-g1x^parg*(-g2x^2+parg*g2x*(log(g1x)+g2x))/hgx)+(1-zkg)*(g2x*(y-uj)/(uj+parg))
  dal <- -zkg*(exp(njl)*g1x^parg*(log(g1x)+g2x)/hgx^2)
  daa <- daa*parg^4
  dab <- dab*parg^2
  if (any(lambdag)==0){
    Iy <- matrix(0, length(y), 1)
  }
  dll <- as.vector(Iy*(exp(njl)*g1x^parg/hgx^2)-exp(njl)/(1+exp(njl))^2)
  dbb <- as.vector(Iy*(-(parg*g1x^parg*g2x/hgx)^2+parg^2*g1x^parg*g2x^2*(1-1/uj)/hgx)-(1-Iy)*(parg*g2x*(1+(y-uj)/(parg+uj))))
  dlb <- as.vector(Iy*(parg*exp(njl)*g1x^parg*g2x/hgx^2))
  I1 <- matrix(1, length(y), 1)
  # print("Dimensões de dbb:")
  # print(dim(dbb))
  # print("Dimensões de x:")
  # print(dim(x))
  # print(dbb)
  # II <- -(t(I1) %*% daa) %*% I1 || -(t(I1) %*% dab) %*% X || -(t(I1) %*% dal) %*% G // (-t(X) %*% (dab * I1) || -t(X) %*% dbb %*% X || -t(X) %*% dlb %*% G) // (-t(G) %*% (dal * I1) || -t(G) %*% (X %*% dlb) || -t(G) %*% t(G) %*% dll %*% G)
  # II=-(I1#daa)`*I1||-(I1#dab)`*X||-(I1#dal)`*G//(-X`*(dab#I1)||-(X#dbb)`*X||-(X#dlb)`*G)//(-G`*(dal#I1)||-G`*(X#dlb)||-(G#dll)`*G);
  II <- rbind(cbind(-(t(I1 * daa)) %*% I1, -(t(I1 * dab)) %*% x, -(t(I1 * dal)) %*% G), 
              cbind(-(t(x) %*% (dab * I1)), -t(x*dbb) %*% x, -t(x * dlb) %*% G), 
              cbind(-(t(G) %*% (dal * I1)), -t(G) %*% (x * dlb), -t(G * dll) %*% G))
  # conferir se II é matricial
  if (all(lambdag)>0 & alphag==E^-6){
    II <- II[2:nrow(II), 2:nrow(II)]
  }
  else if(any(lambdag)==0 & alphag>E^-6){
    II <- II[1:(ncol(x)+1), 1:(ncol(x)+1)]
  }
  else if(any(lambdag)== 0 & alphag==E^-6){
    II <- II[2:(ncol(x)+1), 2:(ncol(x)+1)]
  }
  varabetalambdag <- diag(solve(II))
  stdabetalambdag <- sqrt(abs(varabetalambdag))
  varcovg <- solve(II)
  # *print bg lambdag alphag devg llikeg stdabetalambdag;
  
  # /*****************************************/
  print(c("bandwidth", h))
  long <- unlist(data[, long])
  lat <- unlist(data[, lat])
  COORD <- matrix(c(long, lat), ncol=2)
  if (is.null(grid)){
    POINTS <- matrix(c(long, lat), ncol=2)
  }
  else{
    long2 <- unlist(grid[, long])
    lat2 <- unlist(grid[, lat])
    POINTS <- matrix(c(long2, lat2), ncol=2)
  }
  mm <- nrow(COORD) #substituicao: 'm' por 'mm', pois m representa nosso modelo no R
  bi <- matrix(0, ncol(x)*mm, 4)
  li <- matrix(0, ncol(G)*mm, 4)
  alphai <- matrix(0, mm, 3)
  BB <- matrix(0, ncol(x)*N, N)
  BB1 <- matrix(0, ncol(G)*N, N)
  sumwi <- matrix(0, mm, 1)
  varbi <- matrix(0, ncol(x)*mm, 1)
  varli <- matrix(0, ncol(G)* mm, 1)
  S <- matrix(0, mm, 1) #poderia substituir por rep()
  # S <- rep(0, mm)
  Si <- matrix(0, mm, 1)
  # Si <- rep(0, mm)
  S2 <- matrix(0, mm, 1)
  # S2 <- rep(0, mm)
  biT <- matrix(0, mm, ncol(x)+1)
  ym <- y-mean(y)
  
  # /******** calculating distance **********/;
  #substituicao: _dist_ por dist_
  sequ <- 1:N
  for (i in 1:mm){
    for (j in 1:N){
      seqi <- rep(i, N)
        dx <- sp::spDistsN1(COORD, COORD[i,])
        distan <- cbind(seqi, sequ, dx)
        if (distancekm){
          distan[,3] <- distan[,3]*111
        }  
     }
    u <- nrow(distan)
    w <- rep(0, u)
    for(jj in 1:u){
      if(method=="fixed_g"){
        w[jj] <- exp(-0.5*(dist[jj, 3]/h)^2)
      }
      else if(method=="fixed_bsq"){
        w[jj] <- (1 -(dist[jj, 3]/h)^2)^2
      }
    }
    if(method=="adaptive_bsq"){
      distan <- distan[order(distan[,3]),]
      distan <- cbind(distan, 1:nrow(distan))
      w <- matrix(0,N,2)
      hn <- distan[h,3]
      for (jj in 1:N) {
        if(distan[jj, 4] <= h){
          w[jj, 1] <- (1 -(distan[jj, 3]/hn)^2)^2
        } 
        else{
          w[jj, 1] <- 0
        }
        w[jj, 2] <- distan[jj, 2]
      }
      w <- w[order(w[, 2]), 1]
    }
    # if (method == "adaptive_bsq"){
    #   distan <- distan[order(distan[,3]),]
    #   distan <- cbind(distan, 1:N)
    #   w <- matrix(0,N,2)
    #   hn <- distan[h,3]
    #   
    #   for (jj in 1:N) {
    #     if(dist[jj, 4] <= h){
    #       w[jj, 1] <- (1 -(dist[jj, 3]/hn)^2)^2
    #     }
    #     else{
    #       w[jj, 1] <- 0
    #     }
    #     w[jj, 2] <- dist[jj, 2]
    #   }
    #   w <- w[order(w[, 2]), 1]
    # }
    # obs: essa condicional é igual à anterior
    # /****** MODEL SELECTION *************/
    Iy <- Iy2
    b <- bg
    b2 <- b
    nj <- x%*%b+Offset
    uj <- exp(nj)
    par <- parg
    par2 <- par
    lambda <- lambdag
    njl <- ifelse(njl > maxg, maxg, njl)
    njl <- ifelse(njl < (-maxg), -maxg, njl) 
    #retirada parte comentada do codigo (linhas 1535 a 1604 do SAS)
    if(model != "zip" & model != "zinb" ){
      zk <- 0
    }
    else{
      lambda0 <- (ncol(pos0)-sum((parg/(uj+parg))^parg))/N
      if (lambda0 > 0){
        lambda0 <- log(lambda0/(1-lambda0))
        lambda <- matrix(c(lambda0, rep(0, ncol(G)-1)), ncol=1)
        njl <- G%*%lambda
      }
      zk <- 1/(1+exp(-njl)*(par/(par+uj))^par)
      zk <- ifelse(y>0, 0, zk)
    }
    dllike <- 1
    llike <- 0
    j <- 1
    while (abs(dllike) > 0.00001 & j <= 600) { #fecha na linha 2049 do SAS? Confirmar com professor
      ddpar <- 1
      dpar <- 1
      parold <- par
      aux1 <- 1
      aux2 <- 1
      int <- 1
      if (model == "zip" || model == "poisson") {
        alpha <- E^-6
        par <- 1/alpha
      }
      if(model == "zinb" || model == "negbin"){
        #/*if par<=1E-5 then do;if i>1 then par=1/alphai[i-1,2];end;*/. 
        if(par >= E^6){
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
          zk <- ifelse(y>0, 0, zk)
          #zk[y>0] <- 0
          if (any(lambda) == 0){
            zk <- 0
          }
        }
        while (abs(dpar)>0.000001 & aux2<200){
          par <- ifelse(par < E^-10, E^-10, par)
          gf <- sum(w*wt*(1-zk)*(digamma(par+y)-digamma(par)+log(par)+1-log(par+uj) - (par+y)/(par+uj)))
          hess <- sum(w*wt*(1-zk)*(trigamma(par+y)-trigamma(par)+1/par - 2/(par+uj) + (y + par)/(par+uj))^2)
          hess <- ifelse(hess == 0, E^-23, hess)
          par0 <- par
          par <- as.vector(par0 - solve(hess)%*%gf)
          dpar <- par - par0
          
          if (par >= E^6){
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
            #zk[y > 0] <- 0
            zk <- ifelse(y>0, 0, zk)
            if (any(lambda) == 0){
              zk <- 0
            }
          }
          # print(par, aux2, dpar)
          aux2 <- aux2 + 1
        }
        if (par <= E^-5){
          par <- E^6
          b <- bg
          uj <- exp(x%*%b+Offset)
          lambda <- lambdag
          njl <- G %*% lambda
          njl <- ifelse(njl>maxg, maxg, njl)
          njl <- ifelse(njl<(-maxg),-maxg,njl)
          zk <- 1/(1+exp(-njl)*(parg/(parg+uj))^parg)
          zk <- ifelse(y>0, 0, zk)
          if (any(lambda) == 0){
            zk <- 0
          }
        }
        alpha <- 1/par
        # print(i, j, b, lambda, par, alpha, aux2)
      }
      dev <- 0
      ddev <- 1
      nj <- x%*%b+Offset
      nj <- ifelse(nj>700,700,nj)
      nj <- ifelse(nj<(-700), -700, nj)
      uj <- exp(nj)
      while (abs(ddev)>0.000001 & aux1<100){
        uj <- ifelse(uj>E^100, E^100, uj)
        Ai <- as.vector((1-zk)*((uj/(1+alpha*uj) + (y-uj) * (alpha*uj/(1+2*alpha*uj+alpha^2*uj^2)))))
        Ai <- ifelse(Ai <=0, E^-5, Ai)
        uj <- ifelse(uj < E^-150, E^-150, uj)
        denz <- (((uj/(1+alpha*uj)+(y-uj)*(alpha*uj/(1+2*alpha*uj+alpha^2*uj^2))))*(1+alpha*uj))
        denz <- ifelse(denz == 0, E^-5, denz)
        zj <- (nj+(y-uj)/denz)-Offset
        if(det(t(x) %*% (w*Ai*x*wt)) == 0){
          b <- matrix(0, nvar, 1)
        } 
        else {
          b <- solve(t(x)%*%(w*Ai*x*wt), tol=E^-60)%*%t(x)%*%(w*Ai*wt*zj)
        }
        nj <- x%*%b+Offset
        nj <- ifelse(nj>700,700,nj)
        nj <- ifelse(nj<(-700), -700, nj)
        uj <- exp(nj)
        olddev <- dev
        uj <- ifelse(uj > E^10, E^10, uj) 
        uj <- ifelse(uj == 0, E^10, uj) 
        if (par == E^6) {
          # gamma1=/*(gamma(par+y)/(gamma(y+1)#gamma(par)))#*/(uj/(uj+par))##y#exp(-uj);
          gamma1 <- (uj/(uj+par))^y*exp(-uj)
        } 
        else {
          # gamma1=/*(gamma(par+y)/(gamma(y+1)#gamma(par)))#*/(uj/(uj+par))##y#(par/(uj+par))##par;
          gamma1 <- (uj/(uj+par))^y*(par/(uj+par))^par
        }
        gamma1 <- ifelse(gamma1 <=0, E^-10, gamma1)
        dev <- sum((1-zk)*log(gamma1))
        ddev <- dev - olddev
        # print(b, par, aux1, dev, olddev, ddev)
        aux1 <- aux1 + 1
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
        lambdatemp <- round(lambdatemp, 7)
        if (model=="zinb"){
          condition <- (j>300 & length(alphatemp)>length(unique(alphatemp)) & length(lambdatemp)>length(unique(lambdatemp)))
        }
        else if (model=="zip"){
          #print('i', 'j', 'lambdatemp', '(length(lambdatemp))', '(length(unique(lambdatemp)))')
          #print(i, j, lambdatemp, (length(lambdatemp)), (length(unique(lambdatemp))))
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
          while (abs(ddev)>0.000001 & aux3<100){
            Aii <- as.vector(pi*(1-pi))
            Aii <- ifelse(Aii<=0, E^-5, Aii)	
            zj <- njl+(zk-pi)/Aii
            if (det(t(G*Aii*w*wt)%*%G)==0){ #multiplicador
              lambda <- matrix(0, ncol(G), 1)
            }
            else{
              lambda <- solve(t(G*Aii*w*wt)%*%G, tol=E^-60)%*%t(G*Aii*w*wt)%*%zj
            }
            njl <- G%*%lambda
            njl <- ifelse(njl>maxg, maxg, njl)
            njl <- ifelse(njl<(-maxg), -maxg, njl)
            pi <- exp(njl)/(1+exp(njl))
            olddev <- dev
            dev <- sum(zk*njl-log(1+exp(njl)))
            ddev <- dev-olddev
            #   print(c("lambda", "aux3", "dev", "olddev", "ddev"))
            #   print(c(lambda, aux3, dev, olddev, ddev))
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
      dllike <- llike-oldllike
      #   print(c("i", "j", "b", "alpha", "lambda", "llike", "dllike"))
      #   print(c(i, j, b, alpha, lambda, llike, dllike))
      j <- j+1
    }
    # /**** COMPUTING VARIANCE OF BETA AND LAMBDA ******/
    if(det(t(x)%*%(w*Ai*x*wt))==0){
      C <- matrix(0, ncol(x),nrow(x)) 
    } else{
      C <- solve(t(x)%*%(w*Ai*x*wt), tol=E^-60)%*%t(x)%*%(w*Ai*wt)
    }
    Ci <- matrix(0, ncol(G), 1)
    if(model == "zip" || model == "zinb"){
      if(det(t(G)%*%(w*Aii*G*wt))==0){
        Ci <- matrix(0, ncol(G), nrow(G))
      }
      else{
        Ci <- solve(t(G)%*%(w*Aii*G*wt), tol=E^-60)%*%t(G)%*%(w*Aii*wt)
      }
      if(any(lambda) == 0){
        Ci <- matrix(0, ncol(G), N)
      }
    }
    g1x <- par/(par+uj)
    g2x <- uj/(par+uj)
    hgx <- exp(njl)+g1x^par
    hgx <- ifelse(hgx > E^10, E^10, hgx)
    daa <- w*wt*(zk*((g1x^par*(log(g1x)+g2x))^2*(1 - 1/hgx)/hgx + g1x^par * (g2x^2/par) / hgx) + 
                   (1-zk) * (trigamma(par+y) - trigamma(par) - 2 / (uj+par) + 1 / par + (y+par) / (uj+par))^2)
    dab <- w*wt*(zk*(g1x^(2*par + 1)*uj*(log(g1x) + g2x) / hgx^2 - g1x^par * (-g2x^2+par*g2x*(log(g1x) + g2x)) / hgx) + 
                   (1-zk) * (g2x*(y-uj) / (uj+par)))
    dal <- -w*wt*zk*(exp(njl)*g1x^par*(log(g1x)+g2x) / hgx^2)
    daa <- daa*par^4
    dab <- dab*par^2
    if (any(lambda) == 0) {
      Iy <- matrix(0, nrow(y), 1)
      exphx <- 1 + exp(njl)
      exphx <- ifelse(exphx>E^90, E^90, exphx) 
      dll <- w*wt*(Iy*(exp(njl)*g1x^par/hgx^2) - exp(njl)/exphx)^2
      dbb <- sqrt(w)*wt*(Iy*(-(par*g1x^par*g2x/hgx)^2+par^2*g1x^par*g2x^2*(1 - 1/uj)/hgx) - 
                           (1 - Iy)*(par*g2x*(1 + (y-uj)/(par+uj))))
      dlb <- w*wt*Iy*(par*exp(njl)*g1x^par*g2x/hgx^2)
      dll <- ifelse(is.na(dll), E^100, dll)
      daa <- ifelse(is.na(daa), E^100, daa)
      dab <- ifelse(is.na(dab), E^100, dab)
      dal <- ifelse(is.na(dal), E^100, dal)
      dbb <- ifelse(is.na(dbb), E^100, dbb)
      dlb <- ifelse(is.na(dlb), E^100, dlb)
      I1 <- matrix(1, nrow(y), 1)
    }
    if(any(b) == 0 & any(lambda) == 0){
      II <- matrix(0, ncol(x)+ncol(G)+1, ncol(x)+ncol(G)+1)
    } else if(any(lambda) == 0){
      dbb <- w*wt*(Iy*(-(par*g1x^par*g2x/hgx)^2 + par^2*g1x^par*g2x^2*(1 - 1/uj)/hgx)-(1 - Iy)*(par*g2x*(1 + (y-uj) / (par+uj))))
      if(det(t(x%*%dbb%*%dbb/Ai) %*% x) == 0){
        II <- rbind(
          cbind(-(t(I1*daa))%*%I1, -(t(I1*dab))%*%x, -(t(I1*dal))%*%G),
          cbind(-t(x)%*%(dab*I1), -(t(x%*%dbb)%*%x), -t(x*dlb)%*%G),
          cbind(-t(G)%*%(dal*I1), -t(G)%*%(x*dlb), -(t(G*dll)%*%G))
        )   
      }
      else{
        II <- rbind(
          cbind(-t(I1*daa)%*%I1, -t(I1*dab)%*%x, -t(I1*dal)%*%G),
          cbind(-t(x)%*%(dab*I1), -(t(x*dbb)%*%x)%*%solve(t(x*dbb*dbb/Ai)%*%x)%*%t(x*dbb)%*%x, -(t(x*dlb)%*%G)),
          cbind(-t(G)%*%(dal*I1), -t(G)%*%(x*dlb), -t(G*dll)%*%G)      
        )   
      }
    } else{
      II <- rbind(
        cbind(-t(I1*daa)%*%I1, -t(I1*dab)%*%x, -t(I1*dal)%*%G),
        cbind(-t(x)%*%(dab*I1), -(t(x*dbb)%*%x), -t(x*dlb)%*%G),
        cbind(-t(G)%*%(dal*I1), -t(G)%*%(x*dlb), -t(G*dll)%*%G)
      )     
    }
    if (all(lambda) > 0 & alpha == E^-6) {
      II <- II[2:nrow(II), 2:nrow(II)]
    } else if (any(lambda) == 0 & alpha > E^-6) {
      II <- II[1:(ncol(x)+1), 1:(ncol(x)+1)]
    } else if (any(lambda) == 0 & alpha == E^-6) {
      II <- II[2:(ncol(x)+1), 2:(ncol(x)+1)]
    }
    if (det(II) == 0) {
      if (all(lambda) > 0 & alpha == E^-6) {
        II <- II[1:ncol(x), 1:ncol(x)]
        if (det(II) == 0) {
          varabetalambda <- rbind(matrix(0, nrow(II), 1), matrix(0, ncol(G), 1))
        } 
        else {
          varabetalambda <- rbind(diag(solve(II, tol=E^-60)), matrix(0, ncol(G), 1))
        }
      } 
      else if (any(lambda) == 0 & alpha == E^-6) {
        II <- II[1:ncol(x), 1:ncol(x)]
        if (det(II) == 0) {
          varabetalambda <- rbind(matrix(0, nrow(II), 1), matrix(0, ncol(G), 1))
        } 
        else {
          varabetalambda <- rbind(diag(solve(II,tol=E^-60)), matrix(0, ncol(G), 1))
        }
      } 
      else {
        II <- II[1:(ncol(x)+1), 1:(ncol(x)+1)]
        if (det(II) == 0) {
          varabetalambda <- rbind(matrix(0, nrow(II), 1), matrix(0, ncol(G), 1))
        } 
        else {
          varabetalambda <- rbind(diag(solve(II, tol=E^-60)), matrix(0, ncol(G), 1))
        }
      }
    } 
    else {
      varabetalambda <- diag(1/solve(II, tol=E^-60))
    }
    if (all(lambda) > 0 & alpha > E^-6){
      varb <- varabetalambda[2:(ncol(x)+1)]
      varl <- varabetalambda[(ncol(x)+2):length(varabetalambda)]
      alphai[i, 1] <- i
      alphai[i, 2] <- alpha
      alphai[i, 3] <- sqrt(abs(varabetalambda[1]))
    } 
    else if (all(lambda) > 0 & alpha == E^-6) {
      varb <- varabetalambda[1:ncol(x)]
      varl <- varabetalambda[(ncol(x)+1):length(varabetalambda)]
      alphai[i, 1] <- i
      alphai[i, 2] <- alpha
      alphai[i, 3] <- sqrt(1/abs(-(t(I1*daa))%*%I1))
    } else if (any(lambda) == 0 & alpha > E^-6) {
      varb <- varabetalambda[2:(ncol(x)+1)]
      varl <- matrix(0, ncol(G), 1)
      alphai[i, 1] <- i
      alphai[i, 2] <- alpha
      alphai[i, 3] <- sqrt(abs(varabetalambda[1]))
    } else if (any(lambda) == 0 & alpha == E^-6) {
      varb <- varabetalambda[1:ncol(x)]
      varl <- matrix(0, ncol(G), 1)
      alphai[i, 1] <- i
      alphai[i, 2] <- alpha
      alphai[i, 3] <- sqrt(1/abs(-(t(I1*daa))%*%I1))
    }
    
    # /*******************************/
    m1 <- (i-1)*ncol(x)+1
    m2 <- m1+(ncol(x)-1)
    bi[m1:m2, 1] <- i
    bi[m1:m2, 2] <- b
    bi[m1:m2, 3] <- COORD[i, 1]
    bi[m1:m2, 4] <- COORD[i, 2]
    varbi[m1:m2, 1] <- varb
    
    if (model == "zip" || model == "zinb") {
      m1 <- (i-1)*ncol(G)+1
      m2 <- m1 + (ncol(G)-1)
      li[m1:m2, 1] <- i
      li[m1:m2, 2] <- lambda
      li[m1:m2, 3] <- COORD[i, 1]
      li[m1:m2, 4] <- COORD[i, 2]
      varli[m1:m2, 1] <- varl
    }
    if (is.null(grid)) {
      r <- x[i, ]%*%C 
      S[i] <- r[i]
      S2[i] <- r%*%t(r)
      yhat[i] <- uj[i]
      pihat[i] <- njl[i]
    }
    if (model == "zip" | model == "zinb") {
      ri <- G[i, ]%*%Ci
      Si[i] <- ri[i]
      yhat2[i] <- uj[i]
      yhat[i] <- (uj*(1-exp(njl)/(1+exp(njl))))[i]
    }
    # /** creating non-stationarity matrix **/
    if (method != "adaptive_bsq"){
      CCC <- cbind(x, w, wt)
      m1 <- (i-1)*ncol(x)+1
      m2 <- m1 + (ncol(x)-1)
      
      if(det(t(CCC[, 1:ncol(x)])%*%(CCC[, ncol(CCC)-1]*CCC[, 1:ncol(x)]*CCC[, ncol(CCC)])) == 0){
        BB[m1:m2, ] <- matrix(0, ncol(x), nrow(x))
      } else {
        BB[m1:m2, ] <- solve(t(CCC[, 1:ncol(x)])%*%(CCC[, ncol(CCC)-1]*CCC[, 1:ncol(x)]*CCC[, ncol(CCC)]))%*%t(CCC[, 1:ncol(x)])*(CCC[, ncol(CCC)-1]*CCC[, ncol(CCC)])
      }
      if (model == "zip" || model == "zinb"){
        CCCl <- cbind(G, w, wt)
        m1 <- (i-1)*ncol(G)+1
        m2 <- m1+(ncol(G)-1)
        
        if (det(t(CCCl[, 1:ncol(G)])%*%(CCCl[, ncol(CCCl)-1]*CCCl[, 1:ncol(G)]*CCCl[, ncol(CCCl)])) == 0) {
          BBl[m1:m2, ] <- matrix(0, ncol(G), nrow(G))
        } else {
          BBl[m1:m2, ] <- solve(t(CCCl[, 1:ncol(G)])%*%(CCCl[, ncol(CCCl)-1] * CCCl[, 1:ncol(G)]*CCCl[, ncol(CCCl)]))%*%t(CCCl[, 1:ncol(G)])*(CCCl[, ncol(CCCl)-1]*CCCl[, ncol(CCCl)])
        }
      }
    }
    # /*************************************/
    # substituicao: _w_ <- w_
    w_ <- w
    w_ <- w_[order(w_)]
    sumwi[i] <- sum(w_[1:min(length(w_), length(w_)*0.1)])
    
    if (i==1){
      W_f <- cbind(w, 1:length(w))
    }
    else{
      W_f <- rbind(W_f, cbind(w, 1:length(w)))
      # W_f <- as.data.frame(W_f)
      # View(W_f)
      # verificar se essas demais linhas fazem sentido 
    }
  }
  if (is.null(grid)) {
    v1 <- sum(S) + sum(Si)
    v11 <- sum(S) + sum(Si)
    v2 <- sum(S2)
    nparmodel <- n-v11
    if (v11 < v2) {
      v1 <- v11
    }
    res <- y-yhat
    rsqr1 <- t(res*wt)%*%res
    ym <- t(y*wt)%*%y
    rsqr2 <- ym-sum((y*wt)^2)/sum(wt)
    rsqr <- 1-rsqr1/rsqr2
    rsqradj <- 1-((n-1)/(n-v1))*(1-rsqr)
    sigma2 <- n*rsqr1/((n-v1)*sum(wt))
    root_mse <- sqrt(sigma2)
    print(c('Sigma2e', sigma2))
    print(c('Root MSE', root_mse))
    print(c('#GWR parameters', v1))
    print(c('#GWR parameters (model)', nparmodel))
    print(c('#GWR parameters (variance)', v2))
    #sera que todos esses prints fazem sentido? 
    
    influence <- S
    resstd <- res/(sqrt(sigma2)*sqrt(abs(1-influence)))
    CooksD <- resstd^2*influence/(v1*(1-influence))
    df <- n-(nvar+ncol(G))
    stdbi <- sqrt(abs(varbi))
    tstat <- bi[, 2]/stdbi
    probt <- 2*(1-pt(abs(tstat), df))
    malpha <- 0.05*(nvar/v1)
    #substituicai: _malpha_ <- malpha 
    
    if (model == "zip" || model == "zinb") {
      stdli <- sqrt(abs(varli))
      tstati <- li[, 2]
      
      for (j in 1:nrow(stdli)) {
        if (stdli[j] == 0) {
          tstati[j] <- 0
        } else {
          tstati[j] <- li[j, 2]/stdli[j]
        }
      }
      tstati <- ifelse(is.na(tstati), 0, tstati)
      probti <- 2*(1-pt(abs(tstati), df))
      malpha <- 0.05*((nvar+ncol(G))/v1)
    }
    t_critical <- abs(qt(malpha/2, df))
    par_ <- 1/alphai[, 2] #substituicao: _par_ por par_
    
    if (model == "zinb" || model == "zip") {
      if (any(lambda) == 0) {
        ll <- sum(-log(0+exp(pihat[pos0])) + log(0*exp(pihat[pos0]) + (par_[pos0]/(par_[pos0]+yhat2[pos0]))^par_[pos0])) +
          sum(-log(0+exp(pihat[pos1])) + lgamma(par_[pos1]+y[pos1, ]) - lgamma(y[pos1, ]+1) - lgamma(par_[pos1]) +
                y[pos1]*log(yhat2[pos1]/(par_[pos1]+yhat2[pos1]))+par_[pos1] * log(par_[pos1]/(par_[pos1]+yhat2[pos1])))
        
        llnull1 <- sum(-log(1+zk[pos0]) + log(zk[pos0]+(par_[pos0]/(par_[pos0] + y[pos0, ]))^par_[pos0])) +
          sum(-log(1+zk[pos1])+lgamma(par_[pos1] + y[pos1, ])-lgamma(y[pos1, ]+1) - lgamma(par_[pos1]) +
                y[pos1]*log(y[pos1, ]/(par_[pos1] + y[pos1, ])) + par_[pos1]*log(par_[pos1]/(par_[pos1] + y[pos1, ])))
        
        llnull2 <- sum(-log(1+0) + log(0+(par_/(par_ + mean(y)))^par_)) + #pq log de 1+0? 
          sum(-log(1+0) + lgamma(par_+y) - lgamma(y+1) - lgamma(par_) +
                y*log(mean(y)/(par_ + mean(y))) + par_*log(par_/(par_+mean(y))))
      } else {
        ll <- sum(-log(1+exp(pihat[pos0])) + log(exp(pihat[pos0]) + (par_[pos0]/(par_[pos0]+yhat2[pos0]))^par_[pos0])) +
          sum(-log(1+exp(pihat[pos1])) + lgamma(par_[pos1]+y[pos1]) - lgamma(y[pos1]+1) - lgamma(par_[pos1]) +
                y[pos1]*log(yhat2[pos1]/(par_[pos1]+yhat2[pos1])) + par_[pos1]*log(par_[pos1]/(par_[pos1]+yhat2[pos1])))
        
        llnull1 <- sum(-log(1+zk[pos0]) + log(zk[pos0] + (par_[pos0]/(par_[pos0] + y[pos0]))^par_[pos0])) +
          sum(-log(1 + zk[pos1]) + lgamma(par_[pos1] + y[pos1]) - lgamma(y[pos1]+1) - lgamma(par_[pos1]) +
                y[pos1]*log(y[pos1]/(par_[pos1]+y[pos1])) + par_[pos1]*log(par_[pos1]/(par_[pos1] + y[pos1])))
        
        llnull2 <- sum(-log(1+0) + log(0+(mean(par_)/(mean(par_)+mean(y)))^mean(par_))) +
          sum(-log(1+0) + lgamma(mean(par_)+y) - lgamma(y+1) - lgamma(mean(par_)) +
                y*log(mean(y)/(mean(par_)+mean(y))) + mean(par_)*log(mean(par_)/(mean(par_) + mean(y))))
      }
      dev <- 2*(llnull1-ll)
      pctll <- 1-(llnull1-ll)/(llnull1-llnull2)
      AIC <- 2*v1-2*ll
      AICc <- AIC+2*(v1*(v1+1)/(n-v1-1))
      adjpctll <- 1-(llnull1-ll+v1+0.5)/(llnull1-llnull2)
      
      if(model == "zinb"){
        AIC <- 2*(v1+v1/(ncol(x)+ncol(G)))-2*ll
        AICc <- AIC + 2*((v1+v1/(ncol(x)+ncol(G)))*((v1+v1/(ncol(x)+ncol(G)))+1)/ (n-(v1+v1/(ncol(x)+ncol(G))) - 1))
        adjpctll <- 1 - (llnull1-ll+(v1+v1/(ncol(x) + ncol(G))) + 0.5) / (llnull1 - llnull2)
      }
      print(paste("Deviance:", dev))
      print(paste("Full Log Likelihood:", ll))
      print(paste("Pseudo R-squared:", pctll))
      print(paste("Adjusted Pseudo R-squared:", adjpctll))
      print(paste("AIC:", AIC))
      print(paste("AICc:", AICc))
    }
    if(model == "poisson" || model == "negbin"){
      if (ncol(pos02) == 0){ #acredito que podemos substituir ncol por length, mas prefiro manter assim por enquanto
        pos0 <- pos1
        pos0x <- rep(1, length(pos1))
        pos0xl <- rep(1, length(pos1))
      } else {
        pos0x <- (par_[pos0]/(par_[pos0] + yhat[pos0]))^par_[pos0]
        pos0xl <- (par_[pos0]/(par_[pos0] + y[pos0, ]))^par_[pos0]
      }
      ll <- sum(-log(0+exp(pihat[pos0])) + log(0*exp(pihat[pos0]) + pos0x)) + sum(-log(0+exp(pihat[pos1])) +
                                                                                    lgamma(par_[pos1] + y[pos1, ]) - lgamma(y[pos1, ] + 1) - lgamma(par_[pos1]) + y[pos1]*log(yhat[pos1]/(par_[pos1] + yhat[pos1])) +
                                                                                    par_[pos1]*log(par_[pos1]/(par_[pos1] + yhat[pos1])))
      
      llnull1 <- sum(-log(1+zk) + log(zk+pos0xl)) + sum(-log(1+zk) + 
                                                          lgamma(par_[pos1]+y[pos1, ]) - lgamma(y[pos1, ]+1) - lgamma(par_[pos1])+y[pos1]*log(y[pos1, ]/(par_[pos1]+y[pos1, ])) +
                                                          par_[pos1]*log(par_[pos1]/(par_[pos1] + y[pos1, ])))
      
      pos1xx <- 0+(par_[pos0]/(par_[pos0] + mean(y)))*par_[pos0]
      pos1xx <- ifelse(pos0x <= 0, E^-100, pos0x)
      
      llnull2 <- sum(-log(1+0) + log(pos1xx)) +
        sum(-log(1+0) + lgamma(par_[pos1]+y[pos1]) - lgamma(y[pos1]+1) - lgamma(par_[pos1]) +
              y[pos1]*log(mean(y)/(par_[pos1]+ mean(y))) + par_[pos1]*log(par_[pos1]/(par_[pos1]+ mean(y))))
      
      dev <- 2*(llnull1-ll)
      AIC <- -2*ll+2*v1
      AICc <- AIC+2*(v1*(v1+1)/(n-v1-1))
      pctll <- 1-(llnull1-ll)/(llnull1-llnull2)
      adjpctll <- 1-(llnull1-ll+v1+0.5)/(llnull1-llnull2)
      
      if(model == "negbin"){
        AIC <- 2*(v1+v1/ncol(x))-2*ll
        AICc <- AIC+2*(v1+v1/ncol(x))*(v1+v1/ncol(x)+1)/(n-(v1+v1/ncol(x))-1)
        adjpctll <- 1-(llnull1-ll+v1+v1/ncol(x)+0.5)/(llnull1-llnull2)
      }
      print(paste("Deviance:", dev))
      print(paste("Full Log Likelihood:", ll))
      print(paste("Pctll:", pctll))
      print(paste("AdjPctll:", adjpctll))
      print(paste("AIC:", AIC))
      print(paste("AICc:", AICc))
    }
    # substituicoes: _beta_ por beta_ ; _beta2_ por beta2_
    beta_ <- matrix(bi[, 1:2], nrow = N, byrow=TRUE) #shape: cria matriz a partir de outra matriz. Conferir se resultado eh o mesmo no SAS
    beta2_ <- beta_
    if(model == "negbin" || model == "zinb"){
      alpha_= matrix(alphai[, 1:2], n);
      beta2_= cbind(beta_, alpha_)
    }
    i <- seq(2, ncol(beta_), 2)  
    beta_ <- beta_[, i]         
    i <- seq(2, ncol(beta2_), 2) 
    beta2_ <- beta2_[, i]
    qntl <- apply(beta2_, 2, quantile, c(0.25, 0.5, 0.75))
    IQR <- qntl[3, ]-qntl[1, ]
    qntl <- rbind(qntl, IQR)
    descriptb <- rbind(apply(beta2_, 2, "mean"), apply(beta2_, 2, "min"), apply(beta2_, 2, "max"))
    rownames(descriptb) <- c('Mean', 'Min', 'Max')
    if (model=='negbin' | model=="zinb"){
      colnames(qntl) <- c('Intercept', XVAR, 'alpha')
      colnames(descriptb) <- c('Intercept', XVAR, 'alpha')
    }
    else{
      colnames(qntl) <- c('Intercept', XVAR)
      colnames(descriptb) <- c('Intercept', XVAR)
    }
    print("Quantiles of GWR Parameter Estimates")
    print(qntl)
    print("Descriptive Statistics")
    print(descriptb)
    stdbeta_ <- matrix(stdbi, N, byrow=TRUE)
    stdbeta2_ <- stdbeta_
    if (model=="negbin" | model=="zinb"){
      stdalpha_ <- matrix(alphai[, 3], N, byrow=TRUE)
      stdbeta2_ <- cbind(stdbeta_, stdalpha_)
    }
    qntls <- apply(stdbeta2_, 2, quantile, c(0.25, 0.5, 0.75))
    IQR <- qntls[3, ]-qntls[1, ]
    qntls <- rbind(qntls, IQR)
    descripts <- rbind(apply(stdbeta2_, 2, "mean"), apply(stdbeta2_, 2, "min"), apply(stdbeta2_, 2, "max"))
    rownames(descripts) <- c('Mean', 'Min', 'Max')
    if (model=='negbin' | model=="zinb"){
      colnames(qntls) <- c('Intercept', XVAR, 'alpha')
      colnames(descripts) <- c('Intercept', XVAR, 'alpha')
    }
    else{
      colnames(qntls) <- c('Intercept', XVAR)
      colnames(descripts) <- c('Intercept', XVAR)
    }
    print(c("alpha-level=0.05", malpha_))
    print(c("t-Critical", t_critical))
    print(c("degrees of freedom:", df))
    print("Quantiles of GWR Standard Errors")
    print(qntls)
    print("Descriptive Statistics of Standard Errors")
    print(descripts)
    if (model=="zip" | model=="zinb"){
      lambda_ <- matrix(li[, 1:2], N, byrow=TRUE)
      lambda2_ <- lambda_
      i <- seq(2, ncol(lambda_), 2)
      lambda_ <- lambda_[, i]
      i <- seq(2, ncol(lambda2_), 2)
      lambda2_ <- lambda2_[, i]
      qntl <- apply(lambda2_, 2, quantile, c(0.25, 0.5, 0.75))
      IQR <- qntl[3, ]-qntl[1, ]
      qntl <- rbind(qntl, IQR)
      descriptl <- rbind(apply(lambda2_, 2, "mean"), apply(lambda2_, 2, "min"), apply(lambda2_, 2, "max"))
      rownames(descriptl) <- c('Mean', 'Min', 'Max')
      if (int_inf){
        colnames(qntl) <- c('Intercept', xvarinf)
        colnames(descriptl) <- c('Intercept', xvarinf)
      }
      else{
        colnames(qntl) <- c(xvarinf)
        colnames(descriptl) <- c(xvarinf)
      }
      print("Quantiles of GWR Zero Inflation Parameter Estimates")
      print(qntl)
      print("Descriptive Statistics of GWR Zero Inflation Parameter Estimates")
      print(descriptl)
      stdlambda_ <- matrix(stdli, N, byrow=TRUE)
      stdlambda2_ <- stdlambda_
      qntls <- apply(stdlambda2_, 2, quantile, c(0.25, 0.5, 0.75))
      IQR <- qntls[3, ]-qntls[1, ]
      qntls <- rbind(qntls, IQR)
      descriptls <- rbind(apply(stdlambda2_, 2, "mean"), apply(stdlambda2_, 2, "min"), apply(stdlambda2_, 2, "max"))
      rownames(descriptls) <- c('Mean', 'Min', 'Max')
      if (int_inf){
        colnames(qntls) <- c('Intercept', xvarinf)
        colnames(descriptls) <- c('Intercept', xvarinf)
      }
      else{
        colnames(qntls) <- c(xvarinf)
        colnames(descriptls) <- c(xvarinf)
      }
      ##de novo esses prints?##
      print(c("alpha-level=0.05", malpha_))
      print(c("t-Critical", t_critical))
      print(c("degrees of freedom:", df))
      #####
      print("Quantiles of GWR Zero Inflation Standard Errors")
      print(qntls)
      print("Descriptive Statistics of GWR Zero Inflation Standard Errors")
      print(descriptls)
    } #models
  } #grid
  #### Non-Stationarity Test ####
  if (is.null(grid)){
    if(method!="adaptive_bsq"){
      BBk <- matrix(0, N, N)
      Vk <- matrix(0, ncol(x), 1)
      df1k <- matrix(0, ncol(x), 1)
      df2k <- matrix(0, ncol(x), 1)
      for (k in 1:ncol(x)){
        ek <- matrix(0, ncol(x), 1)
        ek[k] <- 1
        for (i in 1:N){
          m1 <- (i-1)*ncol(x)+1
          m2 <- m1+(ncol(x)-1)
          BBk[i, ] <- t(ek)%*%BB[m1:m2, ]
        }
        Vk[k] <- t(y)%*%(1/N)%*%t(BBk)%*%(diag(N)-(1/N)%*%matrix(1, N, N))%*%BBk%*%y
        df1k[k] <- sum(diag((1/N)%*%t(BBk)%*%(diag(N)-(1/N)%*%matrix(1, N, N))%*%BBk))
        df2k[k] <- sum(diag(((1/N)%*%t(BBk)%*%(diag(N)-(1/N)%*%matrix(1, N, N))%*%BBk)^2))
      }
      Vk <- ifelse(abs(Vk)<=E^-8, 0, Vk)
      Fk <- (Vk/df1k)/sigma2
      ndf <- df1k^2/df2k
      ddf <- n-v1
      ddf <- rep(ddf, ncol(x))
      probf <- 1-pf(Fk, ndf, ddf)
      print("Non-Stationarity Test (Leung et al., 2000)")
      rownames(Vk) <- c('Intercept', XVAR)
      colnames(Vk) <- "V"
      print(Vk)
      colnames(Fk) <- "F"
      print(Fk)
      print(c(ndf, ddf))
      print(c("Pr > F", probf))
      if (model=="zip" | model=="zinb"){
        BBkl <- matrix(0, N, N)
        Vkl <- matrix(0, ncol(G), 1)
        df1kl <- matrix(0, ncol(G), 1)
        df2kl <- matrix(0, ncol(G), 1)
        for (k in 1:ncol(G)){
          ekl <- matrix(0, ncol(G), 1)
          ekl[k] <- 1
          for (i in 1:N){
            m1 <- (i-1)*ncol(G)+1
            m2 <- m1+(ncol(G)-1)
            BBkl[i, ] <- t(ekl)%*%BBl[m1:m2, ]
          }
          Vkl[k] <- t(y)%*%(1/N)%*%t(BBkl)%*%(diag(N)-(1/N)%*%matrix(1, N, N))%*%BBkl%*%y
          df1kl[k] <- sum(diag((1/N)*t(BBkl)%*%(diag(N)-(1/N)%*%matrix(1, N, N))%*%BBkl))
          df2kl[k] <- sum(diag(((1/N)%*%t(BBkl)%*%(diag(N)-(1/N)%*%matrix(1, N, N))%*%BBkl)^2))
        }
        Vkl <- ifelse(abs(Vkl)<=E^-8, 0, Vkl)
        Fkl <- (Vkl/df1kl)/sigma2
        ndfl <- df1kl^2/df2kl
        ddfl <- N-v1
        ddfl <- rep(ddfl, ncol(G))
        probfl <- 1-pf(Fkl, ndfl, ddfl)
        print("Non-Stationarity Test (Leung et al., 2000) - Zero Inflation")
        if (int_inf){
          rownames(Vkl) <- c('Intercept', xvarinf)
        }
        else{
          rownames(Vkl) <- xvarinf
        }
        colnames(Vkl) <- "V"
        colnames(Fkl) <- "F"
        print(Vkl)
        print(Fkl)
        print(c(ndfl, ddfl))
        print(c("Pr > F", probfl))
      }
    }
  }
  #### global estimates ####
  if (is.null(weight)){
    dfg <-N-nvar
  }
  if (model=="zinb"){
    b2 <- rbind(bg, alphag)
    if (alphag==E^-6){
      stdg <- rbind(stdabetalambdag[1:ncol(x)], (sqrt(1/abs(hessg))/(parg^2)))
    }
    else{
      stdg <- rbind(stdabetalambdag[2:ncol(x)+1], stdabetalambdag[1])
    }
    tg <- b2/stdg
    dfg <- nrow(y)-ncol(x)
    probtg <- 2*(1-pt(abs(tg), dfg))
    lambdag <- lambdag
    if (alphag==E^-6){
      stdlambdag <- stdabetalambdag[ncol(x)+1:nrow(stdabetalambdag)]
    }
    else{
      stdlambdag <- stdabetalambdag[ncol(x)+2:nrow(stdabetalambdag)]
    }
    tlambdag <- lambdag/stdlambdag
    dflg <- nrow(y)-ncol(G)
    probtlambdag <- 2*(1-pt(abs(tlambdag), dflg))
    p <- ncol(x)+1+ncol(G)
  }
  if (model=="zip"){
    b2 <- bg
    stdg <- stdabetalambdag[1:ncol(x)]
    tg <- b2/stdg
    dfg <- nrow(y)-ncol(x)
    probtg <- 2*(1-pt(abs(tg), dfg))
    lambdag <- lambdag
    stdlambdag <- stdabetalambdag[ncol(x)+1:nrow(stdabetalambdag)]
    tlambdag <- lambdag/stdlambdag
    dflg <- nrow(y)-ncol(G)
    probtlambdag <- 2*(1-pt(abs(tlambdag), dflg))
    p <- ncol(x)+ncol(G)
  }
  if (model=="negbin"){
    b2 <- rbind(bg, alphag)
    if(alphag==E^-6){
      stdg <- rbind(stdabetalambdag[1:nrow(stdabetalambdag)], (sqrt(1/abs(hessg))/(parg^2)))
    }
    else{
      stdg <- rbind(stdabetalambdag[2:nrow(stdabetalambdag)], stdabetalambdag[1])
    }
    tg <- b2/stdg
    dfg <- nrow(y)-ncol(x)
    probtg <- 2*(1-pt(abs(tg), dfg))
    p <- ncol(x)+1
  }
  if (model=="poisson"){
    b2 <- bg
    stdg <- stdabetalambdag
    tg <- b2/stdg
    dfg <- nrow(y)-ncol(x)
    probtg <- 2*(1-pt(abs(tg), dfg))
    p <- ncol(x)
  }
  bg_stdg <- cbind(b2, stdg)
  print("Global Parameter Estimates")
  if (model=="negbin" | model=="zinb"){
    rownames(bg_stdg) <- c('Intercept', XVAR, 'alpha')
  }
  else{
    rownames(bg_stdg) <- c('Intercept', XVAR)
  }
  colnames(bg_stdg) <- c("Par. Est.", "Std Error")
  print(c("t Value", bg_stdg))
  print(c("Pr > |t|", probtg))
  print("NOTE: The denominator degrees of freedom for the t tests is")
  print(dfg)
  if (model=="zip" | model=="zinb"){
    print("Analysis Of Maximum Likelihood Zero Inflation Parameter Estimate")
    if (!is.null(xvarinf)){
      varnamezip1 <- t(varnamezip)
      if (int_inf){
        varnamezip1 <- rbind("Intercept", t(varnamezip))
      }
    }
    else{
      if(int_inf){
        varnamezip1 <- "Intercept"
      }
    }
    print(c("Parameter", varnamezip1))
    print(c("Estimate", lambdag))
    print(c("Standard Error", stdlambdag))
    print(c("t Value", tlambdag))
    print(c("Pr > |t|", probtlambdag))
    print("NOTE: The denominator degrees of freedom for the t tests is")
    print(dflg)
    ll <- sum(-log(1+exp(G[pos0, ]%*%lambdag))+
                log(exp(G[pos0, ]%*%lambdag)+
                      (parg/(parg+exp(x[pos0, ]%*%bg+Offset[pos0, ])))^parg))+
      sum(-log(1+exp(G[pos1, ]%*%lambdag))+
            lgamma(parg+y[pos1, ])-lgamma(y[pos1, ]+1)-
            lgamma(parg)+y[pos1]*log(exp(x[pos1, ]%*%bg+Offset[pos1, ])/(parg+exp(x[pos1, ]%*%bg+Offset[pos1, ])))+
            parg*log(parg/(parg+exp(x[pos1, ]%*%bg+Offset[pos1, ]))))
    AIC <- 2*p-2*ll
    AICc <- AIC+2*(p*(p+1)/(N-p-1))
    llnull1 <- sum(-log(1+zkg[pos0, ])+log(zkg[pos0, ]+
                                             (parg/(parg+y[pos0, ]))^parg))+
      sum(-log(1+zkg[pos1, ])+lgamma(parg+y[pos1, ])-
            lgamma(y[pos1, ]+1)-lgamma(parg)+
            y[pos1]*log(y[pos1, ]/(parg+y[pos1, ]))+
            parg*log(parg/(parg+y[pos1, ])))
    llnull2 <- sum(-log(1+0)+log(0+(parg/(parg+mean(y)))^parg))+
      sum(-log(1+0)+lgamma(parg+y)-lgamma(y+1)-lgamma(parg)+
            y*log(mean(y)/(parg+mean(y)))+parg*log(parg/(parg+mean(y))))
    devg <- 2*(llnull1-ll)
    pctll <- 1-(llnull1-ll)/(llnull1-llnull2)
    adjpctll <- 1-(llnull1-ll+p+0.5)/(llnull1-llnull2)
    print('Deviance', devg)
    print('Full Log Likelihood', ll)
    print(pctll, adjpctll, AIC, AICc)
  }
  if (model=="poisson" | model=="negbin"){
    yhatg <- exp(x%*%bg+Offset)
    if (ncol(pos02)==0){
      pos0 <- pos1
      pos0x <- 1
      pos0xl <- 1
    }
    else{
      pos0x <- (parg/(parg+exp(x[pos0, ]%*%bg+Offset[pos0, ])))^parg
      pos0xl <- (parg/(parg+y[pos0, ]))^parg
    }
    ll <- sum(-log(0+exp(G[pos0, ]%*%lambdag))+
                log(0*exp(G[pos0, ]%*%lambdag)+pos0x))+
      sum(-log(0+exp(G[pos1, ]%*%lambdag))+
            lgamma(parg+y[pos1, ])-lgamma(y[pos1, ]+1)-
            lgamma(parg)+y[pos1]*log(exp(x[pos1, ]%*%bg+
                                           Offset[pos1, ])/(parg+exp(x[pos1, ]%*%bg+
                                                                       Offset[pos1, ])))+
            parg*log(parg/(parg+exp(x[pos1, ]%*%bg+Offset[pos1, ]))))
    llnull1 <- sum(-log(1+zkg)+log(zkg+pos0xl))+
      sum(-log(1+zkg)+lgamma(parg+y[pos1, ])-
            lgamma(y[pos1, ]+1)-lgamma(parg)+
            y[pos1]*log(y[pos1, ]/(parg+y[pos1, ]))+
            parg*log(parg/(parg+y[pos1, ])))
    devg <- 2*(llnull1-ll)
    AIC <- -2*ll+2*nvar
    AICc <- -2*ll+2*nvar*(N/(N-nvar-1))
    tt2 <- y/mean(y)
    tt2 <- ifelse(tt2==0, E^-10, tt2)
    devnullg <- 2*sum(y*log(tt2)-(y-mean(y)))
    pctdevg <- 1-devg/devnullg
    adjpctdevg <- 1-((N-1)/(N-nvar))*(1-pctdevg)
    if (model=="negbin"){
      AIC <- -2*ll+2*(nvar+1)
      AICc <- -2*ll+2*(nvar+1)*(N/(N-(nvar+1)-1))
      tt2 <- y/mean(y)
      tt2 <- ifelse(tt2==0, E^-10, tt2)
      devnullg <- 2*sum(y*log(tt2)-(y+1/alphag)*log((1+alphag*y)/(1+alphag*mean(y))))
      pctdevg <- 1-devg/devnullg
      adjpctdevg <- 1-((N-1)/(N-(nvar+1)))*(1-pctdevg)
    }
    print(c('Deviance', devg))
    print(c('Full Log Likelihood', ll))
    print(c(pctdevg, adjpctdevg, AIC, AICc))
  }
  print('Variance-Covariance Matrix')
  print(varcovg)
  ##################
  if (is.null(grid)){
    pihat <- exp(pihat)/(1+exp(pihat))
    res_ <- cbind(wt, y, yhat, res, resstd, influence, cooksD, sumwi, pihat)
    print(c("wt", "y", "yhat", "res", "resstd", "influence", "cooksD", "sumwi", "pihat"))
    print(res_)
    beta_ <- bi
    colnames(beta_) <- c("id", "B", "x", "y")
    print(beta_)
    bistdt <- cbind(bi, stdbi, tstat, probt)
    parameters_ <- bistdt
    colnames(parameters_) <- c("id", "B", "x", "y", "stdbi", "tstat", "probt")
    tstat_ <- beta_
    for (j in 1:nrow(stdbeta_)){
      for (k in 1:ncol(stdbeta_)){
        if (stdbeta_[j, k]==0){
          tstat_[j, k] <- 0
        }
        else{
          tstat_[j, k] <- beta_[j, k]/stdbeta_[j, k]
        }
      }
    }
    tstat_ <- ifelse(is.na(tstat_), 0, tstat_)
    probt_ <- 2*(1-pt(abs(tstat_), df))
    sig_ <- matrix("not significant at 90%", N, ncol(x))
    for (i in 1:N){
      for (j in 1:ncol(x)){
        if (probt_[i, j]<0.01*(nvar/v1)){
          sig_[i, j] <- "significant at 99%"
        }
        else if (probt_[i, j]<0.05*(nvar/v1)){
          sig_[i, j] <- "significant at 95%"
        }
        else if (probt_[i, j]<0.1*(nvar/v1)){
          sig_[i, j] <- "significant at 90%"
        }
        else{
          sig_[i, j] <- "not significant at 90%"
        }
        if (is.na(probt_[i, j])){
          sig_[i, j] <- "not significant at 90%"
        }
      }
    }
    bistdt_ <- cbind(COORD, beta_, stdbeta_, tstat_, probt_)
    colname1 <- c("Intercept", XVAR)
    label_ <- rbind(rep("std_", ncol(x)), rep("tstat_", ncol(x)), rep("probt_", ncol(x)))
    colname <- cbind(c("x", "y"), colname1, c(label_, rep(colname1, 3)))
    colname <- sub("_ ", "_", colname)
    colname <- sub("_ ", "_", colname)
    if (model=="zip" | model=="zinb"){
      tstatl_ <- lambda_
      for (j in 1:nrow(stdlambda_)){
        for(k in 1:ncol(stdlambda_)){
          if (stdlambda_[j, k]==0){
            tstatl_[j, k] <- 0
          }
          else{
            tstatl_[j, k] <- lambda_[j, k]/stdlambda_[j, k]
          }
        }
      }
      probtl_ <- 2*(1-pt(abs(tstatl_), df))
      sigl_ <- matrix("not significant at 90%", N, ncol(G))
      for (i in 1:N){
        for (j in 1:ncol(G)){
          if (probtl_[i, j]<0.01*((nvar+ncol(G))/v1)){
            sigl_[i, j] <- "significant at 99%"
          }
          else if (probtl_[i, j]<0.05*((nvar+ncol(G))/v1)){
            sigl_[i, j] <- "significant at 95%"
          }
          else if (probtl_[i, j]<0.1*((nvar+ncol(G))/v1)){
            sigl_[i, j] <- "significant at 90%"
          }
          else{
            sigl_[i, j] <- "not significant at 90%"
          }
          if (is.na(probtl_[i, j])){
            sigl_[i, j] <- "not significant at 90%"
          }
        }
      }
      bistdt_ <- cbind(bistdt_, lambda_, stdlambda_, tstatl_, probtl_)
      colname2 <- xvarinf
      if (int_inf){
        colname2 <- c("Intercept", xvarinf)
      }
      label3_ <- rep("Inf_", ncol(G))
      colname3 <- c(label3_, colname2)
      label2_ <- rbind(rep("Inf_std_", ncol(G)), rep("Inf_tstat_", ncol(G)), rep("Inf_probt_", ncol(G)))
      colname <- cbind(c("x", "y"), colname1, c(label_, rep(colname1, 3)), colname3, c(label2_, rep(colname2, 3)))
      sub("_ ", "_", colname)
      sub("_ ", "_", colname)
      labell_ <- rep("sig_Inf_", ncol(G))
      colnamel <- c(labell_, rep(colname2, 1))
      sig_inf_parameters2 <- sigl_
      colnames(sig_inf_parameters2) <- colnamel
    }
    parameters2_ <- bistdt_
    colnames(parameters2_) <- colname
    label_ <- rep("sig_", ncol(x))
    colname_ <- c(label_, rep(colname1, 1))
    sig_parameters2 <- sig_
    colnames(sig_parameters2) <- colname
    if (model=="negbin" | model=="zip"){
      atstat <- alphai[, 2]
      for(j in 1:nrow(alphai)){
        if (alphai[j, 3]==0){
          atstat[j] <- 0
        }
        else{
          atstat[j] <- alphai[j, 2]/alphai[j, 3]
        }
      }
      atstat <- ifelse(is.na(atstat), 0, atstat)
      aprobtstat <- 2*(1-pnorm(abs(atstat)))
      siga_ <- matrix("not significant at 90%", N, 1)
      for (i in 1:N){
        if (aprobtstat[i]<0.01*(nvar/v1)){
          siga_[i] <- "significant at 99%"
        }
        else if (aprobtstat[i]<0.05*(nvar/v1)){
          siga_[i] <- "significant at 95%"
        }
        else if (aprobtstat[i]<0.1*(nvar/v1)){
          siga_[i] <- "significant at 90%"
        }
        else{
          siga_[i] <- "not significant at 90%"
        }
      }
      alphai <- cbind(alphai, atstat, aprobtstat)
      alpha_ <- alphai
      colnames(alpha_) <- c("id", "alpha", "std", "tstat", "probt")
      sig_alpha_ <- siga_
      colnames(sig_alpha_) <- "sig_alpha"
    }
  }
  else{
    beta_ <- bi
    colnames(beta_) <- c("id", "B", "x", "y")
    stdbi <- sqrt(abs(varbi))
    tstat <- bi[, 2]/stdbi
    for (i in 1:nrow(tstat)){
      if (is.na(tstat[i])){
        tstat[i] <- 0
      }
    }
    bistdt <- cbind(bi, stdbi, tstat)
    parameters_ <- bistdt
    colnames(parameters_) <- c("id", "B", "x", "y", "stdbi", "tstat")
    if (model=="negbin" | model=="zinb"){
      atstat <- alphai[, 2]
      for(j in 1:nrow(alphai)){
        if (alphai[j, 3]==0){
          atstat[j] <- 0
        }
        else
          atstat[j] <- alphai[j, 2]/alphai[j, 3]
      }
      atstat <- ifelse(is.na(atstat), 0, atstat)
      aprobtstat <- 2*(1-pnorm(abs(atstat)))
      alphai <- cbind(POINTS, alphai, atstat, aprobtstat)
      alpha_ <- alphai
      colnames(alpha_) <- c("x", "y", "id", "alpha", "std", "tstat", "probt")
    }
    beta_ <- matrix(bi[, 2], mm, byrow=TRUE)
    stdbeta_ <- matrix(stdbi, mm, byrow=TRUE)
    tstat_ <- beta_/stdbeta_
    bistdt_ <- cbind(POINTS, beta_, stdbeta_, tstat_)
    colname1 <- c("Intercept", XVAR)
    label_ <- rbind(rep("std_", nvar), rep("tstat_", nvar))
    colname <- cbind(c("x", "y"), colname1, c(label_, rep(colname1, 2)))
    colname <- sub("_ ", "_", colname)
    colname <- sub("_ ", "_", colname)
    parameters_grid_ <- bistdt_
    colnames(parameters_grid_) <- colname
  }
  if (is.null(grid)){
    parameters2 <- cbind(parameters2_, sig_parameters2_)
    if (model=="zip" | model=="zinb"){
      parameters2 <- cbind(parameters2, sig_inf_parameters2_)
    }
    if (model=="negbin" | model=="zinb"){
      alpha_ <- cbind(alpha_, sig_alpha_)
    }
  }
} # fecha gwzinbr

# -------------------------- TESTE

# Com h=82 (ou seja, sem rodar a golden)
library(readr)
korea_base_artigo <- read_csv("C:/Users/jehhv/OneDrive/Documentos/UnB/2024/TCC2/korea_base_artigo.csv")
korea_base_artigo <- read_csv("D:/Users/jessica.abreu/Documents/UnB/tcc/korea_base_artigo.csv")

#path Ju
korea_base_artigo <- read_csv("C:/Users/Juliana Rosa/OneDrive/Documents/TCC2/GWZINBR-main/korea_base_artigo.csv")

#obs.: mudar caso de teste
#obs2.: tirar NULL dos defaults

startTime <- Sys.time()
gwzinbr(data = korea_base_artigo, 
        formula = n_covid1~Morbidity+high_sch_p+Healthcare_access+
          diff_sd+Crowding+Migration+Health_behavior,
        xvarinf = c("Morbidity", "high_sch_p", "Healthcare_access", "diff_sd",
                    "Crowding", "Migration", "Health_behavior"),
        lat = "x", long = "y", offset = "ln_total", method = "adaptive_bsq",
        model = "zinb", distancekm = TRUE, h=82)
endTime <- Sys.time()
endTime-startTime

# alterar para vetor: dbb dlb e todos os outros
# linha 744  varabetalambda: trocar matrizes por rep(), se necessario
