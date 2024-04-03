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
    G <<- matrix(1, N, 1)
    # G <<- rep(1, N) #matriz coluna
    lambdag <<- matrix(0, ncol(G), 1) #ncol(G) em vez de length(G)
    # lambdag <<- rep(0, length(G))
  }
  else{
    G <<- unlist(data[, xvarinf])
    G <<- cbind(rep(1, N), G)
  }
  if (int_inf){ #o que e int_inf? 
    G <<- cbind(rep(1, N), G)
  }
  # x <<- cbind(rep(1, N), x)
  yhat <<- rep(0, n)
  yhat2 <<- rep(0, n)
  pihat <<- rep(0, n)
  nvar <<- ncol(x)
  wt <<- rep(1, n)
  if (!is.null(weight)){
    wt <<- unlist(data[, weight])
  }
  Offset <<- rep(0, N)
  if (!is.null(offset)){
    Offset <<- unlist(data[, offset])
  }
  Iy <- ifelse(y>0, 1, y)
  Iy <- 1-Iy
  pos0 <<- which(y==0)
  pos02 <<- which(y==0)
  pos1 <<- which(y>0)
  
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
    while (abs(ddev)>0.000001 & cont2<100){
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
        if (j>0){
          #parg <<- 1/(sum((y-uj)^2/uj)/(N-nvar))
          parg <- 1/(sum((y-uj)^2/uj)/(N-nvar))
        }
        while (abs(dpar)>0.0001 & aux2<200){
          if (parg<0){
            #parg <<- 0.00001
            parg <- 0.00001
          }
          parg <- ifelse(parg<E^-10, E^-10, parg)
          gf <- sum((1-zkg)*(digamma(parg+y)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+y)/(parg+uj)))
          hessg <- sum((1-zkg)*(trigamma(parg+y)-trigamma(parg)+1/parg-2/(parg+uj)+(y+parg)/(parg+uj)^2))
          hessg <- ifelse(hess==0, E^-23, hess)
          par0 <- parg
          parg <- as.vector(par0-solve(hessg)%*%gf)
          if (parg > E^5) {
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
    if (model == "zip" |model == 'zinb'){
      devg <- 0
      ddev <- 1
      njl <- G%*%lambdag
      njl <- ifelse(njl > maxg, maxg, njl)
      njl <- ifelse(njl < (-maxg),-maxg, njl)
      pig <- exp(njl)/(1+exp(njl))
      while (abs(ddev)>0.000001 & aux3<100){
        Ai <- pig*(1-pig)
        Ai <- ifelse(Ai<=0, E^-5, Ai)
        zj <- njl+(zkg-pig)*1/Ai
        if (det(t(G)%*%(Ai*G))==0){ #se der erro, conferir aqui
          lambdag <- matrix(0, ncol(G), 1)
        }  
        else {
          #lambdag <<- solve(t(G*Ai)%*%G)%*%t(G*Ai)%*%zj
          lambdag <- solve(t(G)%*%(Ai%*%G))%*%t(G)*(Ai%*%zj) #se der erro, conferir aqui tbm
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
    Iy <- matrix(0, nrow(y), 1)
  }
  dll <- Iy*(exp(njl)*g1x^parg/hgx^2)-exp(njl)/(1+exp(njl))^2
  dbb <- Iy*(-(parg*g1x^parg*g2x/hgx)^2+parg^2*g1x^parg*g2x^2*(1-1/uj)/hgx)-(1-Iy)*(parg*g2x*(1+(y-uj)/(parg+uj)))
  dlb <- Iy*(parg*exp(njl)*g1x^parg*g2x/hgx^2)
  I1 <- matrix(1, nrow(y), 1)
  # II <- -(t(I1) %*% daa) %*% I1 || -(t(I1) %*% dab) %*% X || -(t(I1) %*% dal) %*% G // (-t(X) %*% (dab * I1) || -t(X) %*% dbb %*% X || -t(X) %*% dlb %*% G) // (-t(G) %*% (dal * I1) || -t(G) %*% (X %*% dlb) || -t(G) %*% t(G) %*% dll %*% G)
  II <- rbind(cbind(-(t(I1) %*% daa) %*% I1, -(t(I1) %*% dab) %*% x, -(t(I1) %*% dal) %*% G), 
              cbind(-(t(x) %*% (dab * I1)), -t(x) %*% dbb %*% x, -t(x) %*% dlb %*% G), 
              cbind(-(t(G) %*% (dal * I1)), -t(G) %*% (x %*% dlb),  -t(G) %*% t(G) %*% dll %*% G)) #revisar
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
  
  long <- unlist(data[, long])
  lat <- unlist(data[, lat])
  COORD <- matrix(c(long, lat), ncol=2)
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
  dist_ <- sp::spDistsN1(COORD, COORD[i,]) 
  #substituicao: _dist_ por dist_ 
  sequ <- 1:N
  for(i in 1:mm){
    for(j in 1:N){
      seqi <- rep(i, N)
      distan <- cbind(seqi, sequ, dist_)
      if (distancekm){
        #substituicao: dist por distan
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
        if(dist[jj, 4] <= h){
          w[jj, 1] <- (1 -(dist[jj, 3]/hn)^2)^2
        } 
        else{
          w[jj, 1] <- 0
        }
        w[jj, 2] <- dist[jj, 2]
      }
    }
    if (method == "adaptive_bsq"){
      distan <- distan[order(distan[, 3]),]
      distan <- cbind(distan, 1:N)
      w <- matrix(0,N,2)
      hn <- dist[h, 3]
      
      for (jj in 1:n) {
        if (dist[jj, 4] <= h) {
          w[jj, 1] <- (1 - (dist[jj, 3] / hn)^2)^2
        } else {
          w[jj, 1] <- 0
        }
        w[jj, 2] <- dist[jj, 2]
      }
      
      w <- w[order(w[, 2]), 1]
    }
    # obs: Uma condicional exatamente igual a anterior (if(method=="adaptive_bsq)) viria logo a seguir
    # so que fora do if(grid), que foi retirado do R
    # ver com o professor Alan se é correto 
  }
  # /****** MODEL SELECTION *************/
  Iy <- Iy2
  b <- bg
  b2 <- b
  nj <- X%*%b+Offset
  uj <- exp(nj)
  par <- parg
  par2 <- par
  lambda <- lambdag
  njl <- G%*%lambda
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
      lambda <- cbind(lambda0, rep(0, ncol(G)-1))
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
        uj <- exp(X%*%b+Offset)
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
        hess <- sum(w*wt*(1-zk)*(trigamma(par+y)-trigamma(par)+1/par - 2/(par + uj) + (y + par)/(par+uj))^2)
        hess <- ifelse(hess == 0, E^-23, hess)
        par0 <- par
        par <- par0 - solve(hess)%*%gf
        dpar <- par - par0
        
        if (par >= E^6){
          par <- E^6
          dpar <- 0
          alpha <- 1/par
          b <- bg
          uj <- exp(X%*%b+Offset)
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
        uj <- exp(X%*%b+Offset)
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
    nj <- X%*%b+Offset
    nj <- ifelse(nj>700,700,nj)
    nj <- ifelse(nj<(-700), -700, nj)
    uj <- exp(nj)
    while (abs(ddev)>0.000001 & aux1<100){
      uj <- ifelse(uj>E^100, E^100, uj)
      Ai <- (1-zk)*((uj/(1+alpha*uj) + (y-uj) * (alpha*uj/(1+2*alpha*uj+alpha^2*uj^2))))
      Ai <- ifelse(Ai <=0, E^-5, Ai)
      uj <- ifelse(uj < E^-150, E^-150, uj)
      denz <- (((uj/(1+alpha*uj)+(y-uj)*(alpha*uj/(1+2*alpha*uj+alpha^2*uj^2))))*(1+alpha*uj))
      denz <- ifelse(denz == 0, E^-5, denz)
      zj <- (nj+(y-uj)/denz)-Offset
      if(det(t(x) %*% (w*Ai*x*wt)) == 0){
        b <- matrix(0, nvar, 1)
      } 
      else {
        b <- solve(t(x) %*% (w*Ai*x*wt)) %*% t(x) %*% (w*Ai*wt*zj)
      }
      nj <- X%*%b+Offset
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
          Aii <- pi*(1-pi)
          Aii <- ifelse(Aii<=0, E^-5, Aii)	
          zj <- njl+(zk-pi)/Aii
          if (det(t(G*Aii*w*wt)%*%G)==0){ #multiplicador
            lambda <- matrix(0, ncol(G), 1)
          }
          else{
            lambda <- solve(t(G*Aii*w*wt)%*%G)%*%t(G*Aii*w*wt)%*%zj
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
}

# /**** COMPUTING VARIANCE OF BETA AND LAMBDA ******/
  # Vou construir essa parte do lado de fora até saber onde fecha o while indicado
  if(det(t(x)%*%(w*Ai*x*wt))==0){
    C <- matrix(0, ncol(x),nrow(x)) 
  } else{
    C <- solve(t(x)%*%(w*Ai*x*wt))%*%t(x)%*%(w*Ai*wt)
  }
  Ci <- matrix(0, ncol(G), 1)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
/**** COMPUTING VARIANCE OF BETA AND LAMBDA ******/
                     %if %upcase(&MODEL)=ZIP or %upcase(&MODEL)=ZINB %then
                     %do;
                     
                     if det(G`*(w#Aii#G#wt))=0 then
                                Ci=j(ncol(G), nrow(G), 0);
                                else
                                  Ci=inv(G`*(w#Aii#G#wt))*G`#(w#Aii#wt)`;
                                         
                                         if any(lambda)=0 then
                                         Ci=j(ncol(G), n, 0);
                                         %END;
                                         g1x=par/(par+uj);
                                         g2x=uj/(par+uj);
                                         hgx=exp(njl)+g1x##par;
                                         hgx=choose(hgx>1E10, 1E10, hgx);
                                         daa=w#wt#(zk#((g1x##par#(log(g1x)+g2x))##2#(1-1/hgx)/hgx+g1x##par#(g2x##2/par)/hgx)+(1-zk)#(trigamma(par+y)-trigamma(par)-2/(uj+par)+1/par+(y+par)/(uj+par)##2));
                                         dab=w#wt#(zk#(g1x##(2*par+1)#uj#(log(g1x)+g2x)/hgx##2-g1x##par#(-g2x##2+par#g2x#(log(g1x)+g2x))/hgx)+(1-zk)#(g2x#(y-uj)/(uj+par)));
                                         dal=-w#wt#zk#(exp(njl)#g1x##par#(log(g1x)+g2x)/hgx##2);
                                         daa=daa*par**4;
                                         dab=dab*par**2;
                                         
                                         if any(lambda)=0 then
                                         Iy=j(nrow(y), 1, 0);
                                         exphx=1+exp(njl);
                                         exphx=choose(exphx>1E90, 1E90, exphx);
                                         dll=w#wt#(Iy#(exp(njl)#g1x##par/hgx##2)-exp(njl)/(exphx)##2);
                                         dbb=sqrt(w)#wt#(Iy#(-(par*g1x##par#g2x/hgx)##2+par**2*g1x##par#g2x##2#(1-1/uj)/hgx)-(1-Iy)#(par*g2x#(1+(y-uj)/(par+uj))));
                                         dlb=w#wt#Iy#(par*exp(njl)#g1x##par#g2x/hgx##2);
                                         dll=choose(dll=., 1E100, dll);
                                         daa=choose(daa=., 1E100, daa);
                                         dab=choose(dab=., 1E100, dab);
                                         dal=choose(dal=., 1E100, dal);
                                         dbb=choose(dbb=., 1E100, dbb);
                                         dlb=choose(dlb=., 1E100, dlb);
                                         I1=j(nrow(y), 1, 1);
                                         
                                         if any(b)=0 & any(lambda)=0 then
                                         do;
                                         II=j(ncol(x)+ncol(G)+1, ncol(x)+ncol(G)+1, 0);
                                         end;
                                         else if any(lambda)=0 then
                                         do;
                                         dbb=w#wt#(Iy#(-(par*g1x##par#g2x/hgx)##2+par**2*g1x##par#g2x##2#(1-1/uj)/hgx)-(1-Iy)#(par*g2x#(1+(y-uj)/(par+uj))));
                                         
                                         if det((X#dbb#dbb/Ai)`*X)=0 then
                                                 II=-(I1#daa)`*I1||-(I1#dab)`*X||-(I1#dal)`*G//(-X`*(dab#I1)||-((X#dbb)`*X)||-(X#dlb)`*G)//(-G`*(dal#I1)||-G`*(X#dlb)||-(G#dll)`*G);
                                                      else
                                                        II=-(I1#daa)`*I1||-(I1#dab)`*X||-(I1#dal)`*G//(-X`*(dab#I1)||-((X#dbb)`*X)*inv((X#dbb#dbb/Ai)`*X)*((X#dbb)`*X)||-(X#dlb)`*G)//(-G`*(dal#I1)||-G`*(X#dlb)||-(G#dll)`*G);
                                                             end;
                                                             else
                                                               II=-(I1#daa)`*I1||-(I1#dab)`*X||-(I1#dal)`*G//(-X`*(dab#I1)||-((X#dbb)`*X)||-(X#dlb)`*G)//(-G`*(dal#I1)||-G`*(X#dlb)||-(G#dll)`*G);
                                                                    
                                                                    if all(lambda)>0 & alpha=1E-6 then
                                                                    II=II[2:nrow(II), 2:nrow(II)];
                                                                    else if any(lambda)=0 & alpha>1E-6 then
                                                                    II=II[1:ncol(x)+1, 1:ncol(x)+1];
                                                                    else if any(lambda)=0 & alpha=1E-6 then
                                                                    II=II[2:ncol(x)+1, 2:ncol(x)+1];
                                                                    
                                                                    if det(II)=0 then
                                                                    do;
                                                                    
                                                                    if all(lambda)>0 & alpha=1E-6 then
                                                                    do;
                                                                    II=II[1:ncol(x), 1:ncol(x)];
                                                                    
                                                                    if det(II)=0 then
                                                                    varabetalambda=j(nrow(II), 1, 0)//j(ncol(G), 1, 0);
                                                                    else
                                                                      varabetalambda=vecdiag(inv(II))//j(ncol(G), 1, 0);
                                                                    end;
                                                                    else if any(lambda)=0 & alpha=1E-6 then
                                                                    do;
                                                                    II=II[1:ncol(x), 1:ncol(x)];
                                                                    
                                                                    if det(II)=0 then
                                                                    varabetalambda=j(nrow(II), 1, 0)//j(ncol(G), 1, 0);
                                                                    else
                                                                      varabetalambda=vecdiag(inv(II))//j(ncol(G), 1, 0);
                                                                    end;
                                                                    else
                                                                      do;
                                                                    II=II[1:ncol(x)+1, 1:ncol(x)+1];
                                                                    
                                                                    if det(II)=0 then
                                                                    varabetalambda=j(nrow(II), 1, 0)//j(ncol(G), 1, 0);
                                                                    ;
                                                                    else
                                                                      varabetalambda=vecdiag(inv(II))//j(ncol(G), 1, 0);
                                                                    end;
                                                                    end;
                                                                    else
                                                                      varabetalambda=vecdiag(inv(II));
                                                                    
                                                                    if all(lambda)>0 & alpha>1E-6 then
                                                                    do;
                                                                    varb=varabetalambda[2:ncol(x)+1];
                                                                    varl=varabetalambda[ncol(x)+2:nrow(varabetalambda)];
                                                                    alphai[i, 1]=i;
                                                                    alphai[i, 2]=alpha;
                                                                    alphai[i, 3]=sqrt(abs(varabetalambda[1]));
                                                                    end;
                                                                    else if all(lambda)>0 & alpha=1E-6 then
                                                                    do;
                                                                    varb=varabetalambda[1:ncol(x)];
                                                                    varl=varabetalambda[ncol(x)+1:nrow(varabetalambda)];
                                                                    alphai[i, 1]=i;
                                                                    alphai[i, 2]=alpha;
                                                                    alphai[i, 3]=sqrt(1/abs(-(I1#daa)`*I1));
                                                                                              end;
                                                                                              else if any(lambda)=0 & alpha>1E-6 then
                                                                                              do;
                                                                                              varb=varabetalambda[2:ncol(x)+1];
                                                                                              varl=j(ncol(G), 1, 0);
                                                                                              alphai[i, 1]=i;
                                                                                              alphai[i, 2]=alpha;
                                                                                              alphai[i, 3]=sqrt(abs(varabetalambda[1]));
                                                                                              end;
                                                                                              else if any(lambda)=0 & alpha=1E-6 then
                                                                                              do;
                                                                                              varb=varabetalambda[1:ncol(x)];
                                                                                              varl=j(ncol(G), 1, 0);
                                                                                              alphai[i, 1]=i;
                                                                                              alphai[i, 2]=alpha;
                                                                                              alphai[i, 3]=sqrt(1/abs(-(I1#daa)`*I1));
                                                                                                                        end;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              /*******************************/
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                m1=(i-1)*ncol(x)+1;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              m2=m1+(ncol(x)-1);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              bi[m1:m2, 1]=i;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              bi[m1:m2, 2]=b;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              bi[m1:m2, 3]=POINTS[i, 1];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              bi[m1:m2, 4]=POINTS[i, 2];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              varbi[m1:m2, 1]=varb;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              %if %upcase(&MODEL)=ZIP or %upcase(&MODEL)=ZINB %then
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              %do;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              m1=(i-1)*ncol(G)+1;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              m2=m1+(ncol(G)-1);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              li[m1:m2, 1]=i;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              li[m1:m2, 2]=lambda;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              li[m1:m2, 3]=POINTS[i, 1];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              li[m1:m2, 4]=POINTS[i, 2];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              varli[m1:m2, 1]=varl;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              %end;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              %if &grid=%then
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              %do;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              r=x[i, ]*C;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              S[i]=r[i];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              S2[i]=r*r`;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              yhat[i]=uj[i];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              pihat[i]=njl[i];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              %if %upcase(&MODEL)=ZIP or %upcase(&MODEL)=ZINB %then
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              %do;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              ri=G[i, ]*Ci;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              Si[i]=ri[i];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              yhat2[i]=uj[i];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              yhat[i]=(uj#(1-exp(njl)/(1+exp(njl))))[i];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       %end;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       /** creating non-stationarity matrix **/
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         %IF %UPCASE(&METHOD) ne ADAPTIVE_BSQ %THEN
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       %DO;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       CCC=x||w||wt;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       m1=(i-1)*ncol(x)+1;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       m2=m1+(ncol(x)-1);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       if det(CCC[, 1:ncol(x)]`*(CCC[, ncol(CCC)-1]#CCC[, 1:ncol(x)]#CCC[, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ncol(CCC)]))=0 then
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       BB[m1:m2, ]=j(ncol(x), nrow(x), 0);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       else
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         BB[m1:m2, ]=inv(CCC[, 1:ncol(x)]`*(CCC[, ncol(CCC)-1]#CCC[, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            1:ncol(x)]#CCC[, ncol(CCC)]))*CCC[, 1:ncol(x)]`#(CCC[, 
ncol(CCC)-1]#CCC[, ncol(CCC)])`;

%if %upcase(&MODEL)=ZIP or %upcase(&MODEL)=ZINB %then
%do;
CCCl=G||w||wt;
m1=(i-1)*ncol(G)+1;
m2=m1+(ncol(G)-1);

if det(CCCl[, 1:ncol(G)]`*(CCCl[, ncol(CCCl)-1]#CCCl[, 
                           1:ncol(G)]#CCCl[, ncol(CCCl)]))=0 then
BBl[m1:m2, ]=j(ncol(G), nrow(G), 0);
else
  BBl[m1:m2, ]=inv(CCCl[, 1:ncol(G)]`*(CCCl[, ncol(CCCl)-1]#CCCl[, 
                                       1:ncol(G)]#CCCl[, ncol(CCCl)]))*CCCl[, 1:ncol(G)]`#(CCCl[, 
ncol(CCCl)-1]#CCCl[, ncol(CCCl)])`;
%end;
%END;

/*************************************/
  _w_=w;
  call sort(_w_, {1});
  sumwi[i]=_w_[1:int(nrow(_w_)*1)][+];
  %end;
  
  if i=1 then
  W_f=W||(1:nrow(w))`;
  else
    W_f=W_f//(W||(1:nrow(w))`);
  end;
  create w_f from W_f;
  append from W_f;
  close w_f;
  free w_f;
  
  /***********************************************/
    
    %if &grid=%then
  %do;
  v1=sum(S)+sum(Si);
  v11=sum(S)+sum(Si);
  v2=sum(S2);
  nparmodel=n-v11;
  
  if v11<v2 then
  v1=v11;
  res=(y-yhat);
  rsqr1=(res#wt)`*res;
         ym=(y#wt)`*y;
             rsqr2=ym-((y#wt)[+]**2)/wt[+];
                        rsqr=1-rsqr1/rsqr2;
                        rsqradj=1-((n-1)/(n-v1))*(1-rsqr);
                        sigma2=n*rsqr1/((n-v1)*wt[+]);
                        root_mse=sqrt(sigma2);
                        print sigma2[label='Sigma2e'] root_mse[label='Root MSE'] 
                        v1[label='#GWR parameters'] nparmodel[label='#GWR parameters (model)'] 
                        v2[label='#GWR parameters (variance)'];
                        influence=S;
                        resstd=res/(sqrt(sigma2)*sqrt(abs(1-influence)));
                        CooksD=resstd#resstd#influence/(v1*(1-influence));
                        df=n-(nvar+ncol(G));
                        stdbi=sqrt(abs(varbi));
                        tstat=bi[, 2]/stdbi;
                        probt=2*(1-probt(abs(tstat), df));
                        _malpha_=0.05*(nvar/v1);
                        
                        %if %upcase(&MODEL)=ZIP or %upcase(&MODEL)=ZINB %then
                        %do;
                        stdli=sqrt(abs(varli));
                        tstati=li[, 2];
                        
                        do j=1 to nrow(stdli);
                        
                        if stdli[j]=0 then
                        tstati[j]=0;
                        else
                          tstati[j]=li[j, 2]/stdli[j];
                        end;
                        tstati=choose(tstati=., 0, tstati);
                        probti=2*(1-probt(abs(tstati), df));
                        _malpha_=0.05*((nvar+ncol(G))/v1);
                        %end;
                        _t_critical_=abs(tinv(_malpha_/2, df));
                        _par_=1/alphai[, 2];
                        
                        %IF %UPCASE(&MODEL)=ZINB or %UPCASE(&MODEL)=ZIP %THEN
                        %DO;
                        
                        if any(lambda)=0 then
                        do;
                        ll=sum(-log(0+exp(pihat[pos0]))+log(0*exp(pihat[pos0])+(_par_[pos0]/(_par_[pos0]+yhat2[pos0]))##_par_[pos0]))+
                                                            sum(-log(0+exp(pihat[pos1]))+lgamma(_par_[pos1]+y[pos1, ])-lgamma(y[pos1, 
                                                            ]+1)-lgamma(_par_[pos1])+
                                                              y[pos1]#log(yhat2[pos1]/(_par_[pos1]+yhat2[pos1]))+_par_[pos1]#log(_par_[pos1]/(_par_[pos1]+yhat2[pos1])));
                                                            llnull1=sum(-log(1+zk[pos0])+log(zk[pos0]+(_par_[pos0]/(_par_[pos0]+y[pos0, 
                                                            ]))##_par_[pos0]))+
                                                            sum(-log(1+zk[pos1])+lgamma(_par_[pos1]+y[pos1, ])-lgamma(y[pos1, 
                                                            ]+1)-lgamma(_par_[pos1])+
                                                              y[pos1]#log(y[pos1, ]/(_par_[pos1]+y[pos1, 
                                                            ]))+_par_[pos1]#log(_par_[pos1]/(_par_[pos1]+y[pos1, ])));
                                                            llnull2=sum(-log(1+0)+log(0+(_par_/(_par_+y[:]))##_par_))+
                                                                                      sum(-log(1+0)+lgamma(_par_+y)-lgamma(y+1)-lgamma(_par_)+
                                                                                            y#log(y[:]/(_par_+y[:]))+_par_#log(_par_/(_par_+y[:])));
                                                                                          end;
                                                                                          else
                                                                                            do;
                                                                                          ll=sum(-log(1+exp(pihat[pos0]))+log(exp(pihat[pos0])+(_par_[pos0]/(_par_[pos0]+yhat2[pos0]))##_par_[pos0]))+
                                                                                                                              sum(-log(1+exp(pihat[pos1]))+lgamma(_par_[pos1]+y[pos1])-lgamma(y[pos1]+1)-lgamma(_par_[pos1])+
                                                                                                                                    y[pos1]#log(yhat2[pos1]/(_par_[pos1]+yhat2[pos1]))+_par_[pos1]#log(_par_[pos1]/(_par_[pos1]+yhat2[pos1])));
                                                                                                                                  llnull1=sum(-log(1+zk[pos0])+log(zk[pos0]+(_par_[pos0]/(_par_[pos0]+y[pos0]))##_par_[pos0]))+
                                                                                                                                                                   sum(-log(1+zk[pos1])+lgamma(_par_[pos1]+y[pos1])-lgamma(y[pos1]+1)-lgamma(_par_[pos1])+
                                                                                                                                                                         y[pos1]#log(y[pos1]/(_par_[pos1]+y[pos1]))+_par_[pos1]#log(_par_[pos1]/(_par_[pos1]+y[pos1])));
                                                                                                                                                                       llnull2=sum(-log(1+0)+log(0+(_par_[:]/(_par_[:]+y[:]))##_par_[:]))+
                                                                                                                                                                                                 sum(-log(1+0)+lgamma(_par_[:]+y)-lgamma(y+1)-lgamma(_par_[:])+
                                                                                                                                                                                                       y#log(y[:]/(_par_[:]+y[:]))+_par_[:]#log(_par_[:]/(_par_[:]+y[:])));
                                                                                                                                                                                                     end;
                                                                                                                                                                                                     dev=2*(llnull1-ll);
                                                                                                                                                                                                     pctll=1-(llnull1-ll)/(llnull1-llnull2);
                                                                                                                                                                                                     AIC=2*v1-2*ll;
                                                                                                                                                                                                     AICc=AIC+2*(v1*(v1+1)/(n-v1-1));
                                                                                                                                                                                                     adjpctll=1-(llnull1-ll+v1+0.5)/(llnull1-llnull2);
                                                                                                                                                                                                     
                                                                                                                                                                                                     %IF %UPCASE(&MODEL)=ZINB %THEN
                                                                                                                                                                                                     %DO;
                                                                                                                                                                                                     AIC=2*(v1+v1/(ncol(x)+ncol(G)))-2*ll;
                                                                                                                                                                                                     AICc=AIC+2*((v1+v1/(ncol(x)+ncol(G)))*((v1+v1/(ncol(x)+ncol(G)))+1)/(n-(v1+v1/(ncol(x)+ncol(G)))-1));
                                                                                                                                                                                                     adjpctll=1-(llnull1-ll+(v1+v1/(ncol(x)+ncol(G)))+0.5)/(llnull1-llnull2);
                                                                                                                                                                                                     %END;
                                                                                                                                                                                                     print dev[label='Deviance'] ll[label='Full Log Likelihood'] pctll 
                                                                                                                                                                                                     adjpctll AIC AICc;
                                                                                                                                                                                                     %END;
                                                                                                                                                                                                     
                                                                                                                                                                                                     %IF %UPCASE(&MODEL)=POISSON or %UPCASE(&MODEL)=NEGBIN %THEN
                                                                                                                                                                                                     %DO;
                                                                                                                                                                                                     
                                                                                                                                                                                                     if ncol(pos02)=0 then
                                                                                                                                                                                                     do;
                                                                                                                                                                                                     pos0=pos1;
                                                                                                                                                                                                     pos0x=1;
                                                                                                                                                                                                     pos0xl=1;
                                                                                                                                                                                                     end;
                                                                                                                                                                                                     else
                                                                                                                                                                                                       do;
                                                                                                                                                                                                     pos0x=(_par_[pos0]/(_par_[pos0]+yhat[pos0]))##_par_[pos0];
                                                                                                                                                                                                     pos0xl=(_par_[pos0]/(_par_[pos0]+y[pos0, ]))##_par_[pos0];
                                                                                                                                                                                                     end;
                                                                                                                                                                                                     ll=sum(-log(0+exp(pihat[pos0]))+log(0*exp(pihat[pos0])+pos0x))+
                                                                                                                                                                                                       sum(-log(0+exp(pihat[pos1]))+lgamma(_par_[pos1]+y[pos1, ])-lgamma(y[pos1, 
                                                                                                                                                                                                       ]+1)-lgamma(_par_[pos1])+
                                                                                                                                                                                                         y[pos1]#log(yhat[pos1]/(_par_[pos1]+yhat[pos1]))+
                                                                                                                                                                                                       _par_[pos1]#log(_par_[pos1]/(_par_[pos1]+yhat[pos1])));
                                                                                                                                                                                                       llnull1=sum(-log(1+zk)+log(zk+pos0xl))+
                                                                                                                                                                                                         sum(-log(1+zk)+lgamma(_par_[pos1]+y[pos1, ])-lgamma(y[pos1, 
                                                                                                                                                                                                         ]+1)-lgamma(_par_[pos1])+
                                                                                                                                                                                                           y[pos1]#log(y[pos1, ]/(_par_[pos1]+y[pos1, 
                                                                                                                                                                                                         ]))+_par_[pos1]#log(_par_[pos1]/(_par_[pos1]+y[pos1, ])));
                                                                                                                                                                                                     pos1xx=0+(_par_[pos0]/(_par_[pos0]+y[:]))##_par_[pos0];
                                                                                                                                                                                                     pos1xx=choose(pos0x<=0, 1E-100, pos0x);
                                                                                                                                                                                                     llnull2=sum(-log(1+0)+log(pos1xx))+
                                                                                                                                                                                                       sum(-log(1+0)+lgamma(_par_[pos1]+y[pos1])-lgamma(y[pos1]+1)-lgamma(_par_[pos1])+
                                                                                                                                                                                                             y[pos1]#log(y[:]/(_par_[pos1]+y[:]))+_par_[pos1]#log(_par_[pos1]/(_par_[pos1]+y[:])));
                                                                                                                                                                                                           dev=2*(llnull1-ll);
                                                                                                                                                                                                           AIC=-2*ll+2*v1;
                                                                                                                                                                                                           AICc=AIC+2*(v1*(v1+1)/(n-v1-1));
                                                                                                                                                                                                           pctll=1-(llnull1-ll)/(llnull1-llnull2);
                                                                                                                                                                                                           adjpctll=1-(llnull1-ll+v1+0.5)/(llnull1-llnull2);
                                                                                                                                                                                                           
                                                                                                                                                                                                           %IF %UPCASE(&MODEL)=NEGBIN %THEN
                                                                                                                                                                                                           %DO;
                                                                                                                                                                                                           AIC=2*(v1+v1/ncol(x))-2*ll;
                                                                                                                                                                                                           AICc=AIC+2*(v1+v1/ncol(x))*(v1+v1/ncol(x)+1)/(n-(v1+v1/ncol(x))-1);
                                                                                                                                                                                                           adjpctll=1-(llnull1-ll+v1+v1/ncol(x)+0.5)/(llnull1-llnull2);
                                                                                                                                                                                                           %END;
                                                                                                                                                                                                           print dev[label='Deviance'] ll[label='Full Log Likelihood'] pctll 
                                                                                                                                                                                                           adjpctll AIC AICc;
                                                                                                                                                                                                           %END;
                                                                                                                                                                                                           _beta_=shape(bi[, 1:2], n);
                                                                                                                                                                                                           _beta2_=_beta_;
                                                                                                                                                                                                           
                                                                                                                                                                                                           %IF %UPCASE(&MODEL)=NEGBIN or %UPCASE(&MODEL)=ZINB %THEN
                                                                                                                                                                                                           %DO;
                                                                                                                                                                                                           _alpha_=shape(alphai[, 1:2], n);
                                                                                                                                                                                                           _beta2_=_beta_||_alpha_;
                                                                                                                                                                                                           %END;
                                                                                                                                                                                                           i=do(2, ncol(_beta_), 2);
                                                                                                                                                                                                           _beta_=_beta_[, i];
                                                                                                                                                                                                           i=do(2, ncol(_beta2_), 2);
                                                                                                                                                                                                           _beta2_=_beta2_[, i];
                                                                                                                                                                                                           call qntl(qntl, _beta2_);
                                                                                                                                                                                                           qntl=qntl//(qntl[3, ]-qntl[1, ]);
                                                                                                                                                                                                           descriptb=_beta2_[:, ]//_beta2_[><, ]//_beta2_[<>, ];
                                                                                                                                                                                                           print qntl[label="Quantiles of GWR Parameter Estimates" rowname={"P25", 
                                                                                                                                                                                                             "P50", "P75", "IQR"} colname={'Intercept' &xvar 
                                                                                                                                                                                                               %IF %UPCASE(&MODEL)=NEGBIN or %UPCASE(&MODEL)=ZINB %THEN
                                                                                                                                                                                                               %DO;
                                                                                                                                                                                                               'alpha' %END;
                                                                                                                                                                                                             }], , descriptb[label="Descriptive Statistics" rowname={"Mean", "Min", 
                                                                                                                                                                                                               "Max"} colname={'Intercept' &xvar 
                                                                                                                                                                                                                 %IF %UPCASE(&MODEL)=NEGBIN or %UPCASE(&MODEL)=ZINB %THEN
                                                                                                                                                                                                                 %DO;
                                                                                                                                                                                                                 'alpha' %END;
                                                                                                                                                                                                               }];
                                                                                                                                                                                                           _stdbeta_=shape(stdbi, n);
                                                                                                                                                                                                           _stdbeta2_=_stdbeta_;
                                                                                                                                                                                                           
                                                                                                                                                                                                           %IF %UPCASE(&MODEL)=NEGBIN or %UPCASE(&MODEL)=ZINB %THEN
                                                                                                                                                                                                           %DO;
                                                                                                                                                                                                           _stdalpha_=shape(alphai[, 3], n);
                                                                                                                                                                                                           _stdbeta2_=_stdbeta_||_stdalpha_;
                                                                                                                                                                                                           %END;
                                                                                                                                                                                                           call qntl(qntls, _stdbeta2_);
                                                                                                                                                                                                           qntls=qntls//(qntls[3, ]-qntls[1, ]);
                                                                                                                                                                                                           descripts=_stdbeta2_[:, ]//_stdbeta2_[><, ]//_stdbeta2_[<>, ];
                                                                                                                                                                                                           print _malpha_[label="alpha-level=0.05"] _t_critical_[format=comma6.2 
                                                                                                                                                                                                                                                                 label="t-Critical"] df;
                                                                                                                                                                                                           print qntls[label="Quantiles of GWR Standard Errors" rowname={"P25", 
                                                                                                                                                                                                             "P50", "P75", "IQR"} colname={'Intercept' &xvar 
                                                                                                                                                                                                               %IF %UPCASE(&MODEL)=NEGBIN or %UPCASE(&MODEL)=ZINB %THEN
                                                                                                                                                                                                               %DO;
                                                                                                                                                                                                               'alpha' %END;
                                                                                                                                                                                                             }], , descripts[label="Descriptive Statistics of Standard Errors" 
                                                                                                                                                                                                                             rowname={"Mean", "Min", "Max"} colname={'Intercept' &xvar 
                                                                                                                                                                                                                               %IF %UPCASE(&MODEL)=NEGBIN or %UPCASE(&MODEL)=ZINB %THEN
                                                                                                                                                                                                                               %DO;
                                                                                                                                                                                                                               'alpha' %END;
                                                                                                                                                                                                                             }];
                                                                                                                                                                                                           
                                                                                                                                                                                                           %if %upcase(&MODEL)=ZIP or %upcase(&MODEL)=ZINB %then
                                                                                                                                                                                                           %do;
                                                                                                                                                                                                           _lambda_=shape(li[, 1:2], n);
                                                                                                                                                                                                           _lambda2_=_lambda_;
                                                                                                                                                                                                           i=do(2, ncol(_lambda_), 2);
                                                                                                                                                                                                           _lambda_=_lambda_[, i];
                                                                                                                                                                                                           i=do(2, ncol(_lambda2_), 2);
                                                                                                                                                                                                           _lambda2_=_lambda2_[, i];
                                                                                                                                                                                                           call qntl(qntl, _lambda2_);
                                                                                                                                                                                                           qntl=qntl//(qntl[3, ]-qntl[1, ]);
                                                                                                                                                                                                           descriptl=_lambda2_[:, ]//_lambda2_[><, ]//_lambda2_[<>, ];
                                                                                                                                                                                                           print qntl[label="Quantiles of GWR Zero Inflation Parameter Estimates" 
                                                                                                                                                                                                                      rowname={"P25", "P50", "P75", "IQR"} 
                                                                                                                                                                                                                      %IF %UPCASE(&INT_INF)=YES %THEN
                                                                                                                                                                                                                      %DO;
                                                                                                                                                                                                                      colname={'Intercept' &XVARINF}], , %END;
                                                                                                                                                                                                           %ELSE
                                                                                                                                                                                                           %DO;
                                                                                                                                                                                                           colname={&XVARINF}], , %END;
descriptl[label="Descriptive Statistics" rowname={"Mean", "Min", "Max"} 
          %IF %UPCASE(&INT_INF)=YES %THEN
          %DO;
          colname={'Intercept' &XVARINF}];
%END;
%ELSE
%DO;
colname={&XVARINF}];
%END;
_stdlambda_=shape(stdli, n);
_stdlambda2_=_stdlambda_;
call qntl(qntls, _stdlambda2_);
qntls=qntls//(qntls[3, ]-qntls[1, ]);
descriptls=_stdlambda2_[:, ]//_stdlambda2_[><, ]//_stdlambda2_[<>, ];
print _malpha_[label="alpha-level=0.05"] _t_critical_[format=comma6.2 
                                                      label="t-Critical"] df;
print qntls[label="Quantiles of GWR Zero Inflation Standard Errors" 
            rowname={"P25", "P50", "P75", "IQR"} 
            %IF %UPCASE(&INT_INF)=YES %THEN
            %DO;
            colname={'Intercept' &XVARINF}], , %END;
%ELSE
%DO;
colname={&XVARINF}], , %END;
descriptls[label="Descriptive Statistics of Zero Inflation Standard Errors" 
           rowname={"Mean", "Min", "Max"} 
           %IF %UPCASE(&INT_INF)=YES %THEN
           %DO;
           colname={'Intercept' &XVARINF}];
%END;
%ELSE
%DO;
colname={&XVARINF}];
%END;
%end;
%end;

/****** Non-Stationarity Test *****************/
  
  %if &grid=%then
%do;

%IF %UPCASE(&METHOD) ne ADAPTIVE_BSQ %THEN
%DO;
BBk=j(n, n, 0);
Vk=j(ncol(x), 1, 0);
df1k=j(ncol(x), 1, 0);
df2k=j(ncol(x), 1, 0);

do k=1 to ncol(x);
ek=j(ncol(x), 1, 0);
ek[k]=1;

do i=1 to n;
m1=(i-1)*ncol(x)+1;
m2=m1+(ncol(x)-1);
BBk[i, ]=ek`*BB[m1:m2, ];
end;
Vk[k]=y`*(1/n)*BBk`*(I(n)-(1/n)*J(n, n, 1))*BBk*y;
df1k[k]=trace((1/n)*BBk`*(I(n)-(1/n)*J(n, n, 1))*BBk);
df2k[k]=trace(((1/n)*BBk`*(I(n)-(1/n)*J(n, n, 1))*BBk)**2);
end;
Vk=choose(abs(Vk)<=1E-8, 0, Vk);
Fk=(Vk/df1k)/sigma2;
ndf=df1k##2/df2k;
ddf=n-v1;
ddf=repeat(ddf, ncol(x));
probf=1-probf(Fk, ndf, ddf);
print , , "Non-Stationarity Test (Leung et al., 2000)", , Vk[label='' 
                                                             rowname={'Intercept' &xvar} colname={"V"}] Fk[label='' colname={"F"}] 
ndf ddf probf[format=pvalue6. label="Pr > F"];

%if %upcase(&MODEL)=ZIP or %upcase(&MODEL)=ZINB %then
%do;
BBkl=j(n, n, 0);
Vkl=j(ncol(G), 1, 0);
df1kl=j(ncol(G), 1, 0);
df2kl=j(ncol(G), 1, 0);

do k=1 to ncol(G);
ekl=j(ncol(G), 1, 0);
ekl[k]=1;

do i=1 to n;
m1=(i-1)*ncol(G)+1;
m2=m1+(ncol(G)-1);
BBkl[i, ]=ekl`*BBl[m1:m2, ];
end;
Vkl[k]=y`*(1/n)*BBkl`*(I(n)-(1/n)*J(n, n, 1))*BBkl*y;
df1kl[k]=trace((1/n)*BBkl`*(I(n)-(1/n)*J(n, n, 1))*BBkl);
df2kl[k]=trace(((1/n)*BBkl`*(I(n)-(1/n)*J(n, n, 1))*BBkl)**2);
end;
Vkl=choose(abs(Vkl)<=1E-8, 0, Vkl);
Fkl=(Vkl/df1kl)/sigma2;
ndfl=df1kl##2/df2kl;
ddfl=n-v1;
ddfl=repeat(ddfl, ncol(G));
probfl=1-probf(Fkl, ndfl, ddfl);
print , , 
"Non-Stationarity Test (Leung et al., 2000) - Zero Inflation", , 
%IF %UPCASE(&INT_INF)=YES %THEN
%DO;
Vkl[label='' rowname={'Intercept' &XVARINF} %END;
    %ELSE
    %DO;
    Vkl[label='' rowname={&XVARINF} %END;
        colname={"V"}] Fkl[label='' colname={"F"}] ndfl ddfl 
    probfl[format=pvalue6. label="Pr > F"];
    %end;
    %END;
    %end;
    
    /***** global estimates ***************/
      
      
      %IF &WEIGHT=%THEN
    %DO;
    dfg=n-nvar;
    %END;
    
    %IF %UPCASE(&MODEL)=ZINB %THEN
    %DO;
    b2=bg//alphag;
    
    if alphag=1E-6 then
    stdg=stdabetalambdag[1:ncol(x)]//(sqrt(1/abs(hessg))/(parg**2));
    else
      stdg=stdabetalambdag[2:ncol(x)+1]//stdabetalambdag[1];
    tg=b2/stdg;
    dfg=nrow(y)-ncol(x);
    probtg=2#(1-probt(abs(tg), dfg));
    lambdag=lambdag;
    
    if alphag=1E-6 then
    stdlambdag=stdabetalambdag[ncol(x)+1:nrow(stdabetalambdag)];
    else
      stdlambdag=stdabetalambdag[ncol(x)+2:nrow(stdabetalambdag)];
    tlambdag=lambdag/stdlambdag;
    dflg=nrow(y)-ncol(G);
    probtlambdag=2#(1-probt(abs(tlambdag), dflg));
    p=ncol(x)+1+ncol(G);
    %END;
    
    %IF %UPCASE(&MODEL)=ZIP %THEN
    %DO;
    b2=bg;
    stdg=stdabetalambdag[1:ncol(x)];
    tg=b2/stdg;
    dfg=nrow(y)-ncol(x);
    probtg=2#(1-probt(abs(tg), dfg));
    lambdag=lambdag;
    stdlambdag=stdabetalambdag[ncol(x)+1:nrow(stdabetalambdag)];
    tlambdag=lambdag/stdlambdag;
    dflg=nrow(y)-ncol(G);
    probtlambdag=2#(1-probt(abs(tlambdag), dflg));
    p=ncol(x)+ncol(G);
    %END;
    
    %IF %UPCASE(&MODEL)=NEGBIN %THEN
    %DO;
    b2=bg//alphag;
    
    if alphag=1E-6 then
    stdg=stdabetalambdag[1:nrow(stdabetalambdag)]//(sqrt(1/abs(hessg))/(parg**2));
    else
      stdg=stdabetalambdag[2:nrow(stdabetalambdag)]//stdabetalambdag[1];
    tg=b2/stdg;
    dfg=nrow(y)-ncol(x);
    probtg=2#(1-probt(abs(tg), dfg));
    p=ncol(x)+1;
    %END;
    
    %IF %UPCASE(&MODEL)=POISSON %THEN
    %DO;
    b2=bg;
    stdg=stdabetalambdag;
    tg=b2/stdg;
    dfg=nrow(y)-ncol(x);
    probtg=2#(1-probt(abs(tg), dfg));
    p=ncol(x);
    %END;
    bg_stdg=b2||stdg;
    print "Global Parameter Estimates", , bg_stdg[label=' ' 
                                                  rowname={'Intercept' &xvar 
                                                    %IF %UPCASE(&MODEL)=NEGBIN or %UPCASE(&MODEL)=ZINB %THEN
                                                    %DO;
                                                    'alpha' %END;
                                                  } colname={"Par. Est." "Std Error"}] tg[format=comma6.2 label="t Value"] 
    probtg[format=pvalue6. label="Pr > |t|"], , 
    "NOTE: The denominator degrees of freedom for the t tests is" 
    dfg[label=' ']".";
    
    %IF %UPCASE(&MODEL)=ZIP or %UPCASE(&MODEL)=ZINB %THEN
    %DO;
    print "Analysis Of Maximum Likelihood Zero Inflation Parameter Estimate";
    
    %if &XVARINF ne %then
    %do;
    varnamezip1=varnamezip`;
    
    %IF %UPCASE(&INT_INF)=YES %THEN
    %DO;
    varnamezip1="Intercept"//varnamezip`;
    %END;
    %end;
    %else
      %do;
    
    %IF %UPCASE(&INT_INF)=YES %THEN
    %DO;
    varnamezip1={"Intercept"};
    %END;
    %end;
    print varnamezip1[label="Parameter"] lambdag[label="Estimate" 
                                                 format=12.6] stdlambdag[label="Standard Error" format=12.6] 
    tlambdag[label="t Value" format=12.2] probtlambdag[label="Pr > |t|" 
                                                       format=pvalue6.4], , 
    "NOTE: The denominator degrees of freedom for the t tests is" 
    dflg[label=' ']".";
    %END;
    
    %IF %UPCASE(&MODEL)=ZINB or %UPCASE(&MODEL)=ZIP %THEN
    %DO;
    ll=sum(-log(1+exp(G[pos0, ]*lambdag))+log(exp(G[pos0, 
    ]*lambdag)+(parg/(parg+exp(X[pos0, ]*bg+offset[pos0, ])))##parg))+
    sum(-log(1+exp(G[pos1, ]*lambdag))+lgamma(parg+y[pos1, ])-lgamma(y[pos1, 
    ]+1)-lgamma(parg)+
      y[pos1]#log(exp(X[pos1, ]*bg+offset[pos1, ])/(parg+exp(X[pos1, 
]*bg+offset[pos1, ])))+parg*log(parg/(parg+exp(X[pos1, ]*bg+offset[pos1, 
]))));
                                                                                                                                                                                                     AIC=2*p-2*ll;
                                                                                                                                                                                                     AICc=AIC+2*(p*(p+1)/(n-p-1));
                                                                                                                                                                                                     llnull1=sum(-log(1+zkg[pos0, ])+log(zkg[pos0, ]+(parg/(parg+y[pos0, 
                                                                                                                                                                                                     ]))##parg))+
                                                                                                                                                                                                     sum(-log(1+zkg[pos1, ])+lgamma(parg+y[pos1, ])-lgamma(y[pos1, 
                                                                                                                                                                                                     ]+1)-lgamma(parg)+
                                                                                                                                                                                                       y[pos1]#log(y[pos1, ]/(parg+y[pos1, ]))+parg*log(parg/(parg+y[pos1, ])));
                                                                                                                                                                                                     llnull2=sum(-log(1+0)+log(0+(parg/(parg+y[:]))##parg))+
                                                                                                                                                                                                                               sum(-log(1+0)+lgamma(parg+y)-lgamma(y+1)-lgamma(parg)+
                                                                                                                                                                                                                                     y#log(y[:]/(parg+y[:]))+parg*log(parg/(parg+y[:])));
                                                                                                                                                                                                                                   devg=2*(llnull1-ll);
                                                                                                                                                                                                                                   pctll=1-(llnull1-ll)/(llnull1-llnull2);
                                                                                                                                                                                                                                   adjpctll=1-(llnull1-ll+p+0.5)/(llnull1-llnull2);
                                                                                                                                                                                                                                   print devg[label='Deviance'] ll[label='Full Log Likelihood'] pctll 
                                                                                                                                                                                                                                   adjpctll AIC AICc;
                                                                                                                                                                                                                                   %END;
                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                   %IF %UPCASE(&MODEL)=POISSON or %UPCASE(&MODEL)=NEGBIN %THEN
                                                                                                                                                                                                                                   %DO;
                                                                                                                                                                                                                                   yhatg=exp(x*bg+offset);
                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                   if ncol(pos02)=0 then
                                                                                                                                                                                                                                   do;
                                                                                                                                                                                                                                   pos0=pos1;
                                                                                                                                                                                                                                   pos0x=1;
                                                                                                                                                                                                                                   pos0xl=1;
                                                                                                                                                                                                                                   end;
                                                                                                                                                                                                                                   else
                                                                                                                                                                                                                                     do;
                                                                                                                                                                                                                                   pos0x=(parg/(parg+exp(X[pos0, ]*bg+offset[pos0, ])))##parg;
                                                                                                                                                                                                                                   pos0xl=(parg/(parg+y[pos0, ]))##parg;
                                                                                                                                                                                                                                   end;
                                                                                                                                                                                                                                   ll=sum(-log(0+exp(G[pos0, ]*lambdag))+log(0*exp(G[pos0, 
                                                                                                                                                                                                                                   ]*lambdag)+pos0x))+
                                                                                                                                                                                                                                     sum(-log(0+exp(G[pos1, ]*lambdag))+lgamma(parg+y[pos1, ])-lgamma(y[pos1, 
                                                                                                                                                                                                                                     ]+1)-lgamma(parg)+
                                                                                                                                                                                                                                       y[pos1]#log(exp(X[pos1, ]*bg+offset[pos1, ])/(parg+exp(X[pos1, 
                                                                                                                                                                                                                                     ]*bg+offset[pos1, ])))+parg*log(parg/(parg+exp(X[pos1, ]*bg+offset[pos1, 
                                                                                                                                                                                                                                     ]))));
                                                                                                                                                                                                     llnull1=sum(-log(1+zkg)+log(zkg+pos0xl))+
                                                                                                                                                                                                       sum(-log(1+zkg)+lgamma(parg+y[pos1, ])-lgamma(y[pos1, ]+1)-lgamma(parg)+
                                                                                                                                                                                                             y[pos1]#log(y[pos1, ]/(parg+y[pos1, ]))+parg*log(parg/(parg+y[pos1, ])));
                                                                                                                                                                                                           devg=2*(llnull1-ll);
                                                                                                                                                                                                           AIC=-2*ll+2*nvar;
                                                                                                                                                                                                           AICc=-2*ll+2*nvar*(n/(n-nvar-1));
                                                                                                                                                                                                           tt2=y/y[:];
                                                                                                                                                                                                           tt2=choose(tt2=0, 1E-10, tt2);
                                                                                                                                                                                                           devnullg=2*sum(y#log(tt2)-(y-y[:]));
                                                                                                                                                                                                                          pctdevg=1-devg/devnullg;
                                                                                                                                                                                                                          adjpctdevg=1-((n-1)/(n-nvar))*(1-pctdevg);
                                                                                                                                                                                                                          
                                                                                                                                                                                                                          %IF %UPCASE(&MODEL)=NEGBIN %THEN
                                                                                                                                                                                                                          %DO;
                                                                                                                                                                                                                          AIC=-2*ll+2*(nvar+1);
                                                                                                                                                                                                                          AICc=-2*ll+2*(nvar+1)*(n/(n-(nvar+1)-1));
                                                                                                                                                                                                                          tt2=y/y[:];
                                                                                                                                                                                                                          tt2=choose(tt2=0, 1E-10, tt2);
                                                                                                                                                                                                                          devnullg=2*sum(y#log(tt2)-(y+1/alphag)#log((1+alphag#y)/(1+alphag#y[:])));
                                                                                                                                                                                                                                         pctdevg=1-devg/devnullg;
                                                                                                                                                                                                                                         adjpctdevg=1-((n-1)/(n-(nvar+1)))*(1-pctdevg);
                                                                                                                                                                                                                                         %END;
                                                                                                                                                                                                                                         print devg[label='Deviance'] ll[label='Full Log Likelihood'] pctdevg 
                                                                                                                                                                                                                                         adjpctdevg AIC AICc;
                                                                                                                                                                                                                                         %END;
                                                                                                                                                                                                                                         print 'Variance-Covariance Matrix', , varcovg[label=''];
                                                                                                                                                                                                                                         
                                                                                                                                                                                                                                         /****************************************/
                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                           %if &grid=%then
                                                                                                                                                                                                                                         %do;
                                                                                                                                                                                                                                         pihat=exp(pihat)/(1+exp(pihat));
                                                                                                                                                                                                                                         create _res_ var{wt y yhat res resstd influence cooksD sumwi pihat};
                                                                                                                                                                                                                                         append;
                                                                                                                                                                                                                                         create _beta_ from bi[colname={"id" "B" "x" "y"}];
                                                                                                                                                                                                                                         append from bi;
                                                                                                                                                                                                                                         bistdt=bi||stdbi||tstat||probt;
                                                                                                                                                                                                                                         create _parameters_ from bistdt[colname={"id" "B" "x" "y" "stdbi" "tstat" 
                                                                                                                                                                                                                                           "probt"}];
                                                                                                                                                                                                                                         append from bistdt;
                                                                                                                                                                                                                                         _tstat_=_beta_;
                                                                                                                                                                                                                                         
                                                                                                                                                                                                                                         do j=1 to nrow(_stdbeta_);
                                                                                                                                                                                                                                         
                                                                                                                                                                                                                                         do k=1 to ncol(_stdbeta_);
                                                                                                                                                                                                                                         
                                                                                                                                                                                                                                         if _stdbeta_[j, k]=0 then
                                                                                                                                                                                                                                         _tstat_[j, k]=0;
                                                                                                                                                                                                                                         else
                                                                                                                                                                                                                                           _tstat_[j, k]=_beta_[j, k]/_stdbeta_[j, k];
                                                                                                                                                                                                                                         end;
                                                                                                                                                                                                                                         end;
                                                                                                                                                                                                                                         _tstat_=choose(_tstat_=., 0, _tstat_);
                                                                                                                                                                                                                                         _probt_=2*(1-probt(abs(_tstat_), df));
                                                                                                                                                                                                                                         _sig_=j(n, ncol(x), "not significant at 90%");
                                                                                                                                                                                                                                         
                                                                                                                                                                                                                                         do i=1 to n;
                                                                                                                                                                                                                                         
                                                                                                                                                                                                                                         do j=1 to ncol(x);
                                                                                                                                                                                                                                         
                                                                                                                                                                                                                                         if _probt_[i, j]<0.01*(nvar/v1) then
                                                                                                                                                                                                                                         _sig_[i, j]="significant at 99%";
                                                                                                                                                                                                                                         else if _probt_[i, j]<0.05*(nvar/v1) then
                                                                                                                                                                                                                                         _sig_[i, j]="significant at 95%";
                                                                                                                                                                                                                                         else if _probt_[i, j]<0.1*(nvar/v1) then
                                                                                                                                                                                                                                         _sig_[i, j]="significant at 90%";
                                                                                                                                                                                                                                         else
                                                                                                                                                                                                                                           _sig_[i, j]="not significant at 90%";
                                                                                                                                                                                                                                         
                                                                                                                                                                                                                                         if _probt_[i, j]=. then
                                                                                                                                                                                                                                         _sig_[i, j]="not significant at 90%";
                                                                                                                                                                                                                                         end;
                                                                                                                                                                                                                                         end;
                                                                                                                                                                                                                                         _bistdt_=COORD||_beta_||_stdbeta_||_tstat_||_probt_;
                                                                                                                                                                                                                                         _colname1_={"Intercept" &xvar};
                                                                                                                                                                                                                                         _label_=repeat("std_", ncol(x))//repeat("tstat_", 
                                                                                                                                                                                                                                                                                 ncol(x))//repeat("probt_", ncol(x));
                                                                                                                                                                                                                                         _colname_={"x" "y"}||_colname1_||concat(_label_, repeat(_colname1_`, 3))`;
                                                                                                                                                                                                                                                                                                 call change(_colname_, "_ ", "_");
                                                                                                                                                                                                                                                                                                 call change(_colname_, "_ ", "_");
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 %if %upcase(&MODEL)=ZIP or %upcase(&MODEL)=ZINB %then
                                                                                                                                                                                                                                                                                                 %do;
                                                                                                                                                                                                                                                                                                 _tstatl_=_lambda_;
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 do j=1 to nrow(_stdlambda_);
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 do k=1 to ncol(_stdlambda_);
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 if _stdlambda_[j, k]=0 then
                                                                                                                                                                                                                                                                                                 _tstatl_[j, k]=0;
                                                                                                                                                                                                                                                                                                 else
                                                                                                                                                                                                                                                                                                   _tstatl_[j, k]=_lambda_[j, k]/_stdlambda_[j, k];
                                                                                                                                                                                                                                                                                                 end;
                                                                                                                                                                                                                                                                                                 end;
                                                                                                                                                                                                                                                                                                 _probtl_=2*(1-probt(abs(_tstatl_), df));
                                                                                                                                                                                                                                                                                                 _sigl_=j(n, ncol(G), "not significant at 90%");
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 do i=1 to n;
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 do j=1 to ncol(G);
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 if _probtl_[i, j]<0.01*((nvar+ncol(G))/v1) then
                                                                                                                                                                                                                                                                                                 _sigl_[i, j]="significant at 99%";
                                                                                                                                                                                                                                                                                                 else if _probtl_[i, j]<0.05*((nvar+ncol(G))/v1) then
                                                                                                                                                                                                                                                                                                 _sigl_[i, j]="significant at 95%";
                                                                                                                                                                                                                                                                                                 else if _probtl_[i, j]<0.1*((nvar+ncol(G))/v1) then
                                                                                                                                                                                                                                                                                                 _sigl_[i, j]="significant at 90%";
                                                                                                                                                                                                                                                                                                 else
                                                                                                                                                                                                                                                                                                   _sigl_[i, j]="not significant at 90%";
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 if _probtl_[i, j]=. then
                                                                                                                                                                                                                                                                                                 _sigl_[i, j]="not significant at 90%";
                                                                                                                                                                                                                                                                                                 end;
                                                                                                                                                                                                                                                                                                 end;
                                                                                                                                                                                                                                                                                                 _bistdt_=_bistdt_||_lambda_||_stdlambda_||_tstatl_||_probtl_;
                                                                                                                                                                                                                                                                                                 _colname2_={&xvarinf};
                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                 %IF %UPCASE(&INT_INF)=YES %THEN
                                                                                                                                                                                                                                                                                                 %DO;
                                                                                                                                                                                                                                                                                                 _colname2_={"Intercept" &xvarinf};
                                                                                                                                                                                                                                                                                                 %END;
                                                                                                                                                                                                                                                                                                 _label3_=repeat("Inf_", ncol(G));
                                                                                                                                                                                                                                                                                                 _colname3_=concat(_label3_, _colname2_`);
                                                                                                                                                                                                                                                                                                 _label2_=repeat("Inf_std_", ncol(G))//repeat("Inf_tstat_", 
                                                                                                                                                                                                                                                                                                                                              ncol(G))//repeat("Inf_probt_", ncol(G));
                                                                                                                                                                                                                                                                                                 _colname_={"x" "y"}||_colname1_||concat(_label_, repeat(_colname1_`, 
                                                                                                                                                                                                                                                                                                                                                         3))`||_colname3_`||concat(_label2_, repeat(_colname2_`, 3))`;
                                                                                                                                                                                                                                                                                                                                                                                                    call change(_colname_, "_ ", "_");
                                                                                                                                                                                                                                                                                                                                                                                                    call change(_colname_, "_ ", "_");
                                                                                                                                                                                                                                                                                                                                                                                                    _labell_=repeat("sig_Inf_", ncol(G));
                                                                                                                                                                                                                                                                                                                                                                                                    _colnamel_=concat(_labell_, repeat(_colname2_`, 1))`;
                                                                                                                                                                                                                                                                                                                                                                                                                                       create _sig_inf_parameters2_ from _sigl_[colname=_colnamel_];
                                                                                                                                                                                                                                                                                                                                                                                                                                       append from _sigl_;
                                                                                                                                                                                                                                                                                                                                                                                                                                       %end;
                                                                                                                                                                                                                                                                                                                                                                                                                                       create _parameters2_ from _bistdt_[colname=_colname_];
                                                                                                                                                                                                                                                                                                                                                                                                                                       append from _bistdt_;
                                                                                                                                                                                                                                                                                                                                                                                                                                       close _parameters2_;
                                                                                                                                                                                                                                                                                                                                                                                                                                       _label_=repeat("sig_", ncol(x));
                                                                                                                                                                                                                                                                                                                                                                                                                                       _colname_=concat(_label_, repeat(_colname1_`, 1))`;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        create _sig_parameters2_ from _sig_[colname=_colname_];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        append from _sig_;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        %IF %UPCASE(&MODEL)=NEGBIN or %UPCASE(&MODEL)=ZINB %THEN
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        %DO;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        atstat=alphai[, 2];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        do j=1 to nrow(alphai);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        if alphai[j, 3]=0 then
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        atstat[j]=0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        else
                                                                                                                                                                                                                                                                                                                                                                                                                                                                          atstat[j]=alphai[j, 2]/alphai[j, 3];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        end;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        atstat=choose(atstat=., 0, atstat);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        aprobtstat=2*(1-probnorm(abs(atstat)));
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        _siga_=j(n, 1, "not significant at 90%");
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        do i=1 to n;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        if aprobtstat[i]<0.01*(nvar/v1) then
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        _siga_[i]="significant at 99%";
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        else if aprobtstat[i]<0.05*(nvar/v1) then
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        _siga_[i]="significant at 95%";
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        else if aprobtstat[i]<0.1*(nvar/v1) then
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        _siga_[i]="significant at 90%";
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        else
                                                                                                                                                                                                                                                                                                                                                                                                                                                                          _siga_[i]="not significant at 90%";
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        end;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        alphai=alphai||atstat||aprobtstat;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        create _alpha_ from alphai[colname={"id" "alpha" "std" "tstat" 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "probt"}];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        append from alphai;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        create _sig_alpha_ from _siga_[colname={"sig_alpha"}];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        append from _siga_;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        %END;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        %end;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        %else
                                                                                                                                                                                                                                                                                                                                                                                                                                                                          %do;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        create _beta_ from bi[colname={"id" "B" "x" "y"}];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        append from bi;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        stdbi=sqrt(abs(varbi));
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        tstat=bi[, 2]/stdbi;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        do i=1 to nrow(tstat);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        if tstat[i]=. then
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        tstat[i]=0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        end;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        bistdt=bi||stdbi||tstat;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        create _parameters_ from bistdt[colname={"id" "B" "x" "y" "stdbi" 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "tstat"}];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        append from bistdt;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        close _parameters_;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        %IF %UPCASE(&MODEL)=NEGBIN or %UPCASE(&MODEL)=ZINB %THEN
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        %DO;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        atstat=alphai[, 2];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        do j=1 to nrow(alphai);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        if alphai[j, 3]=0 then
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        atstat[j]=0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        else
                                                                                                                                                                                                                                                                                                                                                                                                                                                                          atstat[j]=alphai[j, 2]/alphai[j, 3];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        end;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        atstat=choose(atstat=., 0, atstat);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        aprobtstat=2*(1-probnorm(abs(atstat)));
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        alphai=POINTS||alphai||atstat||aprobtstat;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        create _alpha_ from alphai[colname={"x" "y" "id" "alpha" "std" "tstat" 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "probt"}];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        append from alphai;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        %END;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        _beta_=shape(bi[, 2], m);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        _stdbeta_=shape(stdbi, m);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        _tstat_=_beta_/_stdbeta_;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        _bistdt_=POINTS||_beta_||_stdbeta_||_tstat_;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        _colname1_={"Intercept" &xvar};
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        _label_=repeat("std_", nvar)//repeat("tstat_", nvar);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        _colname_={"x" "y"}||_colname1_||concat(_label_, repeat(_colname1_`, 2))`;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                call change(_colname_, "_ ", "_");
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                call change(_colname_, "_ ", "_");
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                create _parameters_grid_ from _bistdt_[colname=_colname_];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                append from _bistdt_;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                close _parameters2_;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                %end;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                quit;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                %if &grid=%then
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                %do;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                data _parameters2_;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                merge _parameters2_ _sig_parameters2_;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                run;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                %if %upcase(&MODEL)=ZIP or %upcase(&MODEL)=ZINB %then
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                %do;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                data _parameters2_;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                merge _parameters2_ _sig_inf_parameters2_;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                run;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                %end;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                %IF %UPCASE(&MODEL)=NEGBIN or %UPCASE(&MODEL)=ZINB %THEN
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                %DO;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                data _alpha_;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                merge _alpha_ _sig_alpha_;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                run;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                %end;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                %END;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                %mend GWZINBR;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                /*****************************************************************************************/
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  /* Macro to separe the datasets of GWR Grid Estimates*/
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  /* REQUIRED PARAMETERS
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                /*     PAR = the number of explicative variables without considers the intercept
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                /*****************************************************************************************/
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  %macro MAP(PAR=);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                proc iml;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                use _parameters_;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                read all into p[colname=names];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                var=1+&par;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                n=nrow(p)/2;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                nvar=ncol(p);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                _beta_=shape(p, n);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                label=names;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                do i=1 to var;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                seq1=(i-1)*nvar+1;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                seq2=seq1+(nvar-1);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                name="B0":"B&par";
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                _betas_=_beta_[, seq1:seq2];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                create b from _betas_[colname=label];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                append from _betas_;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                close b;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                sastables=datasets("work");
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                do j=1 to nrow(sastables);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                if sastables[j]=name[i] then
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                call delete(name[i]);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                end;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                call rename(b, name[i]);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                end;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                quit;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                %mend MAP;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                /*****************************************************************************************/
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  /* Macro for Drawing the Maps in 2D */
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  /* REQUIRED PARAMETERS
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                /*  TABLE = the name of the SAS data set to be used
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                /*    VAR = the name of the variable to be plotted
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                /*    MAP = the name of the SAS data set with the geographic coordinates
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                /*  SAMPLEDATA = the name of the SAS data set with the geographic coordinates of the sample
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                /*               data
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                /* METHOD = there are two choices to define the classes:
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  /*          EQUAL: using the same rule of histograms (the number of classes is defined by
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        Sturges's ruke)
			/*            STD: using the standard deviation for creating 6 classes: media-3std;
			/*                 media-2std;media-1std;media+1std;media+2std;media+3std;
			/*  WHERE = using a conditional to draw the maps (for example, plot only the parameters that
			/*          are significant at 5% level)
			/*****************************************************************************************/
			%macro map2d(table=, var=, map=, id=, idtype=N, sampledata=, method=EQUAL, 
					where=);
			data anno&table;
				set &table;
				length function style $10. color $8.;
				retain line 1 xsys ysys '2' hsys '3' color 'red';
				function='label';
				text='U';
				position='5';
				style='marker';
				size=2;
				&where;
			run;

			proc sql noprint;
				select count(*) into:np from anno&table where &var>=0;
				select count(*) into:nn from anno&table where &var<0;
				select min(&var) into:minp from anno&table where &var>=0;
				select max(&var) into:maxp from anno&table where &var>=0;
				select min(&var) into:minn from anno&table where &var<0;
				select max(&var) into:maxn from anno&table where &var<0;
			quit;

			%put &np &nn;

			%if &np >0 and &nn=0 %then
				%do;

					data _null_;
						np=floor(1+3.3*log10(&np))-1;
						call symput('np', np);
					run;

					%put &np;

					%if %upcase(&method)=EQUAL %then
						%do;
							%put &np &minp &maxp %sysevalf((&maxp-&minp)/&np);
							%colorscale(FFFFFF, , FF3333, &np, clist, no);
							%patt;

							data _null_;
								set clist;
								call symput('color'||trim(left(_n_)), 'cx'||rgb);
							run;

							%macro cl;

								data anno&table;
									set anno&table;

									if &var<=&minp+(&maxp-&minp)/&np then
										color="&color1";

									%do i=2 %to &np;
										else if &var<=&minp+&i*(&maxp-&minp)/&np then
											color="&&color&i";
									%end;
								run;

							%mend cl;

							%cl;
						%end;
					%else %if %upcase(&method)=STD %then
						%do;

							data _null_;
								np=6;
								call symput('np', np);
							run;

							proc sql noprint;
								select mean(&var) into:meanp from anno&table where &var>=0;
								select std(&var) into:stdp from anno&table where &var>=0;
							quit;

							%put &meanp &stdp;
							%colorscale(FFFFFF, , FF3333, &np, clist, no);
							%patt;

							data _null_;
								set clist;
								call symput('color'||trim(left(_n_)), 'cx'||rgb);
							run;

							%macro cl;

								data anno&table;
									set anno&table;

									if &var<=&minp+(&meanp-3*&stdp) then
										color="&color1";
									else if &var<=&minp+&meanp-2*&stdp then
										color="&color2";
									else if &var<=&minp+&meanp-&stdp then
										color="&color3";
									else if &var<=&minp+&meanp+&stdp then
										color="&color4";
									else if &var<=&minp+&meanp+2*&stdp then
										color="&color5";
									else
										color="&color6";
								run;

							%mend cl;

							%cl;
						%end;

					data _null_;
						min=left(trim(putn(round(&minp, 0.1), 'commax10.')));
						max=left(trim(putn(round(&maxp, 0.1), 'commax10.')));
						call symput('min', min);
						call symput('max', max);
					run;

					%put &min &max;
					%bar(FF3333, FFFFFF, &minp, &maxp, vertical, y_i=30, x_i=90);
				%end;
			%else %if &nn > 0 and &np=0 %then
				%do;

					data _null_;
						nn=floor(1+3.3*log10(&nn))-1;
						call symput('nn', nn);
					run;

					%put &nn;
					%put &nn &minn &maxn %sysevalf((&maxn-&minn)/&nn);
					%colorscale(FFFFFF, , 3333FF, &nn, clist, no);
					%patt;

					data _null_;
						set clist;
						call symput('color'||trim(left(_n_)), 'cx'||rgb);
					run;

					%macro cl;

						data anno&table;
							set anno&table;

							if &var<=&minn+(&maxn-&minn)/&nn then
								color="&color1";

							%do i=2 %to &nn;
								else if &var<=&minn+&i*(&maxn-&minn)/&nn then
									color="&&color&i";
							%end;
						run;

					%mend cl;

					%cl;

					data _null_;
						min=left(trim(putn(round(&minn, 0.1), 'commax10.')));
						max=left(trim(putn(round(&maxn, 0.1), 'commax10.')));
						call symput('min', min);
						call symput('max', max);
					run;

					%put &min &max;
					%bar(3333FF, FFFFFF, &minn, &maxn, vertical, y_i=30, x_i=90);
				%end;
			%else
				%do;

					data _null_;
						np=floor(1+3.3*log10(&np))-1;
						nn=floor(1+3.3*log10(&nn))-1;
						call symput('np', np);
						call symput('nn', nn);
					run;

					%put &np &nn;
					%put &np &minp &maxp %sysevalf((&maxp-&minp)/&np);
					%put &nn &minn &maxn %sysevalf((&maxn-&minn)/&nn);
					%colorscale(FFFFFF, , FF3333, &np, clist, no);
					%patt;

					data _null_;
						set clist;
						call symput('color'||trim(left(_n_)), 'cx'||rgb);
					run;

					%macro cl;

						data anno&table.p;
							set anno&table(where=(&var>=0));

							if &var<=&minp+(&maxp-&minp)/&np then
								color="&color1";

							%do i=2 %to &np;
								else if &var<=&minp+&i*(&maxp-&minp)/&np then
									color="&&color&i";
							%end;
						run;

					%mend cl;

					%cl;
					%colorscale(3333FF, , FFFFFF, &nn, clist, no);
					%patt;

					data _null_;
						set clist;
						call symput('color'||trim(left(_n_)), 'cx'||rgb);
					run;

					%macro cl;

						data anno&table.n;
							set anno&table(where=(&var<0));

							if &var<=&minn+(&maxn-&minn)/&nn then
								color="&color1";

							%do i=2 %to &nn;
								else if &var<=&minn+&i*(&maxn-&minn)/&nn then
									color="&&color&i";
							%end;
						run;

					%mend cl;

					%cl;

					data anno&table;
						set anno&table.p anno&table.n;
					run;

					data _null_;
						min=left(trim(putn(round(&minn, 0.1), 'commax10.')));
						max=left(trim(putn(round(&maxp, 0.1), 'commax10.')));
						call symput('min', min);
						call symput('max', max);
					run;

					%put &min &max;
					%bar(FF3333, FFFFFF, 0, &maxp, vertical, y_i=50, x_i=90);

					data anno&table;
						set _a_ anno&table;
					run;

					%bar(FFFFFF, 3333FF, &minn, 0, vertical, y_i=16, x_i=90);
				%end;

			data anno&table;
				set _a_ anno&table;
			run;

			data a;
				%if %upcase(&idtype)=N %then
					%do;
						&id=0;
					%end;
				%else
					%do;
						&id='0';
					%end;
				v=2;
			run;

			%if &sampledata ne %then
				%do;

					data annodata;
						set &sampledata;
						length function style $10. color $8.;
						retain line 1 xsys ysys '2' hsys '3' color 'black' when 'a';
						function='label';
						text='J';
						position='5';
						style='special';
						size=2;
					run;

				%end;
			*goptions reset=all reset=global;

			proc gmap data=a map=&map all 
				%if &sampledata ne %then
					%do;
						anno=annodata%end;
				;
				id &id;
				choro v / anno=anno&table nolegend;
				run;
			quit;

		%mend map2d;
