Golden <- function(data, formula, xvarinf, weight,
                   lat, long, output, globalmin=TRUE,
                   method, model="zinb", bandwidth="cv", offset, 
                   force=TRUE, maxg=100, distancekm=FALSE){
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
    #G <<- matrix(1, N, 1)
    G <<- rep(1, N) #matriz coluna
    #lambdag <<- matrix(0, length(G), 1) #ncol(G) em vez de length(G)
    lambdag <<- rep(0, length(G))
  }
  else{
    #G[colname=varnamezip]
    G <<- data[, xvarinf]
    G <<- cbind(rep(1, N), G)
  }
  #wt <<- matrix(1, N, 1)
  wt <<- rep(1, N)
  if (!is.null(weight)){
    wt <<- data[, weight]
  }
  #Offset <- matrix(0, N, 1)
  Offset <<- rep(0, N)
  if (!is.null(offset)){
    Offset <<- data[, offset]
  }
  x <- cbind(rep(1, N), x)
  nvar <<- ncol(x)
  #yhat <<- matrix(0, N, 1)
  yhat <<- rep(0, N)
  #yhat2 <<- matrix(0, N, 1)
  yhat2 <<- rep(0, N)
  #pihat <<- matrix(0, N, 1)
  pihat <<- rep(0, N)
  #alphai <<- matrix(0, N, 1)
  alphai <<- rep(0, N)
  #S <<- matrix(0, N, 1)
  S <<- rep(0, N)
  #Si <<- matrix(0, N, 1)
  Si <<- rep(0, N)
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
        hess <- sum(trigamma(parg+y)-trigamma(parg)+1/parg-2/(parg+uj)+(y+parg)/(parg+uj)^2)
        hess <- ifelse(hess==0, E^-23, hess)
        par0 <- parg
        #parg <<- par0-solve(hess)%*%gf
        parg <<- par0-as.vector(solve(hess))*gf
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
        bg <<- rep(0,ncol(x))
      } 
      else{
        bg <<- solve(t(x)%*%(Ai*x))%*%t(x)%*%(Ai*zj)
      }  
      nj <- x%*%bg+Offset
      nj <- ifelse(nj>700,700,nj)
      uj <- exp(nj)
      olddev <- devg
      uj <- ifelse(uj<E^-150,E^-150,uj)
      uj <- ifelse(uj>100000,100000,uj)
      tt <- y/uj
      tt <- ifelse(tt=0,E^-10,tt)
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
    lambda0 <- (ncol(pos0)-sum((parg/(uj+parg))^parg))/n
    if (lambda0 <= 0) {
      #lambdag <<- matrix(0, ncol(G),1) 
      lambdag <<- rep(0, ncol(G))
    }
    else {
      lambda0 <<- log(lambda0/(1-lambda0))
      #lambdag <<- rbind(lambda0, matrix(0,ncol(G)-1,1))
      lambdag <<- rbind(lambda0, rep(0,ncol(G)-1))
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
  njl <- G*lambdag
  if (model!="zip" & model!="zinb"){
    zkg <- 0
  }
  else{
    zkg <- 1/(1+exp(-G*lambdag)*(parg/(parg+uj))^parg)
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
        parg <<- 1/alphag
      }
      else{
        if (j>0){
          parg <<- 1/(sum((y-uj)^2/uj)/(N-nvar))
        }
        while (abs(dpar)>0.0001 & aux2<200){
          if (parg<0){
            parg <<- 0.00001
          }
          parg <<- ifelse(parg<E^-10, E^-10, parg)
          gf <- sum((1-zkg)*(digamma(parg+y)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+y)/(parg+uj)))
          hess <- sum((1-zkg)*(trigamma(parg+y)-trigamma(parg)+1/parg-2/(parg+uj)+(y+parg)/(parg+uj)^2))
          hess <- ifelse(hess==0, E^-23, hess)
          par0 <- parg
          parg <<- par0-solve(hess)%*%gf #multiplicador
          if (aux2>50 & parg>E^5){
            dpar <- 0.0001
            cont3 <- cont3+1
            if (cont3==1){
              parg <<- 2
            }
            else if (cont3==2){
              parg <<- E^5
            }
            else if (cont3==3){
              parg <<- 0.0001
            }
          }
          else{
            dpar <- parg-par0
          }
          if (parg>E^6){
            parg <<- E^6
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
        Ai <- (1-zkg)*((uj/(1+alphag*uj)+(y-uj)*(alphag*uj/(1+2*alphag*uj+alphag**2*uj^2))))
        Ai <- ifelse(Ai<=0, E^-5, Ai)
        uj <- ifelse(uj<E^-150, E^-150, uj)
        zj <- (nj+(y-uj)/(((uj/(1+alphag*uj)+(y-uj)*(alphag*uj/(1+2*alphag*uj+alphag^2*uj^2))))*(1+alphag*uj)))-Offset
        if (det(t(x)%*%(Ai*x))==0){
          #bg <<- matrix(0, nvar, 1)
          bg <<- rep(0, nvar)
        }
        else{
          bg <<- solve(t(x)%*%(Ai*x))%*%t(x)%*%(Ai*zj)
        }
        nj <- x%*%bg+Offset
        nj <- ifelse(nj>700, 700, nj)
        nj <- ifelse(nj<-700, -700, nj)
        uj <- exp(nj)
        olddev <- devg
        gamma1 <- (uj/(uj+parg))^y*(parg/(uj+parg))^parg #(gamma(par+y)/(gamma(y+1)#gamma(par)))#
        gamma1 <- ifelse(gamma1<=0, E^-10, gamma1)
        devg <- sum((1-zkg)*(log(gamma1)))
        ddev <- devg-olddev
        #print(c('bg', 'aux1', 'devg', 'olddev', 'ddev'))
        #print(c(bg aux1 devg olddev ddev))
        aux1 <- aux1+1
      }
      ddpar <- parg-parold
      cont <- cont+1
    }
    if (model == "zip" |model == 'zinb'){
      devg <- 0
      ddev <- 1
      njl <- G*lambdag
      njl <- ifelse(njl > maxg, maxg, njl)
      njl <- ifelse(njl < maxg,-maxg, njl)
      pig <- exp(njl)/(1+exp(njl))
      while (abs(ddev)>0.000001 & aux3<100){
        Ai <- pig*(1-pig)
        Ai <- ifelse(Ai<=0, E^-5, Ai)
        zj <- njl+(zkg-pig)*1/Ai
        if (det(t(G*Ai)*G)==0){
          #lambdag <<- matrix(0, ncol(G), 1)
          lambdag <<- rep(0, ncol(G))
        }  
        else {
          lambdag <<- solve(t(G*Ai)%*%G)%*%t(G*Ai)*zj 
        }
        njl <- G*lambdag
        njl <- ifelse(njl > maxg, maxg, njl)
        njl <- ifelse(njl < maxg,-maxg, njl)
        pig <- exp(njl)/(1+exp(njl))
        olddev <- devg
        devg <- sum(zkg*njl-log(1+exp(njl)))
        ddev <- devg-olddev
        # print(c('lambdag', 'devg', 'olddev', 'ddev')) 
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
    # print(c('j', 'bg', 'alphag', 'lambdag', 'llikeg', 'dllike')) 
    j <- j+1
  }
  long <- data[, long]
  lat <- data[, lat]
  COORD <<- matrix(c(long, lat), ncol=2)
  distance <<- dist(COORD,"euclidean") # _dist_ <- distance
  sequ <<- 1:N # seq <- sequ
  #funcao cv aqui
  
  # DEFINING GOLDEN SECTION SEARCH PARAMETERS 
  if(method=="fixed_g" | method=="fixed_bsq"){
    ax <- 0
    bx <- as.integer(max(distance)+1)
    if(distancekm=="yes"){
      bx <- bx*111
    }
  }  
  if(method=="adaptive_bsq"){
    ax <- 5
    bx <- n
  }
  r <- 0.61803399
  tol <- 0.1
  %IF %UPCASE(&globalmin)=NO %THEN %DO;
  lower=ax;
  upper=bx;
  xmin=j(1,2,0);
  do GMY=1 to 1;
  ax1=lower[GMY];
  bx1=upper[GMY];
  %END;
  %ELSE %DO;
  lower=ax||(1-r)*bx||r*bx;
  upper=(1-r)*bx||r*bx||bx;
  xmin=j(3,2,0);
  do GMY=1 to 3;
  ax1=lower[GMY];
  bx1=upper[GMY];
  %END;
  h0=ax1;
  h3=bx1;
  h1=bx1-r*(bx1-ax1);
  h2=ax1+r*(bx1-ax1);
  print h0 h1 h2 h3;
  /***************************************/
} # fecha golden 

                                                                                                                                                                                     
/*****************************************/
 

 

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
/**** DEFINING GOLDEN SECTION SEARCH PARAMETERS *****/
 %IF %UPCASE(&method)=FIXED_G or %UPCASE(&method)=FIXED_BSQ %THEN %DO;
ax=0;
bx=int(max(_dist_)+1);
%IF %UPCASE(&distancekm)=YES %THEN %DO;bx=bx*111;%END;
%END;
%IF %UPCASE(&method)=ADAPTIVE_BSQ %THEN %DO;
ax=5;
bx=n;
%END;
r=0.61803399;
tol=0.1;
%IF %UPCASE(&globalmin)=NO %THEN %DO;
lower=ax;
upper=bx;
xmin=j(1,2,0);
do GMY=1 to 1;
ax1=lower[GMY];
bx1=upper[GMY];
%END;
%ELSE %DO;
lower=ax||(1-r)*bx||r*bx;
upper=(1-r)*bx||r*bx||bx;
xmin=j(3,2,0);
do GMY=1 to 3;
ax1=lower[GMY];
bx1=upper[GMY];
%END;
h0=ax1;
h3=bx1;
h1=bx1-r*(bx1-ax1);
h2=ax1+r*(bx1-ax1);
print h0 h1 h2 h3;
/***************************************/
 
 res1=cv(h1);CV1=res1[1];
res2=cv(h2);CV2=res2[1];

%IF %UPCASE(&method)=FIXED_G or %UPCASE(&method)=FIXED_BSQ or %UPCASE(&method)=ADAPTIVE_BSQ %THEN %DO;
if GMY=1 then create &output var{GMY h1 cv1 h2 cv2};
%END;

int=1;
do while(abs(h3-h0) > tol*(abs(h1)+abs(h2)) & int<200);
if CV2<CV1 then do;
h0=h1;
h1=h3-r*(h3-h0);
h2=h0+r*(h3-h0);
CV1=CV2;
res2=cv(h2);CV2=res2[1];
end;
else do;
h3=h2;
h1=h3-r*(h3-h0);
h2=h0+r*(h3-h0);
CV2=CV1;
res1=cv(h1);CV1=res1[1];
end;
%IF %UPCASE(&method)=FIXED_G or %UPCASE(&method)=FIXED_BSQ or %UPCASE(&method)=ADAPTIVE_BSQ %THEN %DO;
append;
%END;
int=int+1;
end;
if CV1<CV2 then do;
golden=CV1;
xmin[GMY,1]=golden;
xmin[GMY,2]=h1;
npar=res1[2];
%IF %UPCASE(&method)=ADAPTIVE_BSQ %THEN %DO;
xmin[GMY,2]=floor(h1);
%end;
end;
else do;
golden=CV2;
xmin[GMY,1]=golden;
xmin[GMY,2]=h2;
npar=res2[2];
%IF %UPCASE(&method)=ADAPTIVE_BSQ %THEN %DO;
xmin[GMY,2]=floor(h2);
%end;
end;
%IF %UPCASE(&method)=FIXED_G or %UPCASE(&method)=FIXED_BSQ or %UPCASE(&method)=ADAPTIVE_BSQ %THEN %DO;
print golden (xmin[GMY,2])[label='xmin'] %if %UPCASE(&bandwidth)=AIC %THEN %DO;npar%END;;
%END;
%ELSE %IF %UPCASE(&method)=ADAPTIVEN %THEN %DO;
hv[1]=xmin[GMY,2];
append from hv;
end;
%END;
end;
create _min_bandwidth_ from xmin[colname={'golden' 'bandwidth'}];
append from xmin;
%IF %UPCASE(&globalmin)=YES %THEN %DO;
print xmin[colname={'golden' 'bandwidth'}];
xming=xmin[loc(xmin[,1]=min(xmin[,1])),2];
print 'Global Minimum', xming[label='(Da Silva and Mendes, 2018)'];
%END;
quit;

%IF %UPCASE(&method)=FIXED_G or %UPCASE(&method)=FIXED_BSQ or %UPCASE(&method)=ADAPTIVE_BSQ %THEN %DO;
%global _h_;
proc sql noprint;
select bandwidth into:_h_ from _min_bandwidth_
having golden=min(golden);
quit;
%put h=&_h_;
%END;
%mend Golden;
