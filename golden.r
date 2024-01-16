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
  y <- model.extract(mf, "response")
  N <- length(y)
  x <- model.matrix(mt, mf)
  if (is.null(xvarinf)){
    #G <- matrix(1, N, 1)
    G <- rep(1, N) #matriz coluna
    #lambdag <- matrix(0, length(G), 1) #ncol(G) em vez de length(G)
    lambdag <- rep(0, length(G))
  }
  else{
    #G[colname=varnamezip]
    G <- data[, xvarinf]
    G <- cbind(rep(1, N), G)
  }
  #wt <- matrix(1, N, 1)
  wt <- rep(1, N)
  if (!is.null(weight)){
    wt <- data[, weight]
  }
  #Offset <- matrix(0, N, 1)
  Offset <- rep(0, N)
  if (!is.null(offset)){
    Offset <- offset
  }
  x <- cbind(rep(1, N), x)
  nvar <- ncol(x)
  #yhat <- matrix(0, N, 1)
  yhat <- rep(0, N)
  #yhat2 <- matrix(0, N, 1)
  yhat2 <- rep(0, N)
  #pihat <- matrix(0, N, 1)
  pihat <- rep(0, N)
  #alphai <- matrix(0, N, 1)
  alphai <- rep(0, N)
  #S <- matrix(0, N, 1)
  S <- rep(0, N)
  #Si <- matrix(0, N, 1)
  Si <- rep(0, N)
  Iy <- ifelse(y>0, 1, y)
  Iy <- 1-Iy
  pos0 <- which(y==0)
  pos02 <- which(y==0)
  pos1 <- which(y>0)
  
  # /**** global estimates ****/
  uj <- (y+mean(y))/2 
  nj <- log(uj)
  parg <- sum((y-uj)^2/uj)/(N-nvar)
  ddpar <- 1
  cont <- 1
  cont3 <- 0
  while (abs(ddpar)>0.000001 & cont<100){
    dpar <- 1
    parold <- parg
    cont1 <- 1
    if (model == "zip" | model == "poisson"){
      parg <- 1/E^(-6)
      alphag <- 1/parg
    }  
    if (model == "zinb" | model == "negbin"){
      if (cont>1){ 
        parg <- 1/(sum((y-uj)^2/uj)/(N-nvar))
      }  
      while (abs(dpar)>0.0001 & cont1<200){
        if (parg<0){
          parg <- 0.00001
        }
        parg <- ifelse(parg<E^-10, E^-10, parg)
        gf <- sum(digamma(parg+y)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+y)/(parg+uj))
        hess <- sum(trigamma(parg+y)-trigamma(parg)+1/parg-2/(parg+uj)+(y+parg)/(parg+uj)^2)
        hess <- ifelse(hess==0, E^-23, hess)
        par0 <- parg
        #parg <- par0-solve(hess)%*%gf
        parg <- par0-as.vector(solve(hess))*gf
        if (parg>E^5){
          dpar <- 0.0001
          cont3 <- cont3+1
          if (cont3==1){
            parg <- 2
          } 
          else if (cont3==2) {
            parg <- E^5
          }
          else if (cont3==3){
            parg <- 0.0001  
          } 
        }
        else{
          dpar <- parg-par0
          cont1 <- cont1+1
        }
        if (parg>E^6){
          parg <- E^6
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
        # bg <- matrix(0, ncol(x),1)
        bg <- rep(0,ncol(x))
      } 
      else{
        bg <- solve(t(x)%*%(Ai*x))%*%t(x)%*%(Ai*zj)
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
      #lambdag <- matrix(0, ncol(G),1) 
      lambdag <- rep(0, ncol(G))
    }
    else {
      lambda0 <- log(lambda0/(1-lambda0))
      #lambdag <- rbind(lambda0, matrix(0,ncol(G)-1,1))
      lambdag <- rbind(lambda0, rep(0,ncol(G)-1))
    }
    pargg <- parg
    ujg <- uj
    if (is.null(nrow(pos0)) | any(lambdag==0)){ #essa condição sempre vai ser atendida? 
      if (is.null(nrow(pos0))) {
        pos0 <- pos1
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
        parg <- 1/alphag
      }
      else{
        if (j>0){
          parg <- 1/(sum((y-uj)^2/uj)/(N-nvar))
        }
        while (abs(dpar)>0.0001 & aux2<200){
          if (parg<0){
            parg <- 0.00001
          }
          parg <- ifelse(parg<E^-10, E^-10, parg)
          gf <- sum((1-zkg)*(digamma(parg+y)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+y)/(parg+uj)))
          hess <- sum((1-zkg)*(trigamma(parg+y)-trigamma(parg)+1/parg-2/(parg+uj)+(y+parg)/(parg+uj)^2))
          hess <- ifelse(hess==0, E^-23, hess)
          par0 <- parg
          parg <- par0-solve(hess)%*%gf #multiplicador
          if (aux2>50 & parg>E^5){
            dpar <- 0.0001
            cont3 <- cont3+1
            if (cont3==1){
              parg=2
            }
            else if (cont3==2){
              parg <- E^5
            }
            else if (cont3==3){
              parg <- 0.0001
            }
          }
          else{
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
      uj <- exp(nj)
      while (abs(ddev)>0.000001 & aux1<100){
        uj <- ifelse(uj>E^100, E^100, uj)
        Ai <- (1-zkg)*((uj/(1+alphag*uj)+(y-uj)*(alphag*uj/(1+2*alphag*uj+alphag**2*uj^2))))
        Ai <- ifelse(Ai<=0, E^-5, Ai)
        uj <- ifelse(uj<E^-150, E^-150, uj)
        zj <- (nj+(y-uj)/(((uj/(1+alphag*uj)+(y-uj)*(alphag*uj/(1+2*alphag*uj+alphag^2*uj^2))))*(1+alphag*uj)))-Offset
        if (det(t(x)%*%(Ai*x))==0){
          #bg <- matrix(0, nvar, 1)
          bg <- rep(0, nvar)
        }
        else{
          bg <- solve(t(x)%*%(Ai*x))%*%t(x)%*%(Ai*zj)
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
    
    
    
    %if %upcase(&model)=ZIP or %upcase(&model)=ZINB %then %do;
    devg=0;
    ddev=1;
    njl=G*lambdag;
    njl=choose(njl>&maxg,&maxg,njl);
    njl=choose(njl<-&maxg,-&maxg,njl);
    pig=exp(njl)/(1+exp(njl));
    do while (abs(ddev)>0.000001 & aux3<100);
    Ai=pig#(1-pig);
    Ai=choose(Ai<=0,1E-5,Ai);
    zj=njl+(zkg-pig)*1/Ai;
    if det((G#Ai)`*G)=0 then lambdag=j(ncol(G),1,0);
    else lambdag=inv((G#Ai)`*G)*(G#Ai)`*zj;
    njl=G*lambdag;
    njl=choose(njl>&maxg,&maxg,njl);
    njl=choose(njl<-&maxg,-&maxg,njl);
    pig=exp(njl)/(1+exp(njl));
    olddev=devg;
    devg=sum(zkg#njl-log(1+exp(njl)));
    ddev=devg-olddev;
    *print lambdag devg olddev ddev;
    aux3=aux3+1;
    end;
    %end;
    
    zkg=1/(1+exp(-njl)#(parg/(parg+uj))##parg);
    zkg=choose(y>0,0,zkg);
    %if %upcase(&model) ne ZIP and %upcase(&model) ne ZINB %then %do;
    zkg=0;
    %end;
    
    oldllike=llikeg;
    llikeg=sum(zkg#(njl)-log(1+exp(njl))+(1-zkg)#(log(gamma1)));
    dllike=llikeg-oldllike;
    *print j bg alphag lambdag llikeg dllike;
    j=j+1;
  }
} # fecha golden 



do while (abs(dllike)>0.00001 & j<600);

ddpar=1;
cont=1;
do while (abs(ddpar)>0.000001 & cont<100);
dpar=1;
parold=parg;
aux1=1;
aux2=1;
aux3=1;
cont3=0;
int=1;
%if %upcase(&model)=ZIP or %upcase(&model)=POISSON %then %do;
alphag=1E-6;parg=1/alphag;
%end;
%else %do;
if j>0 then parg=1/(sum((y-uj)##2/uj)/(n-nvar));
do while (abs(dpar)>0.0001 & aux2<200);
if parg<0 then parg=0.00001;
parg=choose(parg<1E-10,1E-10,parg);
gf=sum((1-zkg)#(digamma(parg+y)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+y)/(parg+uj)));
hess=sum((1-zkg)#(trigamma(parg+y)-trigamma(parg)+1/parg-2/(parg+uj)+(y+parg)/(parg+uj)##2));
hess=choose(hess=0,1E-23,hess);
par0=parg;
parg=par0-inv(hess)*gf;
if aux2>50 & parg>1E5 then do;
dpar=0.0001;
cont3=cont3+1;
if cont3=1 then parg=2;	
else if cont3=2 then parg=1E5;
else if cont3=3 then parg=0.0001; 
end;
else dpar=parg-par0;
if parg>1E6 then do;parg=1E6;dpar=0;end;
aux2=aux2+1;
end;
alphag=1/parg;
%end;

devg=0;
ddev=1;
nj=x*bg+offset;
uj=exp(nj);
do while (abs(ddev)>0.000001 & aux1<100);
uj=choose(uj>1E100,1E100,uj);
Ai=(1-zkg)#((uj/(1+alphag*uj)+(y-uj)#(alphag*uj/(1+2*alphag*uj+alphag**2*uj##2))));
Ai=choose(Ai<=0,1E-5,Ai);
uj=choose(uj<1E-150,1E-150,uj);
zj=(nj+(y-uj)/(((uj/(1+alphag*uj)+(y-uj)#(alphag*uj/(1+2*alphag*uj+alphag**2*uj##2))))#(1+alphag*uj)))-offset;
if det(x`*(Ai#x))=0 then bg=j(nvar,1,0);
else bg=inv(x`*(Ai#x))*x`*(Ai#zj);
nj=x*bg+offset;
nj=choose(nj>700,700,nj);
nj=choose(nj<-700,-700,nj);
uj=exp(nj);
olddev=devg;
gamma1=/*(gamma(par+y)/(gamma(y+1)#gamma(par)))#*/(uj/(uj+parg))##y#(parg/(uj+parg))##parg;
gamma1=choose(gamma1<=0,1E-10,gamma1);
devg=sum((1-zkg)#(log(gamma1)));
ddev=devg-olddev;
*print bg aux1 devg olddev ddev;
aux1=aux1+1;
end;
ddpar=parg-parold;
cont=cont+1;
end;

%if %upcase(&model)=ZIP or %upcase(&model)=ZINB %then %do;
devg=0;
ddev=1;
njl=G*lambdag;
njl=choose(njl>&maxg,&maxg,njl);
njl=choose(njl<-&maxg,-&maxg,njl);
pig=exp(njl)/(1+exp(njl));
do while (abs(ddev)>0.000001 & aux3<100);
Ai=pig#(1-pig);
Ai=choose(Ai<=0,1E-5,Ai);
zj=njl+(zkg-pig)*1/Ai;
if det((G#Ai)`*G)=0 then lambdag=j(ncol(G),1,0);
else lambdag=inv((G#Ai)`*G)*(G#Ai)`*zj;
njl=G*lambdag;
njl=choose(njl>&maxg,&maxg,njl);
njl=choose(njl<-&maxg,-&maxg,njl);
pig=exp(njl)/(1+exp(njl));
olddev=devg;
devg=sum(zkg#njl-log(1+exp(njl)));
ddev=devg-olddev;
*print lambdag devg olddev ddev;
aux3=aux3+1;
end;
%end;

zkg=1/(1+exp(-njl)#(parg/(parg+uj))##parg);
zkg=choose(y>0,0,zkg);
%if %upcase(&model) ne ZIP and %upcase(&model) ne ZINB %then %do;
zkg=0;
%end;

oldllike=llikeg;
llikeg=sum(zkg#(njl)-log(1+exp(njl))+(1-zkg)#(log(gamma1)));
dllike=llikeg-oldllike;
*print j bg alphag lambdag llikeg dllike;
j=j+1;
end;

/*****************************************/

read all var{&long &lat} into COORD;                                                                                                                      
_dist_ = distance(COORD, "L2");
seq=1:n;
/*create _dist_ from _dist_;append from _dist_;*/

start cv(h) global(n, wt, x, y, g, yhat, yhat2, pihat, hv, coord, _dist_, seq, offset, alphai, S, Si, parg, pargg, ujg, bg, lambdag, pos0, pos1, pos02, nvar);
%IF %UPCASE(&method)=ADAPTIVEN %THEN %DO;
hv=j(1,1,0);
yhat=j(1,1,0);
create &output from hv[colname='h'];
do i=1 to n;
%END;
%IF %UPCASE(&method)=FIXED_G or %UPCASE(&method)=FIXED_BSQ or %UPCASE(&method)=ADAPTIVE_BSQ %THEN %DO;
do i=1 to n;
%END;
do j=1 to n;
seqi=j(n,1,i);
dist=seqi||seq`||_dist_[,i];
%IF %UPCASE(&distancekm)=YES %THEN %DO;
dist[,3]=dist[,3]*111;
%END;
end;
u=nrow(dist);
w=j(u,1,0);
do jj=1 to u;
w[jj]=exp(-0.5*(dist[jj,3]/h)**2);
%IF %UPCASE(&method)=FIXED_BSQ or %UPCASE(&method)=ADAPTIVEN %THEN %DO;
w[jj]=(1-(dist[jj,3]/h)**2)**2;
%END;
%if %UPCASE(&bandwidth)=CV %THEN %DO;w[i]=0;%END;
end;
%IF %UPCASE(&method)=FIXED_BSQ or %UPCASE(&method)=ADAPTIVEN %THEN %DO;
position=loc(dist[,3]<=h);
w[position]=0;
%END;
%IF %UPCASE(&method)=ADAPTIVE_BSQ %THEN %DO;
call sort(dist,{3});
dist=dist||(1:nrow(dist))`;
w=j(n,2,0);	
hn=dist[h,3];
do jj=1 to n;
if dist[jj,4]<=h then
w[jj,1]=(1-(dist[jj,3]/hn)**2)**2;
else w[jj,1]=0;
w[jj,2]=dist[jj,2];
end;
%if %UPCASE(&bandwidth)=CV %THEN %DO;w[loc(w[,2]=i)]=0;%END;
call sort(w,{2});
w=w[,1];
%END;

b=bg;
nj=X*b+offset;
uj=exp(nj);
par=parg;
lambda=lambdag;
njl=G*lambda;
njl=choose(njl>&maxg,&maxg,njl);
njl=choose(njl<-&maxg,-&maxg,njl);

%if %upcase(&model) ne ZIP and %upcase(&model) ne ZINB %then %do;
zk=0;
%end;
%else %do;
lambda0=(ncol(pos0)-sum((parg/(uj+parg))##parg))/n;
if lambda0>0 then do;
lambda0=log(lambda0/(1-lambda0));
lambda=lambda0//j(ncol(G)-1,1,0);
njl=G*lambda;
end;
zk=1/(1+exp(-njl)#(par/(par+uj))##par);
zk=choose(y>0,0,zk);
%end;

dllike=1;
llike=0;
j=1;
do while (abs(dllike)>0.00001 & j<=600);

ddpar=1;
*do while (abs(ddpar)>0.000001);
dpar=1;
parold=par;
aux1=1;
aux2=1;
int=1;
%if %upcase(&model)=ZIP or %upcase(&model)=POISSON %then %do;
alpha=1E-6;par=1/alpha;
%end;
%else %do;
if par<=1E-5 then do;if i>1 then par=1/alphai[i-1,2];end;
if par>=1E6 then do;
par=1E6;dpar=0;
alpha=1/par;
b=bg;
uj=exp(X*b+offset);
lambda=lambdag;
njl=G*lambda;
njl=choose(njl>&maxg,&maxg,njl);
njl=choose(njl<-&maxg,-&maxg,njl);
zk=1/(1+exp(-njl)#(parg/(parg+uj))##parg);
zk=choose(y>0,0,zk);
if any(lambda)=0 then zk=0;
end;
do while (abs(dpar)>0.000001 & aux2<200);
par=choose(par<1E-10,1E-10,par);
gf=sum(w#wt#(1-zk)#(digamma(par+y)-digamma(par)+log(par)+1-log(par+uj)-(par+y)/(par+uj)));
hess=sum(w#wt#(1-zk)#(trigamma(par+y)-trigamma(par)+1/par-2/(par+uj)+(y+par)/(par+uj)##2));
hess=choose(hess=0,1E-23,hess);
par0=par;
par=par0-inv(hess)*gf;
dpar=par-par0;
if par>=1E6 then do;
par=1E6;dpar=0;
alpha=1/par;
b=bg;
uj=exp(X*b+offset);
lambda=lambdag;
njl=G*lambda;
njl=choose(njl>&maxg,&maxg,njl);
njl=choose(njl<-&maxg,-&maxg,njl);
zk=1/(1+exp(-njl)#(parg/(parg+uj))##parg);
zk=choose(y>0,0,zk);
if any(lambda)=0 then zk=0;
end;
aux2=aux2+1;
end;
if par<=1E-5 then do;
par=1E6;
b=bg;
uj=exp(X*b+offset);
lambda=lambdag;
njl=G*lambda;
njl=choose(njl>&maxg,&maxg,njl);
njl=choose(njl<-&maxg,-&maxg,njl);
zk=1/(1+exp(-njl)#(parg/(parg+uj))##parg);
zk=choose(y>0,0,zk);
if any(lambda)=0 then zk=0;
end;
alpha=1/par;
*print i j b lambda par alpha aux2;
%end;

dev=0;
ddev=1;
nj=x*b+offset;
nj=choose(nj>700,700,nj);
nj=choose(nj<-700,-700,nj);
uj=exp(nj);
do while (abs(ddev)>0.000001 & aux1<100);
uj=choose(uj>1E100,1E100,uj);
Ai=(1-zk)#((uj/(1+alpha*uj)+(y-uj)#(alpha*uj/(1+2*alpha*uj+alpha**2*uj##2))));
Ai=choose(Ai<=0,1E-5,Ai);
uj=choose(uj<1E-150,1E-150,uj);
denz=(((uj/(1+alpha*uj)+(y-uj)#(alpha*uj/(1+2*alpha*uj+alpha**2*uj##2))))#(1+alpha*uj));
denz=choose(denz=0,1E-5,denz);
zj=(nj+(y-uj)/denz)-offset;
if det(x`*(w#Ai#x#wt))=0 then b=j(nvar,1,0);
else b=inv(x`*(w#Ai#x#wt))*x`*(w#Ai#wt#zj);
nj=x*b+offset;
nj=choose(nj>700,700,nj);
nj=choose(nj<-700,-700,nj);
uj=exp(nj);
olddev=dev;
uj=choose(uj>1E10,1E10,uj);
uj=choose(uj=0,1E-10,uj);
if par=1E6 then gamma1=/*(gamma(par+y)/(gamma(y+1)#gamma(par)))#*/(uj/(uj+par))##y#exp(-uj);
else gamma1=/*(gamma(par+y)/(gamma(y+1)#gamma(par)))#*/(uj/(uj+par))##y#(par/(uj+par))##par;
gamma1=choose(gamma1<=0,1E-10,gamma1);
dev=sum((1-zk)#(log(gamma1)));
ddev=dev-olddev;
*print b par aux1 dev olddev ddev;
aux1=aux1+1;
end;
ddpar=par-parold;
*end;

%if %upcase(&model)=ZIP or %upcase(&model)=ZINB %then %do;

if j=1 then do;alphatemp=alpha;lambdatemp=lambda[1];end;
else do;alphatemp=alphatemp//alpha;lambdatemp=lambdatemp//lambda[1];end;
alphatemp=round(alphatemp,0.0000001);
lambdatemp=round(lambdatemp,0.0000001);
*print i j alphatemp (nrow(alphatemp)) (ncol(unique(alphatemp)));
%if %upcase(&model)=ZINB %then %do;
if j>300 & nrow(alphatemp)>ncol(unique(alphatemp)) & 
nrow(lambdatemp)>ncol(unique(lambdatemp))then do;
%end;
%if %upcase(&model)=ZIP %then %do;
*print i j lambdatemp (nrow(lambdatemp)) (ncol(unique(lambdatemp)));
if j>300 & nrow(lambdatemp)>ncol(unique(lambdatemp))then do;
%end;
lambda=j(ncol(G),1,0);
njl=G*lambda;
zk=j(n,1,0);
end;
else do;

aux3=1;
dev=0;
ddev=1;
njl=G*lambda;
njl=choose(njl>&maxg,&maxg,njl);
njl=choose(njl<-&maxg,-&maxg,njl);
pi=exp(njl)/(1+exp(njl));
do while (abs(ddev)>0.000001 & aux3<100);
Aii=pi#(1-pi);
Aii=choose(Aii<=0,1E-5,Aii);	
zj=njl+(zk-pi)/Aii;
if det((G#Aii#w#wt)`*G)=0 then lambda=j(ncol(G),1,0);
else lambda=inv((G#Aii#w#wt)`*G)*(G#Aii#w#wt)`*zj;
njl=G*lambda;
njl=choose(njl>&maxg,&maxg,njl);
njl=choose(njl<-&maxg,-&maxg,njl);
pi=exp(njl)/(1+exp(njl));
olddev=dev;
dev=sum(zk#njl-log(1+exp(njl)));
ddev=dev-olddev;
*print lambda aux3 dev olddev ddev;
aux3=aux3+1;
end;

end;
%end;

njl=G*lambda;
njl=choose(njl>&maxg,&maxg,njl);
njl=choose(njl<-&maxg,-&maxg,njl);
zk=1/(1+exp(-njl)#(par/(par+uj))##par);
zk=choose(y>0,0,zk);
if any(lambda)=0 then zk=j(n,1,0);

%if %upcase(&model) ne ZIP and %upcase(&model) ne ZINB %then %do;
zk=0;
%end;

oldllike=llike;
llike=sum(zk#(njl)-log(1+exp(njl))+(1-zk)#(log(gamma1)));
dllike=llike-oldllike;
*print i j b alpha lambda llike dllike;
j=j+1;
end;

%IF %UPCASE(&method)=FIXED_G or %UPCASE(&method)=FIXED_BSQ or %UPCASE(&method)=ADAPTIVE_BSQ %THEN %DO;
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
end;
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
  finish cv;
  
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
