/********************************************************************************/
/* Macro for searching the optimum Bandwidth */
/*
/* REQUIRED PARAMETERS
/*    DATA = the name of the SAS data set to be used
/*    YVAR = the name of the dependent or response variable
/*    XVAR = the name of the independent or explicative variables. A blank space
/*           should separate the names. Note: an intercept variable must not be
/*           created in advance
/*  DCOORD = the name of the SAS data set with the geographic coordinates
/*  OUTPUT = the name of the SAS data set to be used as output results
/*    MINV = the minimum distance between two locations i and k to be consider
/* MIDDLEV = the middle distance between two locations i and k to be consider
/*    MAXV = the maximum distance between two locations i and k to be consider
/*  METHOD = there are three choices:
/*           FIXED_G asks the program to compute the bandwidth as fixed gaussian;
/*           FIXED_BSQ to compute the bandwidth as fixed bi-square;
/*           ADAPTIVEN to compute the bandwidth as adaptive bi-square ($n$ values) and
/*			 ADAPTIVE_BSQ to compute the bandwidth as adaptive bi-square ($one$ value)
/*  DISTANCEKM = if the distances between two locations will be computed in Km
/*               using the Basic formulae for calculating spherical distance. The
/*               default value is NO, so the distance is computed using euclidian
/*               distance.
/********************************************************************************/
%macro Golden(DATA=, YVAR=, XVAR=, XVARINF=, WEIGHT=, LAT=, LONG=, OUTPUT=, 
		GLOBALMIN=YES, METHOD=, MODEL=ZINB, BANDWIDTH=CV, OFFSET=, FORCE=YES, 
		MAXG=100, DISTANCEKM=NO);
	proc iml;
		use &DATA;
		read all var {&YVAR} into y;
		read all var {&XVAR} into x;
		n=nrow(y);

		%if &XVARINF=%then
			%do;
				G=j(n, 1, 1);
				lambdag=j(ncol(G), 1, 0);
			%end;
		%else
			%do;
				read all var {&XVARINF} into G[colname=varnamezip];
				G=j(n, 1, 1)||G; 
			%end;
		wt=j(n, 1, 1);

		%IF &WEIGHT NE %THEN
			%DO;
				read all var {&WEIGHT} into wt;
			%END;
		offset=j(n, 1, 0);

		%IF &OFFSET NE %THEN
			%DO;
				read all var {&OFFSET} into offset;
			%END;
		x=j(n, 1, 1)||x;
		nvar=ncol(x);
		yhat=j(n, 1, 0);
		yhat2=j(n, 1, 0);
		pihat=j(n, 1, 0);
		alphai=j(n, 1, 0);
		S=j(n, 1, 0);
		Si=j(n, 1, 0);
		Iy=choose(y>0, 1, y);
		Iy=1-Iy;
		pos0=loc(y=0);
		pos02=loc(y=0);
		pos1=loc(y>0);

		/**** global estimates ****/
		uj=(y+y[:])/2;
		nj=log(uj);
		parg=sum((y-uj)##2/uj)/(n-nvar);
		ddpar=1;
		cont=1;
		cont3=0;

		do while (abs(ddpar)>0.000001 & cont<100);
			dpar=1;
			parold=parg;
			cont1=1;

			%IF %UPCASE(&MODEL)=ZIP or %UPCASE(&MODEL)=POISSON %THEN
				%DO;
					parg=1/1E-6;
					*if cont=1 then print parg;
					alphag=1/parg;
				%END;

			%IF %UPCASE(&MODEL)=ZINB or %UPCASE(&MODEL)=NEGBIN %THEN
				%DO;
					if cont>1 then
						parg=1/(sum((y-uj)##2/uj)/(n-nvar));
					do while (abs(dpar)>0.0001 & cont1<200);
						if parg<0 then
							parg=0.00001;
						parg=choose(parg<1E-10, 1E-10, parg);
						gf=sum(digamma(parg+y)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+y)/(parg+uj));
						hess=sum(trigamma(parg+y)-trigamma(parg)+1/parg-2/(parg+uj)+(y+parg)/(parg+uj)##2);
						hess=choose(hess=0, 1E-23, hess);
						par0=parg;
						parg=par0-inv(hess)*gf;
						if parg>1E5 then
							do;
								dpar=0.0001;
								cont3=cont3+1;
								if cont3=1 then
									parg=2;
								else if cont3=2 then
									parg=1E5;
								else if cont3=3 then
									parg=0.0001;
							end;
						else
							dpar=parg-par0;
						cont1=cont1+1;

						if parg>1E6 then
							do;
								parg=1E6;
								dpar=0;
							end;
					end;
					alphag=1/parg;
				%END;
			devg=0;
			ddev=1;
		    cont2=0;

			do while (abs(ddev)>0.000001 & cont2<100);
				*if cont=1 & cont2=0 then print alphag;
				Ai=(uj/(1+alphag*uj))+(y-uj)#(alphag*uj/(1+2*alphag*uj+alphag**2*uj#uj));
				Ai=choose(Ai<=0, 1E-5, Ai);
				zj=nj+(y-uj)/(Ai#(1+alphag*uj))-offset;

				if det(x`*(Ai#x))=0 then
					bg=j(ncol(x), 1, 0);
				else do;
					*if cont=1 & cont2=0 then print zj;
					bg=inv(x`*(Ai#x))*x`*(Ai#zj);
				end;
				*if cont=1 & cont2=0 then print bg;
				nj=x*bg+offset;
				nj=choose(nj>700, 700, nj);
				uj=exp(nj);
				olddev=devg;
				uj=choose(uj<1E-150, 1E-150, uj);
				uj=choose(uj>100000, 100000, uj);
				tt=y/uj;
				tt=choose(tt=0, 1E-10, tt);
				devg=2*sum(y#log(tt)-(y+1/alphag)#log((1+alphag*y)/(1+alphag*uj)));

				if cont2>100 then
					ddev=0.0000001;
				else
					ddev=devg-olddev;
				cont2=cont2+1;
			end;
			cont=cont+1;
			ddpar=parg-parold;
		end;
		%if &XVARINF ne %then
			%do;
				lambda0=(ncol(pos0)-sum((parg/(uj+parg))##parg))/n;

				if lambda0<=0 then
					lambdag=j(ncol(G), 1, 0);
				else
					do;
						lambda0=log(lambda0/(1-lambda0));
						lambdag=lambda0//j(ncol(G)-1, 1, 0);
					end;
				pargg=parg;
				ujg=uj;

				if nrow(pos0)=0 | any(lambdag)=0 then
					do;

						if nrow(pos0)=0 then
							do;
								pos0=pos1;

								%if %upcase(&MODEL)=ZINB or %upcase(&MODEL)=ZIP %then
									%do;
										call symputx("MODEL", "NEGBIN");
									%end;
							end;

						%if %upcase(&FORCE) ne YES %then
							%do;
								call symputx("MODEL", "NEGBIN");
							%end;
					end;
			%end;
		njl=G*lambdag;

		%if %upcase(&MODEL) ne ZIP and %upcase(&MODEL) ne ZINB %then
			%do;
				zkg=0;
			%end;
		%else
			%do;
				zkg=1/(1+exp(-G*lambdag)#(parg/(parg+uj))##parg);
				zkg=choose(y>0, 0, zkg);
			%end;
		dllike=1;
		llikeg=0;
		j=0;
		contador=0;
		do while (abs(dllike)>0.00001 & j<600);
			contador=contador+1;
			ddpar=1;
			cont=1;
			contador2=0;
			do while (abs(ddpar)>0.000001 & cont<100);
				contador2=contador2+1;
				dpar=1;
				parold=parg;
				aux1=1;
				aux2=1;
				aux3=1;
				cont3=0;
				int=1;

				%if %upcase(&MODEL)=ZIP or %upcase(&MODEL)=POISSON %then
					%do;
						alphag=1E-6;
						parg=1/alphag;
					%end;
				%else
					%do;

						if j>0 then
							parg=1/(sum((y-uj)##2/uj)/(n-nvar));

						do while (abs(dpar)>0.0001 & aux2<200);

							if parg<0 then
								parg=0.00001;
							parg=choose(parg<1E-10, 1E-10, parg);
							gf=sum((1-zkg)#(digamma(parg+y)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+y)/(parg+uj)));
							hess=sum((1-zkg)#(trigamma(parg+y)-trigamma(parg)+1/parg-2/(parg+uj)+(y+parg)/(parg+uj)##2));
							hess=choose(hess=0, 1E-23, hess);
							par0=parg;
							parg=par0-inv(hess)*gf;

							if aux2>50 & parg>1E5 then
								do;
									dpar=0.0001;
									cont3=cont3+1;

									if cont3=1 then
										parg=2;
									else if cont3=2 then
										parg=1E5;
									else if cont3=3 then
										parg=0.0001;
								end;
							else
								dpar=parg-par0;

							if parg>1E6 then
								do;
									parg=1E6;
									dpar=0;
								end;
							aux2=aux2+1;
						end;
						alphag=1/parg;
					%end;
				devg=0;
				ddev=1;
				nj=x*bg+offset;
				uj=exp(nj);

				do while (abs(ddev)>0.000001 & aux1<100);
					uj=choose(uj>1E100, 1E100, uj);
					Ai=(1-zkg)#((uj/(1+alphag*uj)+(y-uj)#(alphag*uj/(1+2*alphag*uj+alphag**2*uj##2))));
					Ai=choose(Ai<=0, 1E-5, Ai);
					uj=choose(uj<1E-150, 1E-150, uj);
					zj=(nj+(y-uj)/(((uj/(1+alphag*uj)+(y-uj)#(alphag*uj/(1+2*alphag*uj+alphag**2*uj##2))))#(1+alphag*uj)))-offset;
					if det(x`*(Ai#x))=0 then
						bg=j(nvar, 1, 0);
					else
						bg=inv(x`*(Ai#x))*x`*(Ai#zj);
					nj=x*bg+offset;
					nj=choose(nj>700, 700, nj);
					nj=choose(nj<-700, -700, nj);
					uj=exp(nj);
					olddev=devg;
					gamma1=/*(gamma(par+y)/(gamma(y+1)#gamma(par)))#*/(uj/(uj+parg))##y#(parg/(uj+parg))##parg;
					gamma1=choose(gamma1<=0, 1E-10, gamma1);
					devg=sum((1-zkg)#(log(gamma1)));
					ddev=devg-olddev;
					*print bg aux1 devg olddev ddev;
					*print comentado;
					aux1=aux1+1;
				end;
				ddpar=parg-parold;
				cont=cont+1;
			end;

			%if %upcase(&MODEL)=ZIP or %upcase(&MODEL)=ZINB %then
				%do;
					devg=0;
					ddev=1;
					njl=G*lambdag;
					njl=choose(njl>&MAXG, &MAXG, njl);
					njl=choose(njl<-&MAXG, -&MAXG, njl);
					pig=exp(njl)/(1+exp(njl));
					contador3 = 0;
					do while (abs(ddev)>0.000001 & aux3<100);
					contador3 = contador3 + 1;
						Ai=pig#(1-pig);
						Ai=choose(Ai<=0, 1E-5, Ai);
						zj=njl+(zkg-pig)*1/Ai;
						if det((G#Ai)`*G)=0 then do;
							lambdag=j(ncol(G), 1, 0);
						end;
						else do;
							lambdag=inv((G#Ai)`*G)*(G#Ai)`*zj;
						end;
						njl=G*lambdag;
						njl=choose(njl>&MAXG, &MAXG, njl);
						njl=choose(njl<-&MAXG, -&MAXG, njl);
						pig=exp(njl)/(1+exp(njl));
						olddev=devg;
						devg=sum(zkg#njl-log(1+exp(njl)));
						ddev=devg-olddev;
						*print lambdag devg olddev ddev;
						*print comentado;
						aux3=aux3+1;
					end;
				%end;
			zkg=1/(1+exp(-njl)#(parg/(parg+uj))##parg);
			zkg=choose(y>0, 0, zkg);

			%if %upcase(&MODEL) ne ZIP and %upcase(&MODEL) ne ZINB %then
				%do;
					zkg=0;
				%end;
			oldllike=llikeg;
			llikeg=sum(zkg#(njl)-log(1+exp(njl))+(1-zkg)#(log(gamma1)));
			dllike=llikeg-oldllike;
			*print j bg alphag lambdag llikeg dllike;
			*print comentado;
			j=j+1;
		end;

		/*****************************************/
		read all var{&LONG &LAT} into COORD;
		close &DATA;
		_dist_=distance(COORD, "L2");
		seq=1:n;

		/*create _dist_ from _dist_;append from _dist_;*/
		start cv(h) global(n, wt, x, y, g, yhat, yhat2, pihat, hv, coord, _dist_, 
			seq, offset, alphai, S, Si, parg, pargg, ujg, bg, lambdag, pos0, pos1, 
			pos02, nvar);

		%IF %UPCASE(&METHOD)=ADAPTIVEN %THEN
			%DO;
				hv=j(1, 1, 0);
				yhat=j(1, 1, 0);
				create &OUTPUT from hv[colname='h'];

				do i=1 to n;
				%END;

			%IF %UPCASE(&METHOD)=FIXED_G or %UPCASE(&METHOD)=FIXED_BSQ or 
				%UPCASE(&METHOD)=ADAPTIVE_BSQ %THEN
					%DO;

					do i=1 to n;
					%END;

				do j=1 to n;
					seqi=j(n, 1, i);
					dist=seqi||seq`||_dist_[, i];

					%IF %UPCASE(&DISTANCEKM)=YES %THEN
						%DO;
							dist[, 3]=dist[, 3]*111;
						%END;
				end; *fecha loop j;
				u=nrow(dist);
				w=j(u, 1, 0);

				do jj=1 to u;
					w[jj]=exp(-0.5*(dist[jj, 3]/h)**2);

					%IF %UPCASE(&METHOD)=FIXED_BSQ or %UPCASE(&METHOD)=ADAPTIVEN %THEN
						%DO;
							w[jj]=(1-(dist[jj, 3]/h)**2)**2;
						%END;

					%if %UPCASE(&BANDWIDTH)=CV %THEN
						%DO;
							w[i]=0;
						%END;
				end; *fecha loop jj;

				%IF %UPCASE(&METHOD)=FIXED_BSQ or %UPCASE(&METHOD)=ADAPTIVEN %THEN
					%DO;
						position=loc(dist[, 3]<=h);
						w[position]=0;
					%END;

				%IF %UPCASE(&METHOD)=ADAPTIVE_BSQ %THEN
					%DO;
						call sort(dist, {3});
						dist=dist||(1:nrow(dist))`;
						w=j(n, 2, 0);
						hn=dist[h, 3];

						do jj=1 to n;

							if dist[jj, 4]<=h then
								w[jj, 1]=(1-(dist[jj, 3]/hn)**2)**2;
							else
								w[jj, 1]=0;
							w[jj, 2]=dist[jj, 2];
						end; *fecha loop jj;

						%if %UPCASE(&BANDWIDTH)=CV %THEN
							%DO;
								w[loc(w[, 2]=i)]=0;
							%END;
						call sort(w, {2});
						w=w[, 1];
					%END;
				*if i=1 then print w;
				*if i=1 then print bg;
				b=bg;
				nj=X*b+offset;
				uj=exp(nj);
				par=parg;
				lambda=lambdag;
				njl=G*lambda;
				njl=choose(njl>&MAXG, &MAXG, njl);
				njl=choose(njl<-&MAXG, -&MAXG, njl);

				%if %upcase(&MODEL) ne ZIP and %upcase(&MODEL) ne ZINB %then
					%do;
						zk=0;
					%end;
				%else
					%do;
						lambda0=(ncol(pos0)-sum((parg/(uj+parg))##parg))/n;

						if lambda0>0 then
							do;
								lambda0=log(lambda0/(1-lambda0));
								lambda=lambda0//j(ncol(G)-1, 1, 0);
								njl=G*lambda;
							end;
						zk=1/(1+exp(-njl)#(par/(par+uj))##par);
						zk=choose(y>0, 0, zk);
					%end;
				dllike=1;
				llike=0;
				j=1;
				contador4=0;
				do while (abs(dllike)>0.00001 & j<=600);
					contador4=contador4+1;
					ddpar=1;
					*do while (abs(ddpar)>0.000001);
					dpar=1;
					parold=par;
					aux1=1;
					aux2=1;
					int=1;

					%if %upcase(&MODEL)=ZIP or %upcase(&MODEL)=POISSON %then
						%do;
							alpha=1E-6;
							par=1/alpha;
						%end;
					%else
						%do;

							if par<=1E-5 then
								do;

									if i>1 then
										par=1/alphai[i-1, 2];
								end;

							if par>=1E6 then
								do;
									par=1E6;
									dpar=0;
									alpha=1/par;
									b=bg;
									uj=exp(X*b+offset);
									lambda=lambdag;
									njl=G*lambda;
									njl=choose(njl>&MAXG, &MAXG, njl);
									njl=choose(njl<-&MAXG, -&MAXG, njl);
									zk=1/(1+exp(-njl)#(parg/(parg+uj))##parg);
									zk=choose(y>0, 0, zk);

									if any(lambda)=0 then
										zk=0;
								end;

							do while (abs(dpar)>0.000001 & aux2<200);
								par=choose(par<1E-10, 1E-10, par);
								gf=sum(w#wt#(1-zk)#(digamma(par+y)-digamma(par)+log(par)+1-log(par+uj)-(par+y)/(par+uj)));
								hess=sum(w#wt#(1-zk)#(trigamma(par+y)-trigamma(par)+1/par-2/(par+uj)+(y+par)/(par+uj)##2));
								hess=choose(hess=0, 1E-23, hess);
								par0=par;
								par=par0-inv(hess)*gf;
								dpar=par-par0;

								if par>=1E6 then
									do;
										par=1E6;
										dpar=0;
										alpha=1/par;
										b=bg;
										uj=exp(X*b+offset);
										lambda=lambdag;
										njl=G*lambda;
										njl=choose(njl>&MAXG, &MAXG, njl);
										njl=choose(njl<-&MAXG, -&MAXG, njl);
										zk=1/(1+exp(-njl)#(parg/(parg+uj))##parg);
										zk=choose(y>0, 0, zk);

										if any(lambda)=0 then
											zk=0;
									end;
								aux2=aux2+1;
							end;

							if par<=1E-5 then
								do;
									par=1E6;
									b=bg;
									uj=exp(X*b+offset);
									lambda=lambdag;
									njl=G*lambda;
									njl=choose(njl>&MAXG, &MAXG, njl);
									njl=choose(njl<-&MAXG, -&MAXG, njl);
									zk=1/(1+exp(-njl)#(parg/(parg+uj))##parg);
									zk=choose(y>0, 0, zk);

									if any(lambda)=0 then
										zk=0;
								end;
							alpha=1/par;
							*if i=244 then print i j b lambda par alpha aux2;
							*print comentado;
						%end; *fecha pct else;
					dev=0;
					ddev=1;
					*if i=1 & contador4=1 then print b;
					nj=x*b+offset;
					nj=choose(nj>700, 700, nj);
					nj=choose(nj<-700, -700, nj);
					*if i=1 & contador4=1 then print nj;
					uj=exp(nj);

					contador5=0;
					*if i=21 & contador4=1 then print (sum(uj));
					do while (abs(ddev)>0.000001 & aux1<100);
						contador5=contador5+1;
						*if i=21 & contador4=1 & contador5=4 then print (sum(uj));
						uj=choose(uj>1E100, 1E100, uj);
						*if i=21 & contador4=1 & contador5=4 then print (sum(zk)) (sum(uj)) (sum(alpha));
						Ai=(1-zk)#((uj/(1+alpha*uj)+(y-uj)#(alpha*uj/(1+2*alpha*uj+alpha**2*uj##2))));
						*if i=21 & contador4=1 & contador5=4 then print (sum(Ai));
						Ai=choose(Ai<=0, 1E-5, Ai);
						uj=choose(uj<1E-150, 1E-150, uj);
						denz=(((uj/(1+alpha*uj)+(y-uj)#(alpha*uj/(1+2*alpha*uj+alpha**2*uj##2))))#(1+alpha*uj));
						denz=choose(denz=0, 1E-5, denz);
						zj=(nj+(y-uj)/denz)-offset;
						*if i=21 & contador4=1 & contador5=4 then print (sum(x`*(w#Ai#x#wt)));
						if det(x`*(w#Ai#x#wt))=0 then
							b=j(nvar, 1, 0);
						else do;
							b=inv(x`*(w#Ai#x#wt))*x`*(w#Ai#wt#zj);
						end;
						*if i=21 & contador5=1 & contador4=4 then print b;
						nj=x*b+offset;
						*if i=128 & contador5=1 & contador4=1 then print (sum(nj));
						nj=choose(nj>700, 700, nj);
						nj=choose(nj<-700, -700, nj);
						*if i=152 & contador4=1 & contador5=1 then print (sum(nj));
						uj=exp(nj);
						*if i=128 & contador5=1 & contador4=1 then print (sum(uj));
						*if i=128 & contador4=1 then print (sum(uj));
						olddev=dev;
						*if i=21 & contador4=1 & contador5=1 then print (sum(uj));
						uj=choose(uj>1E10, 1E10, uj);
						uj=choose(uj=0, 1E-10, uj);
						*if i=152 & contador4=1 & contador5=1 then print (sum(uj));
						if par=1E6 then do;
							gamma1=/*(gamma(par+y)/(gamma(y+1)#gamma(par)))#*/(uj/(uj+par))##y#exp(-uj);
							*if i=128 & contador4=1 & contador5=6 then print (exp(-uj));
						end;
						else
							gamma1=/*(gamma(par+y)/(gamma(y+1)#gamma(par)))#*/(uj/(uj+par))##y#(par/(uj+par))##par;
						*if i=128 & contador4=1 & contador5=6 then print gamma1;
						gamma1=choose(gamma1<=0, 1E-10, gamma1);
						*if i=152 & contador4=1 & contador5=5 then print dev gamma1;
						dev=sum((1-zk)#(log(gamma1)));
						*if i=128 & contador4=1 & contador5=6 then print olddev dev;
						ddev=dev-olddev;
						*if i=244 then print b par aux1 dev olddev ddev;
						*print comentado;
						aux1=aux1+1;
						*if i=128 then print ("loop interno") contador5;
						*if i=128 & contador4=1 & contador5=6 then print ddev aux1;
					end; *fecha loop while;
					ddpar=par-parold;
					*end;

					%if %upcase(&MODEL)=ZIP or %upcase(&MODEL)=ZINB %then
						%do;

							if j=1 then
								do;
									alphatemp=alpha;
									lambdatemp=lambda[1];
								end;
							else
								do;
									alphatemp=alphatemp//alpha;
									lambdatemp=lambdatemp//lambda[1];
								end;
							alphatemp=round(alphatemp, 0.0000001);
							lambdatemp=round(lambdatemp, 0.0000001);
							*if i=244 then print i j alphatemp (nrow(alphatemp)) (ncol(unique(alphatemp)));
							*print comentado;

							%if %upcase(&MODEL)=ZINB %then
								%do;

									if j>300 & nrow(alphatemp)>ncol(unique(alphatemp)) & 
nrow(lambdatemp)>ncol(unique(lambdatemp))then
										do;
										%end;

									%if %upcase(&MODEL)=ZIP %then
										%do;
											*print i j lambdatemp (nrow(lambdatemp)) (ncol(unique(lambdatemp)));
											*print comentado;

											if j>300 & nrow(lambdatemp)>ncol(unique(lambdatemp))then
												do;
												%end;
											lambda=j(ncol(G), 1, 0);
											njl=G*lambda;
											zk=j(n, 1, 0);
										end;
									else
										do;
											aux3=1;
											dev=0;
											ddev=1;
											njl=G*lambda;
											njl=choose(njl>&MAXG, &MAXG, njl);
											njl=choose(njl<-&MAXG, -&MAXG, njl);
											pi=exp(njl)/(1+exp(njl));

											contador6=0;
											do while (abs(ddev)>0.000001 & aux3<100);
												contador6=contador6+1;
												Aii=pi#(1-pi);
												Aii=choose(Aii<=0, 1E-5, Aii);
												zj=njl+(zk-pi)/Aii;

												if det((G#Aii#w#wt)`*G)=0 then do;
													lambda=j(ncol(G), 1, 0);
												end;
												else do;
													lambda=inv((G#Aii#w#wt)`*G)*(G#Aii#w#wt)`*zj;
												end;
												njl=G*lambda;
												njl=choose(njl>&MAXG, &MAXG, njl);
												njl=choose(njl<-&MAXG, -&MAXG, njl);
												pi=exp(njl)/(1+exp(njl));
												olddev=dev;
												dev=sum(zk#njl-log(1+exp(njl)));
												ddev=dev-olddev;
												*if i=244 then print lambda aux3 dev olddev ddev;
												*print comentado;
												aux3=aux3+1;
											end;
										end; *fecha else;
								%end;
							njl=G*lambda;
							njl=choose(njl>&MAXG, &MAXG, njl);
							njl=choose(njl<-&MAXG, -&MAXG, njl);
							zk=1/(1+exp(-njl)#(par/(par+uj))##par);
							zk=choose(y>0, 0, zk);

							if any(lambda)=0 then
								zk=j(n, 1, 0);

							%if %upcase(&MODEL) ne ZIP and %upcase(&MODEL) ne ZINB %then
								%do;
									zk=0;
								%end;
							oldllike=llike;
							llike=sum(zk#(njl)-log(1+exp(njl))+(1-zk)#(log(gamma1)));
							dllike=llike-oldllike;
							*if i=244 then print i j b alpha lambda llike dllike;
							*print comentado;
							j=j+1;
							*if i=128 then print ("loop externo") contador4;
						end; *fecha loop while (condicoes j e dlike);

					%IF %UPCASE(&METHOD)=FIXED_G or %UPCASE(&METHOD)=FIXED_BSQ or 
						%UPCASE(&METHOD)=ADAPTIVE_BSQ %THEN
							%DO;
							*if i=21 then print (uj[i]);
							yhat[i]=uj[i];
							pihat[i]=njl[i];
							alphai[i]=alpha;

							if det(x`*(w#Ai#x#wt))=0 then do;
								S[i]=0;
							end;
							else do;
								S[i]=(x[i, ]*inv(x`*(w#Ai#x#wt))*(x#w#Ai#wt)`)[i];
							end;

							%if %upcase(&MODEL)=ZIP or %upcase(&MODEL)=ZINB %then
								%do;
									yhat[i]=(uj#(1-exp(njl)/(1+exp(njl))))[i];
									yhat2[i]=uj[i];

									if det(G`*(w#Aii#G#wt))=0 then do;
										Si[i]=0;
									end;
									else do;
										Si[i]=(G[i, ]*inv(G`*(w#Aii#G#wt))*(G#w#Aii#wt)`)[i];
									end;

									if any(lambda)=0 then do;
										Si[i]=0;
									end;
								%end;
						end; *fecha o for loop i;
						*print yhat;
						*print wt;
						CV=((y-yhat)#wt)`*(y-yhat);
						_par_=1/alphai;

						%IF %UPCASE(&MODEL)=ZINB or %UPCASE(&MODEL)=ZIP %THEN
							%DO;
								npar=sum(S)+sum(Si);

								%if %UPCASE(&BANDWIDTH)=AIC %THEN
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
											end;
										dev=2*(llnull1-ll);
										npar=sum(S)+sum(Si);
										AIC=2*npar-2*ll;
										AICc=AIC+2*(npar*(npar+1)/(n-npar-1));

										%IF %UPCASE(&MODEL)=ZINB %THEN
											%DO;
												AIC=2*(npar+npar/(ncol(x)+ncol(G)))-2*ll;
												AICc=AIC+2*((npar+npar/(ncol(x)+ncol(G)))*((npar+npar/(ncol(x)+ncol(G)))+1)/(n-(npar+npar/(ncol(x)+ncol(G)))-1));
												*print AIC AICc npar (ncol(G));
											%END;
								%END;
							%END;

						%IF %UPCASE(&MODEL)=POISSON or %UPCASE(&MODEL)=NEGBIN %THEN
							%DO;
								npar=sum(S);
								
								%if %UPCASE(&BANDWIDTH)=AIC %THEN
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
												pos0x=choose(pos0x=0,1E-10,pos0x);
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
										dev=2*(llnull1-ll);
										npar=sum(S);
										AIC=2*npar-2*ll;
										AICc=AIC+2*(npar*(npar+1)/(n-npar-1));

										%IF %UPCASE(&MODEL)=NEGBIN %THEN
											%DO;
												AIC=2*(npar+npar/ncol(x))-2*ll;
												AICc=AIC+2*(npar+npar/ncol(x))*(npar+npar/ncol(x)+1)/(n-(npar+npar/ncol(x))-1);
											%END;
									%END;
							%END;

						%if %UPCASE(&BANDWIDTH)=AIC %THEN
							%DO;
								CV=AICC;
							%END;
					%END;
				%ELSE %IF %UPCASE(&METHOD)=ADAPTIVEN %THEN
					%DO;
						yhat[1]=x[i, ]*b;
						CV=((y[i]-yhat)#wt)`*(y[i]-yhat);
					%END;
				free dist w;
				*print cv;
				res=cv||npar;
				return (res);
				finish cv;

				/**** DEFINING GOLDEN SECTION SEARCH PARAMETERS *****/

				%IF %UPCASE(&METHOD)=FIXED_G or %UPCASE(&METHOD)=FIXED_BSQ %THEN
					%DO;
						ax=0;
						bx=int(max(_dist_)+1);

						%IF %UPCASE(&DISTANCEKM)=YES %THEN
							%DO;
								bx=bx*111;
							%END;
					%END;

				%IF %UPCASE(&METHOD)=ADAPTIVE_BSQ %THEN
					%DO;
						ax=5;
						bx=n;
					%END;
				r=0.61803399;
				tol=0.1;

				%IF %UPCASE(&GLOBALMIN)=NO %THEN
					%DO;
						lower=ax;
						upper=bx;
						xmin=j(1, 2, 0);

						do GMY=1 to 1;
							ax1=lower[GMY];
							bx1=upper[GMY];
						%END;
					%ELSE
						%DO;
							lower=ax||(1-r)*bx||r*bx;
							upper=(1-r)*bx||r*bx||bx;
							xmin=j(3, 2, 0);

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
						*print ("chama cv");
						res1=cv(h1);
						CV1=res1[1];
						*print ("chama cv");
						res2=cv(h2);
						CV2=res2[1];

						%IF %UPCASE(&METHOD)=FIXED_G or %UPCASE(&METHOD)=FIXED_BSQ or 
							%UPCASE(&METHOD)=ADAPTIVE_BSQ %THEN
								%DO;

								if GMY=1 then
									create &OUTPUT var{GMY h1 cv1 h2 cv2};
									%IF %UPCASE(&METHOD)=FIXED_G or %UPCASE(&METHOD)=FIXED_BSQ or 
								%UPCASE(&METHOD)=ADAPTIVE_BSQ %THEN
									%DO;
									append;
								%END;
							%END;
						int=1;

						do while(abs(h3-h0) > tol*(abs(h1)+abs(h2)) & int<200);
							if CV2<CV1 then
								do;
									h0=h1;
									h1=h3-r*(h3-h0);
									h2=h0+r*(h3-h0);
									CV1=CV2;
									*print ("chama cv");
									res2=cv(h2);
									CV2=res2[1];
								end;
							else
								do;
									h3=h2;
									h1=h3-r*(h3-h0);
									h2=h0+r*(h3-h0);
									CV2=CV1;
									*print ("chama cv");
									res1=cv(h1);
									CV1=res1[1];
								end;

							%IF %UPCASE(&METHOD)=FIXED_G or %UPCASE(&METHOD)=FIXED_BSQ or 
								%UPCASE(&METHOD)=ADAPTIVE_BSQ %THEN
									%DO;
									append;
								%END;
							int=int+1;
						end;

						if CV1<CV2 then
							do;
								golden=CV1;
								xmin[GMY, 1]=golden;
								xmin[GMY, 2]=h1;
								npar=res1[2];

								%IF %UPCASE(&METHOD)=ADAPTIVE_BSQ %THEN
									%DO;
										xmin[GMY, 2]=floor(h1);
									%end;
							end;
						else
							do;
								golden=CV2;
								xmin[GMY, 1]=golden;
								xmin[GMY, 2]=h2;
								npar=res2[2];

								%IF %UPCASE(&METHOD)=ADAPTIVE_BSQ %THEN
									%DO;
										xmin[GMY, 2]=floor(h2);
									%end;
							end;

						%IF %UPCASE(&METHOD)=FIXED_G or %UPCASE(&METHOD)=FIXED_BSQ or 
							%UPCASE(&METHOD)=ADAPTIVE_BSQ %THEN
								%DO;
								print golden (xmin[GMY, 2])[label='xmin'] 
								%if %UPCASE(&BANDWIDTH)=AIC %THEN
									%DO;
										npar%END;
								;
							%END;
						%ELSE %IF %UPCASE(&METHOD)=ADAPTIVEN %THEN
							%DO;
								hv[1]=xmin[GMY, 2];
								append from hv;
							end;
						%END;
				end; *fecha os loops do GMY;
				create _min_bandwidth_ from xmin[colname={'golden' 'bandwidth'}];
				append from xmin;

				%IF %UPCASE(&GLOBALMIN)=YES %THEN
					%DO;
						print xmin[colname={'golden' 'bandwidth'}];
						xming=xmin[loc(xmin[, 1]=min(xmin[, 1])), 2];
						print 'Global Minimum', xming[label='(Da Silva and Mendes, 2018)'];
					%END;
				quit;

				%IF %UPCASE(&METHOD)=FIXED_G or %UPCASE(&METHOD)=FIXED_BSQ or 
					%UPCASE(&METHOD)=ADAPTIVE_BSQ %THEN
						%DO;
						%global _h_;

					proc sql noprint;
						select bandwidth into:_h_ from _min_bandwidth_ having golden=min(golden);
					quit;

					%put h=&_h_;
				%END;
		%mend Golden;

		/*******************************************************************************/
		/* Macro for estimating GWZINBR Model */
		/* REQUIRED PARAMETERS
		/*    DATA = the name of the SAS data set to be used
		/*    YVAR = the name of the dependent or response variable
		/*    XVAR = the name of the independent or explicative variables. A blank space
		/*           should separate the names. Note: an intercept variable must not be
		/*           created in advance
		/*   WEIGHT = the name of the sample weight variable
		/*  DCOORD = the name of the SAS data set with the geographic coordinates
		/*    GRID = the name of the SAS data set with the grid of geographic coordinates
		/*           the standard errors of complex data
		/*     DHV = the name of the SAS data set with the bandwidth adaptive ($n$ values),
		/*           which must have an unique variable
		/*       H = A pre-defined bandwidth value for METHOD equal to FIXED or ADAPTIVE1
		/*    MAXV = the maximum distance between two locations i and k to be consider
		/*  METHOD = there are three choices:
		/*           FIXED_G asks the program to compute the bandwidth as fixed gaussian;
		/*           FIXED_BSQ to compute the bandwidth as fixed bi-square;
		/*           ADAPTIVEN to compute the bandwidth as adaptive bi-square ($n$ values) and
		/*			 ADAPTIVE_BSQ to compute the bandwidth as adaptive bi-square ($one$ value)
		/*  DISTANCEKM = if the distances between two locations will be computed in Km
		/*               using the Basic formulae for calculating spherical distance. The
		/*               default value is NO, so the distance is computed using euclidian
		/*               distance.
		/********************************************************************************/
		%macro GWZINBR(DATA=, YVAR=, XVAR=, XVARINF=, WEIGHT=, LAT=, LONG=, GRID=, 
				DHV=, METHOD=, MODEL=ZINB, OFFSET=, DISTANCEKM=NO, FORCE=YES, INT_INF=YES, 
				MAXG=100, H=);
			proc iml;
				use &DATA;
				read all var {&YVAR} into y;
				read all var {&XVAR} into x;
				n=nrow(y);

				%if &XVARINF=%then
					%do;
						G=j(n, 1, 1);
						lambdag=j(ncol(G), 1, 0);
					%end;
				%else
					%do;
						read all var {&XVARINF} into G[colname=varnamezip];

						%IF %UPCASE(&INT_INF)=YES %THEN
							%DO;
								G=j(n, 1, 1)||G;
							%END;
					%end;
				x=j(n, 1, 1)||x;
				yhat=j(n, 1, 0);
				yhat2=j(n, 1, 0);
				pihat=j(n, 1, 0);
				nvar=ncol(x);
				wt=j(n, 1, 1);

				%IF &WEIGHT NE %THEN
					%DO;
						read all var {&WEIGHT} into wt;
					%END;
				offset=j(n, 1, 0);

				%IF &OFFSET NE %THEN
					%DO;
						read all var {&OFFSET} into offset;
					%END;
				Iy=choose(y>0, 1, y);
				Iy=1-Iy;
				Iy2=Iy;
				pos0=loc(y=0);
				pos02=loc(y=0);
				pos1=loc(y>0);

				/**** global estimates ****/
				uj=(y+y[:])/2;
				nj=log(uj);
				parg=sum((y-uj)##2/uj)/(n-nvar);
				ddpar=1;
				cont=1;
				cont3=0;

				do while (abs(ddpar)>0.000001 & cont<100);
					dpar=1;
					parold=parg;
					cont1=1;

					%IF %UPCASE(&MODEL)=ZIP or %UPCASE(&MODEL)=POISSON %THEN
						%DO;
							parg=1/1E-6;
							alphag=1/parg;
						%END;

					%IF %UPCASE(&MODEL)=ZINB or %UPCASE(&MODEL)=NEGBIN %THEN
						%DO;

							if cont>1 then
								parg=1/(sum((y-uj)##2/uj)/(n-nvar));

							do while (abs(dpar)>0.0001 & cont1<200);

								if parg<0 then
									parg=0.00001;
								parg=choose(parg<1E-10, 1E-10, parg);
								gf=sum(digamma(parg+y)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+y)/(parg+uj));
								hessg=sum(trigamma(parg+y)-trigamma(parg)+1/parg-2/(parg+uj)+(y+parg)/(parg+uj)##2);
								hessg=choose(hessg=0, 1E-23, hessg);
								par0=parg;
								parg=par0-inv(hessg)*gf;

								if parg>1E5 then
									do;
										dpar=0.0001;
										cont3=cont3+1;

										if cont3=1 then
											parg=2;
										else if cont3=2 then
											parg=1E5;
										else if cont3=3 then
											parg=0.0001;
									end;
								else
									dpar=parg-par0;
								cont1=cont1+1;

								if parg>1E6 then
									do;
										parg=1E6;
										dpar=0;
									end;
							end;
							alphag=1/parg;
						%END;
						*print alphag;
					devg=0;
					ddev=1;
					cont2=0;

					do while (abs(ddev)>0.000001 & cont2<200);
						Ai=(uj/(1+alphag*uj))+(y-uj)#(alphag*uj/(1+2*alphag*uj+alphag**2*uj#uj));
						;
						Ai=choose(Ai<=0, 1E-5, Ai);
						zj=nj+(y-uj)/(Ai#(1+alphag*uj))-offset;

						if det(x`*(Ai#x))=0 then
							bg=j(ncol(x), 1, 0);
						else
							bg=inv(x`*(Ai#x))*x`*(Ai#zj);
						nj=x*bg+offset;
						nj=choose(nj>700, 700, nj);
						uj=exp(nj);
						olddev=devg;
						uj=choose(uj<1E-150, 1E-150, uj);
						uj=choose(uj>100000, 100000, uj);
						tt=y/uj;
						tt=choose(tt=0, 1E-10, tt);
						devg=2*sum(y#log(tt)-(y+1/alphag)#log((1+alphag*y)/(1+alphag*uj)));

						if cont2>100 then
							ddev=0.0000001;
						else
							ddev=devg-olddev;
						cont2=cont2+1;
					end;
					cont=cont+1;
					ddpar=parg-parold;
				end;

				%if &XVARINF ne %then
					%do;
						lambda0=(ncol(pos0)-sum((parg/(uj+parg))##parg))/n;

						if lambda0<=0 then
							do;
								lambdag=j(ncol(G), 1, 0);
								print 
									"NOTE: Expected number of zeros ("(char(sum((parg/(uj+parg))##parg), 
									6, 2))") >= number of zeros ("(char(ncol(pos0), 4, 
									0))"). No Need of Zero Model.";
							end;
						else
							do;
								print 
									"NOTE: Expected number of zeros ("(char(sum((parg/(uj+parg))##parg), 
									6, 2))") < number of zeros ("(char(ncol(pos0), 4, 
									0))"). Zero Model Used.";
								lambda0=log(lambda0/(1-lambda0));
								lambdag=lambda0//j(ncol(G)-1, 1, 0);
							end;
						pargg=parg;
						ujg=uj;

						if nrow(pos0)=0 | any(lambdag)=0 then
							do;

								if nrow(pos0)=0 then
									do;
										pos0=pos1;

										%if %upcase(&MODEL)=ZINB or %upcase(&MODEL)=ZIP %then
											%do;
												call symputx("MODEL", "NEGBIN");
											%end;
									end;

								%if %upcase(&FORCE) ne YES %then
									%do;
										call symputx("MODEL", "NEGBIN");
									%end;
							end;
					%end;
				njl=G*lambdag;

				%if %upcase(&MODEL) ne ZIP and %upcase(&MODEL) ne ZINB %then
					%do;
						zkg=0;
					%end;
				%else
					%do;
						zkg=1/(1+exp(-G*lambdag)#(parg/(parg+uj))##parg);
						zkg=choose(y>0, 0, zkg);
					%end;
				dllike=1;
				llikeg=0;
				j=0;

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

						%if %upcase(&MODEL)=ZIP or %upcase(&MODEL)=POISSON %then
							%do;
								alphag=1E-6;
								parg=1/alphag;
							%end;
						%else
							%do;
								*if j>0 then parg=1/(sum((y-uj)##2/uj)/(n-nvar));

								do while (abs(dpar)>0.0001 & aux2<200);

									if parg<0 then
										parg=0.00001;
									parg=choose(parg<1E-10, 1E-10, parg);
									gf=sum((1-zkg)#(digamma(parg+y)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+y)/(parg+uj)));
									hessg=sum((1-zkg)#(trigamma(parg+y)-trigamma(parg)+1/parg-2/(parg+uj)+(y+parg)/(parg+uj)##2));
									hessg=choose(hessg=0, 1E-23, hessg);
									par0=parg;
									parg=par0-inv(hessg)*gf;

									if parg>1E5 then
										do;
											dpar=0.0001;
											cont3=cont3+1;

											if cont3=1 then
												parg=2;
											else if cont3=2 then
												parg=1E5;
											else if cont3=3 then
												parg=0.0001;
										end;
									else
										dpar=parg-par0;

									if parg>1E6 then
										do;
											parg=1E6;
											dpar=0;
										end;
									aux2=aux2+1;
								end;
								alphag=1/parg;
							%end;
						devg=0;
						ddev=1;
						nj=x*bg+offset;
						nj=choose(nj>700, 700, nj);
						nj=choose(nj<-700, -700, nj);
						uj=exp(nj);

						do while (abs(ddev)>0.000001 & aux1<100);
							uj=choose(uj>1E100, 1E100, uj);
							Ai=(1-zkg)#((uj/(1+alphag*uj)+(y-uj)#(alphag*uj/(1+2*alphag*uj+alphag**2*uj##2))));
							Ai=choose(Ai<=0, 1E-5, Ai);
							uj=choose(uj<1E-150, 1E-150, uj);
							zj=(nj+(y-uj)/(((uj/(1+alphag*uj)+(y-uj)#(alphag*uj/(1+2*alphag*uj+alphag**2*uj##2))))#(1+alphag*uj)))-offset;

							if det(x`*(Ai#x))=0 then
								bg=j(nvar, 1, 0);
							else
								bg=inv(x`*(Ai#x))*x`*(Ai#zj);
							nj=x*bg+offset;
							nj=choose(nj>700, 700, nj);
							nj=choose(nj<-700, -700, nj);
							uj=exp(nj);
							olddev=devg;
							gamma1=/*(gamma(par+y)/(gamma(y+1)#gamma(par)))#*/(uj/(uj+parg))##y#(parg/(uj+parg))##parg;
							gamma1=choose(gamma1<=0, 1E-10, gamma1);
							devg=sum((1-zkg)#(log(gamma1)));
							ddev=devg-olddev;
							*print bg aux1 devg olddev ddev;
							aux1=aux1+1;
						end;
						ddpar=parg-parold;
						cont=cont+1;
					end;

					%if %upcase(&MODEL)=ZIP or %upcase(&MODEL)=ZINB %then
						%do;
							devg=0;
							ddev=1;
							njl=G*lambdag;
							njl=choose(njl>&MAXG, &MAXG, njl);
							njl=choose(njl<-&MAXG, -&MAXG, njl);
							pig=exp(njl)/(1+exp(njl));

							do while (abs(ddev)>0.000001 & aux3<100);
								Ai=pig#(1-pig);
								Ai=choose(Ai<=0, 1E-5, Ai);
								zj=njl+(zkg-pig)/Ai;

								if det(G`*(Ai#G))=0 then
									lambdag=j(ncol(G), 1, 0);
								else
									lambdag=inv(G`*(Ai#G))*G`*(Ai#zj);
								njl=G*lambdag;
								njl=choose(njl>&MAXG, &MAXG, njl);
								njl=choose(njl<-&MAXG, -&MAXG, njl);
								pig=exp(njl)/(1+exp(njl));
								olddev=devg;
								devg=sum(zkg#njl-log(1+exp(njl)));
								ddev=devg-olddev;
								*print lambdag devg olddev ddev;
								aux3=aux3+1;
							end;
						%end;
					zkg=1/(1+exp(-njl)#(parg/(parg+uj))##parg);
					zkg=choose(y>0, 0, zkg);

					%if %upcase(&MODEL) ne ZIP and %upcase(&MODEL) ne ZINB %then
						%do;
							zkg=0;
						%end;
					oldllike=llikeg;
					llikeg=sum(zkg#njl-log(1+exp(njl))+(1-zkg)#(log(gamma1)));
					dllike=llikeg-oldllike;
					*print j bg alphag lambdag llikeg dllike;
					j=j+1;
				end;
				g1x=parg/(parg+uj);
				g2x=uj/(parg+uj);
				hgx=exp(njl)+g1x##parg;
				daa=zkg#((g1x##parg#(log(g1x)+g2x))##2#(1-1/hgx)/hgx+g1x##parg#(g2x##2/parg)/hgx)+(1-zkg)#(trigamma(parg+y)-trigamma(parg)-2/(uj+parg)+1/parg+(y+parg)/(uj+parg)##2);
				dab=zkg#(g1x##(2*parg+1)#uj#(log(g1x)+g2x)/hgx##2-g1x##parg#(-g2x##2+parg#g2x#(log(g1x)+g2x))/hgx)+(1-zkg)#(g2x#(y-uj)/(uj+parg));
				dal=-zkg#(exp(njl)#g1x##parg#(log(g1x)+g2x)/hgx##2);
				daa=daa*parg**4;
				dab=dab*parg**2;

				if any(lambdag)=0 then
					Iy=j(nrow(y), 1, 0);
				dll=Iy#(exp(njl)#g1x##parg/hgx##2)-exp(njl)/(1+exp(njl))##2;
				dbb=Iy#(-(parg*g1x##parg#g2x/hgx)##2+parg**2*g1x##parg#g2x##2#(1-1/uj)/hgx)-(1-Iy)#(parg*g2x#(1+(y-uj)/(parg+uj)));
				dlb=Iy#(parg*exp(njl)#g1x##parg#g2x/hgx##2);
				I1=j(nrow(y), 1, 1);
				II=-(I1#daa)`*I1||-(I1#dab)`*X||-(I1#dal)`*G//(-X`*(dab#I1)||-(X#dbb)`*X||-(X#dlb)`*G)//(-G`*(dal#I1)||-G`*(X#dlb)||-(G#dll)`*G);

				if all(lambdag)>0 & alphag=1E-6 then
					II=II[2:nrow(II), 2:nrow(II)];
				else if any(lambdag)=0 & alphag>1E-6 then
					II=II[1:ncol(x)+1, 1:ncol(x)+1];
				else if any(lambdag)=0 & alphag=1E-6 then
					II=II[2:ncol(x)+1, 2:ncol(x)+1];
				varabetalambdag=vecdiag(inv(II));
				stdabetalambdag=sqrt(abs(varabetalambdag));
				varcovg=inv(II);
				*print bg lambdag alphag devg llikeg stdabetalambdag;

				/*****************************************/


				%IF %UPCASE(&METHOD)=ADAPTIVEN %THEN
					%DO;
						use &DHV;
						read all into hv;
					%END;

				%IF &H NE %THEN
					%DO;
						h=&H;
						print h[label="Bandwidth"];
					%END;
				%ELSE
					%DO;

						%IF %UPCASE(&METHOD)=FIXED_G or %UPCASE(&METHOD)=FIXED_BSQ or 
							%UPCASE(&METHOD)=ADAPTIVE_BSQ %THEN
								%DO;
								h=&_h_;
								print h[label="Bandwidth"];
							%END;
					%END;
				read all var{&LONG &LAT} into COORD;
				%if &grid=%then
					%do;
						read all var{&LONG &LAT} into POINTS;
					%end;
				close &DATA;

				%if &grid^=%then
					%do;
						use &grid;
						read all var{&LONG &LAT} into POINTS;
						close &GRID;
					%end;
				m=nrow(POINTS);
				bi=j(ncol(x)*m, 4, 0);
				li=j(ncol(G)*m, 4, 0);
				alphai=j(m, 3, 0);
				BB=j(ncol(x)*n, n, 0);
				BBl=j(ncol(G)*n, n, 0);
				sumwi=j(m, 1, 0);
				varbi=j(ncol(x)*m, 1, 0);
				varli=j(ncol(G)*m, 1, 0);
				S=j(m, 1, 0);
				Si=j(m, 1, 0);
				S2=j(m, 1, 0);
				biT=j(m, ncol(x)+1, 0);
				ym=y-y[:];

				/******** calculating distance **********/;
				_dist_=distance(COORD, POINTS, "L2");
				seq=1:n;
				create _dist_ from _dist_;
				append from _dist_;

				do i=1 to m;

					do j=1 to n;
						seqi=j(n, 1, i);
						dist=seqi||seq`||_dist_[, i];

						%IF %UPCASE(&DISTANCEKM)=YES %THEN
							%DO;
								dist[, 3]=dist[, 3]*111;
							%END;
					end;
					u=nrow(dist);
					w=j(u, 1, 0);

					do jj=1 to u;

						%IF %UPCASE(&METHOD)=FIXED_G %THEN
							%DO;
								w[jj]=exp(-0.5*(dist[jj, 3]/h)**2);
							%END;
						%ELSE %IF %UPCASE(&METHOD)=FIXED_BSQ %THEN
							%DO;
								w[jj]=(1-(dist[jj, 3]/h)**2)**2;
							%END;
						%ELSE %IF %UPCASE(&METHOD)=ADAPTIVEN %THEN
							%DO;
								w[jj]=(1-(dist[jj, 3]/hv[i])**2)**2;
							%END;

						%if &grid=%then
							%do;
							end;

							%IF %UPCASE(&METHOD)=ADAPTIVE_BSQ %THEN
								%DO;
									call sort(dist, {3});
									dist=dist||(1:nrow(dist))`;
									w=j(n, 2, 0);
									hn=dist[h, 3];

									do jj=1 to n;

										if dist[jj, 4]<=h then
											w[jj, 1]=(1-(dist[jj, 3]/hn)**2)**2;
										else
											w[jj, 1]=0;
										w[jj, 2]=dist[jj, 2];
									end;
									call sort(w, {2});
									w=w[, 1];
								%END;
						%end;
					%else
						%do;
						end;

						%IF %UPCASE(&METHOD)=ADAPTIVE_BSQ %THEN
							%DO;
								call sort(dist, {3});
								dist=dist||(1:n)`;
								w=j(n, 2, 0);
								hn=dist[h, 3];

								do jj=1 to n;

									if dist[jj, 4]<=h then
										w[jj, 1]=(1-(dist[jj, 3]/hn)**2)**2;
									else
										w[jj, 1]=0;
									w[jj, 2]=dist[jj, 2];
								end;
								call sort(w, {2});
								w=w[, 1];
							%END;
					%end;

				/****** MODEL SELECTION *************/
				Iy=Iy2;
				b=bg;
				b2=b;
				nj=X*b+offset;
				uj=exp(nj);
				par=parg;
				par2=par;
				lambda=lambdag;
				njl=G*lambda;
				njl=choose(njl>&MAXG, &MAXG, njl);
				njl=choose(njl<-&MAXG, -&MAXG, njl);

				%if %upcase(&MODEL) ne ZIP and %upcase(&MODEL) ne ZINB %then
					%do;
						zk=0;
					%end;
				%else
					%do;

						/*
						locdh=loc(_dist_[,i]<=h);
						y1=y[locdh];
						x1=x[locdh,];
						offset1=offset[locdh];
						pos0dh=loc(y1=0);
						
						uj=(y1+y1[:])/2;
						nj=log(uj);
						parg=sum((y1-uj)##2/uj)/(nrow(y1)-nvar);
						ddpar=1;
						cont=1;
						cont3=0;
						do while (abs(ddpar)>0.000001 & cont<100);
						dpar=1;
						parold=parg;
						cont1=1;
						%IF %UPCASE(&MODEL)=ZIP or %UPCASE(&MODEL)=POISSON %THEN %DO;
						parg=1/1E-6;alphag=1/parg;
						%END;
						%IF %UPCASE(&MODEL)=ZINB or %UPCASE(&MODEL)=NEGBIN %THEN %DO;
						if cont>1 then parg=1/(sum((y1-uj)##2/uj)/(nrow(y1)-nvar));
						do while (abs(dpar)>0.0001 & cont1<200);
						if parg<0 then parg=0.00001;
						parg=choose(parg<1E-10,1E-10,parg);
						gf=sum(digamma(parg+y1)-digamma(parg)+log(parg)+1-log(parg+uj)-(parg+y1)/(parg+uj));
						hessg=sum(trigamma(parg+y1)-trigamma(parg)+1/parg-2/(parg+uj)+(y1+parg)/(parg+uj)##2);
						hessg=choose(hessg=0,1E-23,hessg);
						par0=parg;
						parg=par0-inv(hessg)*gf;
						if parg>1E5 then do;
						dpar= 0.0001;
						cont3=cont3+1;
						if cont3=1 then parg=2 ;
						else if cont3=2 then parg=1E5;
						else if cont3=3 then parg=0.0001;
						end;
						else dpar=parg-par0;
						cont1=cont1+1;
						if parg>1E6 then do;parg=1E6;dpar=0;end;
						end;
						alphag=1/parg;
						
						%END;
						devg=0; ddev=1; cont2=0;
						do while (abs(ddev)>0.000001 & cont2<100);
						Ai=(uj/(1+alphag*uj))+(y1-uj)#(alphag*uj/(1+2*alphag*uj+alphag**2*uj#uj));;
						Ai=choose(Ai<=0,1E-5,Ai);
						zj=nj+(y1-uj)/(Ai#(1+alphag*uj))-offset1;
						if det(x1`*(Ai#x1))=0 then bg=j(ncol(x),1,0);
						else bg=inv(x1`*(Ai#x1))*x1`*(Ai#zj);
						nj=x1*bg+offset1;
						nj=choose(nj>700,700,nj);
						uj=exp(nj);
						olddev=devg;
						uj=choose(uj<1E-150,1E-150,uj);
						uj=choose(uj>100000,100000,uj);
						tt=y1/uj;
						tt=choose(tt=0,1E-10,tt);
						devg=2*sum(y1#log(tt)-(y1+1/alphag)#log((1+alphag*y1)/(1+alphag*uj)));
						if cont2>100 then ddev= 0.0000001;
						else ddev=devg-olddev;
						cont2=cont2+1;
						end;
						cont=cont+1;
						ddpar=parg-parold;
						end;
						lambda0=(ncol(pos0dh)-sum((parg/(uj+parg))##parg))/n;
						print i (ncol(pos0dh)) alphag bg (sum((parg/(uj+parg))##parg)) lambda0;
						*/
						lambda0=(ncol(pos0)-sum((parg/(uj+parg))##parg))/n;

						if lambda0>0 then
							do;
								lambda0=log(lambda0/(1-lambda0));
								lambda=lambda0//j(ncol(G)-1, 1, 0);
								njl=G*lambda;
							end;

						/*nj=X*b2+offset;
						uj=exp(nj);
						par=par2;*/
						zk=1/(1+exp(-njl)#(par/(par+uj))##par);
						zk=choose(y>0, 0, zk);
					%end;
				dllike=1;
				llike=0;
				j=1;
				*if i=1 then print par;
				do while (abs(dllike)>0.00001 & j<=600);
				*if i=1 then print dllike;
					ddpar=1;
					*do while (abs(ddpar)>0.000001);
					dpar=1;
					parold=par;
					aux1=1;
					aux2=1;
					int=1;

					%IF %UPCASE(&MODEL)=ZIP or %UPCASE(&MODEL)=POISSON %THEN
						%DO;
							alpha=1E-6;
							par=1/alpha;
						%END;

					%IF %UPCASE(&MODEL)=ZINB or %UPCASE(&MODEL)=NEGBIN %THEN
						%DO;

							/*if par<=1E-5 then do;if i>1 then par=1/alphai[i-1,2];end;*/
							if par>=1E6 then
								do;
									par=1E6;
									dpar=0;
									alpha=1/par;
									b=bg;
									uj=exp(X*b+offset);
									lambda=lambdag;
									njl=G*lambda;
									njl=choose(njl>&MAXG, &MAXG, njl);
									njl=choose(njl<-&MAXG, -&MAXG, njl);
									zk=1/(1+exp(-njl)#(parg/(parg+uj))##parg);
									zk=choose(y>0, 0, zk);

									if any(lambda)=0 then
										zk=0;
								end;

							do while (abs(dpar)>0.000001 & aux2<200);
								par=choose(par<1E-10, 1E-10, par);
								gf=sum(w#wt#(1-zk)#(digamma(par+y)-digamma(par)+log(par)+1-log(par+uj)-(par+y)/(par+uj)));
								*if i=1 & j=1 & aux2=1 then print (sum(w)) (sum(zk)) (sum(uj));
								hess=sum(w#wt#(1-zk)#(trigamma(par+y)-trigamma(par)+1/par-2/(par+uj)+(y+par)/(par+uj)##2));
								*if i=1 & j=1 & aux2=1 then print hess;
								hess=choose(hess=0, 1E-23, hess);
								par0=par;
								*if i=1 & j=1 & aux2=1 then print hess;
								par=par0-inv(hess)*gf;
								*if i=1 & j=1 & aux2=1 then print par;
								dpar=par-par0;

								if par>=1E6 then
									do;
										par=1E6;
										dpar=0;
										alpha=1/par;
										b=bg;
										uj=exp(X*b+offset);
										lambda=lambdag;
										njl=G*lambda;
										njl=choose(njl>&MAXG, &MAXG, njl);
										njl=choose(njl<-&MAXG, -&MAXG, njl);
										zk=1/(1+exp(-njl)#(parg/(parg+uj))##parg);
										zk=choose(y>0, 0, zk);

										if any(lambda)=0 then
											zk=0;
									end;
								*print par aux2 dpar;
								aux2=aux2+1;
							end;

							if par<=1E-5 then
								do;
									par=1E6;
									b=bg;
									uj=exp(X*b+offset);
									lambda=lambdag;
									njl=G*lambda;
									njl=choose(njl>&MAXG, &MAXG, njl);
									njl=choose(njl<-&MAXG, -&MAXG, njl);
									zk=1/(1+exp(-njl)#(parg/(parg+uj))##parg);
									zk=choose(y>0, 0, zk);

									if any(lambda)=0 then
										zk=0;
								end;
							alpha=1/par;
							*print i j b lambda par alpha aux2;
						%END;
					*if i=244 then print lambda;
					dev=0;
					ddev=1;
					nj=x*b+offset;
					nj=choose(nj>700, 700, nj);
					nj=choose(nj<-700, -700, nj);
					uj=exp(nj);

					do while (abs(ddev)>0.000001 & aux1<100);
						uj=choose(uj>1E100, 1E100, uj);
						Ai=(1-zk)#((uj/(1+alpha*uj)+(y-uj)#(alpha*uj/(1+2*alpha*uj+alpha**2*uj##2))));
						Ai=choose(Ai<=0, 1E-5, Ai);
						uj=choose(uj<1E-150, 1E-150, uj);
						denz=(((uj/(1+alpha*uj)+(y-uj)#(alpha*uj/(1+2*alpha*uj+alpha**2*uj##2))))#(1+alpha*uj));
						denz=choose(denz=0, 1E-5, denz);
						zj=(nj+(y-uj)/denz)-offset;

						if det(x`*(w#Ai#x#wt))=0 then
							b=j(nvar, 1, 0);
						else
							b=inv(x`*(w#Ai#x#wt))*x`*(w#Ai#wt#zj);
						nj=x*b+offset;
						nj=choose(nj>700, 700, nj);
						nj=choose(nj<-700, -700, nj);
						uj=exp(nj);
						olddev=dev;
						uj=choose(uj>1E10, 1E10, uj);
						uj=choose(uj=0, 1E-10, uj);
						*if i=1 & j=1 & aux1=1 then print (sum(uj)) par;
						if par=1E6 then
							gamma1=/*(gamma(par+y)/(gamma(y+1)#gamma(par)))#*/(uj/(uj+par))##y#exp(-uj);
						else
							gamma1=/*(gamma(par+y)/(gamma(y+1)#gamma(par)))#*/(uj/(uj+par))##y#(par/(uj+par))##par;
						gamma1=choose(gamma1<=0, 1E-10, gamma1);
						dev=sum((1-zk)#(log(gamma1)));
						ddev=dev-olddev;
						*print b par aux1 dev olddev ddev;
						aux1=aux1+1;
					end;
					ddpar=par-parold;
					*end;

					%if %upcase(&MODEL)=ZIP or %upcase(&MODEL)=ZINB %then
						%do;

							if j=1 then
								do;
									alphatemp=alpha;
									lambdatemp=lambda[1];
								end;
							else
								do;
									alphatemp=alphatemp//alpha;
									lambdatemp=lambdatemp//lambda[1];
								end;
							alphatemp=round(alphatemp, 0.0000001);
							lambdatemp=round(lambdatemp, 0.0000001);
							*print i j alphatemp (nrow(alphatemp)) (ncol(unique(alphatemp)));
							*if i=244 & j=46 then print (nrow(alphatemp)) (ncol(unique(alphatemp)));
							*if i=244 then print j;

							%if %upcase(&MODEL)=ZINB %then
								%do;

									if j>300 & nrow(alphatemp)>ncol(unique(alphatemp)) & 
nrow(lambdatemp)>ncol(unique(lambdatemp))then
										do;
										%end;

									%if %upcase(&MODEL)=ZIP %then
										%do;
											*print i j lambdatemp (nrow(lambdatemp)) (ncol(unique(lambdatemp)));

											if j>300 & nrow(lambdatemp)>ncol(unique(lambdatemp))then
												do;
												%end;
											lambda=j(ncol(G), 1, 0);
											njl=G*lambda;
											zk=j(n, 1, 0);
										end;
									else
										do;
											aux3=1;
											dev=0;
											ddev=1;
											njl=G*lambda;
											*if i=94 & j=13 then print (sum(lambda));
											njl=choose(njl>&MAXG, &MAXG, njl);
											njl=choose(njl<-&MAXG, -&MAXG, njl);
											pi=exp(njl)/(1+exp(njl));
											*if i=94 & j=14 then print lambda;
											do while (abs(ddev)>0.000001 & aux3<100);
												*if i=26 & j=13 & aux3=12 then print pi;
												Aii=pi#(1-pi);
												Aii=choose(Aii<=0, 1E-5, Aii);
												zj=njl+(zk-pi)/Aii;
												*if i=26 then print j;
												*if i=26 & j=12 then print aux3;
												*if i=26 & j=13 & aux3=12 then print (sum(Aii)) (sum(w));
												if det((G#Aii#w#wt)`*G)=0 then
													lambda=j(ncol(G), 1, 0);
												else
													lambda=inv((G#Aii#w#wt)`*G)*(G#Aii#w#wt)`*zj;
												*if i=94 & j=13 & aux3=21 then print (sum(Aii)) (sum(w)) (sum(zj));
												*if i=94 & j=13 & aux3=23 then print (sum(lambda));
												njl=G*lambda;
												njl=choose(njl>&MAXG, &MAXG, njl);
												njl=choose(njl<-&MAXG, -&MAXG, njl);
												pi=exp(njl)/(1+exp(njl));
												olddev=dev;
												dev=sum(zk#njl-log(1+exp(njl)));
												ddev=dev-olddev;
												*print lambda aux3 dev olddev ddev;
												aux3=aux3+1;
												*if i=26 & j=12 & aux3=14 then print dev olddev;
											end;
											*if i=26 then print (sum(lambda));
										end;
								%end;
							*if i=94 & j=13 then print (sum(lambda));
							njl=G*lambda;
							if i=94 & j=13 then do;
								*print (sum(njl));
							end;
							njl=choose(njl>&MAXG, &MAXG, njl);
							njl=choose(njl<-&MAXG, -&MAXG, njl);
							zk=1/(1+exp(-njl)#(par/(par+uj))##par);
							zk=choose(y>0, 0, zk);

							if any(lambda)=0 then
								zk=j(n, 1, 0);

							%if %upcase(&MODEL) ne ZIP and %upcase(&MODEL) ne ZINB %then
								%do;
									zk=0;
								%end;
							oldllike=llike;
							llike=sum(zk#(njl)-log(1+exp(njl))+(1-zk)#(log(gamma1)));
							*if i=1 & j=1 then print llike;
							dllike=llike-oldllike;
							*print i j b alpha lambda llike dllike;
							j=j+1;
							if i=94 & j=15 then do;
								*print oldllike llike dllike;
							end;
						end;

					/**** COMPUTING VARIANCE OF BETA AND LAMBDA ******/
					if det(x`*(w#Ai#x#wt))=0 then
						C=j(ncol(x), nrow(x), 0);
					else
						C=inv(x`*(w#Ai#x#wt))*x`#(w#Ai#wt)`;
					Ci=j(ncol(G), 1, 0);

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
					*if i=1 then print daa;
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
					*if i=1 then print II;
					if all(lambda)>0 & alpha=1E-6 then
						II=II[2:nrow(II), 2:nrow(II)];
					else if any(lambda)=0 & alpha>1E-6 then
						II=II[1:ncol(x)+1, 1:ncol(x)+1];
					else if any(lambda)=0 & alpha=1E-6 then
						II=II[2:ncol(x)+1, 2:ncol(x)+1];
					*if i=1 then print II;
					if det(II)=0 then
						do;
						*print lambda alpha;

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
					else do;
						varabetalambda=vecdiag(inv(II));
					end;
					*if i=1 then print lambda alpha;
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
					*if i=1 then print varb;
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
							*if i=1 then print lambda;
							li[m1:m2, 1]=i;
							li[m1:m2, 2]=lambda;
							li[m1:m2, 3]=POINTS[i, 1];
							li[m1:m2, 4]=POINTS[i, 2];
							*if i=17 | i=18 then do;
							*	print varabetalambda;
							*	print varl;
							*end;
							varli[m1:m2, 1]=varl;
						%end;
					*if i=1 then print (sum(C));
					%if &grid=%then
						%do;
							r=x[i, ]*C;
							*if i=1 then print r;
							S[i]=r[i];
							S2[i]=r*r`;
							yhat[i]=uj[i];
							pihat[i]=njl[i];

							%if %upcase(&MODEL)=ZIP or %upcase(&MODEL)=ZINB %then
								%do;
									ri=G[i, ]*Ci;
									*if i=1 then print ri;
									Si[i]=ri[i];
									yhat2[i]=uj[i];
									if i=94 then do;
									*print uj;
									*print njl;
									end;
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
									else do;
										BB[m1:m2, ]=inv(CCC[, 1:ncol(x)]`*(CCC[, ncol(CCC)-1]#CCC[, 
											1:ncol(x)]#CCC[, ncol(CCC)]))*CCC[, 1:ncol(x)]`#(CCC[, 
											ncol(CCC)-1]#CCC[, ncol(CCC)])`;
										*if i=1 then print (inv(CCC[, 1:ncol(x)]`*(CCC[, ncol(CCC)-1]#CCC[, 
											1:ncol(x)]#CCC[, ncol(CCC)]))*CCC[, 1:ncol(x)]`);
										end;

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
				end; *fecha for i?;
				create w_f from W_f;
				append from W_f;
				close w_f;
				free w_f;

				/***********************************************/
				*print (sum(S)) (sum(Si));
				%if &grid=%then
					%do;
						v1=sum(S)+sum(Si);
						v11=sum(S)+sum(Si);
						v2=sum(S2);
						*print v1 v11 v2;
						nparmodel=n-v11;

						if v11<v2 then
							v1=v11;
						*print yhat;
						res=(y-yhat);
						rsqr1=(res#wt)`*res;
						ym=(y#wt)`*y;
						rsqr2=ym-((y#wt)[+]**2)/wt[+];
						rsqr=1-rsqr1/rsqr2;
						rsqradj=1-((n-1)/(n-v1))*(1-rsqr);
						*print rsqr1 v1;
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
						*print bi;
						%IF %UPCASE(&MODEL)=NEGBIN or %UPCASE(&MODEL)=ZINB %THEN
							%DO;
								_alpha_=shape(alphai[, 1:2], n);
								_beta2_=_beta_||_alpha_;
							%END;
						i=do(2, ncol(_beta_), 2);
						_beta_=_beta_[, i];
						*print (sum(_beta_));
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
						*print stdbi;
						*print _stdbeta_;
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
					%end; *fecha o grid;

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
										*if i=1 & k=1 then do;
										*	print (sum(BB));
										*end;
										BBk[i, ]=ek`*BB[m1:m2, ];
									end;
									*if k=1 then print BBk;
									Vk[k]=y`*(1/n)*BBk`*(I(n)-(1/n)*J(n, n, 1))*BBk*y;
									*if k=1 then print (Vk[k]);
									df1k[k]=trace((1/n)*BBk`*(I(n)-(1/n)*J(n, n, 1))*BBk);
									df2k[k]=trace(((1/n)*BBk`*(I(n)-(1/n)*J(n, n, 1))*BBk)**2);
									*if k=1 then print (((1/n)*BBk`*(I(n)-(1/n)*J(n, n, 1))*BBk)**2);
								end;
								*print (sum(df2k));
								*print df2k;
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