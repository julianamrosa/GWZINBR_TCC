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
				Ai=(uj/(1+alphag*uj))+(y-uj)#(alphag*uj/(1+2*alphag*uj+alphag**2*uj#uj));
				Ai=choose(Ai<=0, 1E-5, Ai);
				*print (sum(Ai));
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
					*print (sum(Ai));
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
						*print (sum(Ai));
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
					*if i=244 & jj=u then print h;

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
					nj=x*b+offset;
					nj=choose(nj>700, 700, nj);
					nj=choose(nj<-700, -700, nj);
					uj=exp(nj);

					contador5=0;
					do while (abs(ddev)>0.000001 & aux1<100);
						contador5=contador5+1;
						uj=choose(uj>1E100, 1E100, uj);
						Ai=(1-zk)#((uj/(1+alpha*uj)+(y-uj)#(alpha*uj/(1+2*alpha*uj+alpha**2*uj##2))));
						Ai=choose(Ai<=0, 1E-5, Ai);
						if i=58 & contador4=3 & contador5=4 then print (sum(Ai));
						uj=choose(uj<1E-150, 1E-150, uj);
						denz=(((uj/(1+alpha*uj)+(y-uj)#(alpha*uj/(1+2*alpha*uj+alpha**2*uj##2))))#(1+alpha*uj));
						denz=choose(denz=0, 1E-5, denz);
						zj=(nj+(y-uj)/denz)-offset;

						if det(x`*(w#Ai#x#wt))=0 then
							b=j(nvar, 1, 0);
						else do;
							b=inv(x`*(w#Ai#x#wt))*x`*(w#Ai#wt#zj);
						end;
						nj=x*b+offset;
						nj=choose(nj>700, 700, nj);
						nj=choose(nj<-700, -700, nj);
						uj=exp(nj);
						olddev=dev;
						uj=choose(uj>1E10, 1E10, uj);
						uj=choose(uj=0, 1E-10, uj);

						if par=1E6 then
							gamma1=/*(gamma(par+y)/(gamma(y+1)#gamma(par)))#*/(uj/(uj+par))##y#exp(-uj);
						else
							gamma1=/*(gamma(par+y)/(gamma(y+1)#gamma(par)))#*/(uj/(uj+par))##y#(par/(uj+par))##par;
						gamma1=choose(gamma1<=0, 1E-10, gamma1);
						dev=sum((1-zk)#(log(gamma1)));
						ddev=dev-olddev;
						*if i=244 then print b par aux1 dev olddev ddev;
						*print comentado;
						aux1=aux1+1;
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
						end; *fecha loop while (condicoes j e dlike);

					%IF %UPCASE(&METHOD)=FIXED_G or %UPCASE(&METHOD)=FIXED_BSQ or 
						%UPCASE(&METHOD)=ADAPTIVE_BSQ %THEN
							%DO;
							yhat[i]=uj[i];
							pihat[i]=njl[i];
							alphai[i]=alpha;

							if det(x`*(w#Ai#x#wt))=0 then do;
								S[i]=0;
							end;
							else do;
								S[i]=(x[i, ]*inv(x`*(w#Ai#x#wt))*(x#w#Ai#wt)`)[i];
								*if i=59 then print (sum(Ai));
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
						CV=((y-yhat)#wt)`*(y-yhat);
						_par_=1/alphai;

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
									end;
								dev=2*(llnull1-ll);
								npar=sum(S)+sum(Si);
								*print (sum(S));
								AIC=2*npar-2*ll;
								AICc=AIC+2*(npar*(npar+1)/(n-npar-1));

								%IF %UPCASE(&MODEL)=ZINB %THEN
									%DO;
										AIC=2*(npar+npar/(ncol(x)+ncol(G)))-2*ll;
										AICc=AIC+2*((npar+npar/(ncol(x)+ncol(G)))*((npar+npar/(ncol(x)+ncol(G)))+1)/(n-(npar+npar/(ncol(x)+ncol(G)))-1));
									%END;
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
						*print npar;
						*print ll;
						*print G;
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
				res=cv||npar;
				*print res;
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
								GMY=1;
								ax1=lower[GMY];
								bx1=upper[GMY];
							%END;
						h0=ax1;
						h3=bx1;
						h1=bx1-r*(bx1-ax1);
						h2=ax1+r*(bx1-ax1);
						print h0 h1 h2 h3;

						/***************************************/
						res1=cv(h1);
						CV1=res1[1];
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
							*print int;
							if CV2<CV1 then
								do;
									*print ("entrou no if");
									h0=h1;
									h1=h3-r*(h3-h0);
									h2=h0+r*(h3-h0);
									CV1=CV2;
									*print h2;
									res2=cv(h2);
									CV2=res2[1];
								end;
							else
								do;
									*print ("entrou no else");
									h3=h2;
									h1=h3-r*(h3-h0);
									h2=h0+r*(h3-h0);
									CV2=CV1;
									*print h1;
									res1=cv(h1);
									CV1=res1[1];
								end;

							%IF %UPCASE(&METHOD)=FIXED_G or %UPCASE(&METHOD)=FIXED_BSQ or 
								%UPCASE(&METHOD)=ADAPTIVE_BSQ %THEN
									%DO;
									append;
								%END;
							int=int+1;
							*print cv1;
							*print cv2;
							end;
		%mend Golden;