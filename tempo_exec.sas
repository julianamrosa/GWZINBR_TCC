libname tcc '/home/u41131808/TCC2 - GWZINBR';

data Korea_base;
	set tcc.korea_base_artigo;
run;

*POISSON;

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=POISSON,METHOD=ADAPTIVE_BSQ,BANDWIDTH=aic, GLOBALMIN = no,DISTANCEKM=YES);
* 9 segs;

%GWZINBR(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, LONG=x,LAT=y,
OFFSET=ln_total,MODEL=poisson,METHOD=ADAPTIVE_BSQ,DISTANCEKM=YES, H=79);
*1 seg;

*NEGBIN;

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=NEGBIN,METHOD=ADAPTIVE_BSQ,BANDWIDTH=aic, GLOBALMIN = no,DISTANCEKM=YES);
*18 segs;

%GWZINBR(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, LONG=x,LAT=y,
OFFSET=ln_total,MODEL=negbin,METHOD=ADAPTIVE_BSQ,DISTANCEKM=YES, H=82);
*2 segs;

*PIZ;

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access Crowding,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZIP,METHOD=ADAPTIVE_BSQ,BANDWIDTH=aic, GLOBALMIN = no,DISTANCEKM=YES);
*3 mins;

%GWZINBR(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access Crowding,LONG=x,LAT=y,
OFFSET=ln_total,MODEL=zip,METHOD=ADAPTIVE_BSQ,DISTANCEKM=YES, H=56);
*35 segs;

*BNIZ;

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access Crowding,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZINB,METHOD=ADAPTIVE_BSQ,BANDWIDTH=aic, GLOBALMIN = no,DISTANCEKM=YES);
*6 mins;

%GWZINBR(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access Crowding,LONG=x,LAT=y,
OFFSET=ln_total,MODEL=zinb,METHOD=ADAPTIVE_BSQ,DISTANCEKM=YES, H=82);
*35 segs;