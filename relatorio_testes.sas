libname tcc '/home/u41131808/TCC2 - GWZINBR';

data Korea_base;
	set tcc.korea_base_artigo;
run;

/*** GOLDEN ***/

*adaptive;

*zinb, adaptive, aic;
%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access Crowding,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZINB,METHOD=ADAPTIVE_BSQ,BANDWIDTH=aic, GLOBALMIN = no,DISTANCEKM=YES);
*6 mins;

*zip, adaptive, aic;
%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access Crowding,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZIP,METHOD=ADAPTIVE_BSQ,BANDWIDTH=aic, GLOBALMIN = no,DISTANCEKM=YES);
*3 mins;

*zinb, adaptive, cv;
%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access Crowding,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZINB,METHOD=ADAPTIVE_BSQ,BANDWIDTH=cv, GLOBALMIN = no,DISTANCEKM=YES);
*2 mins;

*zip, adaptive, cv;
%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access Crowding,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZIP,METHOD=ADAPTIVE_BSQ,BANDWIDTH=cv, GLOBALMIN = no,DISTANCEKM=YES);
*25 segs;

*fixed;

*zinb, fixed, aic;
%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access Crowding,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZINB,METHOD=FIXED_G,BANDWIDTH=aic, GLOBALMIN = no,DISTANCEKM=YES);
*4 mins;

*zip, fixed, aic;
%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access Crowding,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZIP,METHOD=FIXED_G,BANDWIDTH=aic, GLOBALMIN = no,DISTANCEKM=YES);
*4 mins;

*zinb, fixed, cv;
%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access Crowding,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZINB,METHOD=FIXED_G,BANDWIDTH=cv, GLOBALMIN = no,DISTANCEKM=YES);
*3 mins;

*zip, fixed, cv;
%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access Crowding,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZIP,METHOD=FIXED_G,BANDWIDTH=cv, GLOBALMIN = no,DISTANCEKM=YES);
*24 segs;