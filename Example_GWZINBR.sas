/* Geram os mesmos resultados do artigo */

libname tcc '/home/u41131808/TCC2 - GWZINBR';

data Korea_base;
	set tcc.korea_base_artigo;
run;

*data teste;
*  set Korea_base nobs=__nobs;
*  if _n_ le 100;
*run;

/*** TESTES ***/

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access Crowding,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZINB,METHOD=ADAPTIVE_BSQ,BANDWIDTH=aic, GLOBALMIN = no,DISTANCEKM=YES);
*6 mins;

%GWZINBR(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access Crowding,LONG=x,LAT=y,
OFFSET=ln_total,MODEL=zinb,METHOD=ADAPTIVE_BSQ,DISTANCEKM=YES, H=82);
*36 segs;

*PLANO B;
%GWZINBR(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access,LONG=x,LAT=y,
OFFSET=ln_total,MODEL=ZINB,METHOD=ADAPTIVE_BSQ,DISTANCEKM=YES, H=79);

/** ZINB **/

/* Teste 1: adaptive_bsq - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZINB,METHOD=ADAPTIVE_BSQ,BANDWIDTH=CV, GLOBALMIN = no,DISTANCEKM=YES);
*2 mins;

/* Teste 2: adaptive_bsq - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access, LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZINB,METHOD=ADAPTIVE_BSQ,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);
*2 mins;

/* Teste 3: fixed_g - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZINB,METHOD=fixed_g,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);
*4.5 mins;

/* Teste 4: fixed_g - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,XVARINF=Healthcare_access,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZINB,METHOD=fixed_g,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);
*5.5 mins;

/** ZIP **/

/* Teste 5: adaptive_bsq - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access, LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZIP,METHOD=ADAPTIVE_BSQ,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);
*22 segs;

/* Teste 6: adaptive_bsq - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZIP,METHOD=ADAPTIVE_BSQ,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);
*2 mins;

/* Teste 7: fixed_g - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access, LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZIP,METHOD=fixed_g,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);
*18 segs;

/* Teste 8: fixed_g - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZIP,METHOD=fixed_g,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);
*4 mins --> log grande;

/** NEGBIN **/

/* Teste 9: adaptive_bsq - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=NEGBIN,METHOD=ADAPTIVE_BSQ,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);
*19 segs;

/* Teste 10: adaptive_bsq - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=NEGBIN,METHOD=ADAPTIVE_BSQ,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);
*17 segs;

/* Teste 11: fixed_g - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=NEGBIN,METHOD=fixed_g,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);
*22 segs;

/* Teste 12: fixed_g - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=NEGBIN,METHOD=fixed_g,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);
*14 segs;

/** POISSON **/

/* Teste 13: adaptive_bsq - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=POISSON,METHOD=ADAPTIVE_BSQ,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);
*11 segs;

/* Teste 14: adaptive_bsq - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=POISSON,METHOD=ADAPTIVE_BSQ,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);
*9 segs;

/* Teste 15: fixed_g - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=POISSON,METHOD=fixed_g,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);
*6 segs;

/* Teste 16: fixed_g - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=POISSON,METHOD=fixed_g,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);
*13 segs;

*Alan;
%GWZINBR(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p Healthcare_access diff_sd
crowding Migration Health_behavior,XVARINF=Morbidity high_sch_p Healthcare_access diff_sd
crowding Migration Health_behavior,LAT=y,LONG=x,OFFSET=ln_total,METHOD=ADAPTIVE_BSQ,MODEL=ZINB,DISTANCEKM=YES, h=82);
*5 mins;

*Teste 1 do relatório;
%GWZINBR(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access,LONG=x,LAT=y,
OFFSET=ln_total,MODEL=ZINB,METHOD=fixed_g,DISTANCEKM=YES, H=226.73);
*40 segs;

*Teste 2 do relatório;
%GWZINBR(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access,LONG=x,LAT=y,
OFFSET=ln_total,MODEL=ZIP,METHOD=fixed_g,DISTANCEKM=YES, H=733.70);
*<1 seg;

*Teste 3 do relatório;
%GWZINBR(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, LONG=x,LAT=y,
OFFSET=ln_total,MODEL=NEGBIN,METHOD=fixed_g,DISTANCEKM=YES, H=189.74);
*2 segs;

*Teste 4 do relatório;
%GWZINBR(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,
OFFSET=ln_total,MODEL=POISSON,METHOD=fixed_g,DISTANCEKM=YES, H=733.70);
*1 seg;

*Teste 8;
%GWZINBR(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,
OFFSET=ln_total,MODEL=POISSON,METHOD=adaptive_bsq,DISTANCEKM=YES, H=48);
