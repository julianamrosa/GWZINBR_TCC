/* Geram os mesmos resultados do artigo */

libname tcc '/home/u41131808/TCC2 - GWZINBR';

data Korea_base;
	set tcc.korea_base_artigo;
run;

*proc import file='/home/u41131808/TCC2 - GWZINBR/korea_base_zip.csv'
dbms=csv out=Korea_base_zip replace;*delimiter=',';
*run;

/*** TESTES ***/

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

/* Teste 9: adaptive_bsq - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access, LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZIP,METHOD=ADAPTIVE_BSQ,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);
*22 segs;

/* Teste 10: adaptive_bsq - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZIP,METHOD=ADAPTIVE_BSQ,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);
*2 mins;

/* Teste 11: fixed_g - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access, LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZIP,METHOD=fixed_g,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);
*18 segs;

/* Teste 12: fixed_g - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior, XVARINF=Healthcare_access,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZIP,METHOD=fixed_g,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);
*4 mins --> log grande;

/** NEGBIN **/

/* Teste 17: adaptive_bsq - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=NEGBIN,METHOD=ADAPTIVE_BSQ,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);
*19 segs;

/* Teste 18: adaptive_bsq - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=NEGBIN,METHOD=ADAPTIVE_BSQ,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);
*17 segs;

/* Teste 19: fixed_g - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=NEGBIN,METHOD=fixed_g,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);
*22 segs;

/* Teste 20: fixed_g - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=NEGBIN,METHOD=fixed_g,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);
*14 segs;

/** POISSON **/

/* Teste 25: adaptive_bsq - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=POISSON,METHOD=ADAPTIVE_BSQ,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);
*11 segs;

/* Teste 26: adaptive_bsq - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=POISSON,METHOD=ADAPTIVE_BSQ,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);
*9 segs;

/* Teste 27: fixed_g - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=POISSON,METHOD=fixed_g,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);
*6 segs;

/* Teste 28: fixed_g - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=POISSON,METHOD=fixed_g,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);
*13 segs;

%GWZINBR(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p Healthcare_access diff_sd
crowding Migration Health_behavior,XVARINF=Morbidity high_sch_p Healthcare_access diff_sd
crowding Migration Health_behavior,LAT=y,LONG=x,OFFSET=ln_total,METHOD=ADAPTIVE_BSQ,MODEL=ZINB,DISTANCEKM=YES, h=82);