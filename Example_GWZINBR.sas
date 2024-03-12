/* Geram os mesmos resultados do artigo */

libname tcc '/home/u41131808/TCC2 - GWZINBR';

data Korea_base;
	set tcc.korea_base_artigo;
run;

/*** TESTES ***/

/** ZINB **/

/* Teste 1: adaptive_bsq - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZINB,METHOD=ADAPTIVE_BSQ,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);
*25 segs;

/* Teste 2: adaptive_bsq - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZINB,METHOD=ADAPTIVE_BSQ,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);
*40 segs;

/* Teste 3: fixed_g - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZINB,METHOD=fixed_g,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);
*24 segs;

/* Teste 4: fixed_g - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZINB,METHOD=fixed_g,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);
*31 segs;

/* Teste 5: fixed_bsq - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZINB,METHOD=fixed_bsq,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);

/* Teste 6: fixed_bsq - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZINB,METHOD=fixed_bsq,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);

/* Teste 7: adaptiven - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZINB,METHOD=adaptiven,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);

/* Teste 8: adaptiven - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZINB,METHOD=adaptiven,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);

/** ZIP **/

/* Teste 9: adaptive_bsq - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZIP,METHOD=ADAPTIVE_BSQ,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);

/* Teste 10: adaptive_bsq - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZIP,METHOD=ADAPTIVE_BSQ,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);

/* Teste 11: fixed_g - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZIP,METHOD=fixed_g,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);

/* Teste 12: fixed_g - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZIP,METHOD=fixed_g,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);

/* Teste 13: fixed_bsq - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZIP,METHOD=fixed_bsq,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);

/* Teste 14: fixed_bsq - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZIP,METHOD=fixed_bsq,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);

/* Teste 15: adaptiven - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZIP,METHOD=adaptiven,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);

/* Teste 16: adaptiven - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZIP,METHOD=adaptiven,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);

/** NEGBIN **/

/* Teste 17: adaptive_bsq - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=NEGBIN,METHOD=ADAPTIVE_BSQ,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);

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
*13 segs;

/* Teste 21: fixed_bsq - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=NEGBIN,METHOD=fixed_bsq,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);

/* Teste 22: fixed_bsq - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=NEGBIN,METHOD=fixed_bsq,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);

/* Teste 23: adaptiven - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=NEGBIN,METHOD=adaptiven,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);

/* Teste 24: adaptiven - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=NEGBIN,METHOD=adaptiven,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);

/** POISSON **/

/* Teste 25: adaptive_bsq - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=NEGBIN,METHOD=ADAPTIVE_BSQ,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);

/* Teste 26: adaptive_bsq - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=POISSON,METHOD=ADAPTIVE_BSQ,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);

/* Teste 27: fixed_g - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=POISSON,METHOD=fixed_g,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);
*6 segs;

/* Teste 28: fixed_g - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=POISSON,METHOD=fixed_g,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);

/* Teste 29: fixed_bsq - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=POISSON,METHOD=fixed_bsq,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);

/* Teste 30: fixed_bsq - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=POISSON,METHOD=fixed_bsq,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);

/* Teste 31: adaptiven - cv */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=POISSON,METHOD=adaptiven,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);

/* Teste 32: adaptiven - aic */

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=POISSON,METHOD=adaptiven,BANDWIDTH=aic, GLOBALMIN = NO,DISTANCEKM=YES);

%GWZINBR(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p Healthcare_access diff_sd
crowding Migration Health_behavior,XVARINF=Morbidity high_sch_p Healthcare_access diff_sd
crowding Migration Health_behavior,LAT=y,LONG=x,OFFSET=ln_total,METHOD=ADAPTIVE_BSQ,MODEL=ZINB,DISTANCEKM=YES, h=82);