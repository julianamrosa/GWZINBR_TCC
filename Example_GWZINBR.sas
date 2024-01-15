/* Geram os mesmos resultados do artigo */

libname tcc '/home/u57874055/sasuser.v94/TCC';

data Korea_base;
	set tcc.korea_base_artigo;
run;

%Golden(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p
Healthcare_access diff_sd crowding Migration Health_behavior,LONG=x,LAT=y,OUTPUT=band,
OFFSET=ln_total,MODEL=ZINB,METHOD=ADAPTIVE_BSQ,BANDWIDTH=CV, GLOBALMIN = NO,DISTANCEKM=YES);


%GWZINBR(DATA=Korea_base,YVAR=n_covid1,XVAR=Morbidity high_sch_p Healthcare_access diff_sd
crowding Migration Health_behavior,XVARINF=Morbidity high_sch_p Healthcare_access diff_sd
crowding Migration Health_behavior,LAT=y,LONG=x,OFFSET=ln_total,METHOD=ADAPTIVE_BSQ,MODEL=ZINB,DISTANCEKM=YES, h=82);