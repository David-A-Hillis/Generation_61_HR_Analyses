/*The purpose of this file is to perform the mivque0 analysis with PROC MIXED using the full genome data
This file originates from "permutation3.sas" due to the use of the macro.
However, much of the code is lifted from "Single_locus_multi_model.sas" since this will attempt to perform analyses for all models simultaneously.
Rather than using one long file I will try to use the two files with similar format to the R code.*/

%let _timer_start = %sysfunc(datetime()); /* Start timer, for total runtime of program */

/* The following 4 lines change options in SAS to not save unnecessary data*/
options nonotes; 		/* Turns off blue notes in Log so that it doesn't need to be cleared periodically*/
ods noresults;			/* SAS Support - no results accumulating in the results tree*/
ods _all_ close;		/* SAS Support - no output*/
ods graphics off;		/* SAS Support - no ods graphics*/

%let dir_in=d:\Users\User\Desktop\David\SAS\Allelic_NA;	/*Set working directory for files being imported*/
%let dir_out=d:\Users\User\Desktop\David\SAS\Results\Chr15;	/*Set working directory for files being exported*/
filename ina "&dir_in\MIVQUE_haplotype_metadata.csv";			/*Has mouse data*/
filename inb "&dir_in\chr15_noNA_allelic.csv"; 			/*Has Genotypes and marker columns*/
filename inc "&dir_in\Tests.csv"; 								/*Used to make empty "tests" tables*/
filename ind "&dir_in\AIC.csv"; 								/*Used to make empty "log" tables*/
filename ine "&dir_in\CovParm.csv"; 							/*Used to make empty "parms" tables*/
proc import datafile=ina out=meta dbms=csv replace;		/* Loads mouse data in SAS as "meta"*/
proc import datafile=inb out=geno dbms=csv replace;		/* Loads genetic data in SAS as "geno"*/
proc import datafile=inc out=test dbms=csv replace;		/* Loads test results in SAS as "test"*/
proc import datafile=ind out=aic dbms=csv replace;		/* Loads test results in SAS as "aic"*/
proc import datafile=ine out=cparms dbms=csv replace;	/* Loads test results in SAS as "cparms"*/
run;
/* That should take care of loading and initial prep of data*/

/* Formating Stuff*/
data test1; length Effect $9; set test; format DenDF BEST4. ProbF PVALUE6.;
data cparms1; length CovParm $15; set cparms; format Estimate D8.;
data aic1; length Descr $25; set aic;
run;

proc sql; 		/* These are empty tables to be populated in the macro*/
  create table full_tests like test1;
  create table nogroup_tests like test1;
  create table nogroupLine_tests like test1;
  create table nogroupMouse_tests like test1;
  create table noLineWgroup_tests like test1;
  create table noLine_tests like test1;
  create table full_log like aic1;
  create table nogroup_log like aic1;
  create table nogroupLine_log like aic1;
  create table nogroupMouse_log like aic1;
  create table noLineWgroup_log like aic1;
  create table noLine_log like aic1;
  create table full_parms like cparms1;
  create table proto_nogroup_parms like cparms1;
  create table nogroupLine_parms like cparms1;
  create table nogroupMouse_parms like cparms1;
  create table noLineWgroup_parms like cparms1;
  create table proto_noLine_parms like cparms1;
quit;

/* Remove unneeded Group Variable from Parms tables that don't need it*/
data nogroup_parms (drop=Group); set proto_nogroup_parms; run;
data noLine_parms (drop=Group); set proto_noLine_parms; run;

data geno1 (drop=X__CHROM POS Chromosome); set geno; run; 	/* Drops the first 3 columns*/
proc transpose data=geno1 out=geno2; run;		/* Transposes the newly merge dataset*/

%let n=0; run;

sasfile geno1 open;
sasfile meta open;

/* This is the macro which is necessary for running procedures in do loops*/ 
/* Still work in progress*/
%macro allModels; 											/* "allModels" is the name of the macro*/
%do k = 1 %to 282844;	         	          		 		/* This should always start on 1 and end on however many loci are in the chromosome*/
   %let n=%eval(&n+1);

   data allele; 											/* This identifies the correct locus and drop marker*/
		set geno1; 
		where marker=&n;
   data alleleNoM (drop=marker);
   		set allele;
   proc transpose data=alleleNoM out=allele1;				/* This transposes the one locus to be merged with the mouse data*/
   data locus; merge meta allele1;

   proc mixed data=locus method=mivque0; 					/* Full model*/
   		class pop sub mouse;								/* Identifies population and subpopulation in meta data and is shuffled with mouseID*/
   		model COL1=pop/solution;							/* Describes model*/
   		random sub(pop) /group=pop;							/* Sets random effects, nesting, and groups for analysis parameters*/
		random mouse(sub pop) /group = pop;					/* Sets random effects, nesting, and groups for analysis parameters*/
   		ods output CovParms=parms Tests3=type3 FitStatistics=log;	/* Identifies output table names*/
   proc datasets;
   		append base=full_tests data=type3 force;			/* Appends the p-values*/
		append base=full_log data=log force;				/* Appends the AIC*/
		append base=full_parms data=parms force;			/* Appends the CovParms*/


   proc mixed data=locus method=mivque0; 					/* nogroup model*/
   		class pop sub mouse;								/* Identifies population and subpopulation in meta data and is shuffled with mouseID*/
   		model COL1=pop/solution;							/* Describes model*/
   		random sub(pop);									/* Sets random effects, nesting, and groups for analysis parameters*/
		random mouse(sub pop);								/* Sets random effects, nesting, and groups for analysis parameters*/
   		ods output CovParms=parms Tests3=type3 FitStatistics=log;	/* Identifies output table names*/
   proc datasets;
   		append base=nogroup_tests data=type3 force;			/* Appends the p-values*/
		append base=nogroup_log data=log force;				/* Appends the AIC*/
		append base=nogroup_parms data=parms force;			/* Appends the CovParms*/


   proc mixed data=locus method=mivque0; 					/* nogroupLine model*/
   		class pop sub mouse;								/* Identifies population and subpopulation in meta data and is shuffled with mouseID*/
   		model COL1=pop/solution;							/* Describes model*/
   		random sub(pop);									/* Sets random effects, nesting, and groups for analysis parameters*/
		random mouse(sub pop) /group = pop;					/* Sets random effects, nesting, and groups for analysis parameters*/
   		ods output CovParms=parms Tests3=type3 FitStatistics=log;	/* Identifies output table names*/
   proc datasets;
   		append base=nogroupLine_tests data=type3 force;		/* Appends the p-values*/
		append base=nogroupLine_log data=log force;			/* Appends the AIC*/
		append base=nogroupLine_parms data=parms force;		/* Appends the CovParms*/


   proc mixed data=locus method=mivque0; 					/* nogroupMouse model*/
   		class pop sub mouse;								/* Identifies population and subpopulation in meta data and is shuffled with mouseID*/
   		model COL1=pop/solution;							/* Describes model*/
   		random sub(pop) /group = pop;						/* Sets random effects, nesting, and groups for analysis parameters*/
		random mouse(sub pop);								/* Sets random effects, nesting, and groups for analysis parameters*/
   		ods output CovParms=parms Tests3=type3 FitStatistics=log;	/* Identifies output table names*/
   proc datasets;
   		append base=nogroupMouse_tests data=type3 force;	/* Appends the p-values*/
		append base=nogroupMouse_log data=log force;		/* Appends the AIC*/
		append base=nogroupMouse_parms data=parms force;	/* Appends the CovParms*/


   proc mixed data=locus method=mivque0; 					/* noLineWgroup model*/
   		class pop sub mouse;								/* Identifies population and subpopulation in meta data and is shuffled with mouseID*/
   		model COL1=pop/solution;							/* Describes model*/
   		*random sub(pop) /group = pop;						/* Sets random effects, nesting, and groups for analysis parameters*/
		random mouse(sub pop) /group = pop;					/* Sets random effects, nesting, and groups for analysis parameters*/
   		ods output CovParms=parms Tests3=type3 FitStatistics=log;	/* Identifies output table names*/
   proc datasets;
   		append base=noLineWgroup_tests data=type3 force;	/* Appends the p-values*/
		append base=noLineWgroup_log data=log force;		/* Appends the AIC*/
		append base=noLineWgroup_parms data=parms force;	/* Appends the CovParms*/


   proc mixed data=locus method=mivque0; 					/* noLine model*/
   		class pop sub mouse;								/* Identifies population and subpopulation in meta data and is shuffled with mouseID*/
   		model COL1=pop/solution;							/* Describes model*/
   		*random sub(pop) /group = pop;						/* Sets random effects, nesting, and groups for analysis parameters*/
		random mouse(sub pop);								/* Sets random effects, nesting, and groups for analysis parameters*/
   		ods output CovParms=parms Tests3=type3 FitStatistics=log;	/* Identifies output table names*/
   proc datasets;
   		append base=noLine_tests data=type3 force;			/* Appends the p-values*/
		append base=noLine_log data=log force;				/* Appends the AIC*/
		append base=noLine_parms data=parms force;			/* Appends the CovParms*/

   data _null_; 
   		y = &n; 
		z=mod(y,20);
		if z=0 then do;
			put y;
		end;
%end;
%mend allModels;											/* mend ends the macro, differing notes on whether the name needs to be included*/

/* This actually runs the analysis*/
%allModels;
run;

/* Organize parms has been removed, it would require modified code for each model and we are currenly not using it
	These files can be modified in R after running SAS, if needed*/

/* This changes all p-values to 10 decimal places until stated otherwise*/
data full_tests1; set full_tests; format ProbF pvalue12.10; run;
data nogroup_tests1; set nogroup_tests; format ProbF pvalue12.10; run;
data nogroupLine_tests1; set nogroupLine_tests; format ProbF pvalue12.10; run;
data nogroupMouse_tests1; set nogroupMouse_tests; format ProbF pvalue12.10; run;
data noLineWgroup_tests1; set noLineWgroup_tests; format ProbF pvalue12.10; run;
data noLine_tests1; set noLine_tests; format ProbF pvalue12.10; run;

/* Organize log*/
data log; set full_log; format Value BEST12.; run; /* Changes sig figs to 4 digits before decimal and 8 after*/
proc sort data=log; by Descr; run;
data log2 (rename=(Value=Log2) drop=descr); set log; if _n_ gt 282844  then delete; run;
data logAIC (rename=(Value=AIC) drop=marker descr); set log; if _n_ gt 565688  then delete; if _n_ lt 282845  then delete; run;
data logAICC (rename=(Value=AICC) drop=marker descr); set log; if _n_ gt 848532 then delete; if _n_ lt 565689  then delete; run;
data logBIC (rename=(Value=BIC) drop=marker descr); set log; if _n_ lt 848533 then delete; run;
data full_log_final; merge log2 logAIC logAICC logBIC; run;

data log; set nogroup_log; format Value BEST12.; run; /* Changes sig figs to 4 digits before decimal and 8 after*/
proc sort data=log; by Descr; run;
data log2 (rename=(Value=Log2) drop=descr); set log; if _n_ gt 282844  then delete; run;
data logAIC (rename=(Value=AIC) drop=marker descr); set log; if _n_ gt 565688  then delete; if _n_ lt 282845  then delete; run;
data logAICC (rename=(Value=AICC) drop=marker descr); set log; if _n_ gt 848532 then delete; if _n_ lt 565689  then delete; run;
data logBIC (rename=(Value=BIC) drop=marker descr); set log; if _n_ lt 848533 then delete; run;
data nogroup_log_final; merge log2 logAIC logAICC logBIC; run;

data log; set nogroupLine_log; format Value BEST12. run; /* Changes sig figs to 4 digits before decimal and 8 after*/
proc sort data=log; by Descr; run;
data log2 (rename=(Value=Log2) drop=descr); set log; if _n_ gt 282844  then delete; run;
data logAIC (rename=(Value=AIC) drop=marker descr); set log; if _n_ gt 565688  then delete; if _n_ lt 282845  then delete; run;
data logAICC (rename=(Value=AICC) drop=marker descr); set log; if _n_ gt 848532 then delete; if _n_ lt 565689  then delete; run;
data logBIC (rename=(Value=BIC) drop=marker descr); set log; if _n_ lt 848533 then delete; run;
data nogroupLine_log_final; merge log2 logAIC logAICC logBIC; run;

data log; set nogroupMouse_log; format Value BEST12.; run; /* Changes sig figs to 4 digits before decimal and 8 after*/
proc sort data=log; by Descr; run;
data log2 (rename=(Value=Log2) drop=descr); set log; if _n_ gt 282844  then delete; run;
data logAIC (rename=(Value=AIC) drop=marker descr); set log; if _n_ gt 565688  then delete; if _n_ lt 282845  then delete; run;
data logAICC (rename=(Value=AICC) drop=marker descr); set log; if _n_ gt 848532 then delete; if _n_ lt 565689  then delete; run;
data logBIC (rename=(Value=BIC) drop=marker descr); set log; if _n_ lt 848533 then delete; run;
data nogroupMouse_log_final; merge log2 logAIC logAICC logBIC; run;

data log; set noLineWgroup_log; format Value BEST12.; run; /* Changes sig figs to 4 digits before decimal and 8 after*/
proc sort data=log; by Descr; run;
data log2 (rename=(Value=Log2) drop=descr); set log; if _n_ gt 282844  then delete; run;
data logAIC (rename=(Value=AIC) drop=marker descr); set log; if _n_ gt 565688  then delete; if _n_ lt 282845  then delete; run;
data logAICC (rename=(Value=AICC) drop=marker descr); set log; if _n_ gt 848532 then delete; if _n_ lt 565689  then delete; run;
data logBIC (rename=(Value=BIC) drop=marker descr); set log; if _n_ lt 848533 then delete; run;
data noLineWgroup_log_final; merge log2 logAIC logAICC logBIC; run;

data log; set noLine_log; format Value BEST12.; run; /* Changes sig figs to 4 digits before decimal and 8 after*/
proc sort data=log; by Descr; run;
data log2 (rename=(Value=Log2) drop=descr); set log; if _n_ gt 282844  then delete; run;
data logAIC (rename=(Value=AIC) drop=marker descr); set log; if _n_ gt 565688  then delete; if _n_ lt 282845  then delete; run;
data logAICC (rename=(Value=AICC) drop=marker descr); set log; if _n_ gt 848532 then delete; if _n_ lt 565689  then delete; run;
data logBIC (rename=(Value=BIC) drop=marker descr); set log; if _n_ lt 848533 then delete; run;
data noLine_log_final; merge log2 logAIC logAICC logBIC; run;


filename outa "&dir_out\chr15_NA_allelic_full_tests.csv";
filename outb "&dir_out\chr15_NA_allelic_nogroup_tests.csv";
filename outc "&dir_out\chr15_NA_allelic_nogroupLine_tests.csv";
filename outd "&dir_out\chr15_NA_allelic_nogroupMouse_tests.csv";
filename oute "&dir_out\chr15_NA_allelic_noLineWgroup_tests.csv";
filename outf "&dir_out\chr15_NA_allelic_noLine_tests.csv";
filename outg "&dir_out\chr15_NA_allelic_full_log.csv";
filename outh "&dir_out\chr15_NA_allelic_nogroup_log.csv";
filename outi "&dir_out\chr15_NA_allelic_nogroupLine_log.csv";
filename outj "&dir_out\chr15_NA_allelic_nogroupMouse_log.csv";
filename outk "&dir_out\chr15_NA_allelic_noLineWgroup_log.csv";
filename outl "&dir_out\chr15_NA_allelic_noLine_log.csv";
filename outm "&dir_out\chr15_NA_allelic_full_parms.csv";
filename outn "&dir_out\chr15_NA_allelic_nogroup_parms.csv";
filename outo "&dir_out\chr15_NA_allelic_nogroupLine_parms.csv";
filename outp "&dir_out\chr15_NA_allelic_nogroupMouse_parms.csv";
filename outq "&dir_out\chr15_NA_allelic_noLineWgroup_parms.csv";
filename outr "&dir_out\chr15_NA_allelic_noLine_parms.csv";

proc export data=full_tests1 outfile=outa dbms=csv replace;
proc export data=nogroup_tests1 outfile=outb dbms=csv replace;
proc export data=nogroupLine_tests1 outfile=outc dbms=csv replace;
proc export data=nogroupMouse_tests1 outfile=outd dbms=csv replace;
proc export data=noLineWgroup_tests1 outfile=oute dbms=csv replace;
proc export data=noLine_tests1 outfile=outf dbms=csv replace;
proc export data=full_log_final outfile=outg dbms=csv replace;
proc export data=nogroup_log_final outfile=outh dbms=csv replace;
proc export data=nogroupLine_log_final outfile=outi dbms=csv replace;
proc export data=nogroupMouse_log_final outfile=outj dbms=csv replace;
proc export data=noLineWgroup_log_final outfile=outk dbms=csv replace;
proc export data=noLine_log_final outfile=outl dbms=csv replace;
proc export data=full_parms outfile=outm dbms=csv replace;
proc export data=nogroup_parms outfile=outn dbms=csv replace;
proc export data=nogroupLine_parms outfile=outo dbms=csv replace;
proc export data=nogroupMouse_parms outfile=outp dbms=csv replace;
proc export data=noLineWgroup_parms outfile=outq dbms=csv replace;
proc export data=noLine_parms outfile=outr dbms=csv replace;
Run;


/* Stop timer*/
data _null_;							
  dur = datetime() - &_timer_start;
  put 30*'-' / ' TOTAL DURATION:' dur time13.2 / 30*'-';
run;

