capture log close
clear all
set more off
set mem 500m
set maxvar 10000

*************************************************
*												*
*	THIS DO FILE MAKES THE INCOME DISTRIBUTIONS	*
*	AND ADDS THE MTRS FROM TAXSIM 9.0			*
*		USING THE 2007 PSID DATA				*
*	THIS VERSION Michael Best 23/5/11			*
*												*
*************************************************

/* See IRS' instructions for 1040 form at http://www.irs.gov/pub/irs-prior/i1040gi--2006.pdf 
 and a paper that has done this before at http://www.maxwell.syr.edu/uploadedFiles/cpr/publications/aging_studies/age12.pdf
 PSID questionnaire is at ftp://ftp.isr.umich.edu/pub/src/psid/questionnaires/q2007.pdf*/

 
* To run this on your own computer,
* change THIS FILE PATH to the folder you have the dataset in
*		   	   ||
*			  \||/
*			   \/
capture cd "/Users/aliciazhang/Dropbox (Princeton)/Equality of Opportunity/Simulations/1-Replication/data/new_dta/psid"
use "psid2007.dta"

log using "psidtaxsim.txt", text replace

*================================
*	GENERATE THE INPUT DATA		=
*================================

gen int year = 2006	 /*data in 2007 psid is on income in fy2006*/
gen byte state = ER36003

/* Labour Earnings */
	gen double pwages = ER40921 + ER40900 + ER40898
	replace pwages = ER40921 + ER40900 + (ER40898/2) if ER37298==4
	gen double swages = ER40930 + ER40933
	replace swages = ER40930 + ER40933 + (ER40898/2) if ER37298==4

/* number of dependents */
	gen byte depx = ER36020 /* make # dependents = # kids */
	/* Dependents from outside the hh */
		replace depx = depx+1 if ER37537==1 & ER37550==1 /* People who answered they had 1 dependent outside the hh and that that person relied on the respondent for >1/2 of income */
		replace depx = depx+ ER37549 if ER37537>1 & ER37537<=25 & ER37548==1 & ER37549<99 /* People who had >1 dependent ouside hh */
	/* Dependents in the same hh */
	
	/* Step 1: Tag the fus that could be dependents */
	gen monthsincomefromrels=0 																		/* count the number of months receivedd income from relatives*/
			forvalues i=1/9	{
			replace monthsincomefromrels=monthsincomefromrels+1 if ER3725`i'==1
			}
			replace monthsincomefromrels=monthsincomefromrels+1 if ER37260==1
			replace monthsincomefromrels=monthsincomefromrels+1 if ER37261==1
			replace monthsincomefromrels=monthsincomefromrels+1 if ER37262==1
			
			gen incfromrels = 0 																			/* add up the total income for the year from relatives */
			replace incfromrels=ER37248 if ER37249==6 & ER37248<999998
			replace incfromrels=ER37248*monthsincomefromrels if ER37249==5 & ER37248<999998
			replace incfromrels=ER37248*2*monthsincomefromrels if ER37249==4 & ER37248<999998
			replace incfromrels=ER37248*4*monthsincomefromrels if ER37249==3 & ER37248<999998
			
			gen byte fudep=0 																				/* tag the fus that are dependent on a relative and would qualify as a dependent */
			replace fudep=1 if ER41027<=3300 & incfromrels/ER41027 >= 0.5 & ER36028==5 & ER36023!=1 		/* income<3300, >.5 income from family and pays rent, not married */
			replace fudep=1 if ER41027<=3300 & ER36028!=5 & ER36023!=1 										/* income<3300, doesn't pay rent, not married, */
			
			drop monthsincomefromrels incfromrels															/* Tidy Up */
			
		preserve 																							/*Create 4 datasets to merge with in which fus that would qualify as dependents are tagged */
			keep ER36002 fudep 																				/* Keep interview number and tag for dependency */
			
			rename ER36002 ER41047 																			/* rename the interview number to be the interview number of fu1 that lives in the same hh */
			rename fudep fudep1
			save "otherfu1.dta", replace
			
			rename ER41047 ER41050
			rename fudep1 fudep2
			save "otherfu2.dta", replace
		
			rename ER41050 ER41053
			rename fudep2 fudep3
			save "otherfu3.dta", replace
			
			rename ER41053 ER41056
			rename fudep3 fudep4
			save "otherfu4.dta", replace
			
		restore
	
	/* Step 2: Find the fus that are supporting another dependent fu */	
		merge m:1 ER41047 using "otherfu1.dta"											/* merge on the interview number of the fu living in the hh */
		drop if _merge==2
		drop _merge
		merge m:1 ER41050 using "otherfu2.dta"
		drop if _merge==2
		drop _merge
		merge m:1 ER41053 using "otherfu3.dta"
		drop if _merge==2
		drop _merge
		merge m:1 ER41056 using "otherfu4.dta"
		drop if _merge==2
		drop _merge
		
		forvalues i=1/4	{
			erase "otherfu`i'.dta"
		}
		
		replace depx = depx+1 if ER41048>0 & ER41048<7 & fudep==0 & fudep1==1			/*fu1 is related, fu is not dependent, fu1 is dependent */
		replace depx = depx+1 if ER41051>0 & ER41051<7 & fudep==0 & fudep2==1
		replace depx = depx+1 if ER41054>0 & ER41054<7 & fudep==0 & fudep3==1
		replace depx = depx+1 if ER41057>0 & ER41057<7 & fudep==0 & fudep4==1
	
		drop fudep*																		/* Tidy Up */
	
* NEW: Some edits here
/* Marital status */
	gen byte mstat =1
	replace mstat=2 if ER36023==1 														/* Assume all married file jointly (changed from 3 to 8) */
	replace mstat=8 if ER36023!=1 & depx>0
* Unmarried can't have spousal wages
	replace swages = 0 if mstat != 2

/* Number of 65+ taxpayers */
	gen byte agex=0
	replace agex=1 if ER36017>=65 | ER36019>=65
	replace agex=2 if ER36017>=65 & ER36019>=65

/* Dividend Income */
	
	gen double dividends = ER40925 + ER40937
	
/* Interest and other property income */

	gen double prophead = ER40901 + ER40922 + ER40926 + ER40928 								/* Business Asset + Interest + Rent + Trust Fund/Royalty Income */
	gen double propwife = ER40931 + ER40935 + ER40939 + ER40941									/* Business Asset + Interest + Rent + Trust Fund/Royalty Income */

	gen double otherprop = prophead + propwife
	replace otherprop = otherprop + ER40964 - ER37547 if ER37547<9999998						/* Add in net alimony */
	
	drop prophead propwife

/* Pensions */

	gen double pensions = ER40952 + ER40978

/* Social Security Income */

	gen double gssi = ER41021 + ER41023															/* Note ignoring ER41025 = other FU member SS income */
	
/* Transfer Income */

/* Rent (Note 1 respondent pays daily rent. Ignore) */

	gen double rentpaid = 0
	replace rentpaid = ER36065 if ER36066==6													/* Rent paid for those who reported a yearly amount */
	replace rentpaid = ER36065*12 if ER36066==5													/* Rent paid for those who reported a monthly amount */
	replace rentpaid = ER36065*26 if ER36066==4													/* Rent paid for those who reported a two-weekly amount */
	replace rentpaid = ER36065*52 if ER36066==3													/* Rent paid for those who reported a weekly amount */
	
/* Property Taxes */

	gen double proptax = 0
	replace proptax = ER36036 if ER36036<99998

/* Child Care */

	gen double childcare=0
	replace childcare = ER36633 if ER36633<999998
	
/* Unemployment Compensation */

	gen double ui = ER40958 + ER40980
	
/* Number of Children */

	gen byte depchild = ER36020

/* Mortgage interest */

	/* First Mortgage */
		
		gen firstmortgagerate = 0
		replace firstmortgagerate = ER36046 + (ER36047/1000) if ER36046<98 & ER36047<998		/* Generate the rate paid on the mortgage */
		gen double firstmortgage = 0
		replace firstmortgage = (firstmortgagerate/100)*ER36042 if ER36042<9999998 & ER36045>0	/* Make the mortgage interest payments the rate times the principal if payments are non-zero */

	/* Second Mortgage */
		
		gen secondmortgagerate = 0
		replace secondmortgagerate = ER36058 + (ER36059/1000) if ER36058<98 & ER36059<998
		gen double secondmortgage = 0
		replace secondmortgage = (secondmortgagerate/100)*ER36054 if ER36054<9999998 & ER36056>0

	gen double mortgage = firstmortgage + secondmortgage
	
	drop firstmortgage secondmortgage

save "psid_clean.dta", replace

* NEW: Some edits here to make compatible with TAXSIM32
	drop ER*
	drop firstmortgagerate
	drop secondmortgagerate
	gen dep17 = depx
	gen dep18 = depx
	drop depchild
	
save "psid_clean_taxsim.dta", replace
	
	
*====================
*	SEND TO TAXSIM	=	
*====================

/* Get the rates for the head */

taxsim32

use "taxsim_out.dta", clear

capture erase "taxsim_out_head.dta"
save "taxsim_out_head.dta", replace
capture erase "taxsim_out.dta"

/* Get the rates for the "wife" */

use "psid_clean_taxsim.dta", clear

taxsim32, secondary

capture erase "taxsim_out_wife.dta"
use "taxsim_out.dta", clear

rename fiitax fiitaxwife
rename siitax siitaxwife
rename fica ficawife
rename frate fratewife
rename srate sratewife
rename ficar ficarwife

keep taxsimid year state fratewife sratewife ficarwife

save "taxsim_out_wife.dta", replace
	
*============================================================
*	MERGE WITH THE TAXSIM DATA GENERATE THE INPUT DATA		=
*============================================================

use "psid_clean.dta", clear
	
gen taxsimid = _n

merge 1:1 taxsimid using "taxsim_out_head.dta", nogen	
merge 1:1 taxsimid using "taxsim_out_wife.dta", nogen

save "psidtaxsimmerge.dta", replace

/* Tidy Up */

capture erase "taxsim_out.dta"
capture erase "taxsim_out_head.dta"
capture erase "taxsim_out_wife.dta"
capture erase "psid_clean.dta"

log close
