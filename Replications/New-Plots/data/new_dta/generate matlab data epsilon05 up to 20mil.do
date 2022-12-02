capture log close
clear all
set more off
set mem 500m
set maxvar 10000

*************************************************
*												*
*	THIS DO FILE MAKES THE NECESSARY DATA 		*
*		FOR THE MATLAB SIMULATIONS				*
*		 USING THE 2007 PSID DATA				*
*			THIS VERSION 13/6/11				*
*												*
*************************************************

capture cd "/Users/aliciazhang/Dropbox (Princeton)/Equality of Opportunity/Simulations/1-Replication/data/new_dta/psid"

use "psidtaxsimmerge.dta", clear

*capture cd "C:\Users\Michael\Documents\Z Drive\Research\Tax career concerns\Writing"

capture mkdir "Pictures"
capture mkdir "Pictures\intermed"


*========================
*	GET THE DATA READY	=
*========================


/* Drop everything except the earnings variables, tax variables age and the weight */
keep ER41069 ER36017 ER36019 ER37298 ER40921 ER40900 ER40898 ER40933 ER40930 frate srate ficar fratewife sratewife ficarwife
rename ER36017 agehead
rename ER36019 agewife
rename ER41069 weight


/* Make variables that are total earnings */
gen double earnhead = ER40921 + ER40900 + ER40898
replace earnhead = ER40921 + ER40900 + (ER40898/2) if ER37298==4
gen double earnwife = ER40930 + ER40933
replace earnwife = ER40930 + ER40933 + (ER40898/2) if ER37298==4
drop ER*

/* reshape to get a dataset of individuals not couples */
gen id = _n
rename frate fratehead
rename srate sratehead
rename ficar ficarhead
reshape long age earn frate srate ficar, i(id) j(headorwife) string

/* Generate the total marginal tax rate */
gen mtrate = frate + srate + ficar

/* Drop people with non-positive earnings
Note this is FISHY since there are 252 observations with negative earnings! */
drop if earn<=0
drop if age==999

/* Tag the young and the old */

*Note, those at the median age are young
summ age, detail
gen young = (age<=`r(p50)')

/* Generate the grid of points to evaluate at (equal log spacing of .0106) starting from e^5~$148 to e^16.8125~$20,025,159*/

gen double z = exp( 5 + ( (_n - 1)*.0105) ) if _n<1127												/* Grid for earnings */

/* Do Kernel Density Estimates of the income density and distribution */

	/* The overall distribution */
		akdensity earn [aw=weight], bwidth(2000) gen(redundantz1 hz) at(z) nograph cdf(Hz) stdbands(2)
		
	/* The young distribution */
		akdensity earn [aw=weight] if young==1, bwidth(2000) gen(redundantz2 hzy) at(z) nograph cdf(Hzy) stdbands(2)	
		
	/* The old distribution */
		akdensity earn [aw=weight] if young==0, bwidth(2000) gen(redundantz3 hzo) at(z) nograph cdf(Hzo) stdbands(2)
		
/* Find z^m/z to fit the pareto tail */

	/* The overall distribution */
	
		gen double zhz = z * hz	if _n<1002													/* z * h(z) */
		integ zhz z if _n<1002, gen(intz) replace											/* \int z*h(z)dz */
		local ztot = `r(integral)'															/* \int_0^\infty z*h(z)dz */
		gen double zm = (`ztot'-intz)/(1-Hz) if _n<1002										/* E[z|z>zbar] */
		gen double zmz = zm/z if _n<1002													/* z^m / z */
		gen double impliedalphaz = zmz/(zmz-1) if _n<1002									/* The implied alpha */
		
		reg zmz [aw=hz] if z>100000 & z<150000												/* Use the distribution between 100000 & 150000 to estimate the pareto parameter */
		matrix A=e(b)
		local intermedz = A[1,1]															/* Store the estimated zm/z as `intermedz' */
		local alphaz = `intermedz'/(`intermedz'-1)											/* Store the pareto parameter as `alphaz' */
		
	/* The young distribution */
	
		gen double zyhzy = z * hzy	if _n<1002												/* zy * h(zy) */
		integ zyhzy z if _n<1002, gen(intzy) replace										/* \int zy*h(zy)dzy */
		local zytot = `r(integral)'															/* \int_0^\infty zy*h(zy)dzy */
		gen double zym = (`zytot'-intzy)/(1-Hzy) if _n<1002									/* E[zy|zy>zbar] */
		gen double zymzy = zym/z if _n<1002													/* zy^m / zy */
		gen double impliedalphazy = zymzy/(zymzy-1) if _n<1002								/* The implied alpha */
		
		reg zymzy [aw=hzy] if z>100000 & z<150000													/* Use the distribution between 100000 & 150000 to estimate the pareto parameter */
		matrix B=e(b)
		local intermedzy = B[1,1]															/* Store the estimated zm/z as `intermedzy' */
		local alphazy = `intermedzy'/(`intermedzy'-1)										/* Store the pareto parameter as `alphazy' */

	/* The old Distribution */
	
		gen double zohzo = z * hzo	if _n<1002												/* zo * h(zo) */
		integ zohzo z if _n<1002, gen(intzo) replace										/* \int zo*h(zo)dzo */
		local zotot = `r(integral)'															/* \int_0^\infty zo*h(zo)dzo */
		gen double zom = (`zotot'-intzo)/(1-Hzo) if _n<1002									/* E[zo|zo>zbar] */
		gen double zomzo = zom/z if _n<1002													/* zo^m / zo */
		gen double impliedalphazo = zomzo/(zomzo-1) if _n<1002								/* The implied alpha */
		
		reg zomzo [aw=hzo] if z>100000 & z<150000											/* Use the distribution between 100000 & 150000 to estimate the pareto parameter */
		matrix C=e(b)
		local intermedzo = C[1,1]															/* Store the estimated zm/z as `intermedzo' */
		local alphazo = `intermedzo'/(`intermedzo'-1)										/* Store the pareto parameter as `alphazo' */

/* Put in the Pareto tails */
																							/* Note _n==677 has z==150,578.53 */
	/* the overall distribution */
		local Hz150k = (1 - Hz[660])														/* Grab the proportion of people above 150K */
		*local alphaz = 2.156																/* This is the alpha found using those between 100K and 150K */
		local zminalpha = (150000^`alphaz')*`Hz150k'										/* This is the zmin parameter (to the power of alpha) for the pareto distribution that scales it to have the same proportion of people above 150K given the alpha parameter */
		replace hz = `alphaz'*`zminalpha'*(z^(-(`alphaz'+1))) if z> 150000 & z<21000000		/* The pareto density */
		replace Hz = 1-(`zminalpha'*(z^-`alphaz')) if z>150242 & z<21000000					/* The pareto distribution */
							
	/* the young distribution */
		local Hzy150k = (1 - Hzy[660])														/* Grab the proportion of people above 150K */
		*local alphazy = 2.464																/* This is the alpha found using those between 100K and 150K */
		local zyminalpha = (150000^`alphazy')*`Hzy150k'										/* This is the zmin parameter (to the power of alpha) for the pareto distribution that scales it to have the same proportion of people above 150K given the alpha parameter */
		replace hzy = `alphazy'*`zyminalpha'*(z^(-(`alphazy'+1))) if z> 150000 & z<21000000	/* The pareto density */
		replace Hzy = 1-(`zyminalpha'*(z^-`alphazy')) if z>150000 & z<21000000				/* The pareto distribution */

	/* the old distribution */
		local Hzo150k = (1 - Hzo[660])														/* Grab the proportion of people above 150K */
		*local alphazo = 2.076																/* This is the alpha found using those between 100K and 150K */
		local zominalpha = (150000^`alphazo')*`Hzo150k'										/* This is the zmin parameter (to the power of alpha) for the pareto distribution that scales it to have the same proportion of people above 150K given the alpha parameter */
		replace hzo = `alphazo'*`zominalpha'*(z^(-(`alphazo'+1))) if z> 150000 & z<21000000	/* The pareto density */
		replace Hzo = 1-(`zominalpha'*(z^-`alphazo')) if z>150000 & z<21000000				/* The pareto distribution */

/* Do a local Polynomial regression for the MTR */

	lpoly mtrate earn [aw=weight], bwidth(3000) degree(1) generate(redundatz4 mtr2) at(z)	/* (linear) polynomial local regression of the MTR on earnings */
	ipolate mtr2 z, gen(mtr)																/* Linearly interpolate the missing values */
	replace mtr = mtr2[841] if _n>841 & _n<1127    											/* There's not enough data far out, so take the tax rate out there as the one it seems to converge to */


/*****************************************************/
/***************	THE PARAMETERS 	******************/
/*****************************************************/

global epsilon = 0.5

/*****************************************************/
/*****************************************************/	


*========================================
*	Back out the wage rates of the old	=
*========================================

/* The wages of the old are ( zo / (1-tau)^epsilon )^(1/1+epsilon) */

capture drop wageold
gen double wageold = (earn / ((1-(mtrate/100))^${epsilon}))^(1/(1+${epsilon}))

/* The density and distribution of omega */

gen double omegagrid = exp( (_n - 1)*.0108 ) if _n<1127												/* Grid for ability (equal log spacing of .0112) starting from 1 to e^12.152~189,472.7 */

akdensity wageold [aw=weight] if young==0, bwidth(20) gen(redundantz5 jomegao) at(omegagrid) nograph cdf(Jomegao) stdbands(2)

/* Generate omegam/omega to give the omega distribution a pareto tail also */

drop in 1127/l

gen double omegajomega = omegagrid * jomegao														/* omega * j(omega) */
integ omegajomega omegagrid, gen(intomega) replace													/* \int omega*j(omega)domega */
local omegatot = `r(integral)'																	/* \int_0^\infty omega*j(omega)domega */
gen double omegam = (`omegatot'-intomega)/(1-Jomegao)											/* E[omega|omega>omegagrid] */
gen double omegamomega = omegam/omegagrid															/* omega^m / omega */
gen double impliedalphaomega = omegamomega/(omegamomega-1)										/* the implied alpha */

reg omegamomega [aw=jomegao] if Jomegao>0.87831 & omegagrid<4000									/* Use the highest 10% of the distribution below 4000 to estimate the pareto parameter */
matrix D=e(b)																					
local intermedomega = D[1,1]																	/* Store the estimated omegam/omega as `intermedomega' */
local alphaomega = `intermedomega'/(`intermedomega'-1)											/* Store the pareto parameter as `alphaomega' */

/* Put in the Pareto Tails */

local Jomega4k = (1-Jomegao[769])																/* The proportion of people with omega>4000 Note. _n=769 has omegagrid==4001.41 */
local omegamin = (4000^`alphaomega')*`Jomega4k'													/* The lower bound parameter of the pareto distribution to make the tail integrate to the proportion of people with omega>4000 */
replace jomegao = `alphaomega'*`omegamin'*(omegagrid^(-(`alphaomega'+1))) if omegagrid>4000				/* The pareto density */
replace Jomegao = 1-(`omegamin'*(omegagrid^-`alphaomega')) if omegagrid>4000							/* The pareto distribution */

/* Export the stuff to matlab */

	/* put the global and pareto parameters in the variable "parameters" */
	gen double parameters=.
	replace parameters = ${epsilon} 				in 1
	replace parameters = `zyminalpha'^(1/`alphazy') in 2
	replace parameters = `alphazy' 					in 3
	replace parameters = `zominalpha'^(1/`alphazo') in 4
	replace parameters = `alphazo' 					in 5
	replace parameters = 0 							in 6/l
	
*outsheet z Hzy hzy Hzo hzo mtr omegagrid jomegao Jomegao parameters in 1/1126 using "Z:\Research\Tax career concerns\Simulations\My Matlab Code\statae05020mil.dat", replace nonames

* NEW: Add an export for z variables, not split by y or o
	gen double zparameters=.
	replace zparameters = ${epsilon} 				in 1
	replace zparameters = `zminalpha'^(1/`alphaz')  in 2
	replace zparameters = `alphaz' 					in 3
	replace zparameters = 0 						in 4/l
	
	keep z Hz hz mtr zparameters
	save "psidtaxsim_20mil.dta", replace

