/* ============================================================
Title: Reduced-Form Analysis
Author: Ryan Wang
Project: Organ Tranplant Herding
Last updated: 2025-12-02

Purpose: Construct analysis sample, generate descriptive
         statistics and figures, and run reduced-form
         herding tests.

Inputs:
    - Clean Organ Data Long.dta (not included in repository)

Outputs:
    - data_variables.csv
    - fig1_decision_distribution.png
    - fig2_ability_differential.png
    - tab1_measuring_center_ability.text
    - tab2_first_test.text
    - tab3_second_test.text
	- tab3_second_test_with_risk.text
    - tab4_second_test_restricted.text
	
Notes:
	- This .do file is a modified excerpt of code written as
	  part of my RA work for Kaivan, approved for release
	- Requires estout package (esttab/estadd commands)
============================================================== */


clear all

cd "ra-sample/data" // non-exist folder

use "Clean Organ Data Long.dta", clear // not included in repository

cd "ra-sample/output"

global sfile "ra-sample/output"



****************************************************
* STEP 0: Create Variable Inventory                *
****************************************************

preserve

ds
local vars `r(varlist)'

capture postutil clear
tempfile meta
postfile handle str32 varname str80 varlabel using `meta'

foreach v of local vars {
    local lbl : variable label `v'
    post handle ("`v'") ("`lbl'")
}
postclose handle
use `meta', clear

export delimited using "$sfile/var_inventory.csv", replace

restore



****************************************************
* STEP 1: Generate Appropriate Sub-Sample          *
****************************************************

* Fix mis-labelled year for a specific liver so that its first
* position is in 2006 (otherwise queue is broken)
replace year=2006 if organid=="0818340" & year==2005

* Keep only organs with full queue information up to position 3
keep if fullqueueto3 == 1
label variable fullqueueto3 "no gaps in queue up to position 3"

* Restrict to main analysis period (2006–2015)
keep if year >= 2006

* Donor brain-death indicator, with readable labels
label define dbd 0 "cardiac death" 1 "brain death"
label values donorbraindeath dbd

drop fasttrack // remove original fasttrack flag (re-computed below)

* Indicator for whether an organ ever goes to fast-track at any position
gen fasttrackever = .
replace fasttrackever = 1 if off_method == 4
sort donorid orgsubtype fasttrack
bysort organid: replace fasttrackever = fasttrackever[_n-1] if missing(fasttrackever)
label variable fasttrackever "organ goes to fast-track at some position"
replace fasttrackever = 0 if missing(fasttrackever)

* Indicator for whether a given position is on/after the first fast-track offer
gen fasttrackposition = 0
replace fasttrackposition = 1 if off_method == 4
sort organid position, stable
bysort organid position: ///
    replace fasttrackposition = 1 if fasttrackposition[_n-1] == 1
* once it goes to fast-track, all subsequent offers are also fast-track


* neverused = 1 for organs that are never accepted at any position
gen neverused = 0 if result == 5
sort organid neverused, stable
by organid: replace neverused = neverused[_n-1] if missing(neverused)
sort organid position, stable
replace neverused = 1 if neverused == .
label variable neverused "organs that are never accepted and used at any position"

* discardposition = 1 at the final position of organs that are ultimately discarded
gen discardposition = 0
sort organid position, stable
bysort organid: replace discardposition = 1 if _n == _N & (result == 2 | result == 4)
order discardposition, after(position)
replace discardposition = 0 if neverused == 0
label variable discardposition "organs that are discarded at the current position"

* Limit to first 8 positions and to livers/kidneys only
drop if position >= 9
drop if organtype == 3



****************************************************
* Figure 1: Distribution of Decisions by Position  *
****************************************************

preserve

gen counter = 1
collapse (count) n_decisions = counter, by(organtype position)

* Keep only livers (1) and kidneys (2)
drop if organtype != 1 & organtype != 2

* Compute within-organ-type share of decisions at each position
sort organtype, stable
bysort organtype: egen total_decisions = total(n_decisions)
gen prop = n_decisions / total_decisions

label define organlbl 1 "Liver" 2 "Kidney"
label values organtype organlbl

local barcol  "0 114 189"
local bwidth  0.80

* Liver panel
twoway bar prop position if organtype == 1, ///
    title("Liver") ///
    ytitle("proportion of total decisions") ///
    xtitle("queue position") ///
    xlabel(1(1)8, nogrid) ///
    ylabel(0(0.05)0.45, angle(horizontal) format(%4.1f)) ///
    yscale(range(0 0.45)) ///
    aspectratio(1.2) ///
    barwidth(`bwidth') ///
    fcolor("`barcol'") lcolor("`barcol'") lwidth(vthin) ///
    legend(off) plotregion(fcolor(white)) graphregion(fcolor(white)) ///
    name(fig_liver, replace)

* Kidney panel
twoway bar prop position if organtype == 2, ///
    title("Kidney") ///
    ytitle("proportion of total decisions") ///
    xtitle("queue position") ///
    xlabel(1(1)8, nogrid) ///
    ylabel(0(0.05)0.45, angle(horizontal) format(%4.1f)) ///
    yscale(range(0 0.45)) ///
    aspectratio(1.2) ///
    barwidth(`bwidth') ///
    fcolor("`barcol'") lcolor("`barcol'") lwidth(vthin) ///
    legend(off) plotregion(fcolor(white)) graphregion(fcolor(white)) ///
    name(fig_kidney, replace)


graph combine fig_liver fig_kidney, col(2) ///
    graphregion(color(white) margin(15 15 5 5)) ///
    xsize(9) ysize(4.5)

graph drop fig_liver fig_kidney

graph export "$sfile/fig1_decision_distribution.png", width(1600) replace

restore



****************************************************
* STEP 3: Reduced-Form Herding Tests               *
****************************************************

*-----------------------------------------------*
* 3A. Baseline P(reject) by center × organ      *
*-----------------------------------------------*

* Estimate center fixed effects in first position
xi: regress rejected i.centreid if position == 1 & organtype == 1, robust
predict ihatliver if position == 1 & organtype == 1, xb    // fitted prob
predict errorliver, stdp                                   // std error of prediction

xi: regress rejected i.centreid if position == 1 & organtype == 2, robust
predict ihatkidney if position == 1 & organtype == 2, xb
predict errorkidney, stdp

* Construct center×organ baseline P(reject) using fitted values
gen probreject = ihatliver if position == 1 & organtype == 1
replace probreject = ihatkidney if position == 1 & organtype == 2

* 95% CI bounds for baseline P(reject)
gen lb95 = ihatliver - invnormal(0.975)*errorliver if position == 1 & organtype == 1
replace lb95 = ihatkidney - invnormal(0.975)*errorkidney if position == 1 & organtype == 2

gen ub95 = ihatliver + invnormal(0.975)*errorliver if position == 1 & organtype == 1
replace ub95 = ihatkidney + invnormal(0.975)*errorkidney if position == 1 & organtype == 2

* Fill baseline P(reject) down all positions for a given center×organ type
sort centreid organtype position, stable
bysort centreid organtype: replace probreject = probreject[_n-1] if missing(probreject)

label var lb95 "95% CI lower: P(reject) at pos1"
label var ub95 "95% CI upper: P(reject) at pos1"
label var probreject "Baseline P(reject) (center×organ), filled to all positions"

* Separate baseline probabilities for positions 1–3 of a given organ
gen probrejectfirst  = .
gen probrejectsecond = .
gen probrejectthird  = .
sort organid position, stable
bysort organid: replace probrejectfirst  = probreject[1]
bysort organid: replace probrejectsecond = probreject[2]
bysort organid: replace probrejectthird  = probreject[3]

label var probrejectfirst  "Baseline P(reject) of center in pos1"
label var probrejectsecond "Baseline P(reject) of center in pos2"
label var probrejectthird  "Baseline P(reject) of center in pos3"

*-----------------------------------------------*
* 3B. Proxy measures of center quality          *
*-----------------------------------------------*

* For livers (organtype==1) high P(reject) = "high quality";
* for kidneys (organtype==2) low P(reject) = "high quality"
gen proxy_quality = probreject if organtype == 1
replace proxy_quality = 1 - probreject if organtype == 2

gen proxy_qualityfirst  = probrejectfirst  if organtype == 1
replace proxy_qualityfirst  = 1 - probrejectfirst  if organtype == 2

gen proxy_qualitysecond = probrejectsecond if organtype == 1
replace proxy_qualitysecond = 1 - probrejectsecond if organtype == 2

gen proxy_qualitythird  = probrejectthird  if organtype == 1
replace proxy_qualitythird  = 1 - probrejectthird


****************************************************
* Table 1: Measuring Center Ability                *
****************************************************

* Average risk index (dli for livers, ukdri for kidneys)
sort organtype centreid, stable
bysort organtype centreid: egen averagedli   = mean(dli)   if organtype == 1 & position == 1
bysort organtype centreid: egen averageukdri = mean(ukdri) if organtype == 2 & position == 1

* Fill averages down all positions within center×organ type
sort organtype centreid position, stable
bysort organtype centreid: replace averagedli   = averagedli[_n-1]   if missing(averagedli)
bysort organtype centreid: replace averageukdri = averageukdri[_n-1] if missing(averageukdri)

* For each organ, risk index of the first-position center
bysort organid (position): gen averagerisk_first = averagedli[1]   if organtype == 1
bysort organid (position): replace averagerisk_first = averageukdri[1] if organtype == 2
label variable averagerisk_first "average risk index of center 1 in first position"

* Regressions of rejection on baseline P(reject) and risk index
eststo clear

regress rejected probrejectfirst i.centreid if position == 2 & organtype == 1, robust
quietly summarize rejected if e(sample)
estadd scalar ymean = r(mean)
estadd local Rbar "No"
estadd local fe2  "Yes"
eststo

regress rejected probrejectfirst averagerisk_first i.centreid ///
    if position == 2 & organtype == 1, robust
quietly summarize rejected if e(sample)
estadd scalar ymean = r(mean)
estadd local Rbar "Yes"
estadd local fe2  "Yes"
eststo

regress rejected probrejectfirst i.centreid if position == 2 & organtype == 2, robust
quietly summarize rejected if e(sample)
estadd scalar ymean = r(mean)
estadd local Rbar "No"
estadd local fe2  "Yes"
eststo

regress rejected probrejectfirst averagerisk_first i.centreid ///
    if position == 2 & organtype == 2, robust
quietly summarize rejected if e(sample)
estadd scalar ymean = r(mean)
estadd local Rbar "Yes"
estadd local fe2  "Yes"
eststo

* Export Table 1
esttab * using "$sfile/tab1_measuring_center_ability.text", replace ///
    mtitles(liver liver kidney kidney) nocons label ///
    title("Measuring Center Ability") ///
    keep(probrejectfirst) ///
    varlabels(probrejectfirst "$\overline{p}_1$") ///
    stats(Rbar fe2 ymean N, fmt(%9s %9s 3 0) ///
          labels("$\overline{R}_1$" "Center 2 fixed effects" ///
                 "Mean of dependent variable" "N")) ///
    se(3) starlevels(* 0.10 ** 0.05 *** 0.001)
eststo clear


****************************************************
* Table 2: First Test of Herding (Position 2)      *
****************************************************

* Center ID of first-position center for each organ
bysort organid (position): gen centreid_first = centreid[1]

eststo clear

* Interaction of quality of first two centers, no fixed effects
regress rejected c.proxy_qualityfirst##c.proxy_qualitysecond ///
    if position == 2 & organtype == 1, robust
quietly summarize rejected if e(sample)
estadd scalar ymean = r(mean)
estadd local fe1  "No"
estadd local fe2  "No"
eststo

* With center fixed effects for both first and second centers
regress rejected c.proxy_qualityfirst#c.proxy_qualitysecond ///
    i.centreid i.centreid_first if position == 2 & organtype == 1, robust
quietly summarize rejected if e(sample)
estadd scalar ymean = r(mean)
estadd local fe1  "Yes"
estadd local fe2  "Yes"
eststo

* Repeat for kidneys
regress rejected c.proxy_qualityfirst##c.proxy_qualitysecond ///
    if position == 2 & organtype == 2, robust
quietly summarize rejected if e(sample)
estadd scalar ymean = r(mean)
estadd local fe1  "No"
estadd local fe2  "No"
eststo

regress rejected c.proxy_qualityfirst#c.proxy_qualitysecond ///
    i.centreid i.centreid_first if position == 2 & organtype == 2, robust
quietly summarize rejected if e(sample)
estadd scalar ymean = r(mean)
estadd local fe1  "Yes"
estadd local fe2  "Yes"
eststo

* Export Table 2
esttab * using "$sfile/tab2_first_test.text", replace ///
    label ///
    title("First Test of Herding (based on decisions in second position)") ///
    mtitles(liver liver kidney kidney) ///
    coeflabels(proxy_qualityfirst "Center 1 ability q1" ///
               proxy_qualitysecond "Center 2 ability q2" ///
               c.proxy_qualityfirst#c.proxy_qualitysecond "q1 x q2") ///
    keep(proxy_qualityfirst proxy_qualitysecond ///
         c.proxy_qualityfirst#c.proxy_qualitysecond) ///
    stats(fe1 fe2 ymean N, fmt(%9s %9s 3 0)) ///
    nocons se(3) starlevels(* 0.10 ** 0.05 *** 0.001)

eststo clear


****************************************************
* Table 3: Second Test of Herding (Position 3)     *
****************************************************

eststo clear

* Liver sample
regress rejected c.proxy_qualityfirst c.proxy_qualitysecond ///
    if position == 3 & organtype == 1, robust
quietly summarize proxy_qualityfirst if e(sample)
estadd scalar q1_bar = r(mean)
quietly summarize proxy_qualitysecond if e(sample)
estadd scalar q2_bar = r(mean)

test proxy_qualitysecond - proxy_qualityfirst = 0
estadd scalar F_stat = r(F)
local sign = sign(_b[proxy_qualitysecond] - _b[proxy_qualityfirst])
estadd scalar p_val = normal(`sign'*sqrt(r(F)))
eststo

* Kidney sample
regress rejected c.proxy_qualityfirst c.proxy_qualitysecond ///
    if position == 3 & organtype == 2, robust
quietly summarize proxy_qualityfirst if e(sample)
estadd scalar q1_bar = r(mean)
quietly summarize proxy_qualitysecond if e(sample)
estadd scalar q2_bar = r(mean)

test proxy_qualitysecond - proxy_qualityfirst = 0
estadd scalar F_stat = r(F)
local sign = sign(_b[proxy_qualitysecond] - _b[proxy_qualityfirst])
estadd scalar p_val = normal(`sign'*sqrt(r(F)))
eststo

* Export Table 3
esttab * using "$sfile/tab3_second_test.text", replace ///
    label ///
    title("Second Test of Herding (based on decisions in third position)") ///
    mtitles(liver kidney) ///
    varlabels(proxy_qualityfirst  "Center 1 ability q1" ///
              proxy_qualitysecond "Center 2 ability q2" ///
              _cons "Constant") ///
    stats(F_stat p_val q1_bar q2_bar N, fmt(2 3 2 2 0)) ///
    se(3) starlevels(* 0.10 ** 0.05 *** 0.001)

eststo clear


****************************************************
* Table 3 V2: Add average risk index covariate     *
****************************************************

eststo clear

* Liver with and without risk index
regress rejected c.proxy_qualityfirst c.proxy_qualitysecond ///
    if position == 3 & organtype == 1, robust
quietly summarize proxy_qualityfirst if e(sample)
estadd scalar q1_bar = r(mean)
quietly summarize proxy_qualitysecond if e(sample)
estadd scalar q2_bar = r(mean)

test proxy_qualitysecond - proxy_qualityfirst = 0
estadd scalar F_stat = r(F)
local sign = sign(_b[proxy_qualitysecond] - _b[proxy_qualityfirst])
estadd scalar p_val = normal(`sign'*sqrt(r(F)))
eststo

regress rejected c.proxy_qualityfirst c.proxy_qualitysecond averagerisk_first ///
    if position == 3 & organtype == 1, robust
quietly summarize proxy_qualityfirst if e(sample)
estadd scalar q1_bar = r(mean)
quietly summarize proxy_qualitysecond if e(sample)
estadd scalar q2_bar = r(mean)

test proxy_qualitysecond - proxy_qualityfirst = 0
estadd scalar F_stat = r(F)
local sign = sign(_b[proxy_qualitysecond] - _b[proxy_qualityfirst])
estadd scalar p_val = normal(`sign'*sqrt(r(F)))
eststo

* Kidney with and without risk index
regress rejected c.proxy_qualityfirst c.proxy_qualitysecond ///
    if position == 3 & organtype == 2, robust
quietly summarize proxy_qualityfirst if e(sample)
estadd scalar q1_bar = r(mean)
quietly summarize proxy_qualitysecond if e(sample)
estadd scalar q2_bar = r(mean)

test proxy_qualitysecond - proxy_qualityfirst = 0
estadd scalar F_stat = r(F)
local sign = sign(_b[proxy_qualitysecond] - _b[proxy_qualityfirst])
estadd scalar p_val = normal(`sign'*sqrt(r(F)))
eststo

regress rejected c.proxy_qualityfirst c.proxy_qualitysecond averagerisk_first ///
    if position == 3 & organtype == 2, robust
quietly summarize proxy_qualityfirst if e(sample)
estadd scalar q1_bar = r(mean)
quietly summarize proxy_qualitysecond if e(sample)
estadd scalar q2_bar = r(mean)

test proxy_qualitysecond - proxy_qualityfirst = 0
estadd scalar F_stat = r(F)
local sign = sign(_b[proxy_qualitysecond] - _b[proxy_qualityfirst])
estadd scalar p_val = normal(`sign'*sqrt(r(F)))
eststo

* Export Table 3 with risk covariate
esttab * using "$sfile/tab3_second_test_with_risk.text", replace ///
    label ///
    title("Second Test of Herding (based on decisions in third position)") ///
    mtitles(liver kidney) ///
    varlabels(proxy_qualityfirst  "Center 1 ability q1" ///
              proxy_qualitysecond "Center 2 ability q2" ///
              averagerisk_first   "$\overline{R}_1$" ///
              _cons "Constant") ///
    stats(F_stat p_val q1_bar q2_bar N, fmt(2 3 2 2 0)) ///
    se(3) starlevels(* 0.10 ** 0.05 *** 0.001)

eststo clear


****************************************************
* Figure 4: Ability Differential Distribution      *
****************************************************

* q1minusq2 = difference in quality between center 1 and center 2
gen q1minusq2 = proxy_qualityfirst - proxy_qualitysecond

preserve

* Focus on kidneys, position 3 (main sample)
keep if organtype == 2 & position == 3

* Cutoffs for increasingly restrictive samples
local c1   = 0.30
local c2   = 0.05
local c1p  = 0.275
local c2p  = 0.075

local ytickn  = 0.375
local ytickp  = 0.36
local ytick2n = 2.5
local ytick2p = 2.46
local yarrow  = 0.1

twoway ///
    (kdensity q1minusq2 if organtype == 2 & position == 3, ///
        bwidth(0.05) lwidth(medthick)) ///
    (pci 0  `c1'  `ytickp'  `c1',  lcolor(black) lpattern(dash) lwidth(medthick)) ///
    (pci 0 -`c1'  `ytickn' -`c1',  lcolor(black) lpattern(dash) lwidth(medthick)) ///
    (pci 0  `c2'  `ytick2p' `c2',  lcolor(red)   lpattern(dash) lwidth(medthick)) ///
    (pci 0 -`c2'  `ytick2n' -`c2', lcolor(red)   lpattern(dash) lwidth(medthick)) ///
    (pcarrowi `yarrow' -`c1p' `yarrow' -`c2p', lcolor(black) lwidth(medthick) mcolor(black)) ///
    (pcarrowi `yarrow'  `c1p' `yarrow'  `c2p', lcolor(black) lwidth(medthick) mcolor(black)) ///
    , ///
    xtitle("{it:q}{sub:1} — {it:q}{sub:2}", size(small)) ///
    ytitle("density", size(small)) ///
    xlabel(-0.6(0.2)0.9, labsize(small) format(%4.1f) ///
           tlcolor(gs8) tlength(-1) labcolor(black) labgap(2) nogrid) ///
    ylabel(0(0.5)3, angle(horizontal) labsize(small) nogrid format(%4.1f)) ///
    yscale(range(0 .)) ///
    graphregion(margin(large)) ///
    plotregion(margin(zero)) ///
    aspectratio(1) xsize(9) ysize(4.5) ///
    legend(order(2 "least restrictive sample" 4 "most restrictive sample") ///
           size(vsmall) position(2) ring(0) col(1) region(lstyle(none))) ///
    plotregion(color(white)) graphregion(color(white))

graph export "$sfile/fig2_ability_differential.png", width(1600) replace

restore



****************************************************
* Table 4: Second Test of Herding (Restrictions)   *
****************************************************

eststo clear

* Vary sample restriction |q1 - q2| ≤ c, conditioning on risk index
regress rejected c.proxy_qualityfirst c.proxy_qualitysecond averagerisk_first ///
    if position == 3 & organtype == 2 & abs(q1minusq2) <= 0.30, robust
test proxy_qualitysecond - proxy_qualityfirst = 0
estadd scalar F_stat = r(F)
local sign = sign(_b[proxy_qualitysecond] - _b[proxy_qualityfirst])
estadd scalar p_val = normal(`sign'*sqrt(r(F)))
eststo

regress rejected c.proxy_qualityfirst c.proxy_qualitysecond averagerisk_first ///
    if position == 3 & organtype == 2 & abs(q1minusq2) <= 0.20, robust
test proxy_qualitysecond - proxy_qualityfirst = 0
estadd scalar F_stat = r(F)
local sign = sign(_b[proxy_qualitysecond] - _b[proxy_qualityfirst])
estadd scalar p_val = normal(`sign'*sqrt(r(F)))
eststo

regress rejected c.proxy_qualityfirst c.proxy_qualitysecond averagerisk_first ///
    if position == 3 & organtype == 2 & abs(q1minusq2) <= 0.10, robust
test proxy_qualitysecond - proxy_qualityfirst = 0
estadd scalar F_stat = r(F)
local sign = sign(_b[proxy_qualitysecond] - _b[proxy_qualityfirst])
estadd scalar p_val = normal(`sign'*sqrt(r(F)))
eststo

regress rejected c.proxy_qualityfirst c.proxy_qualitysecond averagerisk_first ///
    if position == 3 & organtype == 2 & abs(q1minusq2) <= 0.05, robust
test proxy_qualitysecond - proxy_qualityfirst = 0
estadd scalar F_stat = r(F)
local sign = sign(_b[proxy_qualitysecond] - _b[proxy_qualityfirst])
estadd scalar p_val = normal(`sign'*sqrt(r(F)))
eststo

* Export Table 4
esttab * using "$sfile/tab4_second_test_restricted.text", replace ///
    label ///
    title("Second Test of Herding (restricted samples)") ///
    drop(averagerisk_first) ///
    varlabels(proxy_qualityfirst "Center 1 ability q1" ///
              proxy_qualitysecond "Center 2 ability q2" ///
              _cons "Constant") ///
    mtitles("[-0.30,0.30]" "[-0.20,0.20]" "[-0.10,0.10]" "[-0.05,0.05]") ///
    stats(F_stat p_val N, fmt(2 3 0)) ///
    se(3) starlevels(* 0.10 ** 0.05 *** 0.001) ///
    note("average quality of organs received by center 1 in first position included as covariate")

eststo clear
