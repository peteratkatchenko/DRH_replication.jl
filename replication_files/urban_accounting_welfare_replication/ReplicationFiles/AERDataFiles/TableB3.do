* This program generates TABLE B3 of "Urban Accounting and Welfare" by Klaus Desmet and Esteban Rossi-Hansberg
* October 11, 2012


clear
set more off
set memory 512m

use FrictionMeasures.dta

log using TableB3.log, replace

log off

* DEFINE VARIABLES

gen oneminustauav=(oneminustau2005+oneminustau2006+oneminustau2007+oneminustau2008)/4
gen jointincconstaxav=(jointincconstax2005+jointincconstax2006+jointincconstax2007+jointincconstax2008)/4

gen laborwedgelev2005=1-oneminustau2005
gen laborwedgelev2006=1-oneminustau2006
gen laborwedgelev2007=1-oneminustau2007
gen laborwedgelev2008=1-oneminustau2008

gen inctaxav=(inctax2005+inctax2006+inctax2007+inctax2008)/4
gen laborwedgelevav=(laborwedgelev2005+laborwedgelev2006+laborwedgelev2007+laborwedgelev2008)/4
gen constaxav=(constax2005+constax2006+constax2007+constax2008)/4

gen fedconsper2005=fedcons2005/pop2005
gen fedconsper2006=fedcons2006/pop2006
gen fedconsper2007=fedcons2007/pop2007
gen fedconsper2008=fedcons2008/pop2008
gen fedconsperav=(fedconsper2005+fedconsper2006+fedconsper2007+fedconsper2008)/4

gen stateconsper2005=statecons2005/pop2005
gen stateconsper2006=statecons2006/pop2006
gen stateconsper2007=statecons2007/pop2007
gen stateconsperav=(stateconsper2005+stateconsper2006+stateconsper2007)/3

gen govconsper2005=fedconsper2005+stateconsper2005
gen govconsper2006=fedconsper2006+stateconsper2006
gen govconsper2007=fedconsper2007+stateconsper2007
gen govconsperav=(govconsper2005+govconsper2006+govconsper2007)/3
gen laborwedgelevavthree=(laborwedgelev2005+laborwedgelev2006+laborwedgelev2007)/3

gen commuteper2008=aggcomm2008/pop2008
gen logcommuteper2008=log(commuteper2008)
gen sharecommute2008=(commuteper2008/60)/(avhours2008+commuteper2008/60)

gen percunionprivav=(percunionpriv2005+percunionpriv2006+percunionpriv2007+percunionpriv2008)/4
gen percunionpubav=(percunionpub2005+percunionpub2006+percunionpub2007+percunionpub2008)/4
gen percuniontotav=(percuniontot2005+percuniontot2006+percuniontot2007+percuniontot2008)/4

log on

* CORRELATIONS TABLE B3

* TAXATION 

* Correlation between laborwedge and inctax (first correlation Table B3)
pwcorr laborwedgelevav inctaxav, sig

* Correlation between laborwedge and constax (second correlation Table B3)
pwcorr laborwedgelevav constaxav, sig


* Correlation between (1-laborwedge) and (1-inctax)/(1+constax) (third correlation Table B3)
pwcorr oneminustauav jointincconstaxav, sig


* GOVERNMENT SPENDING

* Correlation between laborwedge and federal govt spending per person (fourth correlation Table B3)
pwcorr laborwedgelevav fedconsperav, sig


* Correlation between laborwedge and state & local spending per person (fifth correlation Table B3)
pwcorr laborwedgelevavthree stateconsperav, sig 

* Correlation between laborwedge and total govt spending per person (sixth correlation Table B3)
pwcorr laborwedgelevavthree govconsperav, sig 


* COMMUTING COSTS

* Correlation between labor wedge and minutes lost commuting per person (seventh correlation Table B3)
pwcorr laborwedgelev2008 sharecommute2008, sig


* UNIONIZATION

* Correlation between labor wedge and percentage unionization in total economy (eighth correlation Table B3)
pwcorr laborwedgelevav percuniontotav, sig

* Correlation between labor wedge and percentage unionization in private sector (ninth correlation Table B3)
pwcorr laborwedgelevav percunionprivav, sig

* Correlation between labor wedge and percentage unionization in public sector (tenth correlation Table B3)
pwcorr laborwedgelevav percunionpubav, sig


* LAND REGULATION
* Correlation between labor wedge and Wharton Land Regulation Index (eleventh correlation Table B4)
pwcorr laborwedgelev2008 WRLURI, sig

log close

more


