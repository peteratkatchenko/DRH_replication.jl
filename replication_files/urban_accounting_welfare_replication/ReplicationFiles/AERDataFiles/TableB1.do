* This program generates TABLE B1 of "Urban Accounting and Welfare" by Klaus Desmet and Esteban Rossi-Hansberg
* October 11, 2012

clear
set more off
set memory 512m

use amenities.dta

log using TableB1.log, replace

pwcorr amenities lowjanf, sig
pwcorr amenities days90, sig 
pwcorr amenities days32, sig 
pwcorr amenities annualprecip, sig 
pwcorr amenities daysprecip, sig 
pwcorr amenities days02precip, sig 
pwcorr amenities relhumidjuly, sig 
pwcorr amenities julyheat, sig

pwcorr amenities pra_transport, sig 
pwcorr amenities pra_edu, sig 
pwcorr amenities pra_crime, sig 
pwcorr amenities pra_arts, sig 
pwcorr amenities pra_health, sig 
pwcorr amenities pra_recreation, sig 

pwcorr amenities crr_edu, sig
pwcorr amenities crr_health, sig
pwcorr amenities crr_crime, sig
pwcorr amenities crr_transport, sig
pwcorr amenities crr_leisure, sig
pwcorr amenities crr_arts_culture, sig
pwcorr amenities crr_qol, sig

pwcorr amenities proxcoast, sig
pwcorr amenities proxwater, sig

log close

