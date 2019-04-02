READ ME FOR TIMETOGERMINATE REPO
updated April 1 2019 by Dan

Data:
  2 data sheets
  germ_datasheet.xls is the main datafile from the germination experiment
  carex_data.xlsx is a data sheet just for carex grayii which was left in chamber longer than all other species
  
input:
  dat_chill_X.csv- a series of intermediated files that were cleaned from the original data
 germ_data_forDRC.csv- A file formated for DRC event time formal
 dail_dat_nointerval.csv- a data file with both cumlative and daily germination counts without interval censoring.
 
cleaning:
A.germination_cleaning.R- Main cleaning file
 drc_starter_clean- follows A to clean for drc format
 make_cum_daily_no_int- cleans raw data to make a data shee that is cumulative and doesn't ahve interval censoring
 make_survial-formats data for survial analysis
 
 nonlinear- A folder for nonlinear modeling based on cumulative germination data.
  rewind_fake.R Fake data gemerate and stan model prep

stan:
  A. fakeseedmodel.stan-  loglogistic stan model a germination time curve with no treatment levels **runs and returns starter paramenters**
  B. altfakeseed.model.stan- same as above but a more complicted function- ***many divergent transitions**
  C. fakeseed_chillonly.stan A model for germination courses with a 0,1 chilling treatment based on modified version of B, as of 4/2/19 now with alternative parameterization from wikipedia.  ***only throws ~40 divergent transitions***
  D. fakeseedgoodchill.stan A model for germination courses with a 0,1 chilling treatment based on modified version of A **300ish divergent transitions.
  E. fakeseed_forceonly.stan Same as C but for forcing
  F. fakeseed_forcechill_noint.stan A broken model to try to incorperate both forcing and chilling onto fake germinationd ata
  G. twoparamfake.stan incomplete script to attempt F in a different way.