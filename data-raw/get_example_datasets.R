# ======================================================================
# Dinkelman (2011, AER) - Effects of Rural Electrification on Employment
# Based on the data used to construct Table 4 of the paper
# ======================================================================
download.file("https://www.aeaweb.org/aer/data/dec2011/20080791_data.zip",
              "./data-raw/dinkelman.zip")
unzip("./data-raw/dinkelman.zip",
      files = "20080791_dataset/data/matched_censusdata.dta",
      exdir = "./data-raw", junkpaths = TRUE)
dinkelman <- haven::read_dta("./data-raw/matched_censusdata.dta")
system("rm ./data-raw/dinkelman.zip")

# Restrict to areas with at least 100 adults in both years
dinkelman <- subset(dinkelman, largeareas == 1)

keep <- c(# Treatment
          "T", # ESKOM (electrification) Project? (Dummy)
          # Outcomes
          "d_prop_emp_f", # Change in female employment rate
          "d_prop_emp_m", # Change in male employment rate
          # Instrument
          "mean_grad_new", # Land gradient
          # Fixed Effects
          "dccode0", # District index (factor)
          # Controls
          "kms_to_subs0",
          "baseline_hhdens0",
          "base_hhpovrate0",
          "prop_head_f_a0",
          "sexratio0_a",
          "prop_indianwhite0",
          "kms_to_road0",
          "kms_to_town0",
          "prop_matric_m0",
          "prop_matric_f0",
          "d_prop_waterclose",
          "d_prop_flush")
dinkelman <- dinkelman[, keep]
rm(keep)

# Rename "T" to escom so it doesn't clash with "TRUE"
names(dinkelman)[which(names((dinkelman)) == "T")] <- "escom"

# Rescale certain variables to match magnitudes in paper
rescale <- c("mean_grad_new", "kms_to_subs0", "baseline_hhdens0",
             "kms_to_road0", "kms_to_town0")
dinkelman[,rescale] <- dinkelman[,rescale] / 10
rm(rescale)

# Create binary version of instrument
dinkelman$steep <- (dinkelman$mean_grad_new > 1) # use 1 as cutoff ~= median and FAO's def of strongly sloping

# Replace dccode0 (a factor for district) with a group of dummy variables that
# represent the same information, treating F21 as the "baseline"
district <- as.factor(dinkelman$dccode0)
dinkelman$dccode0 <- NULL
dummies <- model.matrix(~ district)
dummies <- dummies[,-1]
dinkelman <- cbind(dinkelman, dummies)
rm(district, dummies)

devtools::use_data(dinkelman, overwrite = TRUE)
rm(dinkelman)


# ======================================================================
# Angrist et al (2002, AER) - Vouchers for Private Schooling
# Data needed to replicate selected results from Table 7
# ======================================================================
download.file("http://economics.mit.edu/files/1395", "./data-raw/tab7.sas7bdat")
angrist <- haven::read_sas("./data-raw/tab7.sas7bdat")

# Convert variable names to lowercase
names(angrist) <- tolower(names(angrist))

# Construct highest grade completed outcome
grades <- angrist[,c("finish6", "finish7", "finish8")]
g <- function(x){
  out <- 5
  if(all(x == c(1, 0, 0)))
    out <- 6
  if(all(x == c(1, 1, 0)))
    out <- 7
  if(all(x == c(1, 1, 1)))
    out <- 8
  return(out)
}
highestGrade <- apply(grades, 1, g)
angrist$highestGrade <- highestGrade
rm(grades, highestGrade, g)

# Construct total repetitions outcome
rept <- angrist[,c("rept6", "rept7", "rept8")]
f <- function(x){
  x[2] <- x[2] - 1
  x[3] <- x[3] - 1
  return(sum(x, na.rm = TRUE))
}
totalRepeats <- apply(rept, 1, f)
angrist$totalRepeats <- totalRepeats
rm(totalRepeats, f, rept)

# In their SAS file, Angrist et al restrict their sample using the following
# condition, on which they do not comment in either the paper or the code.
# It is not innocuous, since it rules out individuals for whom we are not in
# fact missing the variables needed to carry out some of the regressions in
# Table 7. Moreover, some of these variables on which they exclude do not even
# appear to be used in the paper. Nevertheless, we will restrict our sample as
# they do so that our results are comparable.
keepObs <- with(angrist, !is.na(scyfnsh) & !is.na(finish6) & !is.na(rept6) &
                  !is.na(totalRepeats) & !is.na(svy) &
                  !is.na(inschl) & !is.na(finish7) &
                  !is.na(vouch0) & !is.na(prsch_c) &
                  !is.na(finish8) & !is.na(prscha_2) &
                  !is.na(totscyrs) & !is.na(rept) &
                  ((bog95smp == 1) | (bog97smp == 1) | (jam93smp == 1)))

angrist <- angrist[keepObs,]
rm(keepObs)

# Note that we do not include some of the Dummies mentioned by Angrist et al.
# This is because they are collinear and are automatically dropped when running
# OLS or IV. The ones that we excluded are:
#   strata6
#   stratams
#   dbogota
#   d1993
#   d1997
#   dmonth12
# We chose these since they are the ones that are automatically exluded by R's
# lm function when it runs the OLS, reduced form, and first stage regressions.

keepVars <- c(
# Covariates
  "age", # The only covariate that isn't a dummy
  "svy",
  "hsvisit",
  "djamundi",
  "phone",
  "sex2",
  "d1995",
  "strata1",
  "strata2",
  "strata3",
  "strata4",
  "strata5",
  "dmonth1",
  "dmonth2",
  "dmonth3",
  "dmonth4",
  "dmonth5",
  "dmonth6",
  "dmonth7",
  "dmonth8",
  "dmonth9",
  "dmonth10",
  "dmonth11",
  # Regressor
  "usesch",
  # Instrument
  "vouch0",
  # Outcomes
  "highestGrade",
  "totalRepeats",
  # Indicators of different sample-city years
  "bog95smp", # Bogota 1995
  "bog97smp", # Bogota 1997
  "jam93smp", # Jamundi 1993
  "prscha_1", # Started 6th grade in private school
  "prscha_2", # Started 7th grade in private school
  "prsch_c" # Currently in private school
)

angrist <- angrist[, keepVars]
rm(keepVars)

# Note: there is one individual in the dataset whose sex is unknown. Although
# Angrist et al do not make this clear, the rule they use to exclude
# observations in their SAS file already drops this individual. Thus, although
# they include a dummy for sex being missing (sex_miss) in their regressions,
# this is just being dropped by SAS since it's the same for everyone.
# Accordingly, we don't bother to store it.

devtools::use_data(angrist, overwrite = TRUE)
rm(angrist)

# ======================================================================
# Gelbach (2002, AER) - Public Schooling for Young Children
# Based on data used in Tables 6 and 7
# ======================================================================
download.file("http://gelbach.law.upenn.edu/phil/interact.dta",
              "./data-raw/interact.dta")
gelbach <- haven::read_dta("./data-raw/interact.dta")

# Restrict sample to mothers whose youngest child is 5 years old
gelbach <- subset(gelbach, youngest == 5)

# Keep only the variables used in Tables 6 and 7
keep <- c(
          # Controls
          "num612",   # Number of own children in household aged 6-12
          "num1317",  # Number of own children in household aged 13-17
          "numge18",  # Number of own children in household aged >= 18
          "othrlt18", # Number of other household members aged < 18
          "othrge18", # Number of other household members aged >= 18
          "grade",    # Mother's years of education
          "white",    # White? (dummy variable)
          "centcity", # Live in central city? (dummy variable)
          "age",      # Age of mother
          "age2",     # Squared age of mother
          # Factors for Fixed Effects
          "stater", # State of residence
          "stateb", # State of birth
          # Instruments
          "qtr1", # Born in Quarter I? (dummy variable)
          "quarter", # Quarter of birth
          # Regressor
          "public", # Attend public school? (dummy)
          # outcomes
          "hours", # hours worked last week
          "weeksw79", # weeks worked in 1979
          "salary") # wage and salary income in 1979
gelbach <- gelbach[keep]
rm(keep)

devtools::use_data(gelbach, overwrite = TRUE)
rm(gelbach)

# ===========================================================================
# Burde and Linden (2013, AEJ Applied) - Bringing Education to Afghan Girls
# The paper does not report IV. We use the controls from Table 2.
# ===========================================================================
download.file("https://www.aeaweb.org/aej/app/data/2012-0252_data.zip",
              "./data-raw/afghan.zip")
unzip("./data-raw/afghan.zip",
      files = "Data_20120252 2015-06-15/afghanistan_anonymized_data.dta",
      exdir = "./data-raw", junkpaths = TRUE)
afghan <- haven::read_dta("./data-raw/afghanistan_anonymized_data.dta")
system("rm ./data-raw/afghan.zip")

# Remove outliers following the authors' STATA code
outlier <- with(afghan, (f07_num_ppl_hh_cnt > 20 & f07_observed == 1) |
                  (f07_jeribs_cnt > 10 & f07_observed == 1) |
                  (f07_num_sheep_cnt > 50 & f07_observed == 1) |
                  (s08_num_ppl_hh_cnt > 20 & s08_observed == 1) |
                  (s08_jeribs_cnt > 10 & s08_observed == 1) |
                  (s08_num_sheep_cnt > 50 & s08_observed == 1))
afghan <- afghan[!outlier,]

# The variable "treatment" is actually our instrument
# --- living in a village that gets a school

# The variable "f08_formal_school" is our treatment
# --- enrolled in formal school, Spring 2008

# The variable "s08_formal_school" is our outcome
# --- total normalized test score, Spring 2008
afghan <- with(afghan, data.frame(
                               "enrolled" = s08_formal_school,
                               "testscore" = s08_both_norma_total,
                               "buildschool" = treatment,
                               #"c" = as.factor(clustercode), #only for cluster SE
                               "headchild" = s08_heads_child_cnt,
                               "female" = s08_girls_cnt,
                               "age" = s08_age_cnt,
                               "yrsvill" = s08_duration_village_cnt,
                               "farsi" = s08_farsi_cnt,
                               "tajik" = s08_tajik_cnt,
                               "farmers" = s08_farmer_cnt,
                               "agehead" = s08_age_head_cnt,
                               "educhead" = s08_yrs_ed_head_cnt,
                               "nhh" = s08_num_ppl_hh_cnt,
                               "land" = s08_jeribs_cnt,
                               "sheep" = s08_num_sheep_cnt,
                               "distschool" = s08_nearest_scl,
                               "chagcharan" = chagcharan)) # District dummy
                                # Not mentioned in paper!

# remove missing observations
afghan <- na.omit(afghan)

# only look at girls
afghan <- subset(afghan, female == 1)

devtools::use_data(afghan, overwrite = TRUE)
rm(afghan)

