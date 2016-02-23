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

# Rescale certain variables to match magnitudes in paper
rescale <- c("mean_grad_new", "kms_to_subs0", "baseline_hhdens0",
             "kms_to_road0", "kms_to_town0")
dinkelman[,rescale] <- dinkelman[,rescale] / 10
rm(rescale)

# Create binary version of instrument
dinkelman$steep <- (dinkelman$mean_grad_new > 1) # use 1 as cutoff ~= median and FAO's def of strongly sloping

devtools::use_data(dinkelman)
rm(dinkelman)


# ======================================================================
# Angrist et al (2002, AER) - Vouchers for Private Schooling
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

keep <- c(# Covariates
          "age",
          "svy",
          # Dummies
          "strata1",
          "strata2",
          "strata3",
          "strata4",
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
          "dmonth12",
          # Regressor
          "usesch",
          # Instrument
          "vouch0",
          # Outcomes
          "highestGrade",
          "totalRepeats",
          # Region Indicators
          "jam93smp", # Jamundi
          "bog95smp") # Bogota

angrist <- angrist[, keep]
rm(keep)

devtools::use_data(angrist)
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

devtools::use_data(gelbach)
rm(gelbach)



