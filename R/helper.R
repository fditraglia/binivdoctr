# Create factor with groups indicating all possible interactions of a number of
# dummy variables.
dummyInteractions <- function(dummies){
  out <- apply(dummies, 1, paste, collapse = "")
  out <- strtoi(out, base = 2)
  out <- factor(out)
  levels(out) <- 1:length(levels(out))
  return(out)
}

