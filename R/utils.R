# median impute
impute <- function(vec){
  vec[is.na(vec)] <- median(vec, na.rm = TRUE)
  vec
}

# logical or infix function
`%||%` <- function(a, b) if (!is.null(a)) a else b
