# check if a vector has at least two valid responses
invalid_grm <- function(x) max(x, na.rm = TRUE) < 2
