
replace_nas_by_row_mean <- function(met_mat)
{
  met_mat_imp <- apply(met_mat, MARGIN = 1, FUN = function(x) { 
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x
  })
  
}


replace_nas_by_column_mean <- function(met_mat)
{
  met_mat_imp <- apply(met_mat, MARGIN = 2, FUN = function(x) { 
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x
  })
  
}
