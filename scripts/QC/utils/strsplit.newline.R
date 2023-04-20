

strsplit.newline <- function (x, nchar, split = '_', sep = '\n') 
{
  
  lapply(x, function(y) {
    
    arr <-  unlist(strsplit(y, split = ''))
    sequence <- unique(c(seq(0,length(arr), by = nchar), length(arr)))
    #if (length(sequence) > 2) browser()
    
    arr_sep <- unlist(lapply(seq_along(sequence)[-1], function(i) {
      sub_arr <- arr[sequence[(i-1)]:(sequence[i])]
      sub_arr <- ifelse(i != length(sequence) & i != 2, paste0(sub_arr[-1], collapse = ''), paste0(sub_arr, collapse = ''))
      sub_arr <- ifelse(i != length(sequence), paste0(sub_arr, sep, collapse = ''), sub_arr)
      return(sub_arr)
    }))
    
    return(paste(arr_sep, collapse = ''))
  
  })
  
  
}