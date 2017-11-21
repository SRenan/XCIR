# One day I will take the time to have an actual data structure instead
# of a gigantic data.table.


# data.table
# read_count_cutoff
# inactivation_cutoff
# escape_cutoff
## Maybe
# Cross pooled samples bi_allelic cutoff
# skewedSamples
# source file the data was read from
setClass("xcirObj",
         representation = representation(data = "data.table",
                                         info = "data.frame"))

setMethod("show", signature = signature(object = "xcirObj"), function(object){
  print(paste0("An xcirObj. use showInfo(object) to  get the filter values."))
  print(object@data)
  return(NULL)
})

#setMethod("as.data.frame", signature = signature(object = "xcirObj"), function(object){
#  ret <- as.data.frame(object@data)
#  return(ret)
#})
