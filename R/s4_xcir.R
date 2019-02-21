# One day I will take the time to have an actual data structure instead
# of a gigantic data.table.


# data.table
# source file the data was read from
# read_count_cutoff
# het_cutoff
# escape_cutoff
# Genome version
## Maybe
# Cross pooled samples bi_allelic cutoff
# skewedSamples

# FUNCTIONS
# The read functions return xcirObj
# The annotation function takes either data.table or an unannotated xcirObj

################################################################################
# CLASS
################################################################################

#' An S4 class to hold XCIR objects
#'
#' @export
setClass("xcirObj",
         representation = representation(data = "data.table",
                                         info = "data.table"))
################################################################################
# CONSTRUCTOR
################################################################################
setMethod("initialize", "xcirObj", function(.Object, data = data.table(), info = data.table(), ...){
  obj <- callNextMethod(.Object, ...)
  obj@data <- data
  obj@info <- info
  return(obj)
})

################################################################################
# METHODS
################################################################################

setMethod("show", signature = signature(object = "xcirObj"), function(object){
  print(paste0("An xcirObj. use showInfo(object) to  get the filter values."))
  print(object@data)
  return(NULL)
})
setMethod("as.data.frame", signature = signature(x = "xcirObj"), function(x, ...){
  ret <- as.data.frame(x@data, ...)
  return(ret)
})


