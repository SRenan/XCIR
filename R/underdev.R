# This file for functions under development before they get added to the files where they belong

getZ <- function(xcirout){
  Zdt <- copy(xcirout)
  Zdt[, mup := mean(p_value), by = "GENE"]
  Zdt[, sdp := sd(p_value), by = "GENE"]
  Zdt[, Z := (p_value - mup)/sdp]

}
