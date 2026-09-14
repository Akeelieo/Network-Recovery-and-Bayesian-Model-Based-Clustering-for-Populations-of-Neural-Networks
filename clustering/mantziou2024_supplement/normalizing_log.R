#Function for dealing with underflow problem through normalizing (sum to 1) the logarithms
#that the membership function returns

logxpy_3 <- function(lx,ly,lz) {max(lx,ly,lz) + log1p(exp(-abs(lx-ly))+exp(-abs(lx-lz)))}

logxpy_2 <- function(lx,ly) {max(lx,ly) + log1p(exp(-abs(lx-ly)))}

logxpy_C <- function(l){
  clusters<-length(l)
  max_l<-max(l)
  sum_exp<-0
  for (i in 2:clusters){
    sum_exp<-sum_exp+exp(-abs(l[1]-l[i]))
  }
  return(max_l+log1p(sum_exp))
}
