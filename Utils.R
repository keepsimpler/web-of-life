#means = aaply(A, .margins = c(2, 3), mean)
#vars = aaply(A, .margins = c(2), var)

cummean2.rev <- function(A) {
  B = aaply(A, .margins = c(2), function(onecol) {
    cumsum(rev(onecol)) / seq_along(rev(onecol))
  })  
  t(B)
}
cummean2 <- function(A) {
  B = aaply(A, .margins = c(2), function(onecol) {
    cumsum(onecol) / seq_along(onecol)
  })  
  t(B)
}
cummean3.rev <- function(A) {
  B = apply(A, MARGIN = c(2, 3), function(onecol) {
    cumsum(rev(onecol)) / seq_along(rev(onecol))
  })
  #aperm(B, perm = c(3, 1, 2))
}
cummean3 <- function(A) {
  B = apply(A, MARGIN = c(2, 3), function(onecol) {
    cumsum(onecol) / seq_along(onecol)
  })
  #aperm(B, perm = c(3, 1, 2))
}


dim.reduct <- function(A) {
  B = apply(A, MARGIN = c(1), function(oneslice) {
    one = c(diag(oneslice), oneslice[lower.tri(oneslice)])
  })
  t(B)
}

get.vars <- function(mou.out, steps = 10000) {
  vars = aaply(mou.out, .margins = c(1, 3), var)
  range.of.equli = (steps / 2 + 1) : (steps + 1)
  vars = vars[, range.of.equli, , ]
  vars.mean = aaply(vars, .margins = c(1, 3, 4), mean)
  vars.mean
}

get.vars.sum <- function(vars.mean) {
  adply(vars.mean, .margins = c(1), function(vars) {
    c(varsum = sum(vars), selfvarsum = sum(diag(vars)), covsum = sum(vars) - sum(diag(vars)), maxvarsum = sum(sqrt(diag(vars)))^2)
  })
}

get.means.mean <- function(lv2.out) {
  ldply(lv2.out, function(one) {
    mean(one$Nstar)
  })
}
