library(plyr)
library(reshape2)
library(dplyr)
library(ggplot2)
library(grid)  # need viewport() function

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

get.vars.mean <- function(vars, steps = 10000) {
  #vars = aaply(mou.out, .margins = c(1, 3), var)
  range.of.equli = (5000 + 2) : (steps + 1)
  vars = vars[, range.of.equli, , ]
  vars.mean = aaply(vars, .margins = c(1, 3, 4), mean)
  vars.mean
}

get.vars.sum <- function(vars.mean, isLong = TRUE) {
  vars.sum = adply(vars.mean, .margins = c(1), function(vars) {
    c(varsum = sum(vars), selfvarsum = sum(diag(vars)), covsum = sum(vars) - sum(diag(vars)),
      syn =  sum(vars) / sum(sqrt(diag(vars)))^2)   #maxvarsum = sum(sqrt(diag(vars)))^2, 
  })
  if (isLong) {
    vars.sum = melt(vars.sum, id.vars = c('X1'), measure.vars = c('varsum', 'selfvarsum', 'covsum', 'syn'))
  }
  vars.sum
}

get.vars.self <- function(vars.mean) {
  vars.self = adply(vars.mean, .margins = c(1), function(one) {
    data.frame(selfvar = diag(one))
  })
  colnames(vars.self)[1] <- 'gindex'
  vars.self
}

get.Nstars <- function(lv2.out) {
  ldply(1:length(lv2.out), function(i) {
    data.frame(gindex = i, Nstar = lv2.out[[i]]$Nstar, sindex = 1:length(lv2.out[[i]]$Nstar))
  })  
}
get.Nstars.mean <- function(lv2.out) {
  ldply(lv2.out, function(one) {
    c(Nstars.mean = mean(one$Nstar))
  })
}

get.Theta.sum <- function(lv2.out) {
  ldply(lv2.out, function(one) {
    Theta = - diag(1 / one$Nstar) %*% one$Phi
    sum(solve(Theta))
  })
}

get.Phi.lev <- function(lv2.out) {
  ldply(lv2.out, function(one) {
    eigs = eigen(one$Phi)$values
    lev = max(eigs)
  })  
}

display.vars.sum <- function(vars.sum) {
  #vars.sum = vars.sum %.% filter(X1 != 12 & variable != 'syn')
  #pdf("tuzhongtu.pdf",width = 4, height = 4)  # output to pdf file
  p <- ggplot(data = vars.sum %.% filter(X1 != 12 & variable != 'syn'), 
              aes(x = factor(X1), y = value, colour = factor(variable), group = factor(variable))) + 
    geom_line() + geom_point() #+ theme_bw()   # + scale_y_log10()
  print(p)
  subvp <- viewport(width = 0.5, height = 0.4, x = 0.4, y = 0.7, )
  q <- ggplot(data = vars.sum %.% filter(X1 != 12 & variable == 'syn'), 
              aes(x = factor(X1), y = value, colour = factor(variable), group = factor(variable))) + 
    geom_line() + geom_point() #+ theme_bw()
  print(q)  #, vp = subvp
  #dev.off()  # close pdf file as output device
}

display.si <- function(vars.self, Nstars) {
  vars.self$Nstar = Nstars$Nstar
  vars.self$sindex = Nstars$sindex
  p <- ggplot(data = vars.self, aes(x = Nstar, y = selfvar)) + geom_point() + 
    facet_wrap(~ gindex, ncol = 6, scales = 'free') + 
    scale_y_log10('log(species variance)') + scale_x_log10('log(species abundance)') +
    geom_smooth(method = 'lm', size = 1) +
    #  geom_text(data=eq, aes(x=9, y=1,  label=V1),parse = TRUE, inherit.aes=FALSE) +
    #  theme(strip.text.x = element_text(size = 14, angle = 0)) +
    theme_bw()
  print(p)
  
  vars.self.2 = vars.self %.% filter(gindex != 12)
  vars.self.2$si = vars.self.2$Nstar / sqrt(vars.self.2$selfvar)
  p <- ggplot(data = vars.self.2, 
              aes(x = factor(gindex), y = si, colour = factor(sindex), group = factor(sindex))) + 
    geom_line() + geom_point() #+ theme_bw()   # + scale_y_log10()
  print(p)
}

display.st <- function(Nstars.mean, vars.sum) {
  
}

# ldply(ls(), function(objname) {c(objname,object.size(get(objname)))})  # get size of all objects in current environment