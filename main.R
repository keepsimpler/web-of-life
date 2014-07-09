#!/usr/bin/Rscript

## First step, generate random bipartite graphs with different structural properties.
#source('graphs.rewiring.r')  
load(file = 'data/graphs.20.RData')

sigma0 = 0.1
h0 = 0.01
## Second step, simulate the deterministic LV model until the feasible fixed point is approached
source('lv2.r')
graph.perfect.nested = graphs.20[[length(graphs.20)]]$B
gamma0.max = get.gamma0.max(graph.perfect.nested, sigma0 = sigma0)
sub.graphs.num = 24  # samping uniformly from the set of graphs
graphs.index = floor(seq(from = 1, to = length(graphs.20), length.out = sub.graphs.num))
graphs.20 = graphs.20[- graphs.index]

lv2.out = sim.lv2(graphs = graphs.20, h0 = h0, gamma0 = gamma0.max, isout = TRUE)

## Third step, the multivariate OU process is used to simulate the stochastic dyanmics in equilibrium
source('mou.r')

mou.out.withNstar = sim.mou.lv2(lv2.out = lv2.out, sigma0 = sigma0, steps = 20000, withNstar = TRUE)
vars = aaply(mou.out.withNstar, .margins = c(1, 3), var)
save(vars, file = 'data/vars-0.1-0.01.RData')

#mou.out.noNstar = sim.mou.lv2(lv2.out = lv2.out, sigma0 = sigma0, withNstar = FALSE)
#save(mou.out.noNstar, file = 'mou.out.noNstar-0.1.RData')

