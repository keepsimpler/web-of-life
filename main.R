#!/usr/bin/Rscript

## First step, generate random bipartite graphs with different structural properties.
#source('graphs.rewiring.r')  
load(file = 'graphs.20.RData')

sigma0 = 0.1
## Second step, simulate the deterministic LV model until the feasible fixed point is approached
source('lv2.r')
graph.perfect.nested = graphs.20[[length(graphs.20)]]$B
gamma0.max = get.gamma0.max(graph.perfect.nested, sigma0 = sigma0)
sub.graphs.num = 18  # samping uniformly from the set of graphs
graphs.index = floor(seq(from = 1, to = length(graphs.20), length.out = sub.graphs.num))
graphs.20 = graphs.20[- graphs.index]

lv2.out = sim.lv2(graphs = graphs.20, h0 = 0, gamma0 = gamma0.max, isout = TRUE)

## Third step, the multivariate OU process is used to simulate the stochastic dyanmics in equilibrium
source('mou.r')

mou.out.withNstar = sim.mou.lv2(lv2.out = lv2.out, sigma0 = sigma0, withNstar = TRUE)
save(mou.out.withNstar, file = 'mou.out.withNstar-0.1.RData')

mou.out.noNstar = sim.mou.lv2(lv2.out = lv2.out, sigma0 = sigma0, withNstar = FALSE)
save(mou.out.noNstar, file = 'mou.out.noNstar-0.1.RData')

