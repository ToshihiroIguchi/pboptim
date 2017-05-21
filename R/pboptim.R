#Population Base Optimization
#pboptim package

library(DEoptim)
library(psoptim)
library(GA)

pboptim <- function(initialpar = NULL, fn ,lower, upper,
                    method = c("DEoptim", "pso", "GA"),
                    population = 20, generation = 10,
                    trace = TRUE, maxit = FALSE
                    ){
  result <- NULL

  #Differential Evolution Optimization
  result$DEoptim <- DEoptim(fn = fn, lower = lower, upper = upper,
                            control = DEoptim.control(
                              VTR = if(maxit){Inf}else{-Inf},
                              strategy = 2,
                              trace = trace,
                              initialpop = initialpar,
                              itermax = generation,
                              NP = population))

  #Particle Swarm Optimization
  #初期値は定義できない。
  result$pso <- psoptim(FUN = fn, xmin = lower, xmax = upper,
                        w = 0.9, c1 = 0.2, c2 = 0.2,
                        n = population, max.loop = generation,
                        vmax = c(4, 4), anim = FALSE, seed = rnorm(1))

  #Genetic Algorithms
  result$ga <- ga(type = c("real-valued"),
                  fitness = if(maxit){fn}else{function(x){-fn(x)}},
                  min = lower, max = upper,
                  popSize = population, maxiter = generation, run = generation,
                  suggestions = initialpar,
                  optimArgs = list(control = list(fnscale = if(maxit){1}else{-1})))

  return(result)

}

test <- pboptim(function(x){x[1]^2+x[2]^2},
                initialpar = matrix(rep(1,times = 40),ncol = 2),
                lower = c(-2,-2),upper = c(2,2))



#optimx(par, fn, gr=NULL, hess=NULL, lower=-Inf, upper=Inf,
#       method=c("Nelder-Mead","BFGS"), itnmax=NULL, hessian=FALSE,
#       control=list(),

#p1       p2         value fevals gevals niter convcode kkt1 kkt2
#BFGS         9.889487 5.038007  1.190725e+03     14      5    NA        0   NA   NA
#CG           9.890986 5.037998  1.190725e+03    307    101    NA        1   NA   NA
#Nelder-Mead  9.892127 5.038049  1.190725e+03     41     NA    NA        0   NA   NA
#L-BFGS-B     9.889511 5.038000  1.190725e+03     13     13    NA        0   NA   NA


#参考
#http://masaqol.hatenablog.com/entry/2016/07/31/233153
#https://stat.ethz.ch/R-manual/R-devel/library/stats/html/optim.html
#http://qiita.com/HirofumiYashima/items/c2e9ff4e3e9283ec1691
#


#参考
library(optimx)

loglik <- function(x, param) {
  -sum(dnbinom(x, size = param[1], mu = param[2], log = TRUE))
}

set.seed(71)
data <- rnbinom(500, size = 10, mu = 5)
initparam <- c(10, 5)
optimx(par = initparam, fn = loglik,
       x = data, control = list(all.methods = TRUE))
