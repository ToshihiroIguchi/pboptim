#Population Base Optimization
#pboptim package

library(DEoptim)
library(pso)
library(GA)
library(copulaedas)

pboptim <- function(fn ,lower, upper, initialpar = NULL,
                    method = c("deo", "pso", "GA","edm"),
                    population = 20, generation = 10,
                    trace = TRUE, maximize = FALSE
                    ){
  #初期設定
  result <- NULL

  #入力値のチェック
  #上限と下限のベクトルの長さ
  if(length(lower) != length(upper)){
    stop("The lengths of lower and upper vectors are different.")
  }
  vec_len <- length(lower)

  #上限と下限の範囲が逆転していないか。
  if(min(upper - lower) < 0){
    stop("Upper is smaller than lower.")
  }

  #初期値がNULLでないとき、初期値の長さが上下限の長さと一致しているか。
  if(!is.null(initialpar) && length(initialpar) != vec_len){
    stop("The length of the initial value is
         different from the length of the upper and lower limit values.")
  }


  #Differential Evolution Optimization
  cat("Differential Evolution Optimization \n")

  #DEoptim用の初期値を作成。
  if(!is.null(initialpar)){
    deo_par_rand <- rep((upper - lower)*runif(vec_len) + lower, times = (population -1))
    deo_par <- matrix(c(initialpar, deo_par_rand), ncol = vec_len)
  }else{
    deo_par <- NULL
  }

  #DEoptimによる計算
  result$deo <- DEoptim(fn = fn, lower = lower, upper = upper,
                            control = DEoptim.control(
                              VTR = if(maximize){Inf}else{-Inf},
                              strategy = 2,
                              trace = trace,
                              initialpop = deo_par,
                              itermax = generation,
                              NP = population))

  #Particle Swarm Optimization
  cat("Particle Swarm Optimization \n")

  #初期値を作成
  if(!is.null(initialpar)){
    pso_par <- initialpar
  }else{
    pso_par <- rep(NA, times = vec_len)
  }

  #psoptimによる計算
  result$pso <- psoptim(
    par = pso_par,
    fn = fn, lower = lower, upper = upper,
    control = list(maxit = generation, s = population,
                   fnscale = if(maximize){-1}else{1},
                   trace = 1, REPORT = 1, trace.stats = 1)
  )


  #Genetic Algorithms
  cat("Genetic Algorithms \n")
  result$ga <- ga(type = c("real-valued"),
                  fitness = if(maximize){fn}else{function(x){-fn(x)}},
                  min = lower, max = upper,
                  popSize = population, maxiter = generation, run = generation,
                  suggestions = initialpar,
                  optimArgs = list(control = list(fnscale = if(maxit){1}else{-1})),
                  monitor = FALSE)


  #Estimation of Distribution Algorithm
  cat("Estimation of Distribution Algorithm \n")
  setMethod("edaTerminate", "EDA", edaTerminateMaxGen)
  GCEDA <- CEDA(copula = "normal", margin = "norm",
               popSize = population, maxGen = generation)
  GCEDA@name <- "GCEDA"
  result$eda <- edaRun(GCEDA, fn, lower, upper)


  #結果の整理
  #世代と推移

  result$eva <- data.frame(
    DEO = result$deo$member$bestvalit,
    PSO = result$pso$stats$error,
    GA = result$ga@summary[,1]*if(maximize){1}else{-1}

  )



  return(result)

}

test <- pboptim(function(x){x[1]^2+x[2]^2},
                initialpar = c(1,1),
                lower = c(-2,-2),upper = c(2,2))


par(mfrow = c(1,1))









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
