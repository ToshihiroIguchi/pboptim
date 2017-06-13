#Population Base Optimization
#pboptim package

library(DEoptim)
library(pso)
library(GA)
library(copulaedas)


#Population based optimization
#DEO, PSO, GA, EDAのラッパー。
pboptim <- function(fn ,lower, upper, initialpar = NULL,
                    method = c("DEO", "PSO", "GA", "EDA"),
                    population = 20, generation = 10,
                    trace = TRUE, maximize = FALSE
                    ){

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

  #初期設定
  result <- NULL
  gen_list <- list(DEO = NULL, PSO = NULL, GA = NULL, EDA = NULL)
  #gen_df <- data.frame()


  #結果一覧保存用データフレーム
  sum_df <- NULL
  sum_df_names <- c("method", paste0("par", c(1 : vec_len)), "value", "time")

  #Differential Evolution Optimization
  if(!is.na(match("DEO", method))){
    t0 <- proc.time()
    cat("Differential Evolution Optimization \n\n")

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

    gen_list$DEO <- result$deo$member$bestvalit
    t <- proc.time() - t0

    #結果の整理
    sum_df0 <- as.data.frame(matrix(c(result$deo$optim$bestmem,
                                      result$deo$optim$bestval, t[3]), nrow = 1))
    sum_df0 <- data.frame(method = "DEO", sum_df0)
    names(sum_df0) <- sum_df_names
    sum_df <- rbind(sum_df, sum_df0)
  }

  #Particle Swarm Optimization
  if(!is.na(match("PSO", method))){
    t0 <- proc.time()
    cat("Particle Swarm Optimization \n\n")

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
                     trace = 1, REPORT = 1, trace.stats = 1))
    gen_list$PSO <- result$pso$stats$error
    t <- proc.time() - t0

    #結果の整理
    sum_df0 <- as.data.frame(matrix(c(result$pso$par,
                                      result$pso$value, t[3]), nrow = 1))
    sum_df0 <- data.frame(method = "PSO", sum_df0)
    names(sum_df0) <- sum_df_names
    sum_df <- rbind(sum_df, sum_df0)
  }

  #Genetic Algorithms
  if(!is.na(match("GA", method))){
    t0 <- proc.time()
    cat("Genetic Algorithms \n\n")
    result$ga <- ga(type = c("real-valued"),
                    fitness = if(maximize){fn}else{function(x){-fn(x)}},
                    min = lower, max = upper,
                    popSize = population, maxiter = generation, run = generation,
                    suggestions = initialpar,
                    optimArgs = list(control = list(fnscale = if(maxit){1}else{-1})),
                    monitor = FALSE)
    gen_list$GA <- result$ga@summary[,1]*if(maximize){1}else{-1}
    ga_best <- result$ga@fitnessValue*if(maximize){1}else{-1}
    t <- proc.time() - t0

    #結果の整理
    sum_df0 <- as.data.frame(matrix(c(result$ga@solution,
                                      ga_best, t[3]), nrow = 1))
    sum_df0 <- data.frame(method = "GA", sum_df0)
    names(sum_df0) <- sum_df_names
    sum_df <- rbind(sum_df, sum_df0)
  }

  #Estimation of Distribution Algorithm
  if(!is.na(match("EDA", method))){
    t0 <- proc.time()
    cat("Estimation of Distribution Algorithm \n\n")
    setMethod("edaReport", "EDA", edaReportSimple)
    setMethod("edaTerminate", "EDA", edaTerminateMaxGen)
    GCEDA <- CEDA(copula = "normal", margin = "norm",
                  popSize = population, maxGen = generation)
    GCEDA@name <- "GCEDA"
    eda_cap <- capture.output(result$eda <- edaRun(GCEDA, fn, lower, upper))

    #画面出力から推移読み取り
    #最初と最後の位置
    pos1 <- regexpr("Generation", eda_cap[2])[1]+11
    pos2 <- regexpr("Minimum", eda_cap[2])[1]+7

    #切り出し
    best_eda <- NULL
    for(i in 3:length(eda_cap)){
      best_eda <- c(best_eda,
                    as.numeric(substr(eda_cap[i], pos1, pos2)))
    }
    gen_list$EDA <- best_eda
    t <- proc.time() - t0

    #結果の整理
    sum_df0 <- as.data.frame(matrix(c(result$eda@bestSol,
                                      result$eda@bestEval, t[3]), nrow = 1))
    sum_df0 <- data.frame(method = "EDA", sum_df0)
    names(sum_df0) <- sum_df_names
    sum_df <- rbind(sum_df, sum_df0)
  }

  #まとめが存在するか確認
  if(is.null(sum_df)){
    stop("The result does not exist. Please review method.")
  }

  #まとめのmethodを文字型に変更
  sum_df$method <- as.character(sum_df$method)

  #最適な結果
  sum_best_df <- sum_df[order(sum_df$value * if(maximize){-1}else{1}), ][1, ]
  result$bestmethod <- sum_best_df[1, 1]
  result$bestpar <- sum_best_df[1, c(2: (vec_len + 1))]
  result$bestvalue <- sum_best_df[1, (2 + vec_len)]

  #設定値の保管
  result$setting <- list(fn = fn, lower = lower, upper = upper,
                         population = population, generation = generation,
                         maximize = maximize)

  #結果の格納と戻り値
  result$maximize <- maximize
  result$gen_list <- gen_list
  result$summary <- sum_df
  class(result) <- "pboptim"
  return(result)
}

#結果表示
summary.pboptim <- function(obj){

  cat("A function to be ")
  if(obj$setting$maximize){cat("maximized.\n")}else{cat("minimized.\n\n")}
  print(obj$setting$fn)
  cat("\n")

  cat("The population size is ")
  cat(as.integer(obj$setting$population))
  cat("\n")

  cat("The generation size is ")
  cat(as.integer(obj$setting$generation))
  cat("\n\n")

  cat("Results.\n")
  print(obj$summary)
  cat("\n")

  cat("Best method is ")
  cat(obj$bestmethod)
  cat("\n\n")

  cat("Best parameters. \n")
  print(obj$bestpar)
  cat("\n")

  cat("Best value is ")
  cat(obj$bestvalue)
  cat("\n\n")
}


plot.pboptim <- function(obj){
  #plot用のデータ作成
  gl <- obj$gen_list
  plot_df <- data.frame(generation = c(1:obj$setting$generation))
  if(!is.null(gl$DEO)){plot_df <- data.frame(plot_df, DEO = gl$DEO)}
  if(!is.null(gl$PSO)){plot_df <- data.frame(plot_df, PSO = gl$PSO)}
  if(!is.null(gl$GA)){plot_df <- data.frame(plot_df, GA = gl$GA)}
  if(!is.null(gl$EDA)){plot_df <- data.frame(plot_df, EDA = gl$EDA)}

  #線の色と種類定義
  data_n <- c(1:length(plot_df[1,]))
  cols <- c("black", "red", "blue", "green")[data_n]
  ltys <- c(1,2,3,4)[data_n]

  #凡例の位置
  legend_pos <- if(obj$setting$maximize){"bottomright"}else{"topright"}

  #プロット
  matplot(plot_df[,1], plot_df[,-1],
          type = "l", xlab = "Generation", ylab = "Value",
          col = cols, lty = ltys)
  legend(legend_pos, legend = names(plot_df[,-1]), col = cols, lty = ltys)

}


set.seed(108)
test <- pboptim(function(x){x[1]^2+x[2]^2},
                method = c("DEO", "PSO"),
                initialpar = c(1,1),
                lower = c(-2,-2),upper = c(2,2))


plot(test)






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
