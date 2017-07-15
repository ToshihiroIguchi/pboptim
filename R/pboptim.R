#Population based search
#pboptim package


#DEO, PSO, GA, SOMAのラッパー。
pboptim <- function(fn ,lower, upper, initialpar = NULL,
                    method = c("DEO", "PSO", "GA", "RS"),
                    population = 20, generation = 10,
                    trace = TRUE, maximize = FALSE,
                    DEOcontrol = list(
                      strategy = 2, bs = FALSE, CR = 0.5, F = 0.8, p = 0.2, c = 0),
                    PSOcontrol = list(
                      s = floor(10+2*sqrt(length(pso_par))),
                      k = 3, w = 0.7213, c.p = 1.193, type = "SPSO2007"),
                    GAcontrol = list(
                      pcrossover = 0.8,
                      pmutation = 0.1,
                      elitism = base::max(1, round(population * 0.05))
                    ),
                    SOMAcontrol = list(
                      pathLength = 3,
                      stepLength = 0.11,
                      perturbationChance = 0.1
                    )
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
  gen_list <- list(DEO = NULL, PSO = NULL, GA = NULL, SOMA = NULL)
  fns <- if(maximize){-1}else{1}
  fn2 <- if(!maximize){fn}else{function(x){-fn(x)}}

  #結果一覧保存用データフレーム
  sum_df <- NULL
  sum_df_names <- c("method", paste0("par", c(1 : vec_len)), "value", "time")

  #Differential Evolution Optimization
  if(!is.na(match("DEO", method))){
    t0 <- proc.time()
    cat("Differential Evolution Optimization \n\n")

    #DEoptim用の初期値を作成。
    if(!is.null(initialpar)){
      deo_par_rand <- c((upper - lower)*runif(vec_len * (population -1))
                        + rep(lower, times = (population -1)))

      deo_par <- matrix(c(initialpar, deo_par_rand), ncol = vec_len, byrow = TRUE)

    }else{
      deo_par <- NULL
    }

    #DEoptimによる計算
    result$deo <- DEoptim(fn = fn2,
                          lower = lower, upper = upper,
                          control = DEoptim.control(
                            trace = trace,
                            initialpop = deo_par,
                            itermax = generation,
                            NP = population,

                            strategy = DEOcontrol$strategy,
                            bs = DEOcontrol$bs,
                            CR = DEOcontrol$CR,
                            F = DEOcontrol$F,
                            p = DEOcontrol$p,
                            c = DEOcontrol$c
                            ))

    gen_list$DEO <- result$deo$member$bestvalit * fns
    t <- proc.time() - t0

    #結果の整理
    sum_df0 <- as.data.frame(matrix(c(result$deo$optim$bestmem,
                                      result$deo$optim$bestval * fns, t[3]), nrow = 1))
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
                     fnscale = fns,
                     trace = if(trace){1}else{-1},

                     REPORT = 1,
                     trace.stats = trace,

                     s = PSOcontrol$s,
                     k = PSOcontrol$k,
                     w = PSOcontrol$w,
                     c.p = PSOcontrol$c.p,
                     type = PSOcontrol$type

                     #S=12, K=3, p=0.2297, w0=0.7213, w1=0.7213, c.p=1.193, c.g=1.193
                     #v.max=NA, d=5.657, vectorize=FALSE, hybrid=off

                     ))

    gen_list$PSO <- result$pso$stats$error * fns
    t <- proc.time() - t0

    #結果の整理
    sum_df0 <- as.data.frame(matrix(c(result$pso$par,
                                      result$pso$value * fns, t[3]), nrow = 1))
    sum_df0 <- data.frame(method = "PSO", sum_df0)
    names(sum_df0) <- sum_df_names
    sum_df <- rbind(sum_df, sum_df0)
  }

  #Genetic Algorithms
  if(!is.na(match("GA", method))){
    t0 <- proc.time()
    cat("Genetic Algorithms \n\n")
    result$ga <- ga(type = c("real-valued"),
                    fitness = function(x){-fn2(x)},
                    min = lower, max = upper,
                    popSize = population, maxiter = generation, run = generation,
                    suggestions = initialpar,
                    optimArgs = list(control = list(fnscale = if(maxit){1}else{-1})),
                    monitor = FALSE,

                    pcrossover = GAcontrol$pcrossover,
                    pmutation = GAcontrol$pmutation,
                    elitism = GAcontrol$elitism

                    )
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

  #Self-Organising Migrating Algorithm
  if(!is.na(match("SOMA", method))){
    t0 <- proc.time()
    cat("Self-Organising Migrating Algorithm \n\n")
    cat("Caution! The initial value is not reflected. \n\n")

    result$soma <- soma(fn2,
                        bounds = list(min = lower, max = upper),
                        options = list(minRelativeSep = -Inf, minAbsoluteSep = -Inf,
                                       nMigrations = generation,
                                       populationSize = population,

                                       pathLength = SOMAcontrol$pathLength,
                                       stepLength = SOMAcontrol$stepLength,
                                       perturbationChance = SOMAcontrol$perturbationChance
                                       ))

    gen_list$SOMA <- result$soma$history * fns
    soma_best <- result$soma$cost[result$soma$leader] * fns
    t <- proc.time() - t0

    #結果の整理
    sum_df0 <- as.data.frame(matrix(c(result$soma$population[, result$soma$leader],
                                      soma_best, t[3]), nrow = 1))
    sum_df0 <- data.frame(method = "SOMA", sum_df0)
    names(sum_df0) <- sum_df_names
    sum_df <- rbind(sum_df, sum_df0)
  }

  #Random Search
  if(!is.na(match("RS", method))){
    t0 <- proc.time()
    cat("Random Search \n\n")

    #fn ,
    #lower, upper, initialpar = NULL,
    #iterations = 1000,
    #trace = TRUE,
    #maximize = FALSE){

    #randomsearchによる計算
    result$rs <- randomsearch(fn = fn, lower = lower, upper = upper,
                              initialpar = initialpar, maximize = maximize,
                              iterations = population * generation)

    gen_list$RS <- result$rs$result_log[
      seq(population, population * generation, population)]
    t <- proc.time() - t0

    #結果の整理
    sum_df0 <- as.data.frame(matrix(c(result$rs$bestpar,
                                      result$rs$bestvalue, t[3]), nrow = 1))
    sum_df0 <- data.frame(method = "RS", sum_df0)
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
  sum_best_df <- sum_df[order(sum_df$value * fns), ][1, ]
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
  if(!is.null(gl$SOMA)){plot_df <- data.frame(plot_df, SOMA = gl$SOMA)}
  if(!is.null(gl$RS)){plot_df <- data.frame(plot_df, RS = gl$RS)}

  #線の色と種類定義
  data_n <- c(1:length(plot_df[1,]))
  cols <- rainbow(length(plot_df[1,]))

  ltys <- data_n

  #凡例の位置
  legend_pos <- if(obj$setting$maximize){"bottomright"}else{"topright"}

  #プロット
  matplot(plot_df[,1], plot_df[,-1],
          type = "l", xlab = "Generation", ylab = "Value",
          col = cols, lty = ltys,
          main = "Population based search")

  #凡例表示
  legend(legend_pos, legend = names(plot_df)[-1], col = cols, lty = ltys)

}


