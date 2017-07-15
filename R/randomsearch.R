randomsearch <- function(fn ,
                         lower, upper, initialpar = NULL,
                         iterations = 1000,
                         trace = TRUE,
                         maximize = FALSE){
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

  #関数に適用するためのmatrixを作る。
  #一様乱数の配列。
  randmat <- matrix(runif(length(upper) * iterations - length(initialpar)),
                    ncol = length(upper), byrow = TRUE)

  #1の行列作成
  mat1 <- matrix(rep(1, times = (length(upper) * iterations - length(initialpar))),
                 ncol = length(upper), byrow = TRUE)

  #初期値を作成
  initmat <- t(t(randmat) * (upper - lower) + t(mat1) * lower)
  initmat <- rbind(initialpar, initmat)

  #配列演算
  result <- apply(initmat, 1, fn)

  #最大化の場合
  if(maximize){result <- -result}

  #最小値計算
  bestvalue <- Inf
  result_log <- NULL
  for(i in 1 : length(result)){
    if(result[i] < bestvalue){
      bestvalue <- result[i]
      result_log[i] <- bestvalue
      bestpar <- initmat[i,]
    }else{
      result_log[i] <- result_log[i - 1]
    }
  }

  #最大化の場合
  if(maximize){
    result <- -result
    result_log  <- -result_log
    bestvalue <- -bestvalue
  }


  #戻り値
  ret <- list()
  ret$result <- result
  ret$result_log <- result_log
  ret$bestvalue <- bestvalue
  ret$bestpar <- bestpar
  ret$iterations <- iterations
  ret$maximize <- maximize
  class(ret) <- "randomsearch"

  return(ret)
}

summary.randomsearch <- function(data){
  cat("Random search optimization result.\n")
  cat("A function to be ")
  if(data$maximize){cat("maximized.\n\n")}else{cat("minimized.\n\n")}

  cat(paste0("Iterations size is ", data$iterations, "\n\n"))
  cat("Best parameters.\n")
  cat(data$bestpar)
  cat("\n\nBest value is ")
  cat(data$bestvalue)
}

plot.randomsearch <- function(data){
  plot(data$result_log, type = "l",
       xlab = "Iterations", ylab = "Value",
       main ="Random search optimization result")
}




