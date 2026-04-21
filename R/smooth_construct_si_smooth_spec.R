#'
#' Single index effects for mgcv
#' 
#' @name smooth.construct.si.smooth.spec
#' @rdname smooth.construct.si.smooth.spec
#' @importFrom MASS Null
#' @export
#'
smooth.construct.si.smooth.spec <- function(object, data, knots){
  
  # Most information on single index matrix and penalty is inside "si" list.
  si <- object$xt$si
  if( is.null(si) ){ si <- object$xt$si <- list() }
  
  # Inner model matrix (to be projected via single index)
  Xi <- data[[object$term]]
  
  # Need to center Xi and save colMeans because we need to subtract it when predicting using new data
  Xi <- scale(Xi, scale = FALSE)
  si$xm <- attr(Xi, "scaled:center")
  
  di <- ncol( Xi )
  n <- nrow( Xi )
  
  # Dealing with inner penalty
  Si <- si$S
  no_pen <- is.null(Si) && is.null(si$pord)
  if( no_pen ){ # Case [a] no penalisation
    si$X <- Xi
    si$B <- diag(nrow = ncol(Xi))
    si$rank <- 0 
  } else {
    if( is.null(Si) ){ # Case [b] "P-splines" penalty
      Si <- .psp(d = di, ord = si$pord)
      rankSi <- ncol(Xi) - si$pord
    } else { # Case [c] custom penalty Si
      rankSi <- rankMatrix(Si)
    }
    # Reparametrise Xi so that the penalty on the single index vector is diagonal
    si <- append(si, gamFactory:::.diagPen(X = Xi, S = Si, r = rankSi))
  }
  
  # alpha is vector of inner coefficients, si$alpha is a vector of initial values for it.
  # alpha0 is an offset such that the full_alpha = alpha + alpha0
  if( is.null(si$a0) ){
    if( no_pen ){
      si$a0 <- rep(0, di)
    } else {
      si$a0 <- rep(1, di)
    }
  }
  if( is.null(si$alpha) ){ 
    if( is.null(si$a0) || all(si$a0 == 0) ){
      si$alpha <- rep(1, di) 
    } else {
      si$alpha <- rep(0, di) 
    }
  }
  # Reparametrise and then impose that variance should be 1
  si$alpha <- solve(si$B, si$alpha)
  si$a0 <- solve(si$B, si$a0)
  
  positive_si <- isTRUE(si$positive_si) 
  
  if (positive_si) {
    alpha_outer <- pmax((si$alpha + si$a0), 1e-4)
    tmp <- sd(si$X %*% alpha_outer)
    alpha_outer_std <- alpha_outer / tmp
    
    si$alpha <- log(alpha_outer_std)
    si$a0 <- rep(0, di)
    
    ax <- drop( si$X %*% exp(si$alpha + si$a0) )
    
    # # 2. 【求初始方差】：计算在这个安全起点下，单指数向量的标准差
    # tmp <- sd(si$X %*% exp(si$alpha + si$a0))
    # 
    # # 3. 【核心修正】：在对数空间完成标准化
    # # 等价于让 Outer 除以 tmp。经过这一步，si$alpha 变成了 -log(tmp)
    # si$alpha <- si$alpha - log(tmp)
    # 
    # # 4. 【生成基底锚点】：正推算出真实的 ax，用于放置 B-spline 节点
    # # 经过上一步的修正，这里算出来的 ax 方差严格等于 1，完美！
    # ax <- drop( si$X %*% exp(si$alpha + si$a0) )
    
  } else {
    # 无约束逻辑保持原样
    tmp <- sd(si$X %*% (si$alpha + si$a0))
    si$alpha <- si$alpha / tmp
    si$a0 <- si$a0 / tmp
    ax <- drop( si$X %*% (si$alpha + si$a0) )
  }
  # ... 后面接 B-spline 构建
  
  data[[object$term]] <- ax
  
  # Construct the B-splines corresponding to the outer smooth effect 
  out <- .build_nested_bspline_basis(object = object, data = data, knots = knots, si = si)
  
  # Add inner penalty matrix (diagonalised and padded with zeros corresponding to the outer coefficients)
  if( !no_pen ){
    dsmo <- out$bs.dim - di
    si <- out$xt$si
    out$S[[2]] <- rbind(cbind(si$S, matrix(0, di, dsmo)),
                        cbind(matrix(0, dsmo, di), matrix(0, dsmo, dsmo)))
    out$null.space.dim <- out$null.space.dim + (out$bs.dim - si$rank)
    out$rank <- c(out$rank, si$rank)
  }
  
  class(out) <- c("si", "nested")
  return( out )
} 
