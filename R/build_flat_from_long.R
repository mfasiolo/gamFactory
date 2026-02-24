# # ===============================================================
# # # Test case:
# set.seed(1)
# n_day <- 365
# y <- rnorm(n_day, 100, 5)
# n_time <- 3
# 
# # Step 1: Randomly assign each day a time_point ∈ {1, 2, 3}
# time_point <- sample(1:n_time, n_day, replace = TRUE)
# 
# # Step 2: Compute time_seq
# time_seq <- (seq_len(n_day) - 1L) * n_time + time_point
# 
# # Construct: long-format matrices with 1093 rows (2 rows fewer than 365*3 = 1095)
# n_row_long <- 1093
# 
# X_ori <- matrix(rnorm(n_row_long * 4, 15, 5), n_row_long, 4)
# colnames(X_ori) <- c("East", "South", "West", "North")
# 
# w_ori <- cbind(
#   intercept = rep(1, n_row_long),
#   time_diff = rnorm(n_row_long, 8, 2)
# )
# 
# # Apply the function
# res <- build_flat_from_long(y, time_seq, X_ori, w_ori)
# 
# # Check results
# res$n_rep        # should be 3
# res$t_last       # tail(time_seq, 1)
# dim(res$X_flat)  # 365 x (3*4) = 365 x 12
# dim(res$w_flat)  # 365 x (3*2) = 365 x 6
# 
# tail(res$X_flat, 2)
# tail(res$w_flat, 2)
# # ===============================================================


# ===============================================================
# Build flat matrices (per-day slicing) from long-format X_ori / w_ori
# ===============================================================
# Inputs:
#   y        : numeric vector, length = n_day
#   time_seq : integer vector, length = n_day, last element is t_last
#   X_ori    : numeric matrix/data.frame, long format ordered by day:
#              day1: t1..t_nrep, day2: t1..t_nrep, ..., dayD: t1..t_nrep
#   w_ori    : same long-format ordering as X_ori
#
# Outputs:
#   $dat     : data.frame(y = y, xw = I(cbind(X_flat, w_flat, times = time_seq)))
#   $X_flat  : n_day x (n_rep * ncol(X_ori))
#   $w_flat  : n_day x (n_rep * ncol(w_ori))
#   $t_last  : last element of time_seq
#   $n_rep   : repetition count (ceiling)
# ===============================================================
build_flat_from_long <- function(y, time_seq, X_ori, w_ori) {
  # ---------- 0) Basic checks ----------
  stopifnot(is.numeric(y), length(y) >= 1L)
  stopifnot(is.numeric(time_seq), length(time_seq) == length(y))
  n_day <- length(y)
  t_last <- tail(time_seq, 1L)
  
  # Coerce to numeric matrices
  X_ori <- as.matrix(X_ori); storage.mode(X_ori) <- "double"
  w_ori <- as.matrix(w_ori); storage.mode(w_ori) <- "double"
  nx <- nrow(X_ori); px <- ncol(X_ori)
  nw <- nrow(w_ori); pw <- ncol(w_ori)
  if (is.null(colnames(X_ori))) colnames(X_ori) <- paste0("X", seq_len(px))
  if (is.null(colnames(w_ori))) colnames(w_ori) <- paste0("w", seq_len(pw))
  
  # ---------- 1) Constraint on t_last ----------
  max_rows_long <- max(nx, nw)
  if (t_last > max_rows_long) {
    stop(sprintf("t_last (%d) > max(nrow(X_ori), nrow(w_ori)) = %d. Please check time_seq or long-table lengths.",
                 t_last, max_rows_long))
  }
  
  # ---------- 2) Decide n_rep & pad the tail with Inf to full daily slices ----------
  n_rep <- ceiling(max_rows_long / n_day)
  total_needed <- n_day * n_rep
  
  pad_matrix_to <- function(M, total_needed) {
    n_now <- nrow(M); p_now <- ncol(M)
    if (n_now < total_needed) {
      n_pad <- total_needed - n_now
      M <- rbind(M, matrix(Inf, nrow = n_pad, ncol = p_now,
                           dimnames = list(NULL, colnames(M))))
    } else if (n_now > total_needed) {
      # If longer, truncate the tail to n_day * n_rep
      M <- M[seq_len(total_needed), , drop = FALSE]
      warning("Input has more rows than needed; truncated to n_day * n_rep.")
    }
    M
  }
  
  X_pad <- pad_matrix_to(X_ori, total_needed)
  w_pad <- pad_matrix_to(w_ori, total_needed)
  
  # ---------- 3) Per-day slicing & flatten (t1..t_nrep in order) ----------
  #   Assumes row order = day1: t1..t_nrep, day2: t1..t_nrep, ...
  make_flat_per_day <- function(M, n_day, n_rep, base_names) {
    p <- ncol(M)
    flat <- matrix(NA_real_, nrow = n_day, ncol = p * n_rep)
    colnames(flat) <- as.vector(sapply(seq_len(n_rep), function(tt) paste0(base_names, "_t", tt)))
    rownames(flat) <- paste0("day", seq_len(n_day))
    
    for (d in seq_len(n_day)) {
      rows_d <- ((d - 1L) * n_rep + 1L):(d * n_rep)  # t1..t_nrep of current day
      # Concatenate column blocks t1, t2, ...
      flat[d, ] <- unlist(lapply(rows_d, function(rr) M[rr, , drop = FALSE]))
    }
    flat
  }
  
  X_flat <- make_flat_per_day(X_pad, n_day, n_rep, colnames(X_ori))
  w_flat <- make_flat_per_day(w_pad, n_day, n_rep, colnames(w_ori))
  
  # ---------- 4) Assemble xw & dat ----------
  xw <- cbind(X_flat, w_flat, times = time_seq)
  dat <- data.frame(y = as.numeric(y), xw = I(xw))
  
  # ---------- 5) Return ----------
  list(
    dat    = dat,
    X_flat = X_flat,
    w_flat = w_flat,
    t_last = t_last,
    n_rep  = n_rep
  )
}
