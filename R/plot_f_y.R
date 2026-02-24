plot_f_y <- function(x, y, 
                           xlab = "f(s_t)", 
                           ylab = "Fitted values",
                           main = "Estimated s_t vs Fitted values",
                           span = 0.75,
                           point_alpha = 0.2,
                           point_col = "black",
                           line_col = "blue") {
  
  # 1. 散点图
  plot(x, y,
       xlab = xlab,
       ylab = ylab,
       main = main,
       pch  = 16,
       col  = adjustcolor(point_col, alpha.f = point_alpha)
  )
  
  # 2. loess 拟合
  lo_fit <- loess(y ~ x, span = span)
  
  # 3. 预测用于画曲线与置信区间
  x_grid <- seq(min(x, na.rm = TRUE),
                max(x, na.rm = TRUE),
                length.out = 200)
  
  lo_pred <- predict(lo_fit,
                     newdata = data.frame(x = x_grid),
                     se = TRUE)
  
  y_fit  <- lo_pred$fit
  y_se   <- lo_pred$se.fit
  y_upper <- y_fit + 1.96 * y_se
  y_lower <- y_fit - 1.96 * y_se
  
  # 4. 加平滑曲线
  lines(x_grid, y_fit,
        col = line_col,
        lwd = 2)
  
  # 5. 加置信区间带
  polygon(
    x = c(x_grid, rev(x_grid)),
    y = c(y_upper, rev(y_lower)),
    col = adjustcolor(line_col, alpha.f = 0.2),
    border = NA
  )
  
  # 6. 图例
  legend("topleft",
         legend = c("Data points", "Loess smooth", "95% CI"),
         pch    = c(16, NA, NA),
         lty    = c(NA, 1, NA),
         lwd    = c(NA, 2, NA),
         col    = c(adjustcolor(point_col, alpha.f = point_alpha),
                    line_col,
                    adjustcolor(line_col, alpha.f = 0.2)),
         bty    = "n",
         pt.cex = 1.2
  )
}
