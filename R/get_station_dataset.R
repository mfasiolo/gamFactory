#' !!!暂时无法使用，待修改
#' 意大利测试数据专用函数
#' 用于数据准备
#' 用于聚合站点的死亡、污染物与气象为日级 tibble（仅需站点代码）
#'
#' 函数只从包内 \code{data/} 懒加载数据，不允许传入路径。
#' 兼容以下数据布局：
#' \itemize{
#'   \item \code{data/list_dati_hour.rda} 内部对象名为 \code{dat.hour}（或 \code{list_dati_hour}，函数会重命名为 \code{dat.hour}）。
#'   \item \code{data/dat_era5.land.rda} 内部对象为 \code{dat_era5.land}（推荐，命名列表），
#'   或者 \code{data/} 下分散着多个 \code{*.land} 对象（如 \code{t2m.land}），函数会自动批量加载并合并为一个列表。
#' }
#'
#' @param station_code 字符串，站点代码（需与 \code{death_pop}/\code{dat.hour}/ERA5 子表列名一致）
#' @param tz 时区字符串，用于 ERA5 时间（默认 \code{"Europe/Rome"}）
#'
#' @return 一个 \code{tibble}，包含列：\code{time_index}, \code{Date}, \code{deaths}，以及各污染物/气象变量的日均
#' @examples
#' \dontrun{
#'   df <- gamFactory::get_station_dataset("TN_FU_4")
#' }
#' @export
#' @importFrom rlang .data
get_station_dataset <- function(station_code, tz = "Europe/Rome") {
  
  # ---------- 0) 懒加载包内数据 ----------
  env <- environment()
  
  # 0.1 death_pop
  utils::data("death_pop", package = "gamFactory", envir = env)
  if (!exists("death_pop", envir = env, inherits = FALSE))
    stop("Could not load 'death_pop' from package data.")
  
  # 0.2 dat.hour（文件名可能叫 list_dati_hour.rda）
  if (!exists("dat.hour", envir = env, inherits = FALSE)) {
    suppressWarnings(try(utils::data("dat.hour", package = "gamFactory", envir = env), silent = TRUE))
  }
  if (!exists("dat.hour", envir = env, inherits = FALSE)) {
    suppressWarnings(try(utils::data("list_dati_hour", package = "gamFactory", envir = env), silent = TRUE))
    # 若对象名仍为 list_dati_hour，则重命名为 dat.hour
    if (exists("list_dati_hour", envir = env, inherits = FALSE) &&
        !exists("dat.hour", envir = env, inherits = FALSE)) {
      dat.hour <- get("list_dati_hour", envir = env, inherits = FALSE)
      assign("dat.hour", dat.hour, envir = env)
      rm(list = "list_dati_hour", envir = env)
    }
  }
  if (!exists("dat.hour", envir = env, inherits = FALSE))
    stop("Could not load 'dat.hour' from package data (tried 'dat.hour' and 'list_dati_hour').")
  
  # 0.3 dat_era5.land：优先加载聚合后的列表；否则批量加载 *.land 并合并
  if (!exists("dat_era5.land", envir = env, inherits = FALSE)) {
    suppressWarnings(try(utils::data("dat_era5.land", package = "gamFactory", envir = env), silent = TRUE))
  }
  if (!exists("dat_era5.land", envir = env, inherits = FALSE)) {
    # 枚举包内可用数据集名称，抓取 *.land
    avail <- try(utils::data(package = "gamFactory")$results[, "Item"], silent = TRUE)
    if (!inherits(avail, "try-error")) {
      land_items <- avail[grepl("\\.land$", avail)]
      if (length(land_items)) {
        suppressWarnings(try(utils::data(list = land_items, package = "gamFactory", envir = env), silent = TRUE))
        land_objs <- Filter(function(x) exists(x, envir = env, inherits = FALSE), land_items)
        if (length(land_objs)) {
          dat_era5.land <- mget(land_objs, envir = env, inherits = FALSE)
          names(dat_era5.land) <- sub("\\.land$", "", names(dat_era5.land))
          assign("dat_era5.land", dat_era5.land, envir = env)
          # 可选：清理分散对象，避免污染命名空间
          # rm(list = land_objs, envir = env)
        }
      }
    }
  }
  if (!exists("dat_era5.land", envir = env, inherits = FALSE))
    stop("Could not load ERA5 data: neither 'dat_era5.land' nor '*.land' datasets found in package data.")
  
  death_pop     <- get("death_pop",     envir = env, inherits = FALSE)
  dat.hour      <- get("dat.hour",      envir = env, inherits = FALSE)
  dat_era5.land <- get("dat_era5.land", envir = env, inherits = FALSE)
  
  # ---------- 1) 小工具：规范日期列 ----------
  normalize_Date_col <- function(df) {
    if (is.null(df)) return(tibble::tibble(Date = as.Date(character())))
    if (!("Date" %in% names(df))) {
      cand <- intersect(c("Date", "date", "DATE", "Data", "data"), names(df))
      if (length(cand)) {
        df <- dplyr::rename(df, Date = .data[[cand[1]]])
      } else {
        df$Date <- as.Date(character())
      }
    }
    df$Date <- as.Date(df$Date)
    df
  }
  
  # ---------- 2) 死亡数据（日汇总） ----------
  death_one <- death_pop |>
    dplyr::filter(.data$Station == station_code) |>
    dplyr::mutate(Date = as.Date(.data$Date)) |>
    dplyr::group_by(.data$Date) |>
    dplyr::summarise(deaths = sum(.data$Deaths, na.rm = TRUE), .groups = "drop")
  
  # ---------- 3) 污染物小时数据 → 日均 ----------
  pollutant_names <- c("NO", "NO2", "NOx", "PM10", "PM2.5", "O3")
  
  pollutant_list <- purrr::map(pollutant_names, function(pname) {
    df <- dat.hour[[pname]]
    if (is.null(df)) {
      warning(sprintf("Pollutant '%s' not found.", pname))
      return(NULL)
    }
    if (!(station_code %in% names(df))) {
      warning(sprintf("Station '%s' not in dat.hour[['%s']].", station_code, pname))
      return(NULL)
    }
    date_col <- intersect(c("date", "Date", "DATE", "Data", "data"), names(df))
    if (length(date_col) == 0L) {
      warning(sprintf("No date column in dat.hour[['%s']].", pname))
      return(NULL)
    }
    date_col <- date_col[1]
    
    df |>
      dplyr::transmute(
        Date  = as.Date(.data[[date_col]]),
        value = .data[[station_code]]
      ) |>
      dplyr::group_by(.data$Date) |>
      dplyr::summarise("{pname}" := mean(.data$value, na.rm = TRUE), .groups = "drop")
  }) |>
    purrr::discard(is.null)
  
  pol_daily <- if (length(pollutant_list)) {
    purrr::reduce(
      pollutant_list,
      ~ dplyr::full_join(.x, .y, by = "Date"),
      .init = tibble::tibble(Date = as.Date(character()))
    ) |>
      dplyr::arrange(.data$Date)
  } else {
    tibble::tibble(Date = as.Date(character()))
  }
  
  # ---------- 4) ERA5 气象数据（列表每项一个变量）→ 日均 ----------
  era5_list <- purrr::map(names(dat_era5.land), function(vname) {
    df <- dat_era5.land[[vname]]
    if (is.null(df)) return(NULL)
    if (!(station_code %in% names(df))) return(NULL)
    
    # 这里假设 df$date 可被 as.POSIXct 正确解析；如原始即 Date，可改为 as.Date(df$date)
    df |>
      dplyr::mutate(Date = as.Date(lubridate::with_tz(as.POSIXct(.data$date, tz = tz), tzone = tz))) |>
      dplyr::transmute(Date, value = .data[[station_code]]) |>
      dplyr::group_by(.data$Date) |>
      dplyr::summarise("{vname}" := mean(.data$value, na.rm = TRUE), .groups = "drop")
  }) |>
    purrr::discard(is.null)
  
  met_daily <- if (length(era5_list)) {
    purrr::reduce(
      era5_list,
      ~ dplyr::full_join(.x, .y, by = "Date"),
      .init = tibble::tibble(Date = as.Date(character()))
    ) |>
      dplyr::arrange(.data$Date)
  } else {
    tibble::tibble(Date = as.Date(character()))
  }
  
  # ---------- 5) 合并前的日期列检查 ----------
  death_one <- normalize_Date_col(death_one)
  pol_daily <- normalize_Date_col(pol_daily)
  met_daily <- normalize_Date_col(met_daily)
  
  if (nrow(death_one) == 0 && nrow(pol_daily) > 0) {
    message(sprintf("Station '%s' not found in death_pop, deaths set to NA.", station_code))
    death_one <- tibble::tibble(Date = pol_daily$Date, deaths = NA_integer_)
  }
  
  # ---------- 6) 合并所有日级数据 ----------
  full_daily <- list(death_one, pol_daily, met_daily) |>
    purrr::reduce(~ dplyr::full_join(.x, .y, by = "Date")) |>
    dplyr::arrange(.data$Date)
  
  # ---------- 7) 构造时间索引 ----------
  if (nrow(full_daily) == 0) {
    warning(sprintf("No data available for station '%s'.", station_code))
    return(tibble::tibble(time_index = numeric(), Date = as.Date(character())))
  }
  
  last_day <- max(full_daily$Date, na.rm = TRUE)
  full_daily |>
    dplyr::mutate(time_index = as.numeric(difftime(.data$Date, last_day, units = "days"))) |>
    dplyr::select(.data$time_index, .data$Date, .data$deaths, dplyr::everything()) |>
    dplyr::arrange(.data$time_index) |>
    tibble::as_tibble()
}
