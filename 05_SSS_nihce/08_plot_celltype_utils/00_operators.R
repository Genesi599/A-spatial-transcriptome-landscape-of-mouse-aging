# ===================================================================
# 00_operators.R
# 基础操作符定义
# Author: Assistant
# Date: 2025-11-06
# ===================================================================

#' %||% 操作符（空值默认值）
#'
#' @param a 主值
#' @param b 默认值（当a为NULL时使用）
#' @return 返回a（如果非NULL）或b
#'
#' @examples
#' NULL %||% "default"  # 返回 "default"
#' "value" %||% "default"  # 返回 "value"
#'
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}