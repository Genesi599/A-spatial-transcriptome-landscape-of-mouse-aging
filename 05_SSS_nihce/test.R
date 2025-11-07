#!/usr/bin/env Rscript
# ===================================================================
# diagnose.R - 诊断 filter 冲突问题
# Usage: Rscript diagnose.R
# ===================================================================

cat("\n🔍 诊断 filter 冲突问题...\n\n")

# 检查包加载顺序
cat("1. 检查已加载的包:\n")
cat("─────────────────────────────────────────\n")
loaded_packages <- search()
cat(paste(loaded_packages, collapse = "\n"))
cat("\n\n")

# 检查 filter 函数来源
cat("2. 检查 filter 函数来源:\n")
cat("─────────────────────────────────────────\n")

if (exists("filter")) {
  cat(sprintf("filter 存在: %s\n", class(filter)))
  
  # 查找所有 filter 函数
  filter_locations <- find("filter")
  cat("filter 可能的来源:\n")
  for (loc in filter_locations) {
    cat(sprintf("  - %s\n", loc))
  }
} else {
  cat("filter 未定义\n")
}

cat("\n")

# 检查 dplyr 是否加载
cat("3. 检查 dplyr:\n")
cat("─────────────────────────────────────────\n")
if ("package:dplyr" %in% search()) {
  cat("✅ dplyr 已加载\n")
  cat(sprintf("   版本: %s\n", packageVersion("dplyr")))
} else {
  cat("❌ dplyr 未加载\n")
}

cat("\n")

# 检查 MASS 是否加载
cat("4. 检查 MASS:\n")
cat("─────────────────────────────────────────\n")
if ("package:MASS" %in% search()) {
  cat("⚠️  MASS 已加载（可能覆盖 dplyr::select）\n")
  cat(sprintf("   版本: %s\n", packageVersion("MASS")))
  cat("\n   建议：使用 MASS::kde2d 而不是 library(MASS)\n")
} else {
  cat("✅ MASS 未加载\n")
}

cat("\n")

# 检查是否有命名冲突
cat("5. 检查函数冲突:\n")
cat("─────────────────────────────────────────\n")

conflicts <- conflicts(detail = TRUE)
if (length(conflicts) > 0) {
  for (func_name in names(conflicts)) {
    cat(sprintf("⚠️  %s 存在冲突:\n", func_name))
    for (pkg in conflicts[[func_name]]) {
      cat(sprintf("   - %s\n", pkg))
    }
  }
} else {
  cat("✅ 未发现函数冲突\n")
}

cat("\n")

# 测试 filter 是否工作
cat("6. 测试 filter 功能:\n")
cat("─────────────────────────────────────────\n")

test_df <- data.frame(x = 1:5, y = letters[1:5])

# 测试1: 直接使用 filter
test1_result <- tryCatch({
  if (requireNamespace("dplyr", quietly = TRUE)) {
    result <- dplyr::filter(test_df, x > 2)
    sprintf("✅ dplyr::filter() 工作正常 (返回 %d 行)", nrow(result))
  } else {
    "❌ dplyr 未安装"
  }
}, error = function(e) {
  sprintf("❌ dplyr::filter() 失败: %s", e$message)
})

cat(paste0("  ", test1_result, "\n"))

cat("\n")

cat("╔═══════════════════════════════════════════════════════════╗\n")
cat("║  诊断完成                                               ║\n")
cat("╚═══════════════════════════════════════════════════════════╝\n\n")

cat("💡 建议:\n")
cat("   1. 先运行: Rscript fix_filter_issues.R\n")
cat("   2. 确保 main.R 中先加载 dplyr，后使用 MASS::kde2d\n")
cat("   3. 如果还有问题，检查 main.R 的包加载顺序\n\n")