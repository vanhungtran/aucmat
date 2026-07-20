test_that("tableroc returns p-value", {
  set.seed(123)
  X <- data.frame(M1 = c(rnorm(95), rep(NA, 5)), M2 = rnorm(100))
  y <- factor(rbinom(100, 1, 0.5), labels = c("neg", "pos"))
  
  # Run the function
  result <- tableroc(X, y, na_impute = "median", rank = TRUE)
  
  # Check if the 'p_value' column exists
  expect_true("p_value" %in% names(result))
  
  # Check if the p-value is a number
  expect_true(is.numeric(result$p_value))
  
  # Check if the p-value is within the valid range [0, 1]
  expect_true(all(result$p_value >= 0 & result$p_value <= 1), na.rm = TRUE)
})
