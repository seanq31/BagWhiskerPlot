test_that("bag_whisker function plot the bag-whisker plot correctly", {
  set.seed(1)
  n <- 200
  clean <- mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, 0.6, 0.6, 1), 2))

  contam_frac1 <- 0.10
  k1 <- floor(n * contam_frac1)
  contam1 <- cbind(rnorm(k1, 6, 1), rnorm(k1, -6, 1))
  dat_clean <- clean
  dat_cont1 <- clean
  dat_cont1[1:k1, ] <- contam1

  bag_whisker(dat_clean, type1 = 'FDR', q = 0.1, main = "Clean")
  bag_whisker(dat_cont1, type1 = 'FDR', q = 0.1, main = paste0("Contaminated ", contam_frac1 * 100, "%"))
})