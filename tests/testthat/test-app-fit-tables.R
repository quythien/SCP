test_that("K2 explorer includes weak-K1 genes with full-universe BH values", {
  app <- read_app_functions()
  full <- data.frame(
    gene = paste0("g", 1:4), symbol = paste0("G", 1:4),
    p_K1 = c(0.0001, 0.1, 0.9, 0.7),
    p_joint = c(0.02, 0.04, 0.00001, 0.5),
    p_2h = c(0.02, 0.05, 0.00001, 0.7),
    A = c(1, 0.3, 0.1, 0.2), phi = c(3, 9, 15, 21),
    A2 = c(0.2, 0.1, 1.4, 0.1), phi2 = c(1, 4, 7, 10),
    sigma = rep(1, 4), mesor = rep(2, 4))
  full$q_joint <- stats::p.adjust(full$p_joint, "BH")
  full$q_K1 <- stats::p.adjust(full$p_K1, "BH")
  capped <- full[full$p_K1 < 0.2, , drop = FALSE]
  capped$pvalue <- capped$p_K1
  pilot <- list(ngenes = 4L, rhythm_fit = capped, explorer_fit = full)
  table <- app$.app_explorer_table(pilot)
  expect_equal(nrow(table), 4L)
  expect_identical(table$ID[1], "g3")
  expect_true(all(table$K2))
  expect_equal(table$q, full$q_joint[match(table$ID, full$gene)])
  expect_equal(table$q[table$ID == "g3"], 0.00004)
})

test_that("legacy K2 subsets do not invent full-universe BH ranks", {
  app <- read_app_functions()
  rf <- data.frame(gene = c("g1", "g2"), pvalue = c(0.001, 0.1),
                    p_joint = c(0.02, 0.04), A = c(1, 0.2), phi = c(3, 9),
                    A2 = c(0.2, 0.3), phi2 = c(1, 4), sigma = c(1, 1))
  table <- app$.app_explorer_table(list(ngenes = 100L, rhythm_fit = rf))
  expect_true(all(is.na(table$q)))
})

test_that("pilot summaries follow the current threshold and top-K", {
  app <- read_app_functions()
  fixture <- pilot_fixture()
  pilot <- fixture$pilot
  pilot$pilot_top_k <- 1L
  pilot$r_median <- 999
  selected <- app$.app_apply_threshold(pilot, "p", 0.01)
  expect_equal(selected$r_median, 2)
  expect_equal(length(selected$amplitude), 1L)
  empty <- app$.app_apply_threshold(pilot, "p", 1e-8)
  expect_equal(empty$prop_rhythmic, 0)
  expect_true(is.na(empty$r_median))
})
