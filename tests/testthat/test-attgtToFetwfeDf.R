# Tests for attgtToFetwfeDf()
# Run with: testthat::test_file("tests/testthat/test-attgtToFetwfeDf.R")

library(testthat)

make_simple_df <- function(){
  grid <- expand.grid(
    county = 1:3,
    year   = 2000:2002
  )
  
  transform(grid,
      first.treat = ifelse(county == 3, 2001, 0),
      lemp        = rnorm(nrow(grid)),
      lpop        = rnorm(nrow(grid))
    )
}

test_that("basic happy path works", {
  df <- make_simple_df()
  tidy <- attgtToFetwfeDf(
    data   = df,
    yname  = "lemp",
    tname  = "year",
    idname = "county",
    gname  = "first.treat",
    covars = c("lpop")
  )
  expect_s3_class(tidy, "data.frame")
  expect_named(tidy, c("time", "unit", "treatment", "y", "lpop"))
  expect_true(all(tidy$time == as.integer(tidy$time)))
  expect_true(all(tidy$treatment %in% 0:1))
  # treated unit 3 should switch on in 2001 and stay on
  unit3 <- subset(tidy, unit == "3")
  expect_equal(unit3$treatment, c(0L,1L,1L))
})

test_that("duplicate (id,time) rows error", {
  df <- make_simple_df()
  dup <- rbind(df, df[1, ])
  expect_error(
    attgtToFetwfeDf(dup, "lemp", "year", "county", "first.treat"),
    "pair must appear at most once"
  )
})

test_that("non‑constant gname within unit triggers error", {
  df <- make_simple_df()
  df$first.treat[df$county == 2 & df$year == 2002] <- 2002
  expect_error(
    attgtToFetwfeDf(df, "lemp", "year", "county", "first.treat"),
    "irreversible treatment"
  )
})

test_that("dropping first‑period treated units", {
  df <- make_simple_df()
  df$first.treat[df$county == 1] <- 2000  # treated from first period
  tidy <- attgtToFetwfeDf(df, "lemp", "year", "county", "first.treat")
  expect_false(any(tidy$unit == "1"))
})

test_that("missing columns are caught", {
  df <- make_simple_df()[, -which(names(make_simple_df())=="lemp")]
  expect_error(
    attgtToFetwfeDf(df, "lemp", "year", "county", "first.treat"),
    "Column(s) not found", fixed=TRUE
  )
})
