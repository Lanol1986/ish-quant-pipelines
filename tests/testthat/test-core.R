test_that("continuous estimated indices cannot be used as counts", {
  expect_error(ish_assert_integer_endpoint(c(0,1.2,2),"estimated"),"not integer-valued")
  expect_silent(ish_assert_integer_endpoint(c(0,1,2),"count"))
})

test_that("missing-value policy is explicit", {
  expect_error(ish_apply_missing_policy(c("1",NA),"error","x"),"missing")
  z <- ish_apply_missing_policy(c("1",NA),"zero_if_export_blank","x")
  expect_equal(z$values,c(1,0))
})

test_that("panel rows do not change biological n", {
  d <- tibble::tibble(animal_id=rep(c("2.1","2.2"),each=2),
                      panel=rep(c("Panel 1","Panel 2"),2))
  expect_equal(nrow(d),4)
  expect_equal(dplyr::n_distinct(d$animal_id),2)
})
