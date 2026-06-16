test_that("togo_matrix_to_list inverts a membership matrix", {
  m <- matrix(c(1, 0, 1, 0, 1, 1), nrow = 3,
              dimnames = list(c("g1", "g2", "g3"), c("setA", "setB")))
  out <- togo_matrix_to_list(m)
  expect_equal(out$setA, c("g1", "g3"))
  expect_equal(out$setB, c("g2", "g3"))
})

test_that("togo_clean_pathway_names strips prefixes and tidies", {
  expect_equal(togo_clean_pathway_names("REACTOME_CELL_CYCLE"), "Cell Cycle")
  expect_match(togo_clean_pathway_names("HALLMARK_DNA_REPAIR"), "DNA")
})
