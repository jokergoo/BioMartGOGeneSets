
test_that("test BioMartGOGeneSets", {
	expect_error(getBioMartGOGeneSets("foooooooo"), "Wrong dataset")
	expect_error(getBioMartGOGeneSets("hsapiens_gene_ensembl", gene_id_type = "foo"), "does not support gene ID type")
})
