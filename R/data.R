

.onLoad <- function(lib, pkg) {
	ORGANISM_TABLE = readRDS(system.file("extdata", "all_supported_organisms.rds", package = "BioMartGOGeneSets"))
	assign("ORGANISM_TABLE", ORGANISM_TABLE, envir = topenv())
}


# == title (variable:BioMartGOGeneSets)
# Version and source information
#
# == example
# BioMartGOGeneSets
BioMartGOGeneSets = list(
	built_date = "2022-09-20",
	source = "https://www.ensembl.org/info/data/biomart/index.html",
	biomaRt_version = "2.52.0",
	GO.db_version = "3.15.0"
)
class(BioMartGOGeneSets) = "BioMartGOGeneSets_info"

# == title
# Print the BioMartGOGeneSets object
#
# == param
# -x A ``BioMartGOGeneSets_info`` object.
# -... Other arguments
#
# == value
# No value is returned.
#
# == example
# BioMartGOGeneSets
print.BioMartGOGeneSets_info = function(x, ...) {
	cat("BioMart GO gene sets\n")
	cat("  Source:", x$source, "\n")
	cat("  Number of organisms:", nrow(ORGANISM_TABLE), "\n")
	cat("  Marts:", paste(unique(ORGANISM_TABLE$mart), collapse = ", "), "\n")
	cat("  Built date: ", x$built_date, ", with biomaRt (", x$biomaRt_version, "), GO.db (", x$GO.db_version, ")\n", sep = "")
}


# == title
# Get GO gene sets
#
# == param
# -dataset A BioMart dataset. For a proper value, please see `supportedOrganisms`.
# -ontology The value should be "BP", "CC", or "MF".
# -as_table Whether to return the value as a data frame?
# -gene_id_type Since BioMart is from Ensembl database, the default gene ID type is Ensembl gene ID.
#       Depending on different organisms, Entrez ID ("entrez_gene") or gene symbol ("gene_symbol") can also be
#       selected as the gene ID type.
#
# == details
# The gene sets are already compiled and are hosted on https://github.com/jokergoo/BioMartGOGeneSets_data ,
# This function just simply retrieves data from there.
#
# == value
# A list of gene IDs or a data frame.
#
# == example
# lt = getBioMartGOGeneSets("hsapiens_gene_ensembl")
# lt = getBioMartGOGeneSets("hsapiens_gene_ensembl", gene_id_type = "entrez")
# tb = getBioMartGOGeneSets("hsapiens_gene_ensembl", as_table = TRUE)
getBioMartGOGeneSets = function(dataset, ontology = "BP", as_table = FALSE, 
	gene_id_type = "ensembl_gene") {

	dataset = validate_dataset(dataset)

	ontology = tolower(ontology)[1]

	if(!ontology %in% c("bp", "mf", "cc")) {
		stop("`ontology` should be 'BP', 'CC' or 'MF'.")
	}

	url = paste0("https://jokergoo.github.io/BioMartGOGeneSets_data/genesets/", ontology, "_", dataset, "_go_genesets.rds")
	
	gene_id_type = tolower(gene_id_type)

	if(grepl("ensembl", gene_id_type)) {
		url = paste0("https://jokergoo.github.io/BioMartGOGeneSets_data/genesets/", ontology, "_", dataset, "_go_genesets_ensembl_gene_id.rds")
	} else if(grepl("entrez", gene_id_type)) {
		gene = getBioMartGenes(dataset)
		if("entrezgene_id" %in% colnames(mcols(gene))) {

			url = paste0("https://jokergoo.github.io/BioMartGOGeneSets_data/genesets/", ontology, "_", dataset, "_go_genesets_entrezgene_id.rds")

		} else {
			stop("Dataset ", dataset, "does not support gene ID type: ", gene_id_type, ".")
		}
	} else if(grepl("symbol", gene_id_type)) {
		gene = getBioMartGenes(dataset)
		if("external_gene_name" %in% colnames(mcols(gene))) {

			url = paste0("https://jokergoo.github.io/BioMartGOGeneSets_data/genesets/", ontology, "_", dataset, "_go_genesets_external_gene_name.rds")

		} else {
			stop("Dataset ", dataset, "does not support gene ID type: ", gene_id_type, ".")
		}
	} else {
		stop("Dataset ", dataset, "does not support gene ID type: ", gene_id_type, ".")
	}

	res = get_data(url)
	res = lapply(res, as.character)

	if(as_table) {
		res = data.frame(go_geneset = rep(names(res), times = vapply(res, length, 0)), gene = unlist(res))
		colnames(res)[2] = gene_id_type
		rownames(res) = NULL
	}

	return(res)
}

# == title
# Get genes from BioMart
#
# == param
# -dataset A BioMart dataset. For a proper value, please see `supportedOrganisms`.
#
# == value
# A `GenomicRanges::GRanges` object.
#
# == example
# gr = getBioMartGenes("hsapiens_gene_ensembl")
# gr
getBioMartGenes = function(dataset) {

	dataset = validate_dataset(dataset)
	url = paste0("https://jokergoo.github.io/BioMartGOGeneSets_data/genes/granges_", dataset, "_genes.rds")

	get_data(url)
}

validate_dataset = function(dataset) {
	ind = which(ORGANISM_TABLE$dataset %in% dataset)

	if(length(ind)) {
		return(dataset)
	} else {
		ind = which(grepl(dataset, ORGANISM_TABLE$dataset) | grepl(dataset, ORGANISM_TABLE$description, ignore.case = TRUE))
		if(length(ind) == 1) {
			return(ORGANISM_TABLE$dataset[ind])
		} else if(length(ind) > 1) {
			message("Found more than one dataset with query '", dataset, "':")
			for(i in seq_along(ind)) {
				message("  [", i, "] dataset: ", ORGANISM_TABLE$dataset[ind[i]], "; description: ", ORGANISM_TABLE$description[ind[i]])
			}
			
			select = readline(prompt = "Please select which dataset you want to use: ")
			select = as.numeric(select)
			if(is.na(select)) {
				stop("Wrong selection.")
			}

			ind = ind[select]
			if(length(ind) == 1) {
				return(ORGANISM_TABLE$dataset[as.numeric(ind)])
			} else {
				stop("Wrong selection.")
			}
		}
	}

	stop("Wrong dataset.")
}

# == title
# All supported organisms
#
# == param
# -html Whether to open the table in the web browser?
#
# == value
# A data frame of supported organisms.
#
# == example
# if(interactive()) {
#     supportedOrganisms()
# }
supportedOrganisms = function(html = TRUE) {
	if(html) {
		browseURL(system.file("extdata", "global_table.html", package = "BioMartGOGeneSets"))
		invisible(ORGANISM_TABLE)
	} else {
		return(ORGANISM_TABLE)
	}
}


.env = new.env()
get_data = function(url) {
	filename = basename(url)
	if(is.null(.env[[filename]])) {
		tmp = tempfile()
		on.exit(file.remove(tmp))
		download.file(url, destfile = tmp, quiet = TRUE)
		obj = readRDS(tmp)
		.env[[filename]] = obj
	} else {
		obj = .env[[filename]]
	}

	return(obj)
}
