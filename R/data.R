

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
	built_date = "2023-02-11",
	source = "https://www.ensembl.org/info/data/biomart/index.html",
	biomaRt_version = "2.54.0",
	GO.db_version = "3.16.0"
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
	cat("BioMart Gene Ontology gene sets\n")
	cat("  Source:", x$source, "\n")
	cat("  Number of organisms:", nrow(ORGANISM_TABLE), "\n")
	cat("  Marts:", paste(unique(ORGANISM_TABLE$mart), collapse = ", "), "\n")
	cat("  Built date: ", x$built_date, ", with biomaRt (", x$biomaRt_version, "), GO.db (", x$GO.db_version, ")\n", sep = "")
}


# == title
# Get GO gene sets
#
# == param
# -dataset A BioMart dataset or a taxon ID. For a proper value, please see `supportedOrganisms`.
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
getBioMartGOGeneSets = function(dataset, ontology = "BP", 
	as_table = FALSE, gene_id_type = "ensembl_gene") {

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
# -dataset A BioMart dataset or a taxon ID. For a proper value, please see `supportedOrganisms`.
# -add_chr_prefix Whether to add "chr" prefix to chromosome names? If it is ture, it uses ``GenomeInfoDb::seqlevelsStyle(gr) = "UCSC"`` to add the prefix.
#
# == details
# Note ``add_chr_prefix`` is just a helper argument. You can basically do the same as:
#
#     gr = getBioMartGenes("hsapiens_gene_ensembl")
#     seqlevelsStyle(gr) = "UCSC"
#
#
# == value
# A `GenomicRanges::GRanges` object.
#
# == example
# gr = getBioMartGenes("hsapiens_gene_ensembl")
# gr
# gr = getBioMartGenes("hsapiens_gene_ensembl", add_chr_prefix = TRUE)
# gr
getBioMartGenes = function(dataset, add_chr_prefix = FALSE) {

	dataset = validate_dataset(dataset)
	url = paste0("https://jokergoo.github.io/BioMartGOGeneSets_data/genes/granges_", dataset, "_genes.rds")

	gr = get_data(url)

	if(add_chr_prefix) {
		GenomeInfoDb::seqlevelsStyle(gr) = "UCSC"
	}

	gr
}

# == title
# Change sequence names
#
# == param
# -gr The input regions
# -dataset A BioMart dataset or a taxon ID. For a proper value, please see `supportedOrganisms`.
# -seqname_style_from Value should be in ``c("Sequence-Name", "GenBank-Accn", "RefSeq-Accn")``. If you
#        are not sure which seqname style is in ``gr``, use `getBioMartGenomeInfo` to obtain list of examples.
# -seqname_style_to Value should be in ``c("Sequence-Name", "GenBank-Accn", "RefSeq-Accn")``.
# -reformat_from A self-defined function to reformat the seqnames. The internal seqname style can be obtained via ``getBioMartGenomeInfo(dataset)``. This 
#       function converts the internal "from" seqnames to fit the user's input regions.
# -reformat_to A self-defined function to reformat the seqnames.
#
# == details
# Please the conversion is not one to one. For those sequences which cannot be corrected mapped to other styles,
# they are just removed.
#
# == value
# A `GenomicRanges::GRanges` object.
#
# == example
# \dontrun{
# gr = getBioMartGenes("giant panda")
# changeSeqnameStyle(gr, "giant panda", "Sequence-Name", "GenBank-Accn")
# }
changeSeqnameStyle = function(gr, dataset, seqname_style_from, seqname_style_to, 
	reformat_from = NULL, reformat_to = NULL) {
	
	if(!seqname_style_from %in% c("Sequence-Name", "GenBank-Accn", "RefSeq-Accn")) {
		stop("`seqname_style_from` should take value in 'Sequence-Name', 'GenBank-Accn' and 'RefSeq-Accn'.")
	}
	if(!seqname_style_to %in% c("Sequence-Name", "GenBank-Accn", "RefSeq-Accn")) {
		stop("`seqname_style_to` should take value in 'Sequence-Name', 'GenBank-Accn' and 'RefSeq-Accn'.")
	}

	dataset = validate_dataset(dataset)
	ind = which(ORGANISM_TABLE$dataset == dataset)
	if(is.na(ORGANISM_TABLE[ind, "ncbi_genome_link"])) {
		stop(paste0("Cannot find source to convert seqnames for dataset '", dataset, "'."))
	}
	genome_tb = get_data(ORGANISM_TABLE[ind, "ncbi_genome_link"], read.table, sep = "\t")
	colnames(genome_tb) = c("Sequence-Name", "Sequence-Role", "Assigned-Molecule", "Assigned-Molecule-Location/Type", "GenBank-Accn", "Relationship", "RefSeq-Accn", "Assembly-Unit", "Sequence-Length", "UCSC-style-name")

	if(seqname_style_from != "Sequence-Name") {
		if(is.null(reformat_from)) {
			reformat_from = function(x) gsub("\\.\\d+$", "", x)
		}
	}
	if(!is.null(reformat_from)) {
		genome_tb[, seqname_style_from] = reformat_from(genome_tb[, seqname_style_from])
	}

	map = structure(genome_tb[, seqname_style_to], 
		    names = genome_tb[, seqname_style_from])


	if(!is.null(reformat_to)) {
		map = structure(unname(reformat_to(map)), names = names(map))
	}

	if(length(table(table(map))) > 1) {
		stop(paste0("Mapping from '", seqname_style_from, "' to '", seqname_style_to, "' is not one to one."))
	}

	chr = as.character(as.vector(seqnames(gr)))
	if(seqname_style_from == "Sequence-Name") {
		chr2 = map[chr]
	} else {
		chr2 = map[gsub("\\.\\d+$", "", chr)]
	}
	l = !is.na(chr2)
	gr2 = GRanges(seqnames = chr2[l], ranges = ranges(gr)[l], strand = strand(gr)[l])
	mcols(gr2) = mcols(gr)[l, ]
	mcols(gr2)$.original_seqname = chr[l]
	gr2
}

# == title
# Get genome information
#
# == param
# -dataset A BioMart dataset or a taxon ID. For a proper value, please see `supportedOrganisms`.
#
# == value
# A list.
#
# == example
# getBioMartGenomeInfo(9606)
getBioMartGenomeInfo = function(dataset) {
	dataset = validate_dataset(dataset)
	if(grepl("^\\d+$", dataset)) {
		ind = which(ORGANISM_TABLE$taxon_id == dataset)
		if(length(ind)) {
			dataset = ORGANISM_TABLE$dataset[ind]
		}
	}
	ind = which(ORGANISM_TABLE$dataset == dataset)
	if(length(ind)) {
		lt = as.list(ORGANISM_TABLE[ind, ])
		lt = lt[c("dataset", "version", "name", "taxon_id", "genbank_accession", "mart")]

		if(!is.na(ORGANISM_TABLE[ind, "ncbi_genome_link"])) {
			genome_tb = get_data(ORGANISM_TABLE[ind, "ncbi_genome_link"], read.table, sep = "\t")
			colnames(genome_tb) = c("Sequence-Name", "Sequence-Role", "Assigned-Molecule", "Assigned-Molecule-Location/Type", "GenBank-Accn", "Relationship", "RefSeq-Accn", "Assembly-Unit", "Sequence-Length", "UCSC-style-name")

			genome_tb = genome_tb[seq_len(min(5, nrow(genome_tb))), ]
			chr_style = as.list(genome_tb[, c("Sequence-Name", "GenBank-Accn", "RefSeq-Accn")])
			lt$seqname_style = chr_style
		}

		return(lt)
	} else {
		stop(paste0("Cannot find dataset '", dataset, "'."))
	}
}

validate_dataset = function(dataset) {
	if(grepl("^\\d+$", dataset)) {
		ind = which(ORGANISM_TABLE$taxon_id == dataset)
		if(length(ind)) {
			dataset = ORGANISM_TABLE$dataset[ind]
		}
	}
	ind = which(ORGANISM_TABLE$dataset %in% dataset)

	if(length(ind)) {
		return(dataset)
	} else {
		ind = which(grepl(dataset, ORGANISM_TABLE$dataset) | grepl(dataset, ORGANISM_TABLE$name, ignore.case = TRUE))
		if(length(ind) == 1) {
			return(ORGANISM_TABLE$dataset[ind])
		} else if(length(ind) > 1) {
			message("Found more than one dataset with query '", dataset, "':")
			for(i in seq_along(ind)) {
				message("  [", i, "] dataset: ", ORGANISM_TABLE$dataset[ind[i]], "; name: ", ORGANISM_TABLE$name[ind[i]], "; taxon_id: ", ORGANISM_TABLE$taxon_id[ind[i]])
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
get_data = function(url, fun = readRDS, ...) {
	filename = basename(url)
	if(is.null(.env[[filename]])) {
		tmp = tempfile()
		on.exit(file.remove(tmp))
		download.file(url, destfile = tmp, quiet = TRUE)
		obj = fun(tmp, ...)
		.env[[filename]] = obj
	} else {
		obj = .env[[filename]]
	}

	return(obj)
}
