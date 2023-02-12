

library(GetoptLong)

organism_meta_file = c(
	"genes_mart" = "~/project/development/BioMartGOGeneSets/inst/extdata/ensembl_vertebrates.csv",
	"protists_mart" = "~/project/development/BioMartGOGeneSets/inst/extdata/ensembl_protists.csv",
	"fungi_mart" = "~/project/development/BioMartGOGeneSets/inst/extdata/ensembl_fungi.csv",
	"metazoa_mart" = "~/project/development/BioMartGOGeneSets/inst/extdata/ensembl_metazoa.csv",
	"plants_mart" = "~/project/development/BioMartGOGeneSets/inst/extdata/ensembl_plants.csv"
)

library(stringdist)

get_datasets = function(ensembl, label = "", overwrite = TRUE) {
	datasets = listDatasets(ensembl)
	# add more columns to `dataset`
	meta = read.csv(organism_meta_file[[label]])
	if(label == "genes_mart") {
		rownames(meta) = meta[, "Ensembl.Assembly"]
		datasets$name = paste0(meta[ datasets$version, "Scientific.name"], " (", meta[ datasets$version, "Common.name"], ")")
		datasets$taxon_id = meta[ datasets$version, "Taxon.ID"]
		datasets$genbank_accession = meta[ datasets$version, "Accession"]
	} else {
		meta[, "Name"] = gsub(" \\(GCA.*\\)$", "", meta[, "Name"])
		meta[, "Name"] = gsub(" - GCA.*$", "", meta[, "Name"])

		ind = sapply(datasets$description, function(x) {
			d = stringdist(paste0(meta[, "Name"], " (", meta[, "Assembly"], ")"), x)
			which.min(d)
		})

		datasets$name = meta[ ind, "Name"]
		datasets$taxon_id = meta[ ind, "Taxon.ID"]
		datasets$genbank_accession = meta[ ind, "Accession"]

	}

	removed_da = NULL
	
	qqcat("####### @{nrow(datasets)} @{label} datasets #######\n")
	for(da in datasets[, 1]) {
		qqcat("@{da}\n")

		if(da %in% c("elucius_gene_ensembl")) {
			removed_da = c(removed_da, da)
			next
		}

		ensembl = useDataset(dataset = da, mart = ensembl)
		all_at = listAttributes(mart = ensembl)
			
		cat("  gene table...\n")
		at = c("chromosome_name", "start_position", "end_position", "strand", "ensembl_gene_id", "gene_biotype", "entrezgene_id", "external_gene_name")
		
		if(!file.exists(qq("original/@{da}_genes.rds")) || overwrite) {
			
			oe = try({gene <- getBM(attributes = intersect(at, all_at[, 1]), mart = ensembl)})
			if(inherits(oe, "try-error")) {
				at = c("chromosome_name", "start_position", "end_position", "strand", "ensembl_gene_id", "gene_biotype")
				gene <- getBM(attributes = at, mart = ensembl)
			}
			saveRDS(gene, file = qq("original/@{da}_genes.rds"), compress = "xz")
		}

		cat("  go gene sets table...\n")
		at = c("ensembl_gene_id", "go_id", "namespace_1003")

		if(!file.exists(qq("original/@{da}_go_genesets_ensembl_gene_id.rds")) || overwrite) {
			
			if(!all(at %in% all_at[, 1])) {
				removed_da = c(removed_da, da)
				qqcat("  remove @{da} because there is no go_id attribute\n")
				next
			}
			qqcat("    ensembl_gene_id...\n")
			oe = try({go = getBM(attributes = at, mart = ensembl)})
			
			if(inherits(oe, "try-error")) {
				go = download_go("ensembl_gene_id", at, ensembl)
			}
			saveRDS(go, file = qq("original/@{da}_go_genesets_ensembl_gene_id.rds"), compress = "xz")
			
		}

		try({
			for(gene_id_type in c("entrezgene_id", "external_gene_name")) {
				at = c(gene_id_type, "go_id", "namespace_1003")

				if(!all(at %in% all_at[, 1])) {
					next
				}

				if(!file.exists(qq("original/@{da}_go_genesets_@{gene_id_type}.rds")) || overwrite) {
					qqcat("    @{gene_id_type}...\n")
					oe = try({ go = getBM(attributes = at, mart = ensembl) })
					if(inherits(oe, "try-error")) {
						go = download_go(gene_id_type, at, ensembl)
					}
					saveRDS(go, file = qq("original/@{da}_go_genesets_@{gene_id_type}.rds"), compress = "xz")
				}
			}
		})
	}

	datasets = datasets[!datasets[, 1] %in% removed_da, ]

	write.csv(datasets, file = qq("datasets_@{label}.csv"), row.names = FALSE)

}

download_go = function(gene_id_type, at, ensembl) {
	gene_id = unique(getBM(attributes = gene_id_type, mart = ensembl)[, 1])
	gene_id = gene_id[gene_id != ""]
	go = NULL
	for(i in seq_len(floor(length(gene_id)/1000))) {
		ind = 1:1000 +1000*(i-1)
		ind = ind[ind <= length(gene_id)]
		qqcat("    genes: @{ind[1]} ~ @{ind[length(ind)]}\n")
		gi = gene_id[ind]
		go = rbind(go, getBM(attributes = at, mart = ensembl, filter = gene_id_type, value = gi))
	}
	go
}


setwd("~/project/development/BioMartGOGeneSets_data")

###################################################################
#### download data from BioMart in the original format   ##########

library(biomaRt)
ensembl = useEnsembl(biomart = "genes", mirror = "www")
# attributes = listAttributes(ensembl)
get_datasets(ensembl, "genes_mart", overwrite = TRUE)

listEnsemblGenomes()

for(mart in c("protists_mart", "fungi_mart", "metazoa_mart", "plants_mart")) {
	ensembl = useEnsemblGenomes(biomart = mart)
	get_datasets(ensembl, mart, overwrite = TRUE)
}

####################################
##### convert the data format ######

#### genes saved as GRanges objects:
### Ensembl ID is taken as the primary ID type
library(GenomicRanges)
all_files = list.files(path = "original", pattern = "genes.rds$", full.names = TRUE)
for(f in all_files) {
	cat(f, "\n")
	df = readRDS(f)
	df = df[df$start_position <= df$end_position, , drop = FALSE]

	chr = as.character(as.vector(tapply(df$chromosome_name, df$ensembl_gene_id, function(x) x[1])))
	start = as.numeric(as.vector(tapply(df$start_position, df$ensembl_gene_id, function(x) x[1])))
	end = as.numeric(as.vector(tapply(df$end_position, df$ensembl_gene_id, function(x) x[1])))
	strand = as.vector(tapply(df$strand, df$ensembl_gene_id, function(x) x[1]))
	ensembl_gene_id = as.vector(tapply(df$ensembl_gene_id, df$ensembl_gene_id, function(x) x[1]))

	gr = GRanges(seqnames = chr, ranges = IRanges(start, end), strand = ifelse(strand > 0, "+", "-"),
		ensembl_gene_id = ensembl_gene_id)
	names(gr) = ensembl_gene_id

	for(nm in intersect(colnames(df), c("gene_biotype", "entrezgene_id", "external_gene_name"))) {
		df[[nm]] = as.character(df[[nm]])
		l = df[[nm]] == ""; l[is.na(l)] = FALSE
		df[l, nm] = NA
		v = tapply(df[[nm]], df$ensembl_gene_id, unique, simplify = FALSE)

		if(all(sapply(v, length) == 1)) {
			mcols(gr)[[nm]] = unlist(v)
		} else {
			mcols(gr)[[nm]] = CharacterList(v)
		}
	}
	saveRDS(gr, qq("processed/genes/granges_@{basename(f)}"), compress = "xz")
}


#### for each GO term, merge all its offspring terms
library(GO.db)

all_terms = data.frame(go_id = GOID(GOTERM), namespace = Ontology(GOTERM))
bp_terms = all_terms$go_id[all_terms$namespace == "BP"]
cc_terms = all_terms$go_id[all_terms$namespace == "CC"]
mf_terms = all_terms$go_id[all_terms$namespace == "MF"]

GOBPOFFSPRING = as.list(GOBPOFFSPRING)
GOCCOFFSPRING = as.list(GOCCOFFSPRING)
GOMFOFFSPRING = as.list(GOMFOFFSPRING)

all_files = list.files(path = "original", pattern = "go_genesets", full.names = TRUE)
# for(f in all_files) {
parallel::mclapply(all_files, function(f) {
	cat(f, "\n")
	df = readRDS(f)
	df = df[df$go_id != "" & df[, 1] != "", , drop = FALSE]
	df[, 1] = as.character(df[, 1])

	df1 = df[df$namespace_1003 == "biological_process", , drop = FALSE]
	gs = split(df1[, 1], df1$go_id)

	gs2 = lapply(bp_terms, function(nm) {
		go_id = c(nm, GOBPOFFSPRING[[nm]])
		unique(unlist(gs[go_id]))
	})
	names(gs2) = bp_terms
	gs2 = gs2[sapply(gs2, length) > 0]
	saveRDS(gs2, qq("processed/genesets/bp_@{basename(f)}"), compress = "xz")

	df1 = df[df$namespace_1003 == "cellular_component", , drop = FALSE]
	gs = split(df1[, 1], df1$go_id)

	gs2 = lapply(cc_terms, function(nm) {
		go_id = c(nm, GOCCOFFSPRING[[nm]])
		unique(unlist(gs[go_id]))
	})
	names(gs2) = cc_terms
	gs2 = gs2[sapply(gs2, length) > 0]
	saveRDS(gs2, qq("processed/genesets/cc_@{basename(f)}"), compress = "xz")

	df1 = df[df$namespace_1003 == "molecular_function", , drop = FALSE]
	gs = split(df1[, 1], df1$go_id)

	gs2 = lapply(mf_terms, function(nm) {
		go_id = c(nm, GOMFOFFSPRING[[nm]])
		unique(unlist(gs[go_id]))
	})
	names(gs2) = mf_terms
	gs2 = gs2[sapply(gs2, length) > 0]
	saveRDS(gs2, qq("processed/genesets/mf_@{basename(f)}"), compress = "xz")
}, mc.cores = 4)


####### save the global table as a rds file #####

tbl = NULL
for(file in list.files(pattern = "datasets.*.csv")) {
	tb0 = read.csv(file)
	tb0$mart = gsub("^datasets_(.*)\\.csv$", "\\1", file)
	tb0$n_gene = 0
	tb0$n_bp_genesets = 0
	tb0$n_cc_genesets = 0
	tb0$n_mf_genesets = 0

	for(i in seq_len(nrow(tb0))) {
		obj = readRDS(qq("processed/genes/granges_@{tb0$dataset[i]}_genes.rds"))
		tb0$n_gene[i] = length(obj)

		obj = readRDS(qq("processed/genesets/bp_@{tb0$dataset[i]}_go_genesets_ensembl_gene_id.rds"))
		tb0$n_bp_genesets[i] = length(obj)

		obj = readRDS(qq("processed/genesets/cc_@{tb0$dataset[i]}_go_genesets_ensembl_gene_id.rds"))
		tb0$n_cc_genesets[i] = length(obj)

		obj = readRDS(qq("processed/genesets/mf_@{tb0$dataset[i]}_go_genesets_ensembl_gene_id.rds"))
		tb0$n_mf_genesets[i] = length(obj)
	}

	tbl[[file]] = tb0
}

tb = do.call(rbind, tbl)
rownames(tb) = tb[, 1]

saveRDS(tb, file = "all_supported_organisms.rds", compress = "xz")

d = NULL
for(i in seq_len(nrow(tb2))) {
	d[i] = stringdist(gsub("( genes )?\\(.*\\)", "", tb2$description[i]), gsub("\\(.*\\)", "", tb2$name[i]))
}

# manually check
tb["aastaci_eg_gene", "name"] = gsub("( genes )?\\(.*\\)", "", tb["aastaci_eg_gene", "description"])
tb["aastaci_eg_gene", "taxon_id"] = "112090"
tb["aastaci_eg_gene", "genbank_accession"] = "GCA_000520075.1"

tb["ainvadans_eg_gene", "name"] = gsub("( genes )?\\(.*\\)", "", tb["ainvadans_eg_gene", "description"])
tb["ainvadans_eg_gene", "taxon_id"] = "157072"
tb["ainvadans_eg_gene", "genbank_accession"] = "GCA_000520115.1"

tb["gultimum_eg_gene", "name"] = gsub("( genes )?\\(.*\\)", "", tb["gultimum_eg_gene", "description"])
tb["gultimum_eg_gene", "taxon_id"] = "431595"
tb["gultimum_eg_gene", "genbank_accession"] = "GCA_000143045.1"


# tb["pgraminisug99_eg_gene", "genbank_accession"] = "GCA_000149925.1"
# tb["choffmanni_gene_ensembl", "genbank_accession"] = "GCA_000164785.2"
# tb["csavignyi_gene_ensembl", "genbank_accession"] = "GCA_000149265.1"
# tb["eeuropaeus_gene_ensembl", "genbank_accession"] = "GCA_000296755.1"
# tb["etelfairi_gene_ensembl", "genbank_accession"] = 
# tb["gaculeatus_gene_ensembl", "genbank_accession"] = 
# tb["oprinceps_gene_ensembl", "genbank_accession"] = 
# tb["pcapensis_gene_ensembl", "genbank_accession"] = 
# tb["pmarinus_gene_ensembl", "genbank_accession"] = 
# tb["pvampyrus_gene_ensembl", "genbank_accession"] = 
# tb["saraneus_gene_ensembl", "genbank_accession"] = 
# tb["tbelangeri_gene_ensembl", "genbank_accession"] = 
# tb["tnigroviridis_gene_ensembl", "genbank_accession"] = 
# tb["ttruncatus_gene_ensembl", "genbank_accession"] = "GCA_001922835.1"
# tb["vpacos_gene_ensembl", "genbank_accession"] = 
# tb["csonorensis_eg_gene", "genbank_accession"] = 
# tb["lsalmonis_eg_gene", "genbank_accession"] = 
# tb["alaibachii_eg_gene", "genbank_accession"] = 

## add ncbi genome link
library(rvest)
tb$ncbi_genome_link = NA
for(i in seq_len(nrow(tb))) {
	cat(i, "/", nrow(tb), "...\n")
	base_link = "https://ftp.ncbi.nih.gov/genomes/all/GCA"
	code1 = substr(tb[i, "genbank_accession"], 5, 5+2)
	code2 = substr(tb[i, "genbank_accession"], 8, 8+2)
	code3 = substr(tb[i, "genbank_accession"], 11, 11+2)
	base_link = qq("@{base_link}/@{code1}/@{code2}/@{code3}/")

	oe = try(html <- read_html(base_link))
	if(!inherits(oe, "try-error")) {
		assembly = html %>% html_elements("a") %>% html_text()
		assembly = assembly[-c(1, length(assembly))]
		assembly = gsub("/$", "", assembly)
		assembly = assembly[ grepl(tb[i, "genbank_accession"], assembly) ]
		if(length(assembly)) {
			link = qq("@{base_link}/@{assembly}/@{assembly}_assembly_report.txt")
			tb$ncbi_genome_link[i] = link
		}
	}
}

saveRDS(tb, file = "all_supported_organisms.rds", compress = "xz")

