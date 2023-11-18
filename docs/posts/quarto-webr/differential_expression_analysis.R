library(edgeR)
library(here)
library(limma)
library(readr)
library(readxl)

# retrieve and parse raw counts
url <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE158152&",
              "format=file&file=GSE158152%5Fdst150%5Fprocessed%2Exlsx")
temp_file <- tempfile(fileext = ".xlsx")
download.file(url, destfile = temp_file)
raw_counts <- read_excel(temp_file, sheet = "raw_counts")

# retrieve mapping between geo samples & the author's sample identifiers
url <- paste0("https://tomsing1.github.io/blog/posts/",
              "nextflow-core-quantseq-3-xia/GEO_sample_ids.csv")
geo_ids <- read.csv(url)

# retrieve sample annotations
url <- paste0("https://tomsing1.github.io/blog/posts/",
              "nextflow-core-quantseq-3-xia/sample_metadata.csv")
sample_anno <- read.csv(url, row.names = "Experiment")
colnames(sample_anno)<- tolower(colnames(sample_anno))
colnames(sample_anno) <- sub(".", "_", colnames(sample_anno), 
                             fixed = TRUE) 
sample_anno <- sample_anno[, c("sample_name", "animal_id", "genotype", "sex",
                               "batch")]
sample_anno$genotype <- factor(sample_anno$genotype, 
                               levels = c("WT", "Het", "Hom"))
sample_anno$sample_title <- geo_ids[
  match(sample_anno$sample_name, geo_ids$sample_name), "sample_id"]

# create DGEList object
count_data <- as.matrix(raw_counts[, grep("DRN-", colnames(raw_counts))])
row.names(count_data) <- raw_counts$feature_id
colnames(count_data) <- row.names(sample_anno)[
  match(colnames(count_data), sample_anno$sample_title)
]

gene_data <- data.frame(
  gene_id = raw_counts$feature_id,
  gene_name = raw_counts$symbol,
  row.names = raw_counts$feature_id
)

col_data <- data.frame(
  sample_anno[colnames(count_data),
              c("sample_title", "animal_id", "sex", "genotype", "batch")],
  workflow = "geo"
)

dge <- DGEList(
  counts = as.matrix(count_data), 
  samples = col_data[colnames(count_data), ], 
  genes = gene_data[row.names(count_data), ]
)

dge <- normLibSizes(dge, method = "TMM")

# linear modeling
dge$samples$group <- with(dge$samples, paste(genotype, batch, sep = "_"))
design <- model.matrix(~ 0 + group + sex, data = dge$samples)
colnames(design) <- sub("^group", "", colnames(design))
keep <- filterByExpr(dge, design = design, min.count = 25)
fit <- voomLmFit(
  dge[keep, row.names(design)], 
  design = design,
  block = dge$samples$animal_id, 
  sample.weights = TRUE, 
  plot = FALSE
)
contrasts <- makeContrasts(
  day1 = "Hom_Day1-WT_Day1",
  day2 = "Hom_Day2-WT_Day2",
  day3 = "Hom_Day3-WT_Day3",
  levels = design
)
fit <- contrasts.fit(fit, contrasts = contrasts)
fit2 <- eBayes(fit, robust=TRUE)

for (contr in colnames(fit2)) {
  tt <- topTable(fit2, coef = contr, number = Inf)
  tt$logFC <- signif(tt$logFC, 2)
  tt$adj.P.Val <- signif(tt$adj.P.Val, 3)
  write_csv(tt[, c("gene_name", "logFC", "adj.P.Val")], 
            here("posts", "quarto-webr", sprintf("%s_results.csv.gz", contr)))
}

  