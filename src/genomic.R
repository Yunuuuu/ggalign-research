library(ggalign)
library(magrittr)
library(forcats)
library(legendry)
if (!dir.exists("figures/genomic")) dir.create("figures/genomic")

# p3 --------------------------------------------------------
# p4 ------------------------------------------------------
# See src/genomic-primary-vs-metastatic.R

# p5 --------------------------------------------------------
clidata <- readRDS("~/Data/TCGA_data/results/TCGA-BLCA/blca_clidata.rds")
rna_data <- readRDS("~/Data/TCGA_data/results/TCGA-BLCA/blca_rna_expr.rds")
methy <- readRDS(
    "~/Data/TCGA_data/results/TCGA-BLCA/blca_methy450_beta_data.rds"
)

anno450 <- data.table::fread("rawdata/HM450.hg38.manifest.tsv.gz")
unique(vapply(strsplit(names(rna_data), "_"), `[[`, character(1), 1))
counts <- dplyr::select(rna_data, gene_name, starts_with("unstranded_"))
counts <- dplyr::rename_with(counts, ~ sub("unstranded_", "", .x))
tpm <- dplyr::select(rna_data, gene_name, starts_with("tpm_unstranded_"))
tpm <- dplyr::rename_with(tpm, ~ sub("tpm_unstranded_", "", .x))
counts[1:10, 1:10]
tpm[1:10, 1:10]

rna_samples <- yjtools::tcga_remove_duplicated_samples(colnames(counts)[-1])
methy_samples <- yjtools::tcga_remove_duplicated_samples(colnames(methy))
samples <- intersect(
    rna_samples$sample_barcode,
    methy_samples$sample_barcode
)
clidata <- clidata %>%
    dplyr::filter(sample_barcode %in% .env$samples) %>%
    dplyr::select(!sample_vial_barcode) %>%
    dplyr::distinct()
anyDuplicated(clidata$sample_barcode)
clidata <- clidata[match(samples, clidata$sample_barcode), ]
rna_samples <- dplyr::filter(rna_samples, sample_barcode %in% .env$samples)
methy_samples <- dplyr::filter(methy_samples, sample_barcode %in% .env$samples)
class(methy)
methy <- methy[, methy_samples$barcode]
dim(methy)
counts <- dplyr::select(counts, gene_name, all_of(rna_samples$barcode))
tpm <- dplyr::select(tpm, gene_name, all_of(rna_samples$barcode))
index <- tpm %>%
    dplyr::mutate(
        id = seq_len(dplyr::n()),
        means = rowSums(dplyr::across(!c(gene_name, id)))
    ) %>%
    dplyr::slice_max(means, by = gene_name, with_ties = FALSE) %>%
    dplyr::pull(id)
counts <- dplyr::slice(counts, .env$index) %>%
    tibble::column_to_rownames("gene_name") %>%
    as.matrix()
tpm <- dplyr::slice(tpm, .env$index) %>%
    tibble::column_to_rownames("gene_name") %>%
    as.matrix()
all(clidata$sample_barcode == rna_samples$sample_barcode)
all(clidata$sample_barcode == methy_samples$sample_barcode)
clidata$sample_type <- factor(
    clidata$sample_type,
    c("Solid Tissue Normal", "Primary Tumor")
)

## differential expression analysis
rna_diff <- DESeq2::DESeqDataSetFromMatrix(
    counts,
    colData = data.frame(clidata),
    design = ~sample_type
)
rna_diff <- DESeq2::DESeq(rna_diff)
DESeq2::resultsNames(rna_diff)
rna_diff <- DESeq2::results(rna_diff, alpha = 0.01, tidy = TRUE)

## differential methylated CpGs
dim(methy)
# remove probes with NA value
methy <- methy[rowSums(is.na(methy)) == 0L, ]
# remove masked probes and only keep CpG probes
methy <- methy[rownames(methy) %in% anno450$probeID[
    anno450$probeType == "cg" & !anno450$MASK_general
], ]
dim(methy)

# fit the linear model
methy_diff <- limma::eBayes(
    limma::lmFit(
        # Converting beta values to m_values
        log2(methy / (1 - methy)),
        model.matrix(~sample_type, data = clidata)
    ),
    proportion = 0.01
)
methy_diff <- limma::topTable(methy_diff, coef = 2L, number = Inf)
delta_beta <- rowMeans(methy[, clidata$sample_type == "Primary Tumor"]) -
    rowMeans(methy[, clidata$sample_type == "Solid Tissue Normal"])

rna_diff_genes <- rna_diff %>%
    dplyr::filter(padj < 0.01 & abs(log2FoldChange) > 1L) %>%
    dplyr::select(row, padj_rna = padj, logFC = log2FoldChange)

methy_diff_cg <- methy_diff %>%
    tibble::rownames_to_column("probeid") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
        gene_name = anno450$gene_HGNC[match(probeid, anno450$probeID)],
        delta_beta = delta_beta[match(probeid, names(delta_beta))],
        .after = probeid
    ) %>%
    dplyr::filter(adj.P.Val < 0.01, abs(delta_beta) > 0.2) %>%
    dplyr::filter(!is.na(gene_name)) %>%
    dplyr::select(probeid, gene_name, padj_methy = adj.P.Val, delta_beta) %>%
    tidyr::separate_longer_delim(gene_name, delim = ";")

target <- rna_diff_genes %>%
    dplyr::group_by(sign(logFC)) %>%
    dplyr::slice_max(abs(logFC), n = 50) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(methy_diff_cg, by = c(row = "gene_name")) %>%
    dplyr::filter(sign(logFC) != sign(delta_beta))

p5 <- stack_crossh(log2(tpm[unique(target$row), ]), sizes = c(0.2, 1, 1)) +
    ggheatmap() +
    theme_no_axes("x") +
    no_expansion("x") +
    theme(plot.margin = margin()) +
    anno_top() +
    align_group(clidata$sample_type) +
    stack_active() +
    cross_link(
        link_line(
            !!!purrr::imap(split(target$probeid, target$row), function(r, l) {
                l ~ r
            })
        ),
        methy[target$probeid, ],
        size = 0.2
    ) +
    no_expansion("x") +
    theme(plot.margin = margin()) +
    ggheatmap() +
    scale_fill_viridis_c(name = "Beta") +
    scale_y_continuous(position = "right") +
    theme_no_axes("x") +
    theme(plot.margin = margin()) +
    no_expansion("x") +
    anno_top() -
    scheme_align(free_spaces = "l") +
    align_group(clidata$sample_type)

ggsave("figures/genomic/p5.pdf", plot = p5, width = 10, height = 7)

