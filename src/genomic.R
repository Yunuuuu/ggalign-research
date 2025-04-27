library(ggalign)
library(magrittr)
library(forcats)
library(legendry)
if (!dir.exists("figures/genomic")) dir.create("figures/genomic")

# p1 -----------------------------------------------------------
mut <- readRDS("~/Data/TCGA_data/results/TCGA-PRAD/prad_snv_maf.rds")
maf <- maftools::read.maf(mut)
gene_data <- fortify_matrix(maf,
    collapse_vars = FALSE, n_top = 20L,
    missing_genes = "remove"
)

pathway_data <- fortify_matrix(tune(maf))
highlight_tp53 <- stack_crossv(
    rlang::set_names(
        colnames(gene_data)[!is.na(gene_data["TP53", , drop = TRUE])]
    )
) +
    # add plot for TP53 data only
    ggalign(mapping = aes(.column_names, .row_names)) +
    scheme_data(function(d) {
        out <- fortify_data_frame(pathway_data) |>
            dplyr::select(!.column_index) |>
            dplyr::filter(.data$.column_names %in% .env$d$value)
        n <- length(unique(d$value))
        prop <- dplyr::summarise(
            out,
            n = length(unique(.column_names[!is.na(value)])),
            .by = .row_names
        )
        prop$prop <- prop$n / n
        prop <- dplyr::filter(
            prop,
            # prop > 1/2,
            prop > 0.1,
            !.row_names %in% c("Other", "Other signaling")
        ) %>%
            dplyr::arrange(.data$n)
        out <- dplyr::filter(out, .data$.row_names %in% prop$.row_names)
        dplyr::mutate(
            out,
            .row_names = factor(.row_names, rev(prop$.row_names)),
            .column_names = fct_reorder(
                .column_names, value,
                function(x) sum(!is.na(x)),
                .na_rm = FALSE, .desc = TRUE
            )
        )
    }) +
    geom_tile(fill = "grey", height = 0.95) +
    geom_tile(
        fill = "white",
        width = 0.9, height = 0.9,
        linejoin = "round", lineend = "round"
    ) +
    geom_tile(
        fill = "#616060fc",
        width = 0.9, height = 0.8,
        linejoin = "round", lineend = "round",
        data = function(d) dplyr::filter(d, !is.na(value))
    ) +
    geom_rect(
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "#616060fc",
        linejoin = "round", lineend = "round",
        data = function(d) {
            out <- dplyr::filter(d, !is.na(value))
            out$x <- as.integer(factor(out$.column_names))
            out$y <- as.integer(out$.row_names)
            out <- dplyr::mutate(
                dplyr::arrange(out, x, y),
                group = cumsum(diff(c(0, y)) != 1L),
                .by = x
            )
            dplyr::summarise(
                dplyr::filter(out, dplyr::n() > 1L, .by = c(x, group)),
                xmin = unique(x) - 0.5, xmax = unique(x) + 0.5,
                ymin = min(y) - 0.5, ymax = max(y) + 0.5,
                .by = c(x, group)
            )
        },
        inherit.aes = FALSE
    ) +
    geom_tile(aes(fill = .row_names, color = after_scale(fill)),
        width = 0.9, height = 0.7, linewidth = 0.1,
        linejoin = "round", lineend = "round",
        function(d) dplyr::filter(d, !is.na(value))
    ) +
    ggsci::scale_fill_d3(palette = "category20", guide = "none") +
    # scale_fill_brewer(palette = "Paired", guide = "none") +
    theme_no_axes("x") +
    ylab(NULL) +
    theme(
        panel.border = element_blank(),
        axis.ticks.y = element_blank()
    ) +

    # add a link
    cross_link(
        link_tetragon(
            colnames(gene_data)[
                !is.na(gene_data["TP53", , drop = TRUE])
            ] ~ waiver(),
            .element = element_polygon(
                fill = "orange", color = NA, alpha = 0.8
            )
        ),
        data = ggalign:::ggalign_data_restore(t(gene_data), gene_data),
        size = 0.5, on_top = FALSE
    ) +
    ggplot2::expand_limits(x = c(0, 1), y = c(0, 1)) +

    annotate("text",
        x = -Inf, y = 0.5, hjust = 0,
        label = "TP53", fontface = "bold", size = 8
    ) +
    theme_no_axes() +

    # add plot for TMB of all mutations
    ggalign(
        data = function(data) {
            data <- ggalign_attr(data, "sample_summary")
            # matrix input will be automatically melted into a long foramted data
            # frame in `ggalign()` function.
            as.matrix(data[2:(ncol(data) - 1L)])
        },
        size = 0.8
    ) +
    geom_bar(aes(.x, value, fill = .column_names), stat = "identity") +
    scale_y_log10(expand = expansion(), guide = "none") +
    scale_fill_brewer("Mutations", palette = "Set3", na.translate = FALSE) +
    ylab("TMB")

gene_onco <- ggoncoplot(gene_data,
    filling = FALSE,
    reorder_column = FALSE,
    reorder_row = TRUE
) +
    geom_subtile(aes(fill = value), ncol = 1) +
    scale_fill_brewer("Mutations", palette = "Set3", na.translate = FALSE) +
    theme_no_axes("x") +
    # since legends from geom_tile (oncoPrint body) and `geom_bar`
    # is different, though both looks like the same, the internal
    # won't merge the legends. we remove the legends of oncoPrint body
    guides(fill = "none") +
    scale_y_continuous(position = "right") +

    # add top annotation
    anno_top(size = 0.2, initialize = FALSE) +
    highlight_tp53 +
    align_order2(memo_order) -
    scheme_align(free_spaces = "l") +

    # add right annotation
    anno_right(size = 0.25) -
    # remove bottom spaces of the right annotation when aligning
    scheme_align(free_spaces = "b") +
    # add the text percent for the alterated samples in the right annotation
    ggalign(
        data = function(data) {
            # Atomic vector will be put in the `value` column of the data frame.
            ggalign_attr(data, "gene_summary")$AlteredSamples /
                ggalign_attr(data, "n_samples")
        },
        size = 2
    ) +
    geom_text(aes(1, label = scales::label_percent()(value)), hjust = 1) +
    scale_x_continuous(
        expand = expansion(),
        name = NULL, breaks = NULL,
        limits = c(0, 1)
    ) +
    theme(plot.margin = margin()) +
    # add the bar plot in the right annotation
    ggalign(data = function(data) {
        data <- ggalign_attr(data, "gene_summary")
        # matrix input will be automatically melted into a long foramted data
        # frame in `ggalign()` function.
        as.matrix(data[2:(ncol(data) - 3L)])
    }) +
    geom_bar(aes(value, fill = .column_names),
        stat = "identity",
        orientation = "y"
    ) +
    theme_no_axes("x") +
    scale_fill_brewer("Mutations", palette = "Set3", na.translate = FALSE) +
    xlab("No. of samples") +
    # add bottom annotation
    anno_bottom(size = 0.2) -
    scheme_align(free_spaces = "l") +
    # add bar plot in the bottom annotation
    ggalign(data = function(data) {
        data <- ggalign_attr(data, "titv")$fraction.contribution
        # matrix input will be automatically melted into a long foramted data
        # frame in `ggalign()` function.
        as.matrix(data[2:7])
    }) +
    geom_bar(aes(y = value, fill = .column_names), stat = "identity") +
    ylab("Ti/Tv") +
    scale_fill_brewer("Ti/Tv", palette = "Set2") +
    scale_y_continuous(breaks = scales::pretty_breaks(3L))

p1 <- stack_crossh(pathway_data, sizes = c(0.8, 1, 0.2)) - # 0.5 1 0.2
    scheme_theme(
        plot.margin = margin(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 16, colour = "black"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 18, face = "bold")
    ) +
    # ggoncoplot(gene_data, filling = FALSE) +
    # geom_subtile(aes(fill = value), direction = "v") +
    # theme_no_axes("x") +
    # stack_active() +
    ggoncoplot(pathway_data, reorder_row = TRUE, reorder_column = FALSE) +
    scale_fill_manual(values = "red", guide = "none", na.translate = FALSE) +
    theme_no_axes("x") +
    theme(axis.text.y = element_text(size = 14)) +
    ggtitle("Pathways") +
    scheme_align(free_spaces = "t") +
    anno_top() +
    align_order2(memo_order, reverse = TRUE) +

    stack_active() +
    cross_link(
        link_line(
            !!!pair_links(!!!purrr::imap(
                ggalign_attr(pathway_data, "gene_list"),
                function(gene, pathway) pathway ~ gene
            )),
            .handle_missing = "remove"
        ),
        gene_data,
        size = 0.2
    ) +
    gene_onco

ggsave("figures/genomic/p1.pdf",
    plot = p1, width = 15, height = 8,
    family = "Helvetica"
)

# p2 ---------------------------------------------------
# Please download data from:
# https://www.cell.com/immunity/fulltext/S1074-7613(18)30121-3?utm_campaign=STMJ_1522958526_SC&utm_channel=WEB&utm_source=WEB&dgcid=STMJ_1522958526_SC#fig1
immune_subtypes <- readxl::read_xlsx(
    "rawdata/genomic/immune_subtypes.xlsx",
    na = c("", "NA")
)
nrow(immune_subtypes)

immune_scores <- data.table::fread(
    "rawdata/genomic/Scores_160_Signatures.tsv.gz"
)
mat <- cor(t(immune_scores[, !(1:2)]), method = "spearman")
signatures <- immune_scores[[2L]]
subtypes <- c(
    "Lymphocyte", "Macrophage", "TGF-beta", "IFN-gamma",
    "Wound Healing"
)
signatures[match(
    c(
        "LIexpression_score", "CSF1_response", "TGFB_score_21050467",
        "Module3_IFN_score", "CHANG_CORE_SERUM_RESPONSE_UP"
    ),
    signatures
)] <- subtypes

rownames(mat) <- colnames(mat) <- signatures

subtype_scores <- immune_subtypes %>%
    dplyr::select(
        all_of(c(
            "Immune Subtype",
            "Wound Healing", "Macrophage Regulation",
            "Lymphocyte Infiltration Signature Score",
            "IFN-gamma Response", "TGF-beta Response"
        ))
    ) %>%
    dplyr::filter(!is.na(`Immune Subtype`))
names(subtype_scores)[match(
    c(
        "Lymphocyte Infiltration Signature Score",
        "Macrophage Regulation", "TGF-beta Response",
        "IFN-gamma Response", "Wound Healing"
    ),
    names(subtype_scores)
)] <- subtypes

p2 <- stack_crossv(mat, sizes = c(0.1, 1, 1)) -
    scheme_theme(
        plot.margin = margin(),
        axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 10, colour = "black"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14)
    ) -
    scheme_align(NULL) +
    ggheatmap() +

    # we will highlight the subtypes
    geom_vline(
        aes(xintercept = .data$x),
        data = function(x) {
            x <- dplyr::filter(x, .data$.column_names %in% subtypes) %>%
                dplyr::select(.x) %>%
                dplyr::distinct() %>%
                dplyr::pull(.x)
            data.frame(x = c(x - 0.5, x + 0.5))
        },
        linewidth = 0.1
    ) +
    theme_no_axes() +
    no_expansion() +
    labs(fill = "Spearman\ncorrelation") +

    ## add top dendrogram
    anno_top(size = 0.1) -
    scheme_align(free_spaces = "l") +
    align_dendro(method = "ward.D2") +
    no_expansion("y") +
    theme_no_axes("y") +

    ## add left dendrogram
    anno_left() -
    scheme_align(free_spaces = "b") +
    align_dendro(method = "ward.D2") +
    theme_no_axes("x") +
    no_expansion("x") +

    # step into the parent stack
    stack_active() +

    # the original author don't present the full range of the signatures for
    # each hub signature, we just link the single signature instead
    cross_link(
        link_tetragon(
            !!!lapply(subtypes, function(nm) nm ~ nm),
            .element = element_polygon(
                fill = RColorBrewer::brewer.pal(5L, "Set2"),
                color = NA
            )
        ),
        # the data will be used for the layout
        data = t(subtype_scores[, 2:6]),
        inherit_index = TRUE,
        size = 0.1
    ) +
    align_group(factor(
        names(subtype_scores)[-1L],
        c("Wound Healing", "Lymphocyte", "Macrophage", "TGF-beta", "IFN-gamma")
    )) +

    # add second heatmap
    ggheatmap() +
    labs(fill = "scores") +
    scheme_align(NULL) +
    facet_grid(
        labeller = labeller(.panel_y = function(x) {
            dplyr::if_else(x == "C5", "C5  \n", x)
        })
    ) +
    theme(strip.clip = "off") +

    # split by immune subtypes
    anno_left() +
    align_group(forcats::fct_rev(subtype_scores[["Immune Subtype"]])) +

    # we add bottom annotation
    anno_bottom() -
    scheme_align(free_spaces = "l") +
    ggfree() +
    scheme_align(NULL) +
    scheme_data(function(x) {
        dplyr::mutate(x,
            value = scale(value, center = FALSE),
            scaled_scores = mean(value),
            .by = c(.panel, .extra_panel)
        )
    }) +
    geom_density(aes(value, fill = scaled_scores), color = NA) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_fill_gradient2("scaled\nscores", low = "blue", high = "red") +
    facet_grid(
        rows = vars(.extra_panel),
        cols = vars(.panel), as.table = FALSE,
        axes = "all_x"
    ) +
    theme_bw() +
    theme_no_axes() +
    theme(
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 14),
        strip.background = element_blank(),
        panel.border = element_blank(),
        axis.line.x.bottom = element_line(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold")
    )
ggsave("figures/genomic/p2.pdf",
    plot = p2, width = 7, height = 10,
    family = "Helvetica"
)

# p3 --------------------------------------------------------
# Please download data from:
# https://www.cell.com/cancer-cell/fulltext/S1535-6108(18)30111-9?utm_campaign=STMJ_1522958526_SC&utm_channel=WEB&utm_source=WEB&dgcid=STMJ_1522958526_SC#fig1
arm_alterations <- readxl::read_xlsx(
    "rawdata/genomic/arm_alteration.xlsx",
    skip = 1L
) %>%
    dplyr::summarise(
        dplyr::across(`1p`:`22q`, ~ mean(.x, na.rm = TRUE)),
        .by = c(Type)
    ) %>%
    tibble::column_to_rownames("Type") %>%
    as.matrix()

p3 <- ggheatmap(t(arm_alterations), filling = NULL) -
    scheme_theme(
        plot.margin = margin(),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14)
    ) +
    geom_tile(aes(fill = value), width = 0.8) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        data = function(data) {
            dplyr::summarise(data,
                ymin = min(.y) - 0.5, ymax = max(.y) + 0.5,
                .by = ".x"
            ) %>%
                dplyr::mutate(
                    xmin = .x - 0.4,
                    xmax = .x + 0.4
                )
        },
        color = "black", fill = NA,
        inherit.aes = FALSE
    ) +
    labs(fill = "Mean\nArm Alteration") +
    scale_x_continuous(position = "top") +
    guides(y = primitive_labels(n.dodge = 2)) +
    theme(
        axis.text.x.top = element_text(
            angle = -90, hjust = 0.5, vjust = 0.5
        ),
        plot.margin = margin()
    ) +
    anno_top(size = 0.2) +
    align_dendro(aes(color = branch),
        k = 6L, method = "ward.D2", distance = "spearman"
    ) +
    scale_color_brewer(palette = "Dark2", guide = "none") +
    scale_y_continuous(breaks = scales::pretty_breaks(3L)) +
    no_expansion("y") +
    theme(plot.margin = margin()) +
    anno_bottom() +
    ggmark(
        mark_tetragon(
            .element = element_polygon(
                color = RColorBrewer::brewer.pal(6L, "Dark2"),
                fill = NA
            )
        ),
        mapping = aes(.panel, .column_names),
        group1 = TRUE, size = 1, obs_size = 0.8
    ) +
    scheme_data(function(data) {
        data <- dplyr::summarise(data,
            value = mean(value),
            .by = c(.panel, .column_names, .column_index)
        )
        data$.column_names <- forcats::fct_reorder(
            data$.column_names, data$.column_index
        )
        data
    }) +
    geom_tile(aes(fill = value), width = 0.8) +
    scale_fill_gradient2(low = "blue", high = "red", guide = "none") +
    no_expansion("x") +
    ylab(NULL) +
    theme_no_axes("x") +
    guides(y = primitive_labels(n.dodge = 2)) +
    theme(
        strip.text = element_blank(),
        strip.background = element_blank(),
        plot.margin = margin(l = 0.1, r = 0.2, t = 0.1, unit = "npc"),
        panel.spacing = unit(10, "mm")
    ) +
    facet_wrap(vars(), scales = "free_x") &
    theme(
        panel.background = element_blank(),
        plot.background = element_blank()
    )
ggsave("figures/genomic/p3.pdf",
    plot = p3, width = 8, height = 10,
    family = "Helvetica"
)

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


# p6 ----------------------------------------------
library(data.table)
maf_files <- fs::dir_ls(
    "~/Data/TCGA_data/results",
    regexp = "snv_maf\\.rds$", recurse = TRUE
)
clidata <- readxl::read_xlsx(
    "rawdata/genomic/TCGA-CDR-SupplementalTableS1.xlsx"
)

pathway_all <- lapply(maf_files, function(file) {
    maf <- maftools::read.maf(readRDS(file))
    fortify_matrix(tune(maf))
})
names(pathway_all) <- sub("_snv_maf\\.rds", "", basename(maf_files))
rownames(pathway_all[[2L]])
genome_integrity_os <- lapply(pathway_all, function(pathway) {
    group <- !is.na(pathway["Genome integrity", , drop = TRUE])
    data <- dplyr::filter(
        clidata,
        bcr_patient_barcode %in% substr(names(group), 1L, 12)
    )
    data$Genome_integrity <- dplyr::if_else(
        data$bcr_patient_barcode %in% substr(names(group)[group], 1L, 12),
        "Alt", "Wild"
    )
    data$Genome_integrity <- factor(data$Genome_integrity, c("Wild", "Alt"))
    out <- survival::coxph(
        survival::Surv(OS.time, OS) ~ Genome_integrity,
        data = data
    )
    out <- broom::tidy(out, conf.int = TRUE)
    out$Genome_integrity <- sum(data$Genome_integrity == "Alt")
    out$Wild <- sum(data$Genome_integrity == "Wild")
    out
})
forest_data <- dplyr::bind_rows(genome_integrity_os, .id = "project")
dplyr::filter(forest_data, p.value < 0.05)

p6 <- stack_alignh(rep(forest_data$project, each = 2L)) -
    scheme_theme(
        plot.margin = margin(t = 5L),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        strip.text.y = element_text(size = 14)
    ) +

    align_order(function(x) rev(seq_len(nrow(x)))) +

    ggalign(
        vctrs::vec_interleave(forest_data$Wild, forest_data$Genome_integrity)
    ) +
    geom_bar(aes(value, fill = factor(.index %% 2L)),
        stat = "identity", orientation = "y"
    ) +
    scale_x_reverse(name = "Number of samples") +
    scale_fill_brewer(
        palette = "Dark2", direction = -1L,
        name = "Genome integrity\nStatus",
        labels = c("Alt", "Wt")
    ) +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.6, 0.3),
        legend.justification.inside = c(1, 1)
    ) +

    # Proejct text
    ggalign(vctrs::vec_interleave(forest_data$project, NA), size = 0.2) +
    geom_rect(
        aes(
            ymin = .index - 0.5, ymax = .index + 0.5,
            xmin = -Inf, xmax = Inf, fill = fill
        ),
        data = function(d) {
            d$fill <- factor(d$.index %% 2L)
            d
        }
    ) +
    scale_fill_manual(values = c("#F0F0F0", "white"), guide = "none") +
    geom_text(aes(1L, label = toupper(value)), data = function(d) {
        dplyr::filter(d, !is.na(value))
    }) +
    theme_no_axes("y") +
    scale_x_continuous(breaks = NULL) +
    xlab("Project") +

    # Status
    # ggalign(rep(c("Wt", "Alt"), times = nrow(forest_data)), size = 0.2) +
    # geom_rect(
    #     aes(
    #         ymin = .index - 0.5, ymax = .index + 0.5,
    #         xmin = -Inf, xmax = Inf, fill = fill
    #     ),
    #     data = function(d) {
    #         d$fill <- factor(d$.index %% 2L)
    #         d
    #     }
    # ) +
    # scale_fill_manual(values = c("#F0F0F0", "white"), guide = "none") +
    # geom_text(aes(1L, label = value), size = 3, data = function(d) {
    #     dplyr::filter(d, !is.na(value))
    # }) +
    # theme_no_axes("y") +
    # scale_x_continuous(breaks = NULL) +
    # xlab("\nGenome integrity\nStatus") +

    # P value
    ggalign(vctrs::vec_interleave(NA, forest_data$p.value), size = 0.2) +
    geom_rect(
        aes(
            ymin = .index - 0.5, ymax = .index + 0.5,
            xmin = -Inf, xmax = Inf, fill = fill
        ),
        data = function(d) {
            d$fill <- factor(d$.index %% 2L)
            d
        }
    ) +
    scale_fill_manual(values = c("#F0F0F0", "white"), guide = "none") +
    geom_text(aes(1L, label = sprintf("%.2f", value)), data = function(d) {
        dplyr::filter(d, !is.na(value))
    }) +
    theme_no_axes("y") +
    scale_x_continuous(breaks = NULL) +
    xlab("Pvalue") +

    ggalign() +
    geom_rect(
        aes(
            ymin = .index - 0.5, ymax = .index + 0.5,
            xmin = -Inf, xmax = Inf, fill = fill
        ),
        data = function(d) {
            d$fill <- factor(d$.index %% 2L)
            d
        }
    ) +
    scale_fill_manual(values = c("#F0F0F0", "white"), guide = "none") +

    geom_point(aes(0), shape = 15L, data = function(d) {
        dplyr::filter(d, .index %% 2L > 0L)
    }) +
    geom_segment(aes(log2(exp(conf.low)), xend = log2(exp(conf.high))), data = function(d) {
        d <- dplyr::filter(d, .index %% 2L == 0L)
        out <- dplyr::inner_join(
            d, dplyr::select(forest_data, project, conf.low, conf.high),
            by = c(value = "project")
        )
        dplyr::filter(out, value != "pcpg")
    }) +
    geom_point(aes(log2(exp(estimate))), data = function(d) {
        d <- dplyr::filter(d, .index %% 2L == 0L)
        out <- dplyr::inner_join(
            d, dplyr::select(forest_data, project, estimate),
            by = c(value = "project")
        )
        dplyr::filter(out, value != "pcpg")
    }) +
    xlab("log2(HR)") &
    coord_cartesian(clip = "off")

ggsave("figures/genomic/p6.pdf",
    plot = p6, width = 7, height = 7,
    family = "Helvetica"
)
