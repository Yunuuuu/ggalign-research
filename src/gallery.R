# pak::pak("Yunuuuu/biotidy"): we use biotidy to create data
library(ggalign)
# Prepare gallery figures to showcase the general style of each feature.
# Each figure uses minimal data—only enough to illustrate the design.

# https://academic.oup.com/mbe/article/38/9/4039/6294410
alltax <- read.csv("rawdata/circular/all_stain_taxonomy.csv")
trda <- MicrobiotaProcess::convert_to_treedata(alltax)
trda <- ape::as.phylo(trda)
template_mat <- matrix(
    NA,
    nrow = length(trda$tip.label), ncol = length(trda$tip.label),
    dimnames = list(trda$tip.label, trda$tip.label)
)

linktab <- read.csv("rawdata/circular/Interaction_link_tab.csv")
all(linktab$Inhibitor %in% trda$tip.label)
all(linktab$Sensitive %in% trda$tip.label)

weighttab <- read.csv("rawdata/circular/Interaction_weight.csv")
tippoint <- read.csv("rawdata/circular/stain_tippoint.csv")
tippoint$Taxa <- factor(tippoint$Taxa,
    levels = c(
        "Actinobacteria",
        "Bacteroidetes",
        "Firmicutes",
        "Deinococcus-Thermus",
        "Alphaproteobacteria",
        "Betaproteobacteria",
        "Gammaproteobacteria"
    )
)
tippoint$names <- gsub("s__Leaf", "", tippoint$Isolation)

BGCsda <- read.csv("rawdata/circular/BGCs_heatmap.csv")
BGCsda$BGCs <- factor(BGCsda$BGCs,
    levels = c(
        "modular.PKS",
        "modular.PKS.NRPS.hybrid",
        "non_modular.PKS", "NRPS",
        "RiPP",
        "Quorum.sensing",
        "terpene",
        "other"
    )
)
BGCsda$Count <- log2(BGCsda$Count + 1)

circle <- circle_discrete(trda,
    radial = coord_radial(
        rotate.angle = TRUE, expand = FALSE,
        start = pi / 2 + pi / 30, end = pi * 2 + pi / 2,
        reverse = "theta"
    ),
    theme = theme(
        legend.background = element_rect(fill = NA),
        legend.title = element_text(size = 6.5),
        legend.text = element_text(size = 5),
        legend.spacing.y = unit(0, "mm"),
        legend.margin = margin(),
        legend.box.margin = margin(),
        legend.box.spacing = unit(0, "mm"),
        plot.margin = margin(l = -15, t = -25, b = -25, r = 2, unit = "mm")
    )
) +
    # inner links ------------------------------------
    ggalign(rownames, size = 2) + # we use rownames to match the data
    scheme_data(function(data) {
        # for data, value is the rownames, we match it with linktab
        dplyr::mutate(
            linktab,
            x = data$.x[match(Inhibitor, data$.names)],
            xend = data$.x[match(Sensitive, data$.names)]
        )
    }) +
    # Currently, we use `curvelink` from ggtree, but we'll provide a native
    # implementation of links for CircleLayout once
    # https://github.com/tidyverse/ggplot2/pull/6421 is merged into ggplot2.
    ggtree:::geom_curvelink(
        mapping = aes(x, 1, xend = xend, yend = 1, color = Interaction),
        outward = FALSE, key_glyph = draw_key_path,
        linewidth = 0.1
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = NULL) +
    scale_colour_manual(
        values = c("chocolate2", "#3690C0", "#009E73"),
        guide = "none"
    ) +
    theme(panel.spacing.r = unit(0, "mm")) +
    no_expansion("y") +

    # add a point shape for each tip -------------------
    ggalign(size = 0.1) +
    geom_point(aes(.x, 1, shape = Taxa, color = Taxa),
        data = function(data) {
            dplyr::inner_join(data, tippoint, by = c(.names = "Isolation"))
        }
    ) +
    scale_shape_manual(
        values = rep(c(16, 18), c(4, 3)),
        guide = "none"
    ) +
    scale_color_manual(
        values = c(
            "#EF3B2C", "#1D91C0", "#FEB24C", "grey60",
            "#7FBC41", "#4D9221", "#276419"
        ),
        guide = "none"
    ) +
    theme(panel.spacing.r = unit(0, "mm")) +
    theme_no_axes("y") +
    no_expansion("y") +

    # draw the Phylogenetics tree ------------------------
    align_phylo(trda,
        tree_type = "cladogram", size = 0.5,
        ladderize = "left"
    ) +
    scale_y_continuous(breaks = NULL) +
    theme(panel.spacing.r = unit(0, "mm")) +
    no_expansion("y") +

    # draw the bar plot ---------------------------------
    ggalign(rownames) + # we use rownames to match the data
    scheme_data(function(data) {
        # for data, value is the rownames, we match it with weighttab
        dplyr::inner_join(data,
            dplyr::rename(weighttab, Type = Number, Number = value),
            by = c(value = "Strain")
        )
    }) +
    geom_bar(aes(y = Number, fill = Type), stat = "identity") +
    scale_fill_manual(
        name = "Strength of interactions",
        values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
        guide = "none"
    ) +
    no_expansion("y") +
    scale_y_continuous(breaks = scales::breaks_pretty(3)) +
    theme(axis.text.y = element_text(size = 20)) &
    theme(
        panel.background = element_blank(),
        text = element_text(family = "Helvetica")
    )
ggsave("figures/gallery/circle.pdf",
    plot = circle, width = 8.5, height = 7,
    family = "Helvetica"
)

# Circle splitted ---------------------------------------
# https://pmc.ncbi.nlm.nih.gov/articles/PMC9856581/#app1-cancers-15-00342
data <- readxl::read_xlsx(
    "rawdata/circular/Supplementary tables.xlsx",
    skip = 1L
)
set.seed(3L)
selected <- sample(2:18, 5)
logFC <- data[selected, c(1, which(as.character(data[1, ]) == "logFC"))]
adj.P <- data[selected, c(1, which(as.character(data[1, ]) == "adj.P"))]
names(adj.P) <- names(logFC)

logFC <- dplyr::mutate(logFC, dplyr::across(!pathway, as.numeric)) |>
    tibble::column_to_rownames(var = "pathway") |>
    as.matrix() |>
    t()
adj.P <- dplyr::mutate(adj.P, dplyr::across(!pathway, as.numeric)) |>
    tibble::column_to_rownames(var = "pathway") |>
    as.matrix() |>
    t()
circle_splitted <- circle_discrete(logFC,
    radial = coord_radial(start = pi / 2, end = pi * 2, expand = FALSE),
    theme = theme(
        plot.margin = margin(l = 5, t = 5, b = 6, unit = "mm"),
        legend.box.spacing = unit(5, "mm")
    ),
    sector_spacing = 5 * pi / 180
) +
    align_dendro(aes(color = branch), k = 5L, size = 1) +
    theme_no_axes("y") +
    scale_color_brewer(
        palette = "Dark2",
        guide = guide_legend(theme = theme(
            legend.text = element_text(size = 18),
            legend.title = element_text(face = "bold", size = 16)
        ))
    ) +

    # heatmap plot
    ggalign(mapping = aes(y = .column_names, fill = value)) +
    geom_tile(width = 1, height = 1) +
    geom_text(aes(label = "*"), data = function(dd) {
        dd$pvalue <- adj.P[cbind(dd$.row_index, dd$.column_index)]
        dplyr::filter(dd, pvalue < 0.05)
    }) +
    scale_fill_gradient2(
        low = "blue", high = "red",
        name = "logFC",
        # option = "plasma",
        limits = c(-1, 1),
        breaks = c(-0.5, 0, 0.5),
        guide = guide_colorbar(theme = theme(
            legend.text = element_text(size = 18),
            legend.title = element_text(face = "bold", size = 16)
        ))
    ) +
    theme(
        axis.text.r = element_text(family = "Helvetica", size = 22),
        axis.text.theta = element_text(family = "Helvetica", size = 15)
    ) +
    guides(
        r = "none", r.sec = "axis",
        theta = guide_axis_theta(angle = 0)
    ) &
    theme(
        text = element_text(family = "Helvetica"),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold", size = 14)
    )
ggsave("figures/gallery/circle_splitted.pdf",
    plot = circle_splitted, width = 9, height = 8,
    family = "Helvetica"
)

# marginal plot -----------------------------------------
i2 <- iris
i2$Group <- rep(c("Group1", "Group2"), 75)
side_plot <- ggside(i2, aes(Sepal.Width, Sepal.Length, color = Species)) +
    geom_point(size = 2) +
    facet_grid(Species ~ Group) +
    theme(
        strip.background = element_blank(),
        strip.text = element_blank()
    ) +
    anno_top(size = 0.3) +
    ggalign() +
    geom_density(aes(Sepal.Width, y = after_stat(density), colour = Species),
        position = "stack"
    ) +
    facet_grid(cols = vars(Group)) +
    theme(
        strip.text = element_text(size = 20, margin = margin(5, 5, 5, 5)),
        strip.background = element_rect(fill = "grey")
    ) +
    anno_right(size = 0.3) +
    ggalign() +
    geom_density(aes(x = after_stat(density), Sepal.Length, colour = Species),
        position = "stack"
    ) +
    facet_grid(rows = vars(Species)) +
    theme(
        strip.text = element_text(size = 20, margin = margin(5, 5, 5, 5)),
        strip.background = element_rect(fill = "grey")
    ) -
    with_quad(scheme_theme(theme_bw()), "tr", main = TRUE) &
    scale_color_brewer(palette = "Dark2", guide = "none") &
    scale_x_continuous(breaks = scales::pretty_breaks(3L)) &
    scale_y_continuous(breaks = scales::pretty_breaks(3L)) &
    theme(
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20, face = "bold")
    )
ggsave("figures/gallery/side_plot.pdf",
    plot = side_plot, width = 8, height = 8,
    family = "Helvetica"
)

# heatmap -----------------------------------------------
set.seed(123)
small_mat <- log1p(SummarizedExperiment::assay(biotidy::mockSE(9, 9)))
heatmap <- ggheatmap(small_mat) +
    scale_fill_viridis_c(
        option = "plasma",
        breaks = scales::breaks_pretty(3L)
    ) +
    guides(fill = guide_colorbar(theme = theme(
        legend.text = element_text(size = 20, family = "Helvetica")
    ))) +
    labs(fill = NULL) +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank()
    ) +

    # add a dendrogram
    anno_top(size = 0.2) +
    align_dendro(aes(color = branch), k = 3) +
    geom_point(aes(color = branch, y = y)) +
    scale_y_continuous(name = NULL, breaks = scales::breaks_pretty(2)) +
    theme(axis.text.y = element_text(
        family = "Helvetica", size = 20, colour = "black"
    )) +

    # add a dendrogram
    anno_right(size = 0.2) +
    align_dendro(aes(color = branch), k = 3) +
    geom_point(aes(color = branch, y = y)) +
    scale_x_continuous(name = NULL, breaks = scales::breaks_pretty(2)) +
    theme(axis.text.x = element_text(
        family = "Helvetica", size = 20,
        colour = "black"
    )) &
    scale_color_brewer(palette = "Dark2", guide = "none") &
    theme(
        plot.margin = margin(),
        plot.background = element_blank()
    )
ggsave("figures/gallery/heatmap.pdf",
    plot = heatmap, width = 6, height = 6, family = "Helvetica"
)

heatmap3d <- ggheatmap(small_mat,
    filling = FALSE,
    theme = theme(
        panel.spacing.x = unit(4, "mm"),
        legend.box.spacing = unit(5, "mm")
    )
) +
    geom_tile3d(
        aes(fill = value, z = value, width = 0.8, height = 0.8),
        color = "black"
    ) +
    scale_fill_viridis_c(
        name = NULL,
        option = "plasma",
        breaks = scales::breaks_pretty(3L)
    ) +
    guides(fill = guide_colorbar(theme = theme(
        legend.text = element_text(
            size = 18,
            family = "Helvetica",
            face = "plain"
        )
    ))) +
    coord_cartesian(clip = "off") +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(t = 20)
    ) +
    no_expansion() +

    # bottom annotation --------------------------------
    anno_bottom(size = 0.2) -
    scheme_align(free_spaces = "lr") +

    # add a pie
    ggalign(function(d) rowMeans(d > 0)) +
    ggalign::geom_pie(
        aes(
            y = 1L, angle = value * 360, fill = value,
            color = after_scale(fill)
        )
    ) +
    scale_y_continuous(limits = c(0.5, 1.5), breaks = NULL, name = NULL) +
    scale_fill_viridis_c("Percent",
        limits = c(0, 1), breaks = scales::pretty_breaks(3L),
        guide = guide_colorbar(
            theme = theme(
                legend.text = element_text(
                    family = "Helvetica", size = 18,
                    colour = "black"
                ),
                legend.title = element_text(
                    family = "Helvetica", size = 20,
                    colour = "black"
                )
            )
        )
    ) +
    theme(plot.margin = margin()) +

    # add a dendrogram
    align_dendro(aes(color = branch), k = 3) +
    geom_point(aes(color = branch, y = y)) +
    scale_y_continuous(
        name = NULL, breaks = scales::breaks_pretty(2),
        expand = expansion(),
        position = "right"
    ) +
    theme(
        axis.text.y = element_text(
            family = "Helvetica", size = 20, colour = "black"
        ),
        plot.margin = margin()
    ) +

    # left anontation ------------------------------------
    anno_left(size = 0.2) -
    scheme_align(free_spaces = "b") +

    # add a dendrogram
    align_dendro(aes(color = branch), k = 3) +
    geom_point(aes(color = branch, y = y)) +
    scale_x_continuous(
        name = NULL,
        breaks = c(0, 15),
        expand = expansion()
    ) +
    theme(
        axis.text.x = element_text(
            family = "Helvetica", size = 20,
            colour = "black"
        ),
        plot.margin = margin()
    ) +

    # add a pie
    ggalign(function(d) rowMeans(d > 0)) +
    ggalign::geom_pie(aes(
        x = 1L, angle = value * 360, fill = value,
        color = after_scale(fill)
    )) +
    scale_x_continuous(limits = c(0.5, 1.5), breaks = NULL, name = NULL) +
    scale_fill_viridis_c("Percent",
        limits = c(0, 1), breaks = scales::pretty_breaks(3L),
        guide = guide_colorbar(
            theme = theme(
                legend.text = element_text(
                    family = "Helvetica", size = 18,
                    colour = "black"
                ),
                legend.title = element_text(
                    family = "Helvetica", size = 20,
                    colour = "black"
                )
            )
        )
    ) +
    theme(plot.margin = margin()) &

    scale_color_brewer(palette = "Dark2", guide = "none") &
    theme(
        plot.background = element_blank(),
        panel.background = element_blank()
    )

ggsave("figures/gallery/heatmap3d.pdf",
    plot = heatmap3d, width = 7, height = 6.5, family = "Helvetica"
)

# Oncoplot ----------------------------------------------
# load data from `maftools`
laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
# clinical information containing survival information and histology. This is optional
laml.clin <- system.file("extdata", "tcga_laml_annot.tsv", package = "maftools")
laml <- maftools::read.maf(
    maf = laml.maf,
    clinicalData = laml.clin,
    verbose = FALSE
)
set.seed(25) # 2, 11, 20, 25, 29, 32
new_laml <- maftools::subsetMaf(laml,
    tsb = sample(maftools::getClinicalData(laml)$Tumor_Sample_Barcode, 20)
)
oncoplot <- ggoncoplot(
    new_laml,
    n_top = 10, collapse_vars = FALSE, filling = FALSE,
    remove_empty_samples = FALSE
) +
    geom_subtile(aes(fill = value),
        ncol = 1L, height = 0.8,
        show.legend = TRUE
    ) +
    theme_no_axes("x") +

    # add top annotation -------------------------------
    anno_top(size = 0.2) +
    ggalign(
        data = function(data) {
            vars <- ggalign_lvls(data)
            data <- ggalign_attr(data, "sample_summary")
            as.matrix(data[intersect(names(data), vars)])
        }
    ) +
    scheme_data(function(data) {
        data$.column_names <- factor(data$.column_names, ggalign_lvls(data))
        data
    }) +
    geom_bar(aes(.x, value, fill = .column_names),
        stat = "identity"
    ) +
    guides(fill = "none") +
    # labs(fill = "top") +
    ylab("TMB") +
    scale_y_continuous(breaks = scales::breaks_pretty(2L)) +

    # add right annotation -----------------------------
    anno_right(size = 0.2) -
    # remove bottom spaces of the right annotation when aligning
    scheme_align(free_spaces = "b") +

    # add the text percent for the alterated samples in the right annotation
    ggalign(data = function(data) {
        # Atomic vector will be put in the `value` column of the data frame.
        ggalign_attr(data, "gene_summary")$AlteredSamples /
            ggalign_attr(data, "n_samples")
    }) +
    geom_text(aes(1, label = scales::label_percent()(value)),
        hjust = 1.2, size = 5, family = "Helvetica"
    ) +
    scale_x_continuous(
        expand = expansion(),
        name = NULL, breaks = NULL,
        limits = c(0, 1)
    ) +
    theme(plot.margin = margin()) +

    # add the bar plot in the right annotation
    ggalign(data = function(data) {
        vars <- ggalign_lvls(data)
        data <- ggalign_attr(data, "variant_weights")
        as.matrix(data[intersect(names(data), vars)])
    }) +
    scheme_data(function(data) {
        data$.column_names <- factor(data$.column_names, ggalign_lvls(data))
        data
    }) +
    geom_bar(aes(value, fill = .column_names),
        stat = "identity",
        orientation = "y"
    ) +
    guides(fill = "none") +
    # labs(fill = "right") +
    xlab("Variant weights") +
    scale_x_continuous(breaks = scales::breaks_pretty(1)) +
    theme(
        axis.line.y = element_line(),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_blank()
    ) -
    # we apply the scale mapping to
    # - the top and right annotation: `position = "tr"`
    # - the main plot: `main = TRUE`
    with_quad(
        scale_fill_brewer(
            "Mutations",
            palette = "Dark2",
            na.translate = FALSE, drop = FALSE
        ),
        position = "tr",
        main = TRUE
    ) +

    # add bottom annotation ----------------------------
    anno_bottom(size = 0.2) +
    # add bar plot in the bottom annotation
    ggalign(data = function(data) {
        data <- ggalign_attr(data, "titv")$fraction.contribution
        as.matrix(data[-1L])
    }) +
    geom_bar(aes(y = value, fill = .column_names), stat = "identity") +
    ylab("Ti/Tv") +
    scale_y_continuous(breaks = scales::breaks_pretty(3L)) +
    scale_fill_brewer("Ti/Tv", palette = "Set2", na.translate = FALSE) -
    with_quad(
        theme(
            axis.text.y = element_text(size = 18),
            axis.title.y = element_text(size = 20, face = "bold"),
            legend.text = element_text(size = 12, colour = "black"),
            legend.title = element_text(size = 20)
        )
    ) &
    theme(
        plot.margin = margin(b = 5),
        text = element_text(family = "Helvetica"),
        plot.background = element_blank()
    )
ggsave(
    "figures/gallery/oncoplot.pdf",
    plot = oncoplot, family = "Helvetica",
    width = 9, height = 6
)

# upset ---------------------------------------------------------
upset_data <- dplyr::filter(BGCsda, Count > 0)
set.seed(1L)
sampled <- sample(unique(upset_data$Strain), 50)
upset_data <- split(upset_data$Strain, upset_data$BGCs)
upset_data <- lapply(upset_data, intersect, sampled)
upset_data$modular.PKS.NRPS.hybrid <- NULL
upset <- ggupset(tune(upset_data),
    point = list(
        data = function(d) {
            dplyr::mutate(d,
                point_group = dplyr::case_when(
                    !point_group ~ NA_character_,
                    .default = .discrete_y
                ),
                point_group = factor(point_group, levels(.discrete_y))
            )
        },
        size = 3
    )
) +
    scale_fill_manual(values = c("#F0F0F0", "white"), guide = "none") +
    scale_color_brewer(
        palette = "Dark2", guide = "none",
        na.translate = TRUE, na.value = "grey"
    ) +
    theme(axis.text.y = element_text(size = 16)) +
    # scale_color_manual(values = c("grey", "black"), guide = "none") +

    # add top annotation
    anno_top(size = 0.5) +
    ggalign(data = function(d) ggalign_attr(d, "intersection_sizes")) +
    geom_bar(aes(y = .data$value), stat = "identity") +
    ylab("Intersection\nSizes") +
    theme(
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16)
    ) +

    # add right annotation
    anno_right(size = 0.5) +
    ggalign(data = function(d) ggalign_attr(d, "set_sizes")) +
    geom_bar(aes(x = .data$value, fill = .discrete_y),
        stat = "identity",
        orientation = "y"
    ) +
    scale_fill_brewer(
        palette = "Dark2", guide = "none",
        na.translate = TRUE, na.value = "grey"
    ) +
    xlab("Set Sizes") +
    scale_x_continuous(breaks = scales::pretty_breaks(3)) +
    theme(
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16)
    ) &
    theme(plot.background = element_blank())
ggsave("figures/gallery/upset.pdf", plot = upset, width = 8, height = 4)

# Link observations -------------------------------------
# Create some sample data
set.seed(1L)
transcriptomics <- log1p(SummarizedExperiment::assay(biotidy::mockSE(9, 9)))

set.seed(2)
proteomics <- log1p(SummarizedExperiment::assay(biotidy::mockSE(9, 9)))
proteomics[c(1:3, 5, 8:9), ] <- transcriptomics[c(1:3, 5, 8, 4), ]

# tanglegram
link_observations <- stack_crossh() +
    # align_kmeans(data = cbind(transcriptomics, proteomics), centers = 3L) +

    # add heatmap ---------------------------
    ggheatmap(transcriptomics) +
    ggtitle("Transcriptomics") +
    theme_no_axes("x") +
    theme(
        plot.margin = margin(),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 20)
    ) +

    ## add left annotations -----------------
    anno_left(size = 0.25) +
    align_dendro(aes(color = branch), merge_dendrogram = TRUE, k = 3L) +
    theme_no_axes("x") +

    # add Links -----------------------------
    stack_active() +
    cross_link(
        link_line(.element = element_line(
            colour = rep(RColorBrewer::brewer.pal(3, "Dark2"), c(4, 3, 2))[
                order(c(7, 4, 6, 9, 8, 3, 2, 1, 5))
            ]
        )),
        size = 0.5
    ) +
    theme(plot.margin = margin()) +

    # add another heatmap -------------------
    ggheatmap(proteomics) +
    ggtitle("Proteomics") +
    scale_y_continuous(position = "right") +
    theme_no_axes("x") +
    theme(
        plot.margin = margin(),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 20)
    ) +

    ## add right annotations ----------------
    anno_right(size = 0.25) +
    align_dendro(aes(color = branch), merge_dendrogram = TRUE, k = 3L) +
    theme_no_axes("x") &

    # set scales
    scale_fill_viridis_c(
        name = "log2(expr + 1)",
        option = "plasma",
        breaks = scales::breaks_pretty(3L),
        guide = guide_colorbar(theme = theme(
            legend.text = element_text(size = 18),
            legend.title = element_text(face = "bold", size = 16)
        ))
    ) &
    scale_color_brewer(palette = "Dark2", guide = "none")
ggsave(
    "figures/gallery/link_observations.pdf",
    plot = link_observations, family = "Helvetica",
    width = 9, height = 6
)

# Annotate observations -------------------------------------
# Create some sample data
set.seed(1L)
expr1 <- log1p(SummarizedExperiment::assay(biotidy::mockSE(9, 9)))

set.seed(2)
expr2 <- log1p(SummarizedExperiment::assay(biotidy::mockSE(9, 9)))

diff_genes <- dplyr::inner_join(
    fortify_data_frame(expr1),
    fortify_data_frame(expr2),
    by = c(".row_names", ".column_names")
) |>
    dplyr::select(.row_names, .column_names, x = value.x, y = value.y) |>
    dplyr::arrange(.row_names) |>
    dplyr::summarise(
        pvalue = t.test(x, y, paired = TRUE)$p.value,
        .by = .row_names
    ) |>
    dplyr::filter(pvalue < 0.05) |>
    dplyr::pull(.row_names)

annotate_observations <- stack_crossh(expr1) -
    scheme_theme(plot.margin = margin()) +
    align_dendro(aes(color = branch), k = 3L, size = 0.25) +
    theme_no_axes("x") +
    scale_x_reverse() +
    ggheatmap(expr1) +
    ggtitle("Group 1") +
    scale_fill_viridis_c(
        name = "log2(expr + 1)",
        option = "plasma",
        breaks = scales::breaks_pretty(3L),
        guide = guide_colorbar(theme = theme(
            legend.text = element_text(size = 18),
            legend.title = element_text(face = "bold", size = 16)
        ))
    ) +
    theme_no_axes("x") +
    theme(
        plot.margin = margin(),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 20)
    ) +
    stack_active() +

    # initialize the layout data and add mark plot
    cross_mark(
        mark_tetragon(!!!lapply(diff_genes, function(x) x ~ x),
            .reorder = "hand1",
            .element = element_polygon(
                fill = rep(RColorBrewer::brewer.pal(6, "Dark2"), each = 2L)
            )
        ),
        data = expr2,
        inherit_panel = TRUE
    ) +
    scheme_data(function(d) {
        left <- fortify_data_frame(expr1) |>
            dplyr::mutate(.hand = "left")
        right <- fortify_data_frame(expr2) |>
            dplyr::mutate(.hand = "right")
        dplyr::inner_join(d,
            dplyr::bind_rows(left, right),
            by = c(.names = ".row_names", ".hand")
        )
    }) +
    geom_boxplot(aes(.hand, value, fill = .hand),
        outlier.size = 0.2,
        linewidth = 0.1
    ) +
    scale_fill_brewer(palette = "Set1", guide = "none") +
    theme_no_axes("xy") +
    theme(
        plot.margin = margin(30, 20, 30, 20),
        strip.background = element_blank(),
        strip.text = element_blank()
    ) +

    # add a heatmap
    ggheatmap() +
    ggtitle("Group 2") +
    scale_fill_viridis_c(
        name = "log2(expr + 1)",
        option = "plasma",
        breaks = scales::breaks_pretty(3L),
        guide = guide_colorbar(theme = theme(
            legend.text = element_text(size = 18),
            legend.title = element_text(face = "bold", size = 16)
        ))
    ) +
    scale_y_continuous(position = "right") +
    theme_no_axes("x") +
    theme(
        plot.margin = margin(),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 20)
    ) +

    stack_active() +
    align_dendro(aes(color = branch), size = 0.25) +
    theme_no_axes("x") &
    scale_color_brewer(palette = "Dark2", guide = "none")
ggsave(
    "figures/gallery/annotate_observations.pdf",
    plot = annotate_observations, family = "Helvetica",
    width = 9, height = 6
)
