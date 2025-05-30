library(ggalign)
library(mia)
if (!dir.exists("figures/microbiome")) dir.create("figures/microbiome")
# pak::pak("Yunuuuu/ggalign")
# pak::pak("tidyverse/ggplot2")

phylum <- data.table::fread("rawdata/bacteria.sample.relabund.phylum.txt")
metadata <- data.table::fread("rawdata/sample_metadata.txt")
table(metadata$sample_type)
metadata <- dplyr::select(metadata, V1, project, sample_type)
subset <- metadata$sample_type == "Primary Tumor"

phylum <- tibble::column_to_rownames(phylum, var = "name")
phylum <- t(as.matrix(phylum))
all(metadata$V1 == rownames(phylum))
phylum <- phylum[subset, ]
phylum_tree <- TreeSummarizedExperiment(
    assays = SimpleList(relabundance = t(phylum)),
    colData = DataFrame(metadata[subset, ])
)
dim(phylum)
head(sort(colSums(phylum), decreasing = TRUE))
names(sort(colSums(phylum)))[1:10]

genus <- data.table::fread("rawdata/bacteria.sample.relabund.genus.txt")
genus <- tibble::column_to_rownames(genus, var = "name")
genus <- t(as.matrix(genus))
all(metadata$V1 == rownames(genus))
genus <- genus[subset, ]

genus_tree <- TreeSummarizedExperiment(
    assays = SimpleList(relabundance = t(genus)),
    colData = DataFrame(metadata[subset, ])
)
genus_tree <- mia::addDominant(
    genus_tree,
    assay.type = "relabundance", name = "dominant_taxa"
)

metadata <- metadata[subset]

# p1
p1 <- stack_alignv(phylum) -
    scheme_theme(
        plot.margin = margin(),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16, face = "bold")
    ) +
    align_group(metadata$project) +

    # a new dendrogram
    align_dendro(aes(color = branch),
        size = 0.2, reorder_dendrogram = TRUE,
        merge_dendrogram = TRUE
    ) +
    scale_color_brewer(palette = "Dark2", name = "Project") +

    scale_y_continuous(
        expand = expansion(),
        breaks = scales::pretty_breaks(2L)
    ) +

    # add bar plot
    ggalign() +
    geom_bar(aes(y = value, fill = .column_names),
        stat = "identity", position = "stack"
    ) +
    scale_fill_brewer(palette = "Set3", name = "Phylum") +
    scale_y_continuous(
        expand = expansion(), breaks = NULL,
        name = "Relative abundance"
    ) +

    # mark observations
    ggmark(
        mark_tetragon(
            .element = element_polygon(
                fill = RColorBrewer::brewer.pal(6L, "Dark2"),
                alpha = 0.5
            )
        ),
        mapping = aes(.column_names, value),
        size = 0.5
    ) +
    geom_boxplot(aes(fill = .column_names), outlier.size = 0.2) +
    scale_fill_brewer(palette = "Set3", name = "Phylum", guide = "none") +
    facet_wrap(vars(), strip.position = "top") +
    scale_y_continuous(
        expand = expansion(),
        breaks = scales::pretty_breaks(3L)
    ) +
    theme(plot.margin = margin(t = 0.1, r = 0.05, unit = "npc")) +
    theme_no_axes("x") +
    ylab("Relative\nabundance") +
    theme(
        strip.background = element_blank(),
        strip.text = element_blank()
    )
ggsave("figures/microbiome/p1.pdf",
    plot = p1, width = 8, height = 7,
    family = "Helvetica"
)

obs_phylum <- t(phylum)
all(metadata$V1 == colnames(obs_phylum))
p2 <- stack_crossv(obs_phylum) -
    scheme_theme(
        theme_bw(),
        theme(
            plot.margin = margin(),
            axis.title.y = element_text(size = 16, face = "bold"),
            axis.text.y = element_text(size = 16, colour = "black"),
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 16, face = "bold")
        )
    ) +
    purrr::imap(
        split(seq_len(nrow(metadata)), metadata$project),
        function(index, project) {
            force(index)
            list(
                ggalign(~ .x[, index], aes(y = value)) +
                    geom_boxplot(aes(group = .x, fill = .row_names),
                        outlier.size = 0.5
                    ) +
                    ggbeeswarm::geom_beeswarm(size = 0.1) +
                    scale_fill_brewer(palette = "Set3", name = "Phylum") +
                    # ggtitle(project) +
                    scale_y_continuous(
                        name = project,
                        limits = c(0, 1),
                        breaks = c(0, 1)
                    ),
                align_order(~ rowMeans(.x[, index])),
                if (project != "STAD") {
                    ggcross(
                        aes(y = .hand, group = .index, color = .names),
                        size = 0.2
                    ) +
                        geom_line(linewidth = 2) +
                        theme_no_axes("y", line = TRUE) +
                        scale_y_discrete(expand = expansion()) +
                        scale_color_brewer(
                            palette = "Set3", name = "Phylum",
                            guide = "none"
                        ) +
                        theme_no_axes("x") +
                        theme(panel.border = element_blank())
                }
            )
        }
    )

ggsave("figures/microbiome/p2.pdf",
    plot = p2, width = 7, height = 7,
    family = "Helvetica"
)

# genus - PCoA -----------------------------------------
genus_tree <- scater::runMDS(
    genus_tree,
    FUN = getDissimilarity,
    name = "PCoA_BC",
    method = "bray",
    exprs_values = "relabundance"
)
all(rownames(reducedDim(genus_tree, "PCoA_BC")) == genus_tree$V1)
# Retrieving the explained variance
e <- attr(reducedDim(genus_tree, "PCoA_BC"), "eig")
var_explained <- e / sum(e[e > 0]) * 100
pcoa_data <- as.data.frame(reducedDim(genus_tree, "PCoA_BC"))
names(pcoa_data) <- c("PCoA 1", "PCoA 2")
p3 <- ggside(
    cbind(metadata, pcoa_data),
    aes(.data[["PCoA 1"]], .data[["PCoA 2"]])
) -
    # set default theme
    scheme_theme(
        plot.margin = margin(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 16, colour = "black"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16, face = "bold")
    ) +
    geom_point(aes(color = project), alpha = 0.5) +
    scale_color_brewer(palette = "Dark2", guide = "none") +
    ggforce::geom_mark_ellipse(aes(color = project), alpha = 0.5) +
    labs(
        x = paste0("PCoA 1 (", round(var_explained[1], 1), "%)"),
        y = paste0("PCoA 2 (", round(var_explained[2], 1), "%)"),
        color = ""
    ) +
    coord_cartesian(clip = "off") +

    # top annotation
    anno_top(size = 0.2, free_guides = NULL) -
    scheme_align(guides = "r", free_spaces = "r") +
    ## add density plot
    ggfree() +
    geom_density(aes(.data[["PCoA 1"]], fill = project), alpha = 0.5) +
    theme(legend.position = c(1, 0), legend.justification = c(0, 0.5)) +

    ## add box plot
    ggfree(waiver(), aes(.data[["PCoA 1"]], project)) +
    geom_boxplot(aes(fill = project), alpha = 0.5) -
    theme_no_axes() +
    guides(fill = "none") +

    # right annotation
    anno_right(size = 0.2) +
    ## add box plot
    ggfree(waiver(), aes(project, .data[["PCoA 2"]])) +
    geom_boxplot(aes(fill = project), alpha = 0.5) +
    ggfree() +
    geom_density(aes(y = .data[["PCoA 2"]], fill = project), alpha = 0.5) -
    theme_no_axes("xy") -
    guides(fill = "none") &
    scale_fill_brewer(palette = "Dark2", name = NULL)
ggsave("figures/microbiome/p3.pdf",
    plot = p3,
    width = 7, height = 7,
    family = "Helvetica"
)

# P4 --------------------------------------------------------
# compare relative abundance
groups <- split(seq_len(nrow(metadata)), metadata$project)
labels <- rownames(phylum)
# groups <- groups[c("COAD", "READ", "ESCA", "STAD", "HNSC")]
p4 <- stack_alignh(t(phylum)) -
    # set default theme
    scheme_theme(
        plot.margin = margin(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 16, colour = "black"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16, face = "bold")
    ) +
    align_dendro(aes(color = branch), size = 0.1, k = 3L) +
    scale_x_reverse() +
    scale_color_brewer(palette = "Dark2") +
    .mapply(function(index, project, next_index, next_project, size) {
        force(index)
        force(next_index)
        force(next_project)
        if (is.na(next_project)) {
            box <- NULL
        } else {
            box <- ggmark(
                mark_line(
                    # function(panel, link) {
                    #     browser()
                    # },
                    !!!lapply(seq_len(11), function(i) {
                        I(i) ~ waiver()
                    }),
                    .element = element_line(
                        colour = RColorBrewer::brewer.pal(3L, "Dark2")[
                            seq_len(2)
                        ]
                    )
                ),
                ~ .x[, c(index, next_index)],
                size = 0.15
            ) +
                scheme_data(function(data) {
                    print(mean(data$.column_names %in% labels[index]))
                    dplyr::mutate(
                        data,
                        group = dplyr::if_else(
                            .column_names %in% labels[index],
                            project, next_project
                        ),
                        group = factor(group, c(project, next_project))
                    )
                }) +
                geom_boxplot(aes(group, value, fill = group),
                    outlier.size = 0.2,
                    linewidth = 0.1
                ) +
                ggstattest::geom_comparetest(aes(group, value),
                    vjust = 1, size = 7, label_colour = "red"
                ) +
                scale_y_continuous(expand = expansion(mult = 0.3)) +
                theme_no_axes("xy") +
                theme(
                    strip.background = element_blank(),
                    strip.text = element_blank(),
                    plot.margin = margin(l = 0.1, r = 0.1, unit = "npc")
                ) +
                xlab(NULL) +
                scale_fill_brewer(palette = "Dark2", guide = "none") +
                facet_wrap(vars(), scales = "free")
        }
        list(
            ggheatmap(~ .x[, index], width = size) +
                theme_no_axes(if (project == "COAD") "x" else "xy") +
                scale_fill_viridis_c(
                    name = "relabund", option = "A", limits = c(0, 1)
                ) +
                theme(
                    plot.margin = margin(),
                    plot.title = element_text(face = "bold", hjust = 0.5)
                ) +
                ggtitle(project),
            stack_active(),
            box
        )
    }, list(
        groups, names(groups),
        groups[c(names(groups)[2:length(groups)], NA)],
        c(names(groups)[2:length(groups)], NA),
        lengths(groups) / sum(lengths(groups))
    ), NULL) +
    align_dendro(aes(color = branch), size = 0.1, k = 3L) +
    scale_color_brewer(palette = "Dark2")
ggsave("figures/microbiome/p4.pdf",
    plot = p4,
    width = 11, height = 5,
    family = "Helvetica"
)

# p5 --------------------------------------------------
p5 <- stack_discreteh(t(phylum), sizes = c(0.1, 1, 1)) -
    # set default theme
    scheme_theme(
        plot.margin = margin(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 16, colour = "black"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16, face = "bold")
    ) +
    # add heatmap
    ggheatmap(theme = theme(panel.spacing.x = unit(1, "mm"))) +
    scale_fill_viridis_c(
        name = "relabund",
        option = "plasma",
        limits = c(0, 1),
        breaks = scales::pretty_breaks(3L)
    ) +
    theme_no_axes("x") +
    theme(axis.text.y = element_text(colour = "black")) +

    anno_left(size = 0.2) +
    align_dendro() +
    scale_x_continuous(breaks = scales::pretty_breaks(3L)) +

    anno_top() +
    align_group(metadata$project) +

    # a new dendrogram
    align_dendro(aes(color = branch),
        reorder_dendrogram = TRUE,
        merge_dendrogram = TRUE
    ) +
    scale_color_brewer(palette = "Dark2", guide = "none") +
    scale_y_continuous(expand = expansion(), guide = "none") +
    no_expansion("y") +

    # add boxplot to annotate observations
    stack_active() -
    scheme_data(function(data) {
        dplyr::mutate(data,
            project = metadata$project[match(.column_names, metadata$V1)]
        )
    }) +
    ggmark(
        mark_tetragon(
            !!!lapply(seq_len(4), function(i) {
                I(i) ~ waiver()
            }),
            .element = element_vec_rep_each(
                element_polygon(
                    fill = RColorBrewer::brewer.pal(4, "Set3"),
                    color = NA
                ),
                2L
            )
        )
    ) +
    geom_boxplot(aes(value, project, fill = project), orientation = "y") +
    ggstattest::geom_comparetest(aes(value, project),
        vjust = 0.5, size = 7, label_colour = "black",
        step_increase = 0.25,
        height = 0.5,
        tip_length = 0.05,
        p.adjust = "BH"
    ) +
    scale_fill_brewer(palette = "Dark2", guide = "none") +
    scale_x_continuous(expand = expansion(add = c(0, 0.2), mult = c(0.05, 0))) +
    facet_grid(scales = "free_x") +
    xlab(NULL) +
    theme_no_axes("y") +
    theme(
        strip.text = element_blank(),
        strip.text.y.right = element_blank(),
        strip.background = element_blank(),
        # axis.text.x = element_text(angle = -60, hjust = 0),
        panel.spacing.y = unit(3, "mm")
    ) +
    theme(plot.margin = margin(l = 0.1, r = 0.1, unit = "npc")) +

    # add barplot
    ggalign(size = 0.5) +
    geom_bar(aes(value, fill = project),
        orientation = "y", stat = "identity",
        width = 1, position = "fill"
    ) +
    scale_fill_brewer(name = "Project", palette = "Dark2") +
    xlab(NULL) +
    no_expansion("x") +
    theme(
        panel.border = element_rect(fill = NA),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )

ggsave("figures/microbiome/p5.pdf",
    plot = p5, width = 11, height = 7,
    family = "Helvetica"
)

# COAD vs STAD ----------------------------------------
library(magrittr)
dim(phylum)
phylum_test <- apply(combn(unique(metadata$project), 2L), 2L, function(pair) {
    index <- metadata$project %in% pair
    data <- phylum[index, ]
    ans <- apply(data, 2L, function(relabundance) {
        ans <- wilcox.test(relabundance ~ group, data = data.frame(
            relabundance = relabundance,
            group = factor(metadata$project[index], pair)
        ))
        broom::tidy(ans)
    }, simplify = FALSE)
    ans <- dplyr::bind_rows(ans, .id = "phylum")
    ans$group1 <- pair[1L]
    ans$group2 <- pair[2L]
    ans$padj <- stats::p.adjust(ans$p.value, "BH")
    ans
}, simplify = FALSE)
phylum_sig <- dplyr::bind_rows(phylum_test) %>%
    dplyr::filter(padj < 0.05) %>%
    dplyr::arrange(group1, group2) %>%
    dplyr::mutate(n = dplyr::n(), .by = c(group1, group2))
coad_vs_stad_phylum <- phylum_sig %>%
    dplyr::filter(group1 == "COAD", group2 == "STAD")

genus_test <- apply(combn(unique(metadata$project), 2L), 2L, function(pair) {
    index <- metadata$project %in% pair
    data <- genus[index, ]
    ans <- apply(data, 2L, function(relabundance) {
        ans <- wilcox.test(relabundance ~ group, data = data.frame(
            relabundance = relabundance,
            group = factor(metadata$project[index], pair)
        ))
        broom::tidy(ans)
    }, simplify = FALSE)
    ans <- dplyr::bind_rows(ans, .id = "genus")
    ans$group1 <- pair[1L]
    ans$group2 <- pair[2L]
    ans$padj <- stats::p.adjust(ans$p.value, "BH")
    ans
}, simplify = FALSE)

genus_sig <- dplyr::bind_rows(genus_test) %>%
    dplyr::filter(padj < 0.05) %>%
    dplyr::arrange(group1, group2) %>%
    dplyr::mutate(n = dplyr::n(), .by = c(group1, group2))
coad_vs_stad_genus <- dplyr::filter(
    genus_sig, group1 == "COAD", group2 == "STAD"
)

taxonomy <- data.table::fread("rawdata/taxonomy.txt")
genus_taxonomy <- taxonomy[type == "genus"]
phylum_taxonomy <- taxonomy[type == "phylum"]
genus_to_phylum <- lapply(phylum_taxonomy$taxonomy, function(prefix) {
    genus_taxonomy$name[startsWith(genus_taxonomy$taxonomy, prefix)]
})
genus_to_phylum <- data.frame(
    phylum = rep(phylum_taxonomy$name, times = lengths(genus_to_phylum)),
    genus = unlist(genus_to_phylum, FALSE, FALSE)
)
genus_to_phylum <- genus_to_phylum %>%
    dplyr::filter(phylum %in% coad_vs_stad_phylum$phylum) %>%
    dplyr::filter(genus %in% coad_vs_stad_genus$genus) %>%
    dplyr::distinct()
anyDuplicated(genus_to_phylum$genus)

p6 <- stack_crossh(t(phylum), sizes = c(0.1, 1, 1)) -
    # set default theme
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

    # add COAD heatmap
    ggheatmap(~ .x[, metadata$project == "COAD"]) +
    scale_fill_viridis_c(
        name = "relabund",
        option = "plasma", 
        limits = c(0, 1),
        breaks = scales::pretty_breaks(3L)
    ) +
    theme_no_axes("x") +

    ## for the heatmap, add top annotation
    anno_top() -
    scheme_align(free_spaces = "lr") +
    align_dendro(aes(color = branch), k = 3L) +
    scale_color_brewer(palette = "Set2", guide = "none") +
    theme_no_axes("y") +
    ggtitle("phylum") +
    no_expansion() +

    # add a box to compare the relabundance
    stack_active() +
    ggmark(mark_line(
        !!!lapply(seq_len(11), function(i) I(i) ~ waiver()),
        .element = element_line(
            colour = RColorBrewer::brewer.pal(3L, "Set1")[seq_len(2)]
        )
    ), data = NULL, size = 0.4) +
    scheme_data(function(data) {
        ans <- phylum[
            rownames(phylum) %in% metadata$V1[
                metadata$project %in% c("COAD", "STAD")
            ],
        ]
        ans <- fortify_data_frame(ans)
        ans$.hand <- dplyr::if_else(
            ans$.row_names %in% metadata$V1[metadata$project == "COAD"],
            "COAD", "STAD"
        )
        data$.hand <- dplyr::if_else(data$.hand == "left", "COAD", "STAD")
        ans <- dplyr::full_join(ans, data,
            by = c(.column_names = ".names", ".hand")
        )
        ans$.hand <- factor(ans$.hand, c("COAD", "STAD"))
        ans
    }) +
    geom_boxplot(aes(.hand, value, fill = .hand),
        outlier.size = 0.2, linewidth = 0.1
    ) +
    ggstattest::geom_comparetest(aes(.hand, value),
        vjust = 0.5, size = 7, label_colour = "red"
    ) +
    scale_y_continuous(expand = expansion(add = c(0, 0.5), mult = c(0.1, 0))) +
    theme_no_axes("xy") +
    theme(
        strip.background = element_blank(),
        strip.text = element_blank(),
        strip.text.y = element_blank(),
        plot.margin = margin(t = 0.01, l = 0.1, r = 0.1, unit = "npc")
    ) +
    xlab(NULL) +
    scale_fill_brewer(palette = "Set1", guide = "none") +
    facet_grid(scales = "free_y") +

    # add STAD heatmap
    ggheatmap(~ .x[, metadata$project == "STAD"]) +
    scale_fill_viridis_c(
        name = "relabund",
        option = "plasma", limits = c(0, 1),
        breaks = scales::pretty_breaks(3L)
    ) +
    theme_no_axes() +
    theme(plot.margin = margin()) +
    anno_top() -
    scheme_align(free_spaces = "lr") +
    align_dendro(aes(color = branch), k = 3L) +
    scale_color_brewer(palette = "Set2", guide = "none") +
    no_expansion("y") +
    theme_no_axes("y") +
    stack_active() +

    # add differential expressed genus
    cross_link(
        link_tetragon(
            !!!purrr::imap(
                split(genus_to_phylum$genus, genus_to_phylum$phylum),
                function(genus, phylum) phylum ~ genus
            ),
            .element = element_polygon(
                fill = RColorBrewer::brewer.pal(5, "Dark2"),
                color = NA
            )
        ),
        t(genus)[
            genus_to_phylum$genus, metadata$project %in% c("COAD", "STAD")
        ],
        on_top = TRUE,
        size = 0.2
    ) +
    theme(plot.margin = margin()) +

    #
    align_group(
        factor(
            genus_to_phylum$phylum,
            levels = na.omit(genus_to_phylum$phylum[
                match(colnames(phylum), genus_to_phylum$phylum)
            ])
        )
    ) +

    ggheatmap() +
    scale_fill_viridis_c(
        name = "relabund",
        option = "plasma", limits = c(0, 1),
        breaks = scales::pretty_breaks(3L)
    ) +
    scale_y_continuous(position = "right") +
    theme_no_axes("x") +
    theme(panel.border = element_rect(fill = NA),
        axis.text.y = element_text(size = 14)
    ) +

    anno_top() -
    scheme_align(free_spaces = "l") +
    align_group(metadata$project[metadata$project %in% c("COAD", "STAD")]) +
    align_dendro(aes(color = branch), key_glyph = draw_key_rect) +
    ggtitle("genus") +
    scale_color_brewer(palette = "Set1", name = "Project") +
    # scale_y_continuous(breaks = scales::breaks_pretty(3)) +
    no_expansion("y") +
    theme_no_axes("y")

ggsave("figures/microbiome/p6.pdf",
    plot = p6, width = 11, height = 7,
    family = "Helvetica"
)
