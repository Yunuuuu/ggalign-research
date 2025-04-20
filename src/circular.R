# Please use ggplot2 > 3.5.1 (or development version 3.5.1.9000)
library(ggalign)
library(magrittr)

if (!dir.exists("figures/circular")) dir.create("figures/circular")

# p1 ----------------------------------------------
# https://pmc.ncbi.nlm.nih.gov/articles/PMC9856581/#app1-cancers-15-00342
data <- readxl::read_xlsx(
    "rawdata/circular/Supplementary tables.xlsx",
    skip = 1L
)
logFC <- data[2:18, c(1, which(as.character(data[1, ]) == "logFC"))]
adj.P <- data[2:18, c(1, which(as.character(data[1, ]) == "adj.P"))]
names(adj.P) <- names(logFC)

logFC <- dplyr::mutate(logFC, dplyr::across(!pathway, as.numeric)) %>%
    tibble::column_to_rownames(var = "pathway") %>%
    as.matrix() %>%
    t()
adj.P <- dplyr::mutate(adj.P, dplyr::across(!pathway, as.numeric)) %>%
    tibble::column_to_rownames(var = "pathway") %>%
    as.matrix() %>%
    t()

p1 <- circle_discrete(logFC,
    radial = coord_radial(
        start = pi / 2, end = pi * 2, expand = FALSE
    ),
    theme = theme(plot.margin = margin(l = 5, unit = "mm"))
) +
    align_dendro(aes(color = branch), k = 4L, size = 0.5) +
    theme_no_axes("y") +
    scale_color_brewer(palette = "Dark2") +
    ggalign(mapping = aes(y = .column_names, fill = value)) +
    geom_tile(width = 1, height = 1) +
    scale_fill_gradient2(low = "blue", high = "red", name = "logFC") +
    geom_text(aes(label = "*"), data = function(dd) {
        dd$pvalue <- adj.P[cbind(dd$.row_index, dd$.column_index)]
        dplyr::filter(dd, pvalue < 0.05)
    }) +
    no_expansion("y") +
    theme(
        panel.spacing.r = unit(0, "mm"),
        axis.text.r = element_text(),
        axis.text.theta = element_text()
    ) +
    guides(
        r = "none", r.sec = "axis",
        theta = guide_axis_theta(angle = 0)
    )
ggsave("figures/circular/p1.pdf", plot = p1, width = 8, height = 8)

# p2 -------------------------------------------
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

p2 <- circle_discrete(trda,
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
        plot.margin = margin(l = -10, t = -20, b = -20, unit = "mm")
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
    ggtree:::geom_curvelink(
        mapping = aes(x, 1, xend = xend, yend = 1, color = Interaction),
        outward = FALSE, key_glyph = draw_key_path,
        linewidth = 0.1
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = NULL) +
    scale_colour_manual(
        values = c("chocolate2", "#3690C0", "#009E73"),
        guide = guide_legend(
            keywidth = 0.8, keyheight = 0.5,
            order = 1, override.aes = list(
                alpha = 1, linewidth = 1
            )
        )
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
        guide = guide_legend(keywidth = 0.5, keyheight = 0.5, order = 2)
    ) +
    scale_color_manual(
        values = c(
            "#EF3B2C", "#1D91C0", "#FEB24C", "grey60",
            "#7FBC41", "#4D9221", "#276419"
        ),
        guide = guide_legend(keywidth = 0.5, keyheight = 0.5, order = 2)
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

    # add BGC data --------------------------------------
    ggalign() +
    scheme_data(function(d) {
        BGCsda$x <- d$.x[match(BGCsda$Strain, d$.names)]
        dplyr::filter(BGCsda, !is.na(x), Count > 0L)
    }) +
    geom_point(aes(x, BGCs, color = BGCs, size = log2(Count + 1)), shape = 15) +
    scale_color_manual(
        values = c(
            "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
            "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"
        ),
        guide = guide_legend(keywidth = 0.4, keyheight = 0.4, order = 3)
    ) +
    theme_no_axes() +
    theme(panel.spacing.r = unit(0, "mm")) +
    scale_size_continuous(
        range = c(1, 2),
        breaks = c(1, 2, 3),
        guide = guide_legend(
            keywidth = 0.4, keyheight = 0.4, order = 4,
            override.aes = list(shape = 0)
        )
    ) +

    # add labels -----------------------------------------
    ggalign(rownames, size = 0.4) + # we use rownames to match the data
    scheme_data(function(data) {
        data$.names <- gsub("s__Leaf", "", data$.names)
        data
    }) +
    geom_text(aes(y = 0, label = .names), angle = 90, size = 2) +
    no_expansion("y") +
    theme_no_axes() +
    theme(panel.spacing.r = unit(0, "mm")) +

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
        name = "Number of interactions",
        values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
        guide = guide_legend(keywidth = 0.5, keyheight = 0.5, order = 5)
    ) +
    no_expansion("y") &
    theme(panel.background = element_blank())

ggsave("figures/circular/p2.pdf", plot = p2, width = 10, height = 8)

# p4 ---------------------------------------------------------
upset_data <- dplyr::filter(BGCsda, Count > 0)
upset_data <- split(upset_data$Strain, upset_data$BGCs)
p4 <- ggupset(tune(upset_data),
    point = list(
        data = function(d) {
            dplyr::mutate(d, point_group = dplyr::case_when(
                !point_group ~ NA_character_,
                .default = .discrete_y
            ))
        }
    )
) -
    scheme_theme(plot.margin = margin()) +
    scale_fill_manual(values = c("#F0F0F0", "white"), guide = "none") +
    scale_color_brewer(
        palette = "Set1", guide = "none",
        na.translate = TRUE, na.value = "grey"
    ) +
    # scale_color_manual(values = c("grey", "black"), guide = "none") +
    anno_top(size = 0.5) +
    ggalign(data = function(d) ggalign_attr(d, "intersection_sizes")) +
    geom_bar(aes(y = .data$value), stat = "identity") +
    ylab("Intersection Sizes") +
    anno_right(size = 0.5) +
    ggalign(data = function(d) ggalign_attr(d, "set_sizes")) +
    geom_bar(aes(x = .data$value),
        stat = "identity",
        orientation = "y"
    ) +
    xlab("Set Sizes") &
    theme(plot.background = element_blank())
ggsave("figures/circular/p4.pdf", plot = p4, width = 8, height = 4)
