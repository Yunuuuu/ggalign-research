library(magrittr)
data <- readxl::read_xlsx("rawdata/genomic/mut_load.xlsx")

# copied from https://github.com/UMCUGenetics/primary-met-wgs-comparison/blob/master/overview_figure/source_script.R
cancer_type_order_fig1 <- c(
    "Pan-cancer", "Glioblastoma multiforme",
    "Upper respiratory tract carcinoma", "Thyroid carcinoma",
    "Lung adenocarcinoma", "Lung squamous cell carcinoma",
    "Diffuse large B-cell lymphoma", "Breast carcinoma",
    "Cholangiocarcinoma", "Hepatocellular carcinoma", "Pancreas carcinoma",
    "Pancreas neuroendocrine", "Colorectal carcinoma",
    "Esophageal carcinoma", "Stomach carcinoma",
    "Kidney renal clear cell carcinoma", "Cervical carcinoma", "Ovarian serous adenocarcinoma", "Uterus carcinoma",
    "Bladder urothelial carcinoma", "Prostate carcinoma", "Leiomyosarcoma", "Liposarcoma", "Skin melanoma"
)

# summarize the biopsy data by binning it into 4 groups: local, lymph, distant,
# and unknown based on the specific sites annotated in the clinical data
data$simplified_biopsiy_site <- NA
data$biopsySite <- data$biopsy_site

for (i in 2:24) {
    cc_type <- cancer_type_order_fig1[i]

    if (i == 2) {
        data$simplified_biopsiy_site[
            data$cancer_type == cc_type & data$biopsySite != "ovarium"
        ] <- "Local"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "ovarium"] <- "Distant"
    }

    if (i == 3) {
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "Primary"] <- "Local"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "null"] <- "Unknown"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "Lymph node"] <- "Lymph"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & is.na(data$simplified_biopsiy_site)] <- "Distant"
    }


    if (i == 4) {
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "null"] <- "Unknown"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "Lymph node"] <- "Lymph"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & is.na(data$simplified_biopsiy_site)] <- "Distant"
    }

    if (i == 5) {
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite %in% c("Primary", "Lung")] <- "Local"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "null"] <- "Unknown"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "Lymph node"] <- "Lymph"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & is.na(data$simplified_biopsiy_site)] <- "Distant"
    }

    if (i == 6) {
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite %in% c("Primary", "Lung")] <- "Local"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "null"] <- "Unknown"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "Lymph node"] <- "Lymph"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & is.na(data$simplified_biopsiy_site)] <- "Distant"
    }

    if (i == 7) {
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "null"] <- "Unknown"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "Lymph node"] <- "Lymph"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & is.na(data$simplified_biopsiy_site)] <- "Distant"
    }

    if (i == 8) {
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite %in% c("Primary", "Breast")] <- "Local"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "null"] <- "Unknown"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "Lymph node"] <- "Lymph"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & is.na(data$simplified_biopsiy_site)] <- "Distant"
    }

    if (i == 9) {
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "Primary"] <- "Local"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite %in% c("null", "other")] <- "Unknown"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "Lymph node"] <- "Lymph"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & is.na(data$simplified_biopsiy_site)] <- "Distant"
    }

    if (i == 10) {
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite %in% c("Primary", "Liver")] <- "Local"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & is.na(data$simplified_biopsiy_site)] <- "Distant"
    }

    if (i == 11) {
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite %in% c("Primary", "pancreas")] <- "Local"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "null"] <- "Unknown"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "Lymph node"] <- "Lymph"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & is.na(data$simplified_biopsiy_site)] <- "Distant"
    }

    if (i == 12) {
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "null"] <- "Unknown"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & is.na(data$simplified_biopsiy_site)] <- "Distant"
    }

    if (i == 13) {
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "Colorectum"] <- "Local"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite %in% c("Unknown", "Other site")] <- "Unknown"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "Lymph node"] <- "Lymph"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & is.na(data$simplified_biopsiy_site)] <- "Distant"
    }

    if (i == 14) {
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite %in% c("Primary", "oesophagus")] <- "Local"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "null"] <- "Unknown"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "Lymph node"] <- "Lymph"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & is.na(data$simplified_biopsiy_site)] <- "Distant"
    }

    if (i == 15) {
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite %in% c("Primary", "stomach", "Stomach")] <- "Local"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "null"] <- "Unknown"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "Lymph node"] <- "Lymph"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & is.na(data$simplified_biopsiy_site)] <- "Distant"
    }

    if (i == 16) {
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite %in% c("Primary", "Kidney")] <- "Local"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "null"] <- "Unknown"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "Lymph node"] <- "Lymph"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & is.na(data$simplified_biopsiy_site)] <- "Distant"
    }

    if (i == 17) {
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite %in% c("Primary", "cervix")] <- "Local"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "null"] <- "Unknown"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "Lymph node"] <- "Lymph"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & is.na(data$simplified_biopsiy_site)] <- "Distant"
    }

    if (i == 18) {
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite %in% c("Primary", "Ovarium CA", "ovary")] <- "Local"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "null"] <- "Unknown"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "Lymph node"] <- "Lymph"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & is.na(data$simplified_biopsiy_site)] <- "Distant"
    }

    if (i == 19) {
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite %in% c("Primary")] <- "Local"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "null"] <- "Unknown"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "Lymph node"] <- "Lymph"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & is.na(data$simplified_biopsiy_site)] <- "Distant"
    }

    if (i == 20) {
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite %in% c("Primary")] <- "Local"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite %in% c("null", "unknown")] <- "Unknown"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "Lymph node"] <- "Lymph"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & is.na(data$simplified_biopsiy_site)] <- "Distant"
    }

    if (i == 21) {
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite %in% c("Primary", "prostate")] <- "Local"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "null"] <- "Unknown"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "Lymph node"] <- "Lymph"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & is.na(data$simplified_biopsiy_site)] <- "Distant"
    }

    if (i == 22) {
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "null"] <- "Unknown"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "Lymph node"] <- "Lymph"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & is.na(data$simplified_biopsiy_site)] <- "Distant"
    }

    if (i == 23) {
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "null"] <- "Unknown"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & is.na(data$simplified_biopsiy_site)] <- "Distant"
    }

    if (i == 24) {
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "Skin/Subcutaneous"] <- "Local"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "Unknown"] <- "Unknown"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & data$biopsySite == "Lymph node"] <- "Lymph"
        data$simplified_biopsiy_site[data$cancer_type == cc_type & is.na(data$simplified_biopsiy_site)] <- "Distant"
    }
}

burden <- data %>%
    dplyr::filter(!is.na(simplified_biopsiy_site)) %>%
    dplyr::filter(!is_blacklisted) %>%
    dplyr::mutate(
        group = dplyr::if_else(
            !is.na(metastatic_location) & metastatic_location == "Primary",
            "Primary", "Metastatic"
        ),
        group = factor(group, c("Primary", "Metastatic"))
    )

names(data)

burden_compare <- readxl::read_xlsx(
    "rawdata/genomic/sig_contribs.xlsx",
    sheet = "smnv_burden_enrichment",
    skip = 1L
)

etio_enrich <- readxl::read_xlsx(
    "rawdata/genomic/sig_contribs.xlsx",
    sheet = "etio_enr.cancer_stage", skip = 1L
)

etio_enrich_metastatic <- etio_enrich %>%
    dplyr::filter(
        mut_type %in% "SBS",
        sig_group %in% c(
            "SBS5/40 | Age",
            "SBS1 | Age, 5mC deamination",
            "SBS2/13 | APOBEC",
            "SBS12/16 | Liver associated",
            "SBS6/15/21/26/44 | MMRd",
            "SBS31/35 | Platinum treatment",
            "SBS17a/17b | ROS, 5FU treatment",
            "SBS18/36 | ROS, ROS repair deficiency",
            "SBS33 | Unknow"
        )
    ) %>%
    tidyr::pivot_wider(
        id_cols = sig_group,
        names_from = cancer_type,
        values_from = median_metastatic
    ) %>%
    tibble::column_to_rownames(var = "sig_group") %>%
    as.matrix()
etio_enrich_primary <- etio_enrich %>%
    dplyr::filter(
        mut_type %in% "SBS",
        sig_group %in% c(
            "SBS5/40 | Age",
            "SBS1 | Age, 5mC deamination",
            "SBS2/13 | APOBEC",
            "SBS12/16 | Liver associated",
            "SBS6/15/21/26/44 | MMRd",
            "SBS31/35 | Platinum treatment",
            "SBS17a/17b | ROS, 5FU treatment",
            "SBS18/36 | ROS, ROS repair deficiency",
            "SBS33 | Unknow"
        )
    ) %>%
    tidyr::pivot_wider(
        id_cols = sig_group,
        names_from = cancer_type,
        values_from = median_primary
    ) %>%
    tibble::column_to_rownames(var = "sig_group") %>%
    as.matrix()
etio_enrich_difference <- etio_enrich %>%
    dplyr::filter(
        mut_type %in% "SBS",
        sig_group %in% c(
            "SBS5/40 | Age",
            "SBS1 | Age, 5mC deamination",
            "SBS2/13 | APOBEC",
            "SBS12/16 | Liver associated",
            "SBS6/15/21/26/44 | MMRd",
            "SBS31/35 | Platinum treatment",
            "SBS17a/17b | ROS, 5FU treatment",
            "SBS18/36 | ROS, ROS repair deficiency",
            "SBS33 | Unknow"
        )
    ) %>%
    tidyr::pivot_wider(
        id_cols = sig_group,
        names_from = cancer_type,
        values_from = median_diff
    ) %>%
    tibble::column_to_rownames(var = "sig_group") %>%
    as.matrix()
all(rownames(etio_enrich_primary) == rownames(etio_enrich_metastatic))
all(colnames(etio_enrich_primary) == colnames(etio_enrich_metastatic))

# visualization with ggalign
library(ggalign)
library(forcats)
names(data)

p <- stack_alignv(data.frame(row.names = colnames(etio_enrich_primary))) +
    align_group(factor(
        colnames(etio_enrich_primary),
        colnames(etio_enrich_primary)
    )) +

    # point plot for the burden
    ggfree() +
    scheme_data(function(data) {
        ans <- burden[c(
            "cancer_type", "sample_id", "group",
            "sbs_load", "dbs_load", "indel_load"
        )]
        ans <- dplyr::mutate(
            ans,
            dplyr::across(ends_with("_load"), as.numeric)
        )
        ans$.panel <- data$.panel[
            match(ans$cancer_type, sub("\\s*\\(.+\\)$", "", data$.names))
        ]
        tidyr::pivot_longer(
            ans,
            all_of(c("sbs_load", "dbs_load", "indel_load")),
            names_to = "Mut", values_to = "burden"
        ) %>%
            dplyr::mutate(
                x = scales::rescale(
                    as.integer(fct_reorder(sample_id, burden)),
                    c(0, 1)
                ),
                .by = c(.panel, group, Mut)
            ) %>%
            dplyr::mutate(
                drawing_group = paste(group, Mut, sep = "-"),
                drawing_group = factor(
                    drawing_group,
                    unique(drawing_group[order(
                        match(Mut, c("sbs_load", "dbs_load", "indel_load")),
                        as.integer(group)
                    )])
                )
            )
    }) +
    geom_point(
        aes(.data$x, .data$burden, color = drawing_group),
        size = 1, shape = 16
    ) +
    geom_line(
        aes(x = x, y = y, color = drawing_group),
        data = function(dd) {
            which.median <- function(x) {
                ordering <- order(x)
                if ((len <- length(x)) == 0L) {
                    integer()
                } else if (len %% 2L == 0) {
                    ordering[len / 2L + 0:1]
                } else {
                    ordering[(len + 1L) / 2L]
                }
            }
            dlist <- vctrs::vec_split(dd, dd[c("drawing_group", ".panel")])
            coord <- .mapply(function(key, val) {
                index <- which.median(.subset2(val, "burden"))
                d <- vctrs::vec_slice(val, index)
                x <- mean(d$x)
                y <- mean(d$burden)
                data.frame(x = x + c(-0.4, 0.4), y = y)
            }, dlist, NULL)
            vctrs::vec_cbind(
                vctrs::vec_rep_each(
                    .subset2(dlist, "key"),
                    vctrs::list_sizes(coord)
                ),
                vctrs::vec_rbind(!!!coord)
            )
        }
    ) +
    scale_y_continuous(name = "Number of mutations", transform = "log10") +
    scale_color_brewer(
        name = NULL,
        palette = "Dark2",
        guide = guide_legend(override.aes = list(size = 5)),
        labels = function(x) {
            vapply(strsplit(sub("_load", "", tolower(x)), "-"),
                function(label) {
                    paste(
                        switch(label[2L],
                            sbs = "SBS",
                            dbs = "DBS",
                            indel = "Indel"
                        ),
                        label[1L]
                    )
                }, character(1L),
                USE.NAMES = FALSE
            )
        }
    ) +
    theme_no_axes("x") +
    theme(
        strip.clip = "off",
        strip.text.x = element_text(angle = -90, hjust = 0),
        strip.background = element_blank(),
        legend.position = "bottom",
        # panel.border = element_rect(fill = NA, colour = "black")
        axis.line.x = element_line(),
        strip.placement = "outside"
    ) +
    no_expansion("x") +
    facet_grid(cols = vars(.panel), switch = "x")

ggsave("figures/genomic/p4.pdf", plot = p, width = 10, height = 7)

# moon plot
# ggheatmap(etio_enrich_metastatic, filling = NULL) +
# scheme_data(function(d) {
#     d$primary <- etio_enrich_primary[
#         cbind(d$.row_index, d$.column_index)
#     ]
#     d$metastatic <- d$value
#     d$difference <- etio_enrich_difference[
#         cbind(d$.row_index, d$.column_index)
#     ]
#     d <- tidyr::pivot_longer(d, c(primary, metastatic),
#         names_to = "group",
#         values_to = "SBS"
#     )
#     d$point_size <- dplyr::case_when(
#         d$difference == 0L ~ 0.1,
#         d$difference < 100L ~ 0.1,
#         d$difference == 0L ~ 0.1,
#     )
#     d$group <- factor(d$group, c("primary", "metastatic"))
#     d
# }) +
# ggstattest::geom_draw(
#     function(data, panel_params, coord) {
#         coords <- lapply(seq_len(nrow(data)), function(i) {
#             d <- vctrs::vec_slice(data, i)
#             if (d$group == 1L) {
#                 # right circle
#                 rs <- seq(-pi / 2, pi / 2, length.out = 100L)
#                 xc <- d$x + d$size * cos(rs) / 2L
#                 yc <- d$y + d$size * sin(rs) / 2L
#             } else {
#                 # left circle
#                 rs <- seq(pi / 2, pi * 3 / 2, length.out = 100L)
#                 xc <- d$x + d$size * cos(rs) / 2L
#                 yc <- d$y + d$size * sin(rs) / 2L
#             }
#             data.frame(x = xc, y = yc)
#         })
#         browser()
#         id_lengths <- vctrs::list_sizes(coords)
#         coords <- coord$transform(vctrs::vec_rbind(!!!coords), panel_params)
#         grid::polygonGrob(
#             x = coords$x,
#             y = coords$y,
#             id.lengths = id_lengths,
#             default.units = "native",
#             gp = grid::gpar(fill = data$fill, col = NA)
#         )
#     }, aes(fill = SBS, group = group, size = difference)
# ) +
# theme_no_axes("x") +
# scale_size_binned(
#     breaks = c(0, 100, 500, 1000),
#     range = c(0, 1),
#     label = c("<0.01", "0.02", "0.04", ">0.04")
# ) +
# scale_fill_gradientn(
#     colours = c(
#         "#3288BD", "#66C2A5", "#ABDDA4", "#E6F598",
#         "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43",
#         "#D53E4F", "#9E0142"
#     )
# )
