library(SingleCellExperiment)
if (!dir.exists("munging")) dir.create("munging")

# read data -------------------------------------------------
# Please download data from:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210347
data <- readRDS("rawdata/GSE210347_counts.Rds")
metadata <- data.table::fread("rawdata/GSE210347_meta.txt.gz")
sample_data <- data.table::fread("rawdata/GSE210347_sample_information.txt.gz")
study <- readxl::read_xls("rawdata/GSE210347_study_metadata.xls")

all(metadata$cellname == colnames(data))
sce_obj <- SingleCellExperiment::SingleCellExperiment(
    list(counts = data),
    colData = metadata
)
dim(sce_obj)

# filter data
set.seed(1234)
sce_obj <- scend::logNormCounts(sce_obj)
sce_obj <- sce_obj[
    ,
    sce_obj$celltype != "undefined" & sce_obj$cluster != "undefined"
]

# MNN
dec_combined <- scend::modelGeneVar(sce_obj, block = sce_obj$SampleID)
SingleCellExperiment::rowSubset(sce_obj) <- scend::getTopHVGs(
    dec_combined,
    n = 2000L
)
set.seed(1234L)
sce_obj <- scend::quickMNN(
    sce_obj,
    block = sce_obj$SampleID,
    d = 50L,
    subset_row = SingleCellExperiment::rowSubset(sce_obj)
)
SummarizedExperiment::assay(sce_obj, "multiBatchNorm") <- NULL

# UMAP
set.seed(1234L)
sce_obj <- scend::runUMAP(sce_obj, dimred = "corrected", min_dist = 0.3)
qs::qsave(sce_obj, "munging/sce_obj.qs")

marker_list <- list(
    Epithelium = c("EPCAM", "KRT19"),
    Lymphocyte = c("CD3D", "CD3E"),
    Myeloid = c("CD14", "CD68"),
    Fibroblast = c("DCN", "COL1A2", "COL1A1"),
    Endothelium = c("PECAM1", "VWF"),
    Plasma = c("IGHG1", "JCHAIN")
)
umap_subset <- lapply(sort(names(marker_list)), function(x) {
    index <- which(sce_obj$celltype == x)
    sce_obj <- sce_obj[, index]
    set.seed(1234)
    umap <- reducedDim(
        scend::runUMAP(sce_obj, dimred = "corrected", min_dist = 0.3),
        "UMAP"
    )
    umap <- setNames(data.frame(umap), c("umap1", "umap2"))
    umap$.panel <- x
    umap$.index <- index
    umap$cluster <- sce_obj$cluster[index]
    umap
})
names(umap_subset) <- sort(names(marker_list))
qs::qsave(umap_subset, "munging/umap_subset.qs")
