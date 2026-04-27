
# Project:
# Single-cell transcriptomic analysis of differentially
# abundant immune cell populations between PTCB and TCB
# in GSE271413


# 1. Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)
library(cowplot)
library(pheatmap)
library(Azimuth)
library(clusterProfiler)
library(org.Hs.eg.db)

theme_set(theme_cowplot())

# 2. Set paths

base_dir <- "D:/Research/GSE271413_RAW"
out_dir  <- "D:/Research/GSE271413_analysis"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# 3. Define sample information

sample_info <- data.frame(
  sample_id = c(
    "PTCB_001", "PTCB_002", "PTCB_003", "PTCB_004", "PTCB_005",
    "TCB_006", "TCB_007", "TCB_008", "TCB_009", "TCB_010"
  ),
  folder = c(
    "GSM8376465_1394_PTCB",
    "GSM8376495_1446_PTCB",
    "GSM8376498_1445_PTCB",
    "GSM8376501_1429_PTCB",
    "GSM8376504_1400_PTCB",
    "GSM8376458_1379_TCB",
    "GSM8376459_1390_TCB",
    "GSM8376477_1330_TCB",
    "GSM8376480_1324_TCB",
    "GSM8376483_1321_TCB"
  ),
  group = c(
    rep("PTCB", 5),
    rep("TCB", 5)
  ),
  stringsAsFactors = FALSE
)

sample_info$patient <- sample_info$sample_id

# 4. Read 10X data and create Seurat objects

seurat_list <- list()

for (i in seq_len(nrow(sample_info))) {
  sid   <- sample_info$sample_id[i]
  sdir  <- file.path(base_dir, sample_info$folder[i])
  grp   <- sample_info$group[i]
  pat   <- sample_info$patient[i]
  
  message("Reading: ", sid)
  
  mat <- Read10X(data.dir = sdir)
  
  obj <- CreateSeuratObject(
    counts = mat,
    project = sid,
    min.cells = 3,
    min.features = 200
  )
  
  obj$sample_id <- sid
  obj$patient   <- pat
  obj$group     <- grp
  
  seurat_list[[sid]] <- obj
}


# 5. Merge all samples
merged_obj <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "GSE271413_PTCB_vs_TCB"
)

saveRDS(merged_obj, file = file.path(out_dir, "01_merged_raw.rds"))


# 6. QC metrics

DefaultAssay(merged_obj) <- "RNA"

merged_obj[["percent.mt"]]   <- PercentageFeatureSet(merged_obj, pattern = "^MT-")
merged_obj[["percent.ribo"]] <- PercentageFeatureSet(merged_obj, pattern = "^RP[SL]")

qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")

pdf(file.path(out_dir, "02_QC_before_filtering_by_patient.pdf"), width = 14, height = 8)
print(
  VlnPlot(
    merged_obj,
    features = qc_features,
    group.by = "patient",
    ncol = 4,
    pt.size = 0,
    raster = FALSE
  ) +
    plot_annotation(title = "QC Metrics Before Filtering")
)
dev.off()

# 7. QC filtering

filtered_obj <- subset(
  merged_obj,
  subset =
    nCount_RNA > 200 &
    nCount_RNA < 20000 &
    nFeature_RNA > 200 &
    nFeature_RNA < 5000 &
    percent.mt < 10
)

filtered_obj[["percent.mt"]]   <- PercentageFeatureSet(filtered_obj, pattern = "^MT-")
filtered_obj[["percent.ribo"]] <- PercentageFeatureSet(filtered_obj, pattern = "^RP[SL]")

pdf(file.path(out_dir, "03_QC_after_filtering_by_patient.pdf"), width = 14, height = 8)
print(
  VlnPlot(
    filtered_obj,
    features = qc_features,
    group.by = "patient",
    ncol = 4,
    pt.size = 0,
    raster = FALSE
  ) +
    plot_annotation(title = "QC Metrics After Filtering")
)
dev.off()

saveRDS(filtered_obj, file = file.path(out_dir, "04_filtered_obj.rds"))


# 8. Standard preprocessing

filtered_obj <- NormalizeData(filtered_obj)
filtered_obj <- FindVariableFeatures(filtered_obj, selection.method = "vst", nfeatures = 3000)
filtered_obj <- ScaleData(filtered_obj)
filtered_obj <- RunPCA(filtered_obj, verbose = FALSE)

# Unintegrated clustering
filtered_obj <- FindNeighbors(filtered_obj, dims = 1:30)
filtered_obj <- FindClusters(filtered_obj, resolution = 0.8)
filtered_obj <- RunUMAP(filtered_obj, dims = 1:30, reduction.name = "umap")

pdf(file.path(out_dir, "05_preintegration_umap.pdf"), width = 14, height = 6)
p1 <- DimPlot(filtered_obj, reduction = "umap", group.by = "patient", raster = FALSE) +
  ggtitle("Pre-integration UMAP by Patient")
p2 <- DimPlot(filtered_obj, reduction = "umap", group.by = "group", raster = FALSE) +
  ggtitle("Pre-integration UMAP by Group")
print(p1 | p2)
dev.off()

# 9. Integration (Seurat v5)

integrated_obj <- IntegrateLayers(
  object = filtered_obj,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)

integrated_obj <- JoinLayers(integrated_obj)

integrated_obj <- FindNeighbors(integrated_obj, reduction = "integrated.cca", dims = 1:20)
integrated_obj <- FindClusters(integrated_obj, resolution = 0.6, cluster.name = "cca_clusters")
integrated_obj <- RunUMAP(
  integrated_obj,
  reduction = "integrated.cca",
  dims = 1:20,
  reduction.name = "umap.cca"
)

saveRDS(integrated_obj, file = file.path(out_dir, "06_integrated_obj.rds"))

pdf(file.path(out_dir, "06_integrated_umap.pdf"), width = 14, height = 10)
p3 <- DimPlot(integrated_obj, reduction = "umap.cca", group.by = "group", raster = FALSE) +
  ggtitle("Integrated UMAP by Group")
p4 <- DimPlot(integrated_obj, reduction = "umap.cca", group.by = "patient", raster = FALSE) +
  ggtitle("Integrated UMAP by Patient")
p5 <- DimPlot(integrated_obj, reduction = "umap.cca", group.by = "cca_clusters",
              label = TRUE, repel = TRUE, raster = FALSE) +
  ggtitle("Integrated UMAP by Cluster")
print((p3 | p4) / p5)
dev.off()


# 10. Azimuth annotation

# Replace this with your real Azimuth reference file
azi_ref <- readRDS("D:/Research/ref.Rds")

# Ensure reference loading names are valid if needed
if ("refDR" %in% names(azi_ref@reductions)) {
  colnames(azi_ref[["refDR"]]@feature.loadings) <- paste0(
    "refdr_",
    seq_len(ncol(azi_ref[["refDR"]]@feature.loadings))
  )
}

query_obj <- integrated_obj

# Match common genes
azi_features <- intersect(rownames(azi_ref), rownames(query_obj))

anchors <- FindTransferAnchors(
  reference = azi_ref,
  query = query_obj,
  normalization.method = "SCT",
  reference.reduction = "refDR",
  dims = 1:50,
  features = azi_features
)

query_obj <- MapQuery(
  anchorset = anchors,
  query = query_obj,
  reference = azi_ref,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    celltype.l3 = "celltype.l3"
  ),
  reference.reduction = "refDR",
  reduction.model = "refUMAP"
)

# Save level-2 predicted label as final working annotation
query_obj$cell_type <- query_obj$predicted.celltype.l2

saveRDS(query_obj, file = file.path(out_dir, "07_annotated_obj.rds"))

pdf(file.path(out_dir, "07_celltype_umap.pdf"), width = 12, height = 8)
print(
  DimPlot(
    query_obj,
    reduction = "umap.cca",
    group.by = "cell_type",
    label = TRUE,
    repel = TRUE,
    raster = FALSE
  ) + ggtitle("Integrated UMAP with Azimuth Cell-type Labels")
)
dev.off()


# 11. Marker check 

DefaultAssay(query_obj) <- "RNA"

marker_genes <- c(
  "CD3D", "IL7R", "LTB",       # T cells
  "NKG7", "GNLY",              # NK
  "MS4A1", "CD79A",            # B cells
  "LYZ", "S100A8", "S100A9",   # myeloid
  "FCER1A", "CST3",            # dendritic
  "HBB", "HBA1"                # erythroid contamination check
)

marker_genes <- marker_genes[marker_genes %in% rownames(query_obj)]

if (length(marker_genes) > 0) {
  pdf(file.path(out_dir, "08_marker_dotplot.pdf"), width = 12, height = 6)
  print(
    DotPlot(query_obj, features = marker_genes, group.by = "cell_type") +
      RotatedAxis() +
      ggtitle("Marker Gene Validation")
  )
  dev.off()
}

# 12. Patient-level cell-type proportions

meta_df <- query_obj@meta.data

cell_counts <- meta_df %>%
  filter(!is.na(cell_type), !is.na(group), !is.na(patient)) %>%
  group_by(group, patient, cell_type) %>%
  summarise(n = n(), .groups = "drop")

total_counts <- meta_df %>%
  filter(!is.na(group), !is.na(patient)) %>%
  group_by(group, patient) %>%
  summarise(total = n(), .groups = "drop")

prop_df <- left_join(cell_counts, total_counts, by = c("group", "patient")) %>%
  mutate(proportion = n / total)

write.csv(prop_df,
          file = file.path(out_dir, "09_celltype_proportions_by_patient.csv"),
          row.names = FALSE)


# 13. Heatmap of cell-type proportions

heatmap_df <- prop_df %>%
  group_by(patient, cell_type) %>%
  summarise(proportion = sum(proportion), .groups = "drop") %>%
  pivot_wider(
    names_from = cell_type,
    values_from = proportion,
    values_fill = 0
  )

heatmap_mat <- as.data.frame(heatmap_df)
rownames(heatmap_mat) <- heatmap_mat$patient
heatmap_mat$patient <- NULL
heatmap_mat <- as.matrix(heatmap_mat)

annotation_row <- data.frame(
  Group = ifelse(grepl("^PTCB", rownames(heatmap_mat)), "PTCB", "TCB")
)
rownames(annotation_row) <- rownames(heatmap_mat)

pdf(file.path(out_dir, "10_celltype_proportion_heatmap.pdf"), width = 10, height = 8)
pheatmap(
  heatmap_mat,
  annotation_row = annotation_row,
  scale = "none",
  main = "Cell-type Proportions by Patient"
)
dev.off()

# 14. Keep robust cell types
celltype_presence <- prop_df %>%
  group_by(cell_type, group) %>%
  summarise(n_patients_present = sum(proportion > 0), .groups = "drop") %>%
  pivot_wider(
    names_from = group,
    values_from = n_patients_present,
    values_fill = 0
  )

keep_celltypes <- celltype_presence %>%
  filter(PTCB >= 3, TCB >= 3) %>%
  pull(cell_type)

prop_df_filt <- prop_df %>%
  filter(cell_type %in% keep_celltypes)


# 15. Differential abundance testing

stat_results <- prop_df_filt %>%
  group_by(cell_type) %>%
  summarise(
    n_PTCB = sum(group == "PTCB"),
    n_TCB = sum(group == "TCB"),
    mean_PTCB = mean(proportion[group == "PTCB"], na.rm = TRUE),
    mean_TCB = mean(proportion[group == "TCB"], na.rm = TRUE),
    p_value = if (n_PTCB >= 2 & n_TCB >= 2) {
      wilcox.test(proportion ~ group)$p.value
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(
    adj_p = p.adjust(p_value, method = "BH"),
    diff = mean_PTCB - mean_TCB,
    abs_diff = abs(diff)
  ) %>%
  filter(mean_PTCB >= 0.02 | mean_TCB >= 0.02) %>%
  arrange(adj_p, desc(abs_diff))

write.csv(stat_results,
          file = file.path(out_dir, "11_celltype_differential_abundance_stats.csv"),
          row.names = FALSE)


# 16. Select top differential cell populations

top_celltypes <- stat_results %>%
  filter(!is.na(adj_p)) %>%
  arrange(adj_p, desc(abs_diff)) %>%
  slice_head(n = 3)

write.csv(top_celltypes,
          file = file.path(out_dir, "12_top3_differential_celltypes.csv"),
          row.names = FALSE)

top_names <- top_celltypes$cell_type


# 17. Plot top differential cell populations

pdf(file.path(out_dir, "13_top_celltype_boxplots.pdf"), width = 10, height = 6)
print(
  ggplot(
    filter(prop_df_filt, cell_type %in% top_names),
    aes(x = group, y = proportion, fill = group)
  ) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.12, size = 2) +
    facet_wrap(~cell_type, scales = "free_y") +
    theme_bw() +
    labs(
      title = "Top Differential Cell Populations: PTCB vs TCB",
      x = "",
      y = "Cell-type proportion"
    )
)
dev.off()

pdf(file.path(out_dir, "14_top_celltypes_highlighted_on_umap.pdf"), width = 8, height = 6)
for (ct in top_names) {
  print(
    DimPlot(
      query_obj,
      reduction = "umap.cca",
      cells.highlight = WhichCells(query_obj, expression = cell_type == ct),
      cols = c("grey90", "blue"),
      raster = FALSE
    ) + ggtitle(ct)
  )
}
dev.off()

# 18. Differential expression in top cell types
# NOTE:
# This is Seurat cell-level DE and is acceptable for
# exploratory analysis. For a stronger final analysis,
# pseudobulk by patient is preferred.

DefaultAssay(query_obj) <- "RNA"
Idents(query_obj) <- "group"

deg_list <- list()

for (ct in top_names) {
  message("Running DE for: ", ct)
  
  sub_obj <- subset(query_obj, subset = cell_type == ct)
  
  # ensure both groups are present
  if (length(unique(sub_obj$group)) == 2) {
    deg_res <- FindMarkers(
      object = sub_obj,
      ident.1 = "PTCB",
      ident.2 = "TCB",
      logfc.threshold = 0.25,
      min.pct = 0.1,
      test.use = "wilcox"
    )
    
    deg_res$gene <- rownames(deg_res)
    deg_res$cell_type <- ct
    deg_list[[ct]] <- deg_res
    
    out_name <- paste0("15_DEG_", gsub("[ /]", "_", ct), "_PTCB_vs_TCB.csv")
    write.csv(deg_res, file = file.path(out_dir, out_name), row.names = FALSE)
  }
}

deg_all <- bind_rows(deg_list)
write.csv(deg_all,
          file = file.path(out_dir, "15_DEG_all_top_celltypes.csv"),
          row.names = FALSE)


# 19. Keep significant DEGs

deg_sig_list <- list()

if (length(deg_list) > 0) {
  fc_col <- if ("avg_log2FC" %in% colnames(deg_list[[1]])) "avg_log2FC" else "avg_logFC"
  
  for (ct in names(deg_list)) {
    deg_sig <- deg_list[[ct]] %>%
      filter(p_val_adj < 0.05, abs(.data[[fc_col]]) > 0.25)
    
    deg_sig_list[[ct]] <- deg_sig
  }
}


# 20. Print/save top genes


top_gene_tables <- list()

for (ct in names(deg_sig_list)) {
  
  # Get top 10 genes sorted by adjusted p-value
  top_tab <- deg_sig_list[[ct]] %>%
    dplyr::arrange(p_val_adj) %>%
    dplyr::select(gene, dplyr::all_of(fc_col), p_val_adj) %>%
    head(10)
  
  # Store in list
  top_gene_tables[[ct]] <- top_tab
  
  # Save as CSV file
  write.csv(
    top_tab,
    file = file.path(out_dir,
                     paste0("16_top10_genes_", gsub("[ /]", "_", ct), ".csv")),
    row.names = FALSE
  )
}




# 21. DotPlots for top DEGs

pdf(file.path(out_dir, "17_top_DEG_dotplots.pdf"), width = 12, height = 7)
for (ct in names(deg_sig_list)) {
  top_genes <- deg_sig_list[[ct]] %>%
    arrange(p_val_adj) %>%
    slice_head(n = 10) %>%
    pull(gene)
  
  top_genes <- intersect(top_genes, rownames(query_obj))
  
  if (length(top_genes) > 0) {
    sub_obj <- subset(query_obj, subset = cell_type == ct)
    
    print(
      DotPlot(sub_obj, features = top_genes, group.by = "group") +
        RotatedAxis() +
        ggtitle(paste("Top DEGs:", ct))
    )
  }
}
dev.off()


# 22. Split DEGs into PTCB-up and TCB-up

deg_split <- list()

for (ct in names(deg_sig_list)) {
  deg_ct <- deg_sig_list[[ct]]
  
  genes_ptcb <- deg_ct %>%
    filter(.data[[fc_col]] > 0) %>%
    pull(gene)
  
  genes_tcb <- deg_ct %>%
    filter(.data[[fc_col]] < 0) %>%
    pull(gene)
  
  deg_split[[ct]] <- list(
    ptcb = unique(genes_ptcb),
    tcb  = unique(genes_tcb)
  )
}


# 23. GO enrichment

go_results <- list()

for (ct in names(deg_split)) {
  genes_ptcb <- deg_split[[ct]]$ptcb
  genes_tcb  <- deg_split[[ct]]$tcb
  
  ego_ptcb <- NULL
  ego_tcb  <- NULL
  
  if (length(genes_ptcb) > 0) {
    ego_ptcb <- enrichGO(
      gene = genes_ptcb,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = "BH"
    )
  }
  
  if (length(genes_tcb) > 0) {
    ego_tcb <- enrichGO(
      gene = genes_tcb,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = "BH"
    )
  }
  
  go_results[[ct]] <- list(ptcb = ego_ptcb, tcb = ego_tcb)
}


# 24. Save GO results tables

for (ct in names(go_results)) {
  if (!is.null(go_results[[ct]]$ptcb) && nrow(as.data.frame(go_results[[ct]]$ptcb)) > 0) {
    write.csv(
      as.data.frame(go_results[[ct]]$ptcb),
      file = file.path(out_dir, paste0("18_GO_PTCB_up_", gsub("[ /]", "_", ct), ".csv")),
      row.names = FALSE
    )
  }
  
  if (!is.null(go_results[[ct]]$tcb) && nrow(as.data.frame(go_results[[ct]]$tcb)) > 0) {
    write.csv(
      as.data.frame(go_results[[ct]]$tcb),
      file = file.path(out_dir, paste0("18_GO_TCB_up_", gsub("[ /]", "_", ct), ".csv")),
      row.names = FALSE
    )
  }
}


# 25. Plot GO enrichment

pdf(file.path(out_dir, "19_GO_dotplots.pdf"), width = 10, height = 7)
for (ct in names(go_results)) {
  if (!is.null(go_results[[ct]]$ptcb) && nrow(as.data.frame(go_results[[ct]]$ptcb)) > 0) {
    print(
      dotplot(go_results[[ct]]$ptcb, showCategory = 10) +
        ggtitle(paste(ct, "- Enriched in PTCB"))
    )
  }
  
  if (!is.null(go_results[[ct]]$tcb) && nrow(as.data.frame(go_results[[ct]]$tcb)) > 0) {
    print(
      dotplot(go_results[[ct]]$tcb, showCategory = 10) +
        ggtitle(paste(ct, "- Enriched in TCB"))
    )
  }
}
dev.off()


























