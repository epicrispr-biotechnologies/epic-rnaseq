#' @title Myo vs Testes Volcano Plot
#' @description This script generates a volcano plot comparing myo and testes samples highlighting specific genes
#' @author Tyler Borrman
#' @date 02/22/2025
#' 

library(EnhancedVolcano)

df <- read.table("/home/tylerborrman/g_drive/EPI-321/off_target/pipelines/epic-rnaseq/epic-rnaseq/Human_Myoblast_D7_vs_Testes/deseq2/Human_Myoblast_D7_vs_Testes_degs.tsv",
    sep="\t",
    header=TRUE
)

l2fc_thresh <- 1
padj_thresh <- 0.05
contrast <- "Myotubes vs Leydig Cells"

png("/home/tylerborrman/g_drive/EPI-321/off_target/pipelines/epic-rnaseq/epic-rnaseq/Human_Myoblast_D7_vs_Testes/deseq2/myo_vs_testes_volcano_labs_myo.png",
width = 3000, height = 3000, res = 300)
EnhancedVolcano(df,
    lab = df$gene,
    title = contrast,
    x = "log2FoldChange",
    y = "padj",
    selectLab = c("MYOD1", "ACTA1", 'MYH2'),
    boxedLabels = TRUE,
    caption = paste(" Total genes:", nrow(df),
      " log2FC > ", l2fc_thresh, ": ", nrow(subset(df,log2FoldChange > l2fc_thresh)),
      " log2FC < -", l2fc_thresh, ": ", nrow(subset(df,log2FoldChange < -l2fc_thresh)),
      " significant log2FC > ", l2fc_thresh, ": ", nrow(subset(df, log2FoldChange > l2fc_thresh & padj < padj_thresh)),
      " significant log2FC < -", l2fc_thresh, ": ", nrow(subset(df, log2FoldChange < -l2fc_thresh & padj < padj_thresh))
      ),
    captionLabSize = 11,
    FCcutoff = l2fc_thresh,
    pCutoff = padj_thresh,
    drawConnectors = TRUE,
    typeConnectors = "open"
)
dev.off()

png("/home/tylerborrman/g_drive/EPI-321/off_target/pipelines/epic-rnaseq/epic-rnaseq/Human_Myoblast_D7_vs_Testes/deseq2/myo_vs_testes_volcano_plain.png",
width = 3000, height = 3000, res = 300)
EnhancedVolcano(df,
    lab = "",
    title = contrast,
    x = "log2FoldChange",
    y = "padj",
    caption = paste(" Total genes:", nrow(df),
      " log2FC > ", l2fc_thresh, ": ", nrow(subset(df,log2FoldChange > l2fc_thresh)),
      " log2FC < -", l2fc_thresh, ": ", nrow(subset(df,log2FoldChange < -l2fc_thresh)),
      " significant log2FC > ", l2fc_thresh, ": ", nrow(subset(df, log2FoldChange > l2fc_thresh & padj < padj_thresh)),
      " significant log2FC < -", l2fc_thresh, ": ", nrow(subset(df, log2FoldChange < -l2fc_thresh & padj < padj_thresh))
      ),
    captionLabSize = 11,
    FCcutoff = l2fc_thresh,
    pCutoff = padj_thresh,
    drawConnectors = TRUE,
    typeConnectors = "open"
)
dev.off()

png("/home/tylerborrman/g_drive/EPI-321/off_target/pipelines/epic-rnaseq/epic-rnaseq/Human_Myoblast_D7_vs_Testes/deseq2/myo_vs_testes_volcano_labs_testes.png",
width = 3000, height = 3000, res = 300)
EnhancedVolcano(df,
    lab = df$gene,
    title = contrast,
    x = "log2FoldChange",
    y = "padj",
    selectLab = c("CYP11A1", "HSD3B2", 'INSL3'),
    boxedLabels = TRUE,
    caption = paste(" Total genes:", nrow(df),
      " log2FC > ", l2fc_thresh, ": ", nrow(subset(df,log2FoldChange > l2fc_thresh)),
      " log2FC < -", l2fc_thresh, ": ", nrow(subset(df,log2FoldChange < -l2fc_thresh)),
      " significant log2FC > ", l2fc_thresh, ": ", nrow(subset(df, log2FoldChange > l2fc_thresh & padj < padj_thresh)),
      " significant log2FC < -", l2fc_thresh, ": ", nrow(subset(df, log2FoldChange < -l2fc_thresh & padj < padj_thresh))
      ),
    captionLabSize = 11,
    FCcutoff = l2fc_thresh,
    pCutoff = padj_thresh,
    drawConnectors = TRUE,
    typeConnectors = "open"
)
dev.off()
