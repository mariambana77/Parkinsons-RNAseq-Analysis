if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
 BiocManager::install("DESeq2")

 library(DESeq2)
 
 # 2. Create the Metadata
 sample_names <- c("Control_1", "Control_2", "Control_3", "PD_1", "PD_2", "PD_3")
 condition <- factor(c("Control", "Control", "Control", "PD", "PD", "PD"))
 metadata <- data.frame(row.names = sample_names, condition = condition)
 
 # 3. Create the Count Matrix (FIXED & ACCURATE)
 # هنا وضعنا بيانات كل عينة في عمود منفصل لتجنب أي تداخل
 genes <- c("TH", "SNCA", "PINK1", "PRKN", "GAPDH", "ACTB", "LRRK2", "BDNF", "GFAP", "IL6")
 
 counts_df <- data.frame(
   Control_1 = c(1500, 200, 800, 900, 2000, 2500, 300, 1200, 400, 100),
   Control_2 = c(1450, 210, 820, 880, 2050, 2400, 310, 1150, 420, 110),
   Control_3 = c(1600, 195, 790, 910, 1980, 2600, 290, 1250, 380, 95),
   PD_1      = c(200, 1800, 400, 450, 1950, 2450, 900, 300, 1200, 800),
   PD_2      = c(150, 1900, 450, 500, 2100, 2550, 950, 250, 1300, 850),
   PD_3      = c(180, 1750, 380, 420, 2000, 2500, 880, 320, 1250, 780),
   row.names = genes
 )
 
 # تحويل البيانات إلى مصفوفة (Matrix) كما يطلب DESeq2
 counts <- as.matrix(counts_df)
 
 # 4. Construct the DESeqDataSet object
 dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = metadata,
                               design = ~ condition)
 
 # 5. Set the reference level (Control is the baseline)
 dds$condition <- relevel(dds$condition, ref = "Control")
 
 # 6. Run the analysis
 dds <- DESeq(dds)
 
 # --- الآن يمكنك استخراج النتائج ورسم الـ Volcano Plot والـ Heatmap بنس الأكواد السابقة ---
 
 # 1. استخراج النتائج وتحويلها إلى res_df
 res <- results(dds)
 res_ordered <- res[order(res$padj), ]
 res_df <- as.data.frame(res_ordered)
 res_df$gene <- rownames(res_df)
 
 # 2. تحديد الجينات المهمة (Upregulated vs Downregulated)
 res_df$Expression <- "Not Significant"
 res_df$Expression[res_df$log2FoldChange > 1 & res_df$padj < 0.05] <- "Upregulated (e.g. SNCA)"
 res_df$Expression[res_df$log2FoldChange < -1 & res_df$padj < 0.05] <- "Downregulated (e.g. TH)"
 
 # 3. استدعاء مكتبات الرسم
 library(ggplot2)
 library(ggrepel)
 
 # 4. رسم الـ Volcano Plot بالشكل الصحيح
 ggplot(data = res_df, aes(x = log2FoldChange, y = -log10(padj), col = Expression, label = gene)) +
   geom_point(size = 3.5, alpha = 0.8) +
   geom_text_repel(max.overlaps = 15, size = 4, show.legend = FALSE) +
   scale_color_manual(values = c("Downregulated (e.g. TH)" = "#2c7bb6", 
                                 "Not Significant" = "grey", 
                                 "Upregulated (e.g. SNCA)" = "#d7191c")) +
   geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed") +
   geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed") +
   theme_minimal(base_size = 14) +
   labs(title = "Volcano Plot: Parkinson's Disease vs Healthy Controls",
        subtitle = "RNA-seq Differential Gene Expression Analysis",
        x = "Log2 Fold Change (Biological Effect Size)",
        y = "-log10 Adjusted P-value (Statistical Significance)") +
   theme(legend.position = "bottom")
 
 # 1. Install pheatmap if you don't have it
 # install.packages("pheatmap")
 library(pheatmap)
 
 # 2. Extract Normalized Counts
 # We CANNOT use raw counts for a heatmap. We must use the normalized data from DESeq2
 normalized_counts <- counts(dds, normalized = TRUE)
 
 # 3. Select genes to plot
 # In a real dataset of 20,000 genes, we select the top 20-50 significant genes. 
 # Here, we will plot all our 10 targeted genes ordered by significance.
 top_genes <- rownames(res_ordered)
 heatmap_data <- normalized_counts[top_genes, ]
 
 # 4. Create an annotation data frame
 # This adds a color bar at the top to distinguish PD patients from Controls easily
 annotation_col <- data.frame(Condition = metadata$condition)
 rownames(annotation_col) <- colnames(heatmap_data)
 
 # 5. Generate the Heatmap
 pheatmap(heatmap_data,
          cluster_rows = TRUE,     # Group genes that go up/down together
          cluster_cols = TRUE,     # Group samples that have similar expression profiles
          scale = "row",           # CRUCIAL: Z-score scaling (explained below)
          annotation_col = annotation_col,
          main = "Expression Heatmap of Targeted Genes in Parkinson's",
          color = colorRampPalette(c("#2c7bb6", "white", "#d7191c"))(50), # Blue -> White -> Red
          fontsize = 12)