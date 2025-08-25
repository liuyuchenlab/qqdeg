####qqdeg
#自定义函数#
#' Title
#'
#' @param file file in document
#' @param object_type gene data or te data
#' @param group1 treat group
#' @param group2 control group
#' @param fc_threshold foldchange value
#' @param species mouse or human
#' @return deg results and enrich results
#' @export
#'
#' @examples
#' example_file <- system.file("extdata", "mouse.xlsx", package = "qqdeg")
#' result <- qqdeg(example_file, "gene", "male-ko", "male-wt", fc_threshold = 1.5,species = "mouse")
###测试###
qqdeg <- function(file, object_type, group1, group2, fc_threshold = 1.5,species = "mouse") {
  # 加载必要的包，并抑制启动消息
  suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(ggpubr)
    library(openxlsx)
    library(clusterProfiler)
    library(ggthemes)
    library(ggrepel)
    library(gmodels)
    library(matrixStats)
    library(DOSE)
    library(topGO)
    library(pathview)
    library(msigdbr)
    if (species == "mouse") {
      library(org.Mm.eg.db)
      OrgDb <- org.Mm.eg.db
      kegg_org <- "mmu"
      gsea_speice <- "mouse"
    } else if (species == "human") {
      library(org.Hs.eg.db)
      OrgDb <- org.Hs.eg.db
      kegg_org <- "hsa"
      gsea_speice <- "human"
    } else {
      stop("species 必须是 'mouse' 或 'human'")
    }
  })



  # 创建输出文件夹
  output_dir <- paste0(group1, "_vs_", group2, "_", object_type)
  dir.create(output_dir, showWarnings = FALSE)

  # 读取数据（自动判断文件类型）
  file_ext <- tools::file_ext(file)  # 获取文件后缀[1](@ref)
  data <- switch(file_ext,
                 "xlsx" = openxlsx::read.xlsx(file, colNames = TRUE, rowNames = FALSE),
                 "csv" = read.csv(file, header = TRUE, stringsAsFactors = FALSE),
                 "txt" = read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE,check.names=F),
                 stop("Unsupported file format. Only .xlsx, .csv and .txt are supported")
  )

  # 根据对象类型提取数据
  if (object_type == "gene") {
    obj_data <- data[grep(data[,1], pattern = ":", invert = TRUE), ]
  } else if (object_type == "te") {
    obj_data <- data[grep(data[,1], pattern = ":", invert = FALSE), ]
  } else {
    stop("Invalid object_type. Must be 'gene' or 'te'.")
  }

  row.names(obj_data) <- obj_data$gene_id
  obj_data <- obj_data[, -1]

  # 提取分组信息（支持多种分隔符：- _ .）
  coldata <- data.frame(
    condition = sub("[-_.]\\d+$", "", colnames(obj_data)),  # 移除末尾的数字及分隔符
    row.names = colnames(obj_data)
  )

  #提取分组矩阵
  target_samples <- rownames(coldata)[coldata$condition %in% c(group1, group2)]
  sub_counts <- obj_data[, target_samples]  # 提取对应样本的表达矩阵

  # 提取分组信息（支持多种分隔符：- _ .）
  coldata <- data.frame(
    condition = sub("[-_.]\\d+$", "", colnames(sub_counts)),  # 移除末尾的数字及分隔符
    row.names = colnames(sub_counts)
  )
  coldata$condition <- as.factor(coldata$condition)

  # 创建DESeq2数据集
  dds <- DESeqDataSetFromMatrix(countData = sub_counts, colData = coldata, design = ~condition)
  keep <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep, ]

  # 标准化和差异表达分析
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", group1, group2))

  #pca
  #VST标准化
  all_vsd <- vst(dds)

  #pca美化
  # 自定义颜色
  pca_group_colors <- c("#1f77b4",
                        "#d62728")
  #美化图片
  pca <- plotPCA(all_vsd, intgroup = "condition") +
    geom_point(aes(color = condition), size = 4) +  # 将颜色映射到 condition
    scale_color_manual(values = pca_group_colors) +  # 应用自定义颜色
    theme_few() +
    theme(legend.position = "top",
          # 调整图例标题的文字大小
          legend.title = element_text(size = 13),
          # 调整图例项目的文字大小
          legend.text = element_text(size = 12),
          # 调整图例符号（小色块或点）的大小
          legend.key.size = unit(1.2, 'cm')
    ) +
    theme(aspect.ratio = 1)  # 设置纵横比为1:1

  # 显示pca 图
  print(pca)

  # 保存pca
  ggsave(filename = file.path(output_dir, paste0("pca_plot_", object_type, "_", group1, "_vs_", group2, ".pdf")),
         plot = pca, width = 8,          # 增大画布宽度
         height = 6,
         dpi = 600,          # 提高分辨率
         device = cairo_pdf)

  # 保存差异分析结果
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized = TRUE)), by = "row.names", sort = FALSE)
  write.xlsx(resdata, file = file.path(output_dir, paste0(object_type, "_", group1, "_vs_", group2, ".xlsx")))

  # 提取显著差异基因
  diff_genes <- subset(res, padj < 0.05 & (log2FoldChange > log2(fc_threshold) | log2FoldChange < -log2(fc_threshold)))
  write.csv(diff_genes, file = file.path(output_dir, paste0("diff_", object_type, "_", group1, "_vs_", group2, ".csv")))

  # 提取上调和下调基因
  up_genes <- subset(res, padj < 0.05 & log2FoldChange > log2(fc_threshold))
  down_genes <- subset(res, padj < 0.05 & log2FoldChange < -log2(fc_threshold))
  write.csv(up_genes, file = file.path(output_dir, paste0(group1, "_vs_", group2, "_up-", object_type, ".csv")))
  write.csv(down_genes, file = file.path(output_dir, paste0(group1, "_vs_", group2, "_down-", object_type, ".csv")))

  # 绘制火山图
  resdata$logP <- -log10(resdata$padj)
  resdata$GROUP <- "not significant"
  resdata$GROUP[(resdata$padj < 0.05) & (resdata$log2FoldChange > log2(fc_threshold))] <- "up-regulated"
  resdata$GROUP[(resdata$padj < 0.05) & (resdata$log2FoldChange < -log2(fc_threshold))] <- "down-regulated"
  resdata$Lable <- ""
  # >>> 新增统计代码 <<<
  up_count <- sum(resdata$GROUP == "up-regulated")
  down_count <- sum(resdata$GROUP == "down-regulated")
  # 去除缺失值
  resdata <- na.omit(resdata)
  resdata <- resdata[order(resdata$padj), ]
  top_up_genes <- head(resdata$Row.names[resdata$GROUP == "up-regulated"], 10)
  top_down_genes <- head(resdata$Row.names[resdata$GROUP == "down-regulated"], 10)
  top_genes <- c(as.character(top_up_genes), as.character(top_down_genes))
  resdata$Lable[match(top_genes, resdata$Row.names)] <- top_genes

  palette <- c("#293890", "#BBBBBB", "#BF1D2D")
  group_colors <- c("down-regulated" = palette[1], "not significant" = palette[2], "up-regulated" = palette[3])

  volcano_plot <- ggplot(resdata, aes(x = log2FoldChange, y = logP)) +
    geom_point(aes(color = GROUP), size = 1) +
    geom_text_repel(aes(label = Lable), size = 3, box.padding = 0.5, max.overlaps = Inf) +
    scale_color_manual(values = group_colors,
                       labels = c(
                         "up-regulated" = paste0("Up regulated (", up_count, ")"),
                         "down-regulated" = paste0("Down regulated (", down_count, ")"),
                         "not significant" = "Not significant"
                       )) +
    theme_classic() +
    geom_hline(yintercept = 1.30, linetype = "dashed") +
    geom_vline(xintercept = c(-log2(fc_threshold), log2(fc_threshold)), linetype = "dashed") +
    expand_limits(y = 0) +
    xlab(bquote(log[2]*FoldChange * "(" ~ .(group1) ~ "/" ~ .(group2) ~ ")")) +
    ylab(expression(-log[10](Adjusted ~ P-value))) +
    theme(aspect.ratio = 1)  # 设置纵横比为1:1

  # 显示火山图
  print(volcano_plot)

  # 保存火山图
  ggsave(filename = file.path(output_dir, paste0("Volcano_plot_", object_type, "_", group1, "_vs_", group2, ".pdf")),
         plot = volcano_plot, width = 8,          # 增大画布宽度
         height = 6,
         dpi = 600,          # 提高分辨率
         device = cairo_pdf)

  # GO富集分析
  go <- NULL
  kk <- NULL
  if (object_type == "gene") {
    up <- read.csv(file.path(output_dir, paste0(group1, "_vs_", group2, "_up-", object_type, ".csv")))
    down <- read.csv(file.path(output_dir, paste0(group1, "_vs_", group2, "_down-", object_type, ".csv")))

    up_genes <- bitr(up$X, fromType = "SYMBOL",
                     toType = "ENTREZID", OrgDb = OrgDb, drop = FALSE)
    down_genes <- bitr(down$X, fromType = "SYMBOL",
                       toType = "ENTREZID", OrgDb = OrgDb, drop = FALSE)

    if (nrow(up_genes) == 0 | nrow(down_genes) == 0) {
      cat("Warning: No genes mapped to ENTREZ IDs. Skipping GO and KEGG analysis.\n")
    } else {
      up.go_all <- enrichGO(gene = up_genes$ENTREZID,
                            OrgDb = OrgDb,
                            keyType = 'ENTREZID', ont = "BP",
                            pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
      down.go_all <- enrichGO(gene = down_genes$ENTREZID,
                              OrgDb = OrgDb,
                              keyType = 'ENTREZID', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

      # 处理GO结果
      up.go_all_readable <- setReadable(up.go_all,
                                        OrgDb = OrgDb, # 确保与您enrichGO时使用的OrgDb一致
                                        keyType = "ENTREZID") # keyType 也需一致
      down.go_all_readable <- setReadable(down.go_all,
                                          OrgDb = OrgDb, # 确保与您enrichGO时使用的OrgDb一致
                                          keyType = "ENTREZID") # keyType 也需一致
      # 现在提取结果
      go.up <- up.go_all_readable@result
      go.up <- go.up[order(go.up$p.adjust), ]
      write.xlsx(go.up, file = file.path(output_dir, paste0("go_up", "_", group1, "_vs_", group2, ".xlsx")))
      go.up$value <- -log10(go.up$p.adjust)
      go.down <- down.go_all_readable@result
      go.down <- go.down[order(go.down$p.adjust), ]
      write.xlsx(go.down, file = file.path(output_dir, paste0("go_down", "_", group1, "_vs_", group2, ".xlsx")))
      go.up$value <- -log10(go.up$p.adjust)
      go.down$new <- -log10(go.down$p.adjust)
      go.down$value <- -(go.down$new)
      go.down <- go.down[ , -which(colnames(go.down) %in% "new")]


      # 合并显著的上下调
      go.up <- head(go.up, 5)
      go.down <- head(go.down, 5)
      go <- rbind(go.up, go.down)

      ## 调整因子水平
      go <- go[order(go$value),]
      # 计算每个 Description 的频率
      desc_counts <- table(go$Description)
      # 筛选出只出现一次的 Description
      unique_desc <- names(desc_counts[desc_counts == 1])
      # 根据唯一的 Description 过滤 go 数据框
      go <- go[go$Description %in% unique_desc, ]
      name <- go$Description
      go$Description <- factor(go$Description, levels = name)

      # 绘制GO富集图
      go_plot <- ggplot(go, aes(x = value, y = Description)) +
        xlab('Enrich Value') +
        ylab("") +
        geom_bar(stat = "identity", aes(fill = ifelse(value < 0, "Down", "Up"))) +
        scale_fill_manual(name = "", values = c("Down" = "#293890", "Up" = "#BF1D2D")) +
        theme_few() +
        ggtitle(paste(group1, "vs", group2, "GO BP Enrichment Analysis")) +
        theme(aspect.ratio = 1)  # 设置纵横比为1:1

      print(go_plot)
      ggsave(filename = file.path(output_dir, paste0("GO_BP_Enrichment_", object_type, "_", group1, "_vs_", group2, ".pdf")), plot = go_plot,width = 8,          # 增大画布宽度
             height = 6,
             dpi = 600,          # 提高分辨率
             device = cairo_pdf)

      #######kegg
      #上调kegg
      up.kk <- enrichKEGG(gene = up_genes$ENTREZID,organism =  kegg_org,  pvalueCutoff = 0.05)
      #下调kegg
      down.kk <- enrichKEGG(gene = down_genes$ENTREZID,organism =  kegg_org,  pvalueCutoff = 0.05)
      #画图
      up.kk_readable <- setReadable(up.kk,
                                    OrgDb = OrgDb, # 确保与您enrichGO时使用的OrgDb一致
                                    keyType = "ENTREZID") # keyType 也需一致
      down.kk_readable <- setReadable(down.kk,
                                      OrgDb = OrgDb, # 确保与您enrichGO时使用的OrgDb一致
                                      keyType = "ENTREZID") # keyType 也需一致
      #提取数据
      upkegg <- up.kk_readable@result
      up.kegg <- upkegg[order(upkegg$p.adjust),]
      write.xlsx(upkegg, file = file.path(output_dir, paste0("kegg_up", "_", group1, "_vs_", group2, ".xlsx")))
      up.kegg$value <- -log10(up.kegg$p.adjust)

      downkegg <- down.kk_readable@result
      down.kegg <- downkegg[order(downkegg$p.adjust),]
      write.xlsx(downkegg, file = file.path(output_dir, paste0("kegg_down", "_", group1, "_vs_", group2, ".xlsx")))
      down.kegg$value <- log10(down.kegg$p.adjust)


      ####合并显著的上下调
      kegg.up <-  head(up.kegg,5)
      kegg.down <- head(down.kegg,5)
      kegg <- rbind(kegg.up,kegg.down)


      kegg$Description <- gsub("-.*", "", kegg$Description)
      ####BAR
      ## 调整因子水平
      kegg <- kegg[order(kegg$value),]
      # 计算每个 Description 的频率
      desc_counts <- table(kegg$Description)
      # 筛选出只出现一次的 Description
      unique_desc <- names(desc_counts[desc_counts == 1])
      # 根据唯一的 Description 过滤 kegg 数据框
      kegg <- kegg[kegg$Description %in% unique_desc, ]
      name <- kegg$Description
      kegg$Description <- factor(kegg$Description, levels = name)
      #绘制富集图
      kegg_plot <- ggplot(kegg, aes(x = value, y = Description),col = col) +
        xlab('Enrich Value') +
        ylab("") +
        geom_bar(stat = "identity", aes(fill = ifelse(value < 0, "Down", "Up"))) +
        scale_fill_manual(name = "", values = c("Down" ="#293890", "Up" = "#BF1D2D")) +
        theme_few() +
        ggtitle(paste(group1, "vs", group2, "KEGG Enrichment Analysis")) +
        theme(aspect.ratio = 1)  # 设置纵横比为1:1

      print(kegg_plot)

      ggsave(filename = file.path(output_dir, paste0("KEGG_Enrichment_", object_type, "_", group1, "_vs_", group2, ".pdf")), plot = kegg_plot,width = 8,          # 增大画布宽度
             height = 6,
             dpi = 600,          # 提高分辨率
             device = cairo_pdf)
      # 生成prerank文件（不筛选）
      if (object_type == "gene") {
        prerank_data <- resdata[, c("Row.names", "log2FoldChange")]
        prerank_data <- prerank_data[order(prerank_data$log2FoldChange, decreasing = TRUE), ]
        write.csv(prerank_data,
                  file.path(output_dir, paste0("prerank_", group1, "_vs_", group2, ".csv")),
                  row.names = FALSE)
      }

      # HALLMARK GSEA分析
      if (object_type == "gene") {
        # 基因列表构建
        genelist <- prerank_data$log2FoldChange
        names(genelist) <- prerank_data$Row.names
        genelist <- sort(genelist, decreasing = TRUE)

        # 获取HALLMARK基因集
        hallmark_gene_sets <- msigdbr(species = gsea_speice, category = "H")
        geneset <- data.frame(
          term = gsub("HALLMARK_", "", hallmark_gene_sets$gs_name),
          gene = hallmark_gene_sets$gene_symbol
        )

        # 执行GSEA分析
        egmt <- GSEA(
          genelist,
          TERM2GENE = geneset,
          pvalueCutoff = 1,
          minGSSize = 1,
          maxGSSize = 500000  # 保持原始参数设置
        )

        # 结果可视化
        if (nrow(egmt@result) > 0) {
          data <- egmt@result[, c("ID", "NES", "setSize", "pvalue","p.adjust","core_enrichment")]
          data <- data[order(data$NES, decreasing = TRUE), ]

          write.xlsx(data, file = file.path(output_dir, paste0("hallmark_gsea", "_", group1, "_vs_", group2, ".xlsx")))

          data$ID <- factor(data$ID, levels = data$ID)
          data$xlab <- 1:nrow(data)
          # Calculate y-axis limits outside of ggplot chain
          y_min <- floor(min(data$NES))
          y_max <- ceiling(max(data$NES))# Select pathways to label - only top and bottom 3 by NES score
          top_pathways <- head(as.character(data$ID), 3)  # Top 3 pathways with highest NES
          bottom_pathways <- tail(as.character(data$ID), 3)  # Bottom 3 pathways with lowest NES
          label <- c(top_pathways, bottom_pathways)
          data_label <- data[data$ID %in% label, ]
          data_label$col <- rep(c("#F6631C", "#2C91E0", "#F3A332", "#018A67"), length.out = nrow(data_label))
          data_label
          fill_color <- "#6D65A3"
          # 可视化参数设置
          top_n <- 6
          label_data <- rbind(
            head(data[order(data$NES, decreasing = TRUE), ], top_n/2),
            tail(data[order(data$NES, decreasing = TRUE), ], top_n/2)
          )

          # 生成气泡图
          p_gsea <- ggplot(data, aes(x = xlab, y = NES, color = NES, size = setSize)) +
            geom_point(aes(alpha = -log10(p.adjust)),
                       shape = 21, stroke = 0.7,
                       fill = fill_color, colour = "black") +
            scale_alpha_continuous(range = c(0.1, 0.9), name = "Significance\n(-log10 p-value)") +
            labs(
              title = paste(group1, "vs", group2, "Pathway Enrichment Analysis"),
              x = "Gene Set Rank",
              y = "Normalized Enrichment Score (NES)",
              size = "Gene Set Size",
              color = "NES"
            ) +
            theme_classic(base_size = 15) +
            # Adjust axes
            scale_x_continuous(breaks = seq(0, 50, by = 10), labels = seq(0, 50, by = 10)) +
            # Dynamic Y-axis scaling based on pre-calculated limits
            scale_y_continuous(
              breaks = seq(y_min, y_max, by = 1),
              labels = seq(y_min, y_max, by = 1),
              limits = c(y_min, y_max),
              expand = c(0.1, 0.1)  # Add a bit of padding
            ) +
            # Customize legend
            guides(
              size = guide_legend(title = "Gene Set Size", order = 1),
              alpha = guide_legend(title = "Significance\n(-log10 p.adjust)", order = 2)
            ) +
            # Beautify theme
            theme(
              plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
              plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
              axis.line = element_line(color = "black", size = 0.6),  # Bold axis lines
              axis.text = element_text(face = "bold"),               # Bold axis text
              axis.title = element_text(size = 14, face = "bold"),   # Set axis title
              legend.title = element_text(face = "bold"),            # Bold legend title
              legend.position = "right",                             # Place legend on the right
              legend.box = "vertical",                               # Vertical legend arrangement
              panel.grid = element_blank(),                          # Remove all grid lines
              panel.background = element_rect(fill = "white")        # Ensure white background
            ) +
            geom_text_repel(data = data_label,
                            aes(x = xlab, y = NES, label = ID),
                            size = 3,
                            color = data_label$col,  # Using the custom colors defined above
                            force = 20,                # Label repulsion force (increased as in reference)
                            point.padding = 0.5,       # Minimum distance between labels and points
                            min.segment.length = 0,    # Minimum guide line length
                            hjust = 1.2,               # Horizontal alignment (from reference code)
                            segment.color = "grey20",  # Guide line color
                            segment.size = 0.3,        # Guide line thickness (as in reference)
                            segment.alpha = 0.8,       # Guide line transparency
                            nudge_y = -0.1) +
            theme(aspect.ratio = 1)  # 设置纵横比为1:1

          print(p_gsea)
          # 保存结果
          ggsave(file.path(output_dir,
                           paste0("HALLMARK_GSEA_", group1, "_vs_", group2, ".pdf")),
                 plot = p_gsea, width = 8, height = 6, dpi = 600, device = cairo_pdf)
        }
      }
    }
  }

  # 将所有结果保存到列表
  result <- list(
    pca = pca,
    volcano_plot = volcano_plot,
    go_plot = if (exists("go_plot")) go_plot else NULL,
    kegg_plot = if (exists("kegg_plot")) kegg_plot else NULL,
    diff_genes = diff_genes,
    up_genes = up_genes,
    down_genes = down_genes,
    obj_data = obj_data,
    coldata = coldata,
    res = res,
    resdata = resdata,
    go = if (exists("go")) go else NULL,
    go.up = if (exists("go.up")) go.up else NULL,
    go.down = if (exists("go.down")) go.down else NULL,
    kegg = if (exists("kegg")) kegg else NULL,
    kegg.up = if (exists("kegg.up")) kegg.up else NULL,
    kegg.down = if (exists("kegg.down")) kegg.down else NULL,
    gesa = if (exists("data")) data else NULL,
    p_gsea = if (exists("p_gsea")) p_gsea else NULL
  )
  # 打印完毕消息
  # 创建虚线
  line <- strrep("-", 50)

  # 输出虚线和文本
  cat(line, "\n")
  cat("数据分析完毕，已保存！\n")
  cat(line, "\n")
  return(result)

}
#


