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
qqdeg <- function(file, object_type, group1, group2, fc_threshold = 1.5, species = "mouse") {

  ###########################################
  # 1. 加载依赖包
  ###########################################
  suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(ggpubr)
    library(openxlsx)
    library(clusterProfiler)
    library(ggthemes)
    library(ggrepel)
    library(DOSE)
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
      stop("物种参数错误！仅支持 'mouse' 或 'human'")
    }
  })

  ###########################################
  # 2. 初始化输出目录
  ###########################################
  file_name <- sub("\\.[^.]*$", "", basename(file))
  output_dir <- paste0(file_name, "_", group1, "_vs_", group2, "_", object_type)
  dir.create(output_dir, showWarnings = FALSE)
  cat("结果将保存至：", output_dir, "\n")

  ###########################################
  # 3. 读取输入数据
  ###########################################
  file_ext <- tools::file_ext(file)
  data <- switch(
    file_ext,
    "xlsx" = openxlsx::read.xlsx(file, colNames = TRUE),
    "csv" = read.csv(file, header = TRUE, check.names = FALSE),
    "txt" = read.table(file, header = TRUE, sep = "\t", check.names = FALSE),
    stop("不支持的文件格式！仅支持 .xlsx, .csv, .txt")
  )

  ###########################################
  # 4. 提取基因/TE数据
  ###########################################
  if (object_type == "gene") {
    obj_data <- data[grep(data[, 1], pattern = ":", invert = TRUE), ]
  } else if (object_type == "te") {
    obj_data <- data[grep(data[, 1], pattern = ":", invert = FALSE), ]
  } else {
    stop("object_type错误！仅支持 'gene' 或 'te'")
  }
  row.names(obj_data) <- obj_data$gene_id
  obj_data <- obj_data[, -1]

  ###########################################
  # 5. 分组信息处理
  ###########################################
  coldata <- data.frame(
    condition = sub("[-_.]\\d+$", "", colnames(obj_data)),
    row.names = colnames(obj_data)
  )
  target_samples <- rownames(coldata)[coldata$condition %in% c(group1, group2)]
  sub_counts <- obj_data[, target_samples]
  coldata <- data.frame(
    condition = sub("[-_.]\\d+$", "", colnames(sub_counts)),
    row.names = colnames(sub_counts)
  )
  coldata$condition <- as.factor(coldata$condition)

  ###########################################
  # 6. DESeq2差异表达分析
  ###########################################
  dds <- DESeqDataSetFromMatrix(
    countData = sub_counts,
    colData = coldata,
    design = ~condition
  )
  keep <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep, ]
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", group1, group2))
  cat("过滤后检测到的基因/TE数量：", nrow(dds), "\n")

  ###########################################
  # 7. 数据标准化（PCA用）
  ###########################################
  threshold <- 2000
  if (nrow(dds) >= threshold) {
    tryCatch({
      all_norm <- vst(dds)
      cat("使用VST标准化\n")
    }, error = function(e) {
      cat("VST失败，切换到RLOG\n")
      all_norm <- rld(dds)
    })
  } else {
    cat("使用RLog标准化\n")
    all_norm <- rlog(dds)
  }

  ###########################################
  # 8. 绘制PCA图
  ###########################################
  pca <- plotPCA(all_norm, intgroup = "condition") +
    geom_point(aes(color = condition), size = 4) +
    scale_color_manual(values = c("#1f77b4", "#d62728")) +
    theme_few() +
    theme(legend.position = "top") +
    theme(aspect.ratio = 1)
  print(pca)
  ggsave(
    file.path(output_dir, paste0("pca_plot_", object_type, "_", group1, "_vs_", group2, ".pdf")),
    plot = pca, width = 8, height = 6, units = "in", dpi = 600, device = cairo_pdf
  )

  ###########################################
  # 9. 整理差异分析结果并保存
  ###########################################
  resdata <- merge(
    as.data.frame(res),
    as.data.frame(counts(dds, normalized = TRUE)),
    by = "row.names", sort = FALSE
  )
  write.xlsx(
    resdata,
    file.path(output_dir, paste0(object_type, "_", group1, "_vs_", group2, ".xlsx"))
  )

  diff_genes <- subset(res, padj < 0.05 & (log2FoldChange > log2(fc_threshold) | log2FoldChange < -log2(fc_threshold)))
  up_genes <- subset(res, padj < 0.05 & log2FoldChange > log2(fc_threshold))
  down_genes <- subset(res, padj < 0.05 & log2FoldChange < -log2(fc_threshold))

  write.csv(diff_genes, file.path(output_dir, paste0("diff_", object_type, "_", group1, "_vs_", group2, ".csv")))
  write.csv(up_genes, file.path(output_dir, paste0(group1, "_vs_", group2, "_up-", object_type, ".csv")))
  write.csv(down_genes, file.path(output_dir, paste0(group1, "_vs_", group2, "_down-", object_type, ".csv")))

  ###########################################
  # 10. 绘制火山图
  ###########################################
  resdata$logP <- -log10(resdata$padj)
  resdata$GROUP <- "not significant"
  resdata$GROUP[(resdata$padj < 0.05) & (resdata$log2FoldChange > log2(fc_threshold))] <- "up-regulated"
  resdata$GROUP[(resdata$padj < 0.05) & (resdata$log2FoldChange < -log2(fc_threshold))] <- "down-regulated"
  resdata$Lable <- ""

  up_count <- sum(resdata$GROUP == "up-regulated")
  down_count <- sum(resdata$GROUP == "down-regulated")
  resdata <- na.omit(resdata)
  resdata <- resdata[order(resdata$padj), ]
  top_up_genes <- head(resdata$Row.names[resdata$GROUP == "up-regulated"], 10)
  top_down_genes <- head(resdata$Row.names[resdata$GROUP == "down-regulated"], 10)
  top_genes <- c(as.character(top_up_genes), as.character(top_down_genes))
  resdata$Lable[match(top_genes, resdata$Row.names)] <- top_genes

  palette <- c("#293890", "#BBBBBB", "#BF1D2D")
  group_colors <- c(
    "down-regulated" = palette[1],
    "not significant" = palette[2],
    "up-regulated" = palette[3]
  )
  volcano_plot <- ggplot(resdata, aes(x = log2FoldChange, y = logP)) +
    geom_point(aes(color = GROUP), size = 1) +
    geom_text_repel(aes(label = Lable), size = 3, max.overlaps = Inf) +
    scale_color_manual(
      values = group_colors,
      labels = c(
        "up-regulated" = paste0("Up (", up_count, ")"),
        "down-regulated" = paste0("Down (", down_count, ")"),
        "not significant" = "Not significant"
      )
    ) +
    theme_classic() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-log2(fc_threshold), log2(fc_threshold)), linetype = "dashed") +
    xlab(bquote(log[2]*"FoldChange (" ~ .(group1) ~ "/" ~ .(group2) ~ ")")) +
    ylab(expression(-log[10]("Adjusted P-value"))) +
    theme(aspect.ratio = 1)
  print(volcano_plot)
  ggsave(
    file.path(output_dir, paste0("Volcano_plot_", object_type, "_", group1, "_vs_", group2, ".pdf")),
    plot = volcano_plot, width = 8, height = 6, units = "in", dpi = 600, device = cairo_pdf
  )

  ###########################################
  # 11. GO和KEGG富集分析（严格卡p.adjust < 0.05）
  ###########################################
  go <- NULL; go_plot <- NULL; go.up <- NULL; go.down <- NULL
  go_up_plot <- NULL; go_down_plot <- NULL
  kegg <- NULL; kegg_plot <- NULL; kegg.up <- NULL; kegg.down <- NULL
  kegg_up_plot <- NULL; kegg_down_plot <- NULL

  if (object_type == "gene") {
    # 转换为数据框并保留基因ID（原up_genes的行名是基因ID，转为"X"列，与原逻辑一致）
    up <- as.data.frame(up_genes)
    up$X <- rownames(up)
    down <- as.data.frame(down_genes)
    down$X <- rownames(down)
    ###########################################
    # 11.1 基因ID映射
    ###########################################
    up_genes_entrez <- bitr(up$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb, drop = FALSE)
    down_genes_entrez <- bitr(down$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb, drop = FALSE)

    if (nrow(up_genes_entrez) == 0 | nrow(down_genes_entrez) == 0) {
      cat("警告：基因ID映射失败，跳过富集分析！\n")
    } else {

      ###########################################
      # 11.2 GO富集分析（BP）- 严格卡p.adjust < 0.05
      ###########################################
      # 上调基因GO富集（函数内先按pvalueCutoff=0.05筛选）
      up.go_all <- enrichGO(
        gene = up_genes_entrez$ENTREZID,
        OrgDb = OrgDb,
        keyType = 'ENTREZID',
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,  # 第一步筛选：p值阈值
        qvalueCutoff = 0.2    # 放宽q值，避免漏筛，后续再用p.adjust过滤
      )

      # 下调基因GO富集
      down.go_all <- enrichGO(
        gene = down_genes_entrez$ENTREZID,
        OrgDb = OrgDb,
        keyType = 'ENTREZID',
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      )

      # 整理上调GO结果：**二次过滤p.adjust < 0.05**
      if (!is.null(up.go_all) && nrow(up.go_all@result) > 0) {
        up.go_all_readable <- setReadable(up.go_all, OrgDb = OrgDb, keyType = "ENTREZID")
        go.up <- up.go_all_readable@result
        go.up <- go.up[go.up$p.adjust < 0.05, ]  # 严格过滤：仅保留p.adjust < 0.05
        if (nrow(go.up) == 0) {  # 若过滤后无结果，标记为NULL
          go.up <- NULL
          cat("上调基因GO富集无显著条目（p.adjust < 0.05），跳过...\n")
        } else {
          go.up <- go.up[order(go.up$p.adjust), ]  # 按p.adjust升序（更显著在前）
          write.xlsx(go.up, file.path(output_dir, paste0("go_up", "_", group1, "_vs_", group2, ".xlsx")))
          go.up$value <- -log10(go.up$p.adjust)  # 富集值（-log10(p.adjust)）
        }
      } else {
        go.up <- NULL
        cat("上调基因无GO富集结果，跳过...\n")
      }

      # 整理下调GO结果：**二次过滤p.adjust < 0.05**
      if (!is.null(down.go_all) && nrow(down.go_all@result) > 0) {
        down.go_all_readable <- setReadable(down.go_all, OrgDb = OrgDb, keyType = "ENTREZID")
        go.down <- down.go_all_readable@result
        go.down <- go.down[go.down$p.adjust < 0.05, ]  # 严格过滤
        if (nrow(go.down) == 0) {
          go.down <- NULL
          cat("下调基因GO富集无显著条目（p.adjust < 0.05），跳过...\n")
        } else {
          go.down <- go.down[order(go.down$p.adjust), ]
          write.xlsx(go.down, file.path(output_dir, paste0("go_down", "_", group1, "_vs_", group2, ".xlsx")))
          go.down$value <- -log10(go.down$p.adjust)
          go.down$value <- -go.down$value  # 下调用负值标记
        }
      } else {
        go.down <- NULL
        cat("下调基因无GO富集结果，跳过...\n")
      }

      ###########################################
      # 11.3 绘制单独的GO上调/下调图（仅用显著条目）
      ###########################################
      # 上调GO图（仅用p.adjust < 0.05的条目）
      if (!is.null(go.up) && nrow(go.up) > 0) {
        go_up_top <- head(go.up, 10)  # 取前10个最显著（p.adjust最小）的条目
        go_up_top <- go_up_top[order(go_up_top$value, decreasing = TRUE), ]  # 按富集值排序
        go_up_top$Description <- factor(go_up_top$Description, levels = rev(go_up_top$Description))

        go_up_plot <- ggplot(go_up_top, aes(x = value, y = Description)) +
          geom_bar(stat = "identity", fill = "#BF1D2D") +
          xlab('Enrich Value (-log10(adjusted p-value))') +
          ylab("") +
          ggtitle(paste(group1, "vs", group2, "GO BP (Up-regulated)")) +
          theme_few() +
          theme(aspect.ratio = 1)
        print(go_up_plot)
        ggsave(
          file.path(output_dir, paste0("GO_BP_Up_", object_type, "_", group1, "_vs_", group2, ".pdf")),
          plot = go_up_plot, width = 8, height = 6, units = "in", dpi = 600, device = cairo_pdf
        )
      } else {
        cat("无显著上调GO条目（p.adjust < 0.05），跳过单独图...\n")
      }

      # 下调GO图（仅用p.adjust < 0.05的条目）
      if (!is.null(go.down) && nrow(go.down) > 0) {
        go_down_top <- head(go.down, 10)
        go_down_top <- go_down_top[order(abs(go_down_top$value), decreasing = TRUE), ]
        go_down_top$Description <- factor(go_down_top$Description, levels = rev(go_down_top$Description))

        go_down_plot <- ggplot(go_down_top, aes(x = -value, y = Description)) +
          geom_bar(stat = "identity", fill = "#293890") +
          xlab('Enrich Value (-log10(adjusted p-value))') +
          ylab("") +
          ggtitle(paste(group1, "vs", group2, "GO BP (Down-regulated)")) +
          theme_few() +
          theme(aspect.ratio = 1)
        print(go_down_plot)
        ggsave(
          file.path(output_dir, paste0("GO_BP_Down_", object_type, "_", group1, "_vs_", group2, ".pdf")),
          plot = go_down_plot, width = 8, height = 6, units = "in", dpi = 600, device = cairo_pdf
        )
      } else {
        cat("无显著下调GO条目（p.adjust < 0.05），跳过单独图...\n")
      }

      ###########################################
      # 11.4 GO合并图（仅用p.adjust < 0.05的条目）
      ###########################################
      if (!is.null(go.up) && !is.null(go.down) && nrow(go.up) > 0 && nrow(go.down) > 0) {
        # 从显著条目中各取前5个
        go.up_sub <- head(go.up, 5)
        go.down_sub <- head(go.down, 5)
        go <- rbind(go.up_sub, go.down_sub)
        # 去重（同一通路只保留一次）
        desc_counts <- table(go$Description)
        unique_desc <- names(desc_counts[desc_counts == 1])
        go <- go[go$Description %in% unique_desc, ]
        # 排序并绘图
        go <- go[order(go$value), ]
        go$Description <- factor(go$Description, levels = go$Description)
        go_plot <- ggplot(go, aes(x = value, y = Description, fill = ifelse(value < 0, "Down", "Up"))) +
          xlab('Enrich Value') + ylab("") +
          geom_bar(stat = "identity") +
          scale_fill_manual(name = "", values = c("Down" = "#293890", "Up" = "#BF1D2D")) +
          theme_few() +
          ggtitle(paste(group1, "vs", group2, "GO BP Enrichment")) +
          theme(aspect.ratio = 1)
        print(go_plot)
        ggsave(
          file.path(output_dir, paste0("GO_BP_Combined_", object_type, "_", group1, "_vs_", group2, ".pdf")),
          plot = go_plot, width = 8, height = 6, units = "in", dpi = 600, device = cairo_pdf
        )
      } else {
        cat("无足够显著GO条目（p.adjust < 0.05），跳过合并图...\n")
      }

      ###########################################
      # 11.5 KEGG富集分析 - 严格卡p.adjust < 0.05
      ###########################################
      # 上调基因KEGG富集
      up.kk <- enrichKEGG(
        gene = up_genes_entrez$ENTREZID,
        organism = kegg_org,
        pvalueCutoff = 0.05  # 第一步筛选
      )

      # 下调基因KEGG富集
      down.kk <- enrichKEGG(
        gene = down_genes_entrez$ENTREZID,
        organism = kegg_org,
        pvalueCutoff = 0.05
      )

      # 整理上调KEGG结果：**二次过滤p.adjust < 0.05**
      if (!is.null(up.kk) && nrow(up.kk@result) > 0) {
        up.kk_readable <- setReadable(up.kk, OrgDb = OrgDb, keyType = "ENTREZID")
        upkegg <- up.kk_readable@result
        up.kegg <- upkegg[upkegg$p.adjust < 0.05, ]  # 严格过滤
        if (nrow(up.kegg) == 0) {
          up.kegg <- NULL
          cat("上调基因KEGG富集无显著条目（p.adjust < 0.05），跳过...\n")
        } else {
          up.kegg <- up.kegg[order(up.kegg$p.adjust), ]
          write.xlsx(upkegg, file.path(output_dir, paste0("kegg_up", "_", group1, "_vs_", group2, ".xlsx")))
          up.kegg$value <- -log10(up.kegg$p.adjust)
        }
      } else {
        up.kegg <- NULL
        cat("上调基因无KEGG富集结果，跳过...\n")
      }

      # 整理下调KEGG结果：**二次过滤p.adjust < 0.05**
      if (!is.null(down.kk) && nrow(down.kk@result) > 0) {
        down.kk_readable <- setReadable(down.kk, OrgDb = OrgDb, keyType = "ENTREZID")
        downkegg <- down.kk_readable@result
        down.kegg <- downkegg[downkegg$p.adjust < 0.05, ]  # 严格过滤
        if (nrow(down.kegg) == 0) {
          down.kegg <- NULL
          cat("下调基因KEGG富集无显著条目（p.adjust < 0.05），跳过...\n")
        } else {
          down.kegg <- down.kegg[order(down.kegg$p.adjust), ]
          write.xlsx(downkegg, file.path(output_dir, paste0("kegg_down", "_", group1, "_vs_", group2, ".xlsx")))
          down.kegg$value <- -log10(down.kegg$p.adjust)
          down.kegg$value <- -down.kegg$value  # 下调用负值标记
        }
      } else {
        down.kegg <- NULL
        cat("下调基因无KEGG富集结果，跳过...\n")
      }

      ###########################################
      # 11.6 绘制单独的KEGG上调/下调图（仅用显著条目）
      ###########################################
      if (!is.null(up.kegg) && nrow(up.kegg) > 0) {
        kegg_up_top <- head(up.kegg, 10)
        kegg_up_top <- kegg_up_top[order(kegg_up_top$value, decreasing = TRUE), ]
        kegg_up_top$Description <- factor(kegg_up_top$Description, levels = rev(kegg_up_top$Description))

        kegg_up_plot <- ggplot(kegg_up_top, aes(x = value, y = Description)) +
          geom_bar(stat = "identity", fill = "#BF1D2D") +
          xlab('Enrich Value (-log10(adjusted p-value))') +
          ylab("") +
          ggtitle(paste(group1, "vs", group2, "KEGG (Up-regulated)")) +
          theme_few() +
          theme(aspect.ratio = 1)
        print(kegg_up_plot)
        ggsave(
          file.path(output_dir, paste0("KEGG_Up_", object_type, "_", group1, "_vs_", group2, ".pdf")),
          plot = kegg_up_plot, width = 8, height = 6, units = "in", dpi = 600, device = cairo_pdf
        )
      } else {
        cat("无显著上调KEGG条目（p.adjust < 0.05），跳过单独图...\n")
      }

      if (!is.null(down.kegg) && nrow(down.kegg) > 0) {
        kegg_down_top <- head(down.kegg, 10)
        kegg_down_top <- kegg_down_top[order(abs(kegg_down_top$value), decreasing = TRUE), ]
        kegg_down_top$Description <- factor(kegg_down_top$Description, levels = rev(kegg_down_top$Description))

        kegg_down_plot <- ggplot(kegg_down_top, aes(x = -value, y = Description)) +
          geom_bar(stat = "identity", fill = "#293890") +
          xlab('Enrich Value (-log10(adjusted p-value))') +
          ylab("") +
          ggtitle(paste(group1, "vs", group2, "KEGG (Down-regulated)")) +
          theme_few() +
          theme(aspect.ratio = 1)
        print(kegg_down_plot)
        ggsave(
          file.path(output_dir, paste0("KEGG_Down_", object_type, "_", group1, "_vs_", group2, ".pdf")),
          plot = kegg_down_plot, width = 8, height = 6, units = "in", dpi = 600, device = cairo_pdf
        )
      } else {
        cat("无显著下调KEGG条目（p.adjust < 0.05），跳过单独图...\n")
      }

      ###########################################
      # 11.7 KEGG合并图（仅用p.adjust < 0.05的条目）
      ###########################################
      if (!is.null(up.kegg) && !is.null(down.kegg) && nrow(up.kegg) > 0 && nrow(down.kegg) > 0) {
        kegg.up_sub <- head(up.kegg, 5)
        kegg.down_sub <- head(down.kegg, 5)
        kegg <- rbind(kegg.up_sub, kegg.down_sub)
        kegg$Description <- gsub("-.*", "", kegg$Description)
        # 去重
        desc_counts <- table(kegg$Description)
        unique_desc <- names(desc_counts[desc_counts == 1])
        kegg <- kegg[kegg$Description %in% unique_desc, ]
        # 排序并绘图
        kegg <- kegg[order(kegg$value), ]
        kegg$Description <- factor(kegg$Description, levels = kegg$Description)
        kegg_plot <- ggplot(kegg, aes(x = value, y = Description, fill = ifelse(value < 0, "Down", "Up"))) +
          xlab('Enrich Value') + ylab("") +
          geom_bar(stat = "identity") +
          scale_fill_manual(name = "", values = c("Down" = "#293890", "Up" = "#BF1D2D")) +
          theme_few() +
          ggtitle(paste(group1, "vs", group2, "KEGG Enrichment")) +
          theme(aspect.ratio = 1)
        print(kegg_plot)
        ggsave(
          file.path(output_dir, paste0("KEGG_Combined_", object_type, "_", group1, "_vs_", group2, ".pdf")),
          plot = kegg_plot, width = 8, height = 6, units = "in", dpi = 600, device = cairo_pdf
        )
      } else {
        cat("无足够显著KEGG条目（p.adjust < 0.05），跳过合并图...\n")
      }

      ###########################################
      # 12. GSEA分析（卡p.adjust < 0.25）
      ###########################################
      gsea_data <- NULL; p_gsea <- NULL
      tryCatch({
        prerank_data <- resdata[, c("Row.names", "log2FoldChange")]
        prerank_data <- prerank_data[order(prerank_data$log2FoldChange, decreasing = TRUE), ]
        write.csv(prerank_data, file.path(output_dir, paste0("prerank_", group1, "_vs_", group2, ".csv")), row.names = FALSE)

        genelist <- prerank_data$log2FoldChange
        names(genelist) <- prerank_data$Row.names
        genelist <- sort(genelist, decreasing = TRUE)

        hallmark_gene_sets <- msigdbr(species = gsea_speice, category = "H")
        geneset <- data.frame(
          term = gsub("HALLMARK_", "", hallmark_gene_sets$gs_name),
          gene = hallmark_gene_sets$gene_symbol
        )

        egmt <- GSEA(
          genelist,
          TERM2GENE = geneset,
          pvalueCutoff = 0.25,  # GSEA常用阈值
          minGSSize = 10,
          maxGSSize = 500000
        )

        if (!is.null(egmt) && nrow(egmt@result) > 0) {
          # GSEA结果过滤：p.adjust < 0.25
          gsea_data <- egmt@result[egmt@result$p.adjust < 0.25, ]
          if (nrow(gsea_data) == 0) {
            cat("无显著GSEA结果（p.adjust < 0.25），跳过...\n")
          } else {
            gsea_data <- gsea_data[order(gsea_data$NES, decreasing = TRUE), ]
            write.xlsx(gsea_data, file.path(output_dir, paste0("hallmark_gsea", "_", group1, "_vs_", group2, ".xlsx")))

            # 绘制GSEA图
            gsea_data$ID <- factor(gsea_data$ID, levels = gsea_data$ID)
            gsea_data$xlab <- 1:nrow(gsea_data)
            y_min <- floor(min(gsea_data$NES))
            y_max <- ceiling(max(gsea_data$NES))
            top_pathways <- head(as.character(gsea_data$ID), 3)
            bottom_pathways <- tail(as.character(gsea_data$ID), 3)
            label <- c(top_pathways, bottom_pathways)
            data_label <- gsea_data[gsea_data$ID %in% label, ]
            data_label$col <- rep(c("#F6631C", "#2C91E0", "#F3A332", "#018A67"), length.out = nrow(data_label))

            p_gsea <- ggplot(gsea_data, aes(x = xlab, y = NES, color = NES, size = setSize)) +
              geom_point(aes(alpha = -log10(p.adjust)), shape = 21, stroke = 0.7, fill = "#6D65A3", colour = "black") +
              scale_alpha_continuous(range = c(0.1, 0.9), name = "Significance\n(-log10 p-value)") +
              labs(
                title = paste(group1, "vs", group2, "HALLMARK GSEA"),
                x = "Gene Set Rank",
                y = "Normalized Enrichment Score (NES)",
                size = "Gene Set Size",
                color = "NES"
              ) +
              theme_classic(base_size = 15) +
              theme(aspect.ratio = 1) +
              geom_text_repel(
                data = data_label,
                aes(x = xlab, y = NES, label = ID),
                size = 3,
                color = data_label$col,
                max.overlaps = Inf
              )
            print(p_gsea)
            ggsave(
              file.path(output_dir, paste0("HALLMARK_GSEA_", group1, "_vs_", group2, ".pdf")),
              plot = p_gsea, width = 8, height = 6, units = "in", dpi = 600, device = cairo_pdf
            )
          }
        } else {
          cat("无GSEA结果，跳过...\n")
        }
      }, error = function(e) {
        cat("GSEA分析失败:", e$message, "跳过...\n")
      })
    }
  }

  ###########################################
  # 13. 整理结果并返回
  ###########################################
  result <- list(
    pca = pca,
    volcano_plot = volcano_plot,
    diff_genes = diff_genes,
    up_genes = up_genes,
    down_genes = down_genes,
    go_plot = go_plot,
    go_up_plot = go_up_plot,
    go_down_plot = go_down_plot,
    go = go,
    go.up = go.up,
    go.down = go.down,
    kegg_plot = kegg_plot,
    kegg_up_plot = kegg_up_plot,
    kegg_down_plot = kegg_down_plot,
    kegg = kegg,
    kegg.up = kegg.up,
    kegg.down = kegg.down,
    gsea_data = gsea_data,
    p_gsea = p_gsea
  )

  cat("\n", strrep("-", 50), "\n分析完成！结果已保存至指定目录。\n", strrep("-", 50), "\n", sep = "")
  return(result)
}
#


