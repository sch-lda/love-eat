# 1. 清空环境并加载包
rm(list = ls())
library(readxl)
library(vegan)
library(ggplot2)
library(dplyr)

# 2. 读取数据 (使用 file.choose() 避免路径错误)
cat("请选择 侯1.xlsx...\n")
meta_data <- read_excel(file.choose()) 

cat("请选择 侯2.xlsx...\n")
otu_data  <- read_excel(file.choose())

# ==========================================
# 3. 关键修复：查看并修正列名
# ==========================================

# --- 修复 A：处理 OTU 表 (侯2.xlsx) ---
# 打印前几列的名字，让我们看看第一列到底叫什么
print("=== 请查看控制台：OTU表的第一列实际名称 ===")
print(head(colnames(otu_data), 5))

# 无论第一列叫什么（Taxonomy, X_Taxonomy, 或其他），我们都把它取出来作为行名
# 假设第一列是分类学信息，其余列是样本数据
first_col_name <- colnames(otu_data)[1] # 获取第一列的名称
otu_matrix <- otu_data[, -1] # 去掉第一列，只保留数字数据

# 强制转换为矩阵并确保是数字
otu_matrix <- as.matrix(otu_matrix)
storage.mode(otu_matrix) <- "numeric" # 关键：强制把所有内容转为数字 (非数字会变NA)

# 设置行名 (物种名)
rownames(otu_matrix) <- otu_data[[first_col_name]]

# --- 修复 B：处理样本信息表 (侯1.xlsx) ---
print("=== 请查看控制台：样本表的第一列实际名称 ===")
print(head(colnames(meta_data), 5))

# 获取样本名列名
first_col_meta <- colnames(meta_data)[1]
sample_ids <- meta_data[[first_col_meta]]

# 提取分组信息 (假设第二列是 Group)
if(ncol(meta_data) >= 2) {
  group_ids <- meta_data[[2]]
} else {
  # 如果只有1列，就给所有样本一个统一的组名
  group_ids <- rep("Group1", nrow(meta_data))
}

# ==========================================
# 4. 执行分析
# ==========================================

# 4.1 样本对齐
# 确保 OTU 表的列名和样本表的行名一致
common_samples <- intersect(sample_ids, colnames(otu_matrix))

if(length(common_samples) == 0) {
  stop("错误：找不到匹配的样本名。请检查 Excel 文件。")
}

# 重新排序
otu_sorted <- otu_matrix[, common_samples, drop = FALSE]
meta_sorted <- meta_data[match(common_samples, sample_ids), ]

# 4.2 计算距离矩阵 (Bray-Curtis)
cat("正在计算距离...\n")
dist_matrix <- vegdist(t(otu_sorted), method = "bray", na.rm = TRUE)

# 4.3 执行 PCoA
cat("正在执行 PCoA...\n")
pcoa_result <- cmdscale(dist_matrix, k = 2, eig = TRUE)
eig <- pcoa_result$eig
pct <- round(eig / sum(eig[eig > 0]) * 100, 2)

# 提取坐标
coords <- as.data.frame(pcoa_result$points)
colnames(coords) <- c("PCoA1", "PCoA2")
coords$Sample <- rownames(coords)

# 4.4 合并数据绘图
final_df <- coords %>%
  left_join(data.frame(Sample = common_samples, 
                       Group = group_ids[match(common_samples, sample_ids)]), 
            by = "Sample")

# ==========================================
# 5. 绘图
# ==========================================
p <- ggplot(final_df, aes(x = PCoA1, y = PCoA2, color = Group, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.5, size = 3) +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  labs(
    title = "PCoA Plot",
    x = paste0("PCoA 1 (", pct[1], "%)"),
    y = paste0("PCoA 2 (", pct[2], "%)")
  ) +
  theme(plot.title = element_text(hjust = 0.5))

print(p)
ggsave("PCoA_Result.png", plot = p, width = 8, height = 6, dpi = 300)