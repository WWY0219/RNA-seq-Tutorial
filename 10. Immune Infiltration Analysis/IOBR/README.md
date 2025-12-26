# IOBR
**IOBR(Immuno-Oncology Biological Research**，免疫肿瘤学生物研究)这个R包整合了8种已发表的方法（CIBERSORT、TIMER、xCell、MCPcounter、ESTIMATE、EPIC、IPS和quanTIseq)，用于解析肿瘤微环境的背景信息<br>
开发者还收集了264个已发表的基因集签名，涵盖了肿瘤微环境、肿瘤代谢、m6A、外泌体、微卫星不稳定性和三级淋巴结构。<br>
通过运行`signature_collection_citation`函数，可以获取这些基因集签名的来源文献，而`signature_collection`函数则返回所有给定签名的详细基因。<br>
IOBR使用了三种计算方法来计算签名得分，包括**PCA**、**z-score**和**ssGSEA**。值得注意的是，IOBR还收集并使用了多种方法进行变量转换、可视化、批量生存分析、特征选择和统计分析。该包支持批量分析和结果可视化
## DESCRIPTION
### Version
IOBR V2.0.0
### Learning
* Github :<https://github.com/IOBR/IOBR>
* ZhiHu-1:<https://blog.csdn.net/zfyyzhys/article/details/141048663>
* Zhihu-2:<https://zhuanlan.zhihu.com/p/715185869>
* Weichat:<https://mp.weixin.qq.com/s/nMeDmyO4M09z7vGnYU2K-A>
## Usage
### 01. 加载包并准备环境
```R
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("workspace")
getwd()
library(IOBR)
library(patchwork)
library(ComplexHeatmap)
library(svglite)
library(qs)
library(dplyr)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(gplots)
set.seed(1234)
# 查看工作路径下的文件
list.files()
dir.create("./CellChat/)
```
加载表达数据
```R
load("~/Desktop/data/combat.Rdata")
head(exp)[1:4,1:4]
```
临床数据
```R
head(pd)
```
确认临床信息与表达数据是否符合
```R
identical(rownames(pd),colnames(exp))
```
内置免疫分析方法
```R
tme_deconvolution_methods
```
> mcpcounter，epic，xcell，cibersort，cibersort_abs，ips，ESTIMATE，SVR，lsei，TIMER，quanTIseq<br>

### Part I: cibersort
```R
ciber_res <- deconvo_tme(eset = exp, 
                         method = "cibersort", 
                         arrays = TRUE, 
                         perm = 100 )
ciber_res
ciber_plot <- cell_bar_plot(input = ciber_res, 
                            features = colnames(ciber_res)[2:23], 
                            title = "CIBERSORT Cell Fraction")
```
### Part II: epic
```R
epic_res <- deconvo_tme(eset = exp,
                        method = "epic",
                        arrays = TRUE)
epic_res
epic_plot <- cell_bar_plot(input = epic_res,
                           features = colnames(epic_res)[2:9], 
                           title = "EPIC Cell Fraction")
```
### Part III: quantiseq
```R
quantiseq_res <- deconvo_tme(eset = exp, tumor = F, 
                             arrays = TRUE, scale_mrna = TRUE, method = "quantiseq")
quantiseq_res
quantiseq_plot <- cell_bar_plot(input = quantiseq_res, features = colnames(quantiseq_res)[2:12], 
                            title = "Quantiseq Cell Fraction")
```
### Part IV: MCPcounter
```R
mcpcounter_res <- deconvo_tme(eset = exp, method = "mcpcounter")
mcpcounter_res
```

### Part V: xCELL
```R
xCELL_res <- deconvo_tme(eset = exp, method = "xcell",arrays = TRUE)
```

### Part V: TIMER
```R
TIMER_res <- deconvo_tme(eset = exp, method = "timer",group_list = group_list)
```
> 这里导入的group_list必须是TIMER能够识别的，比如TCGA中33肿瘤类型<br>


### Part VI: estimate
```R
estimate_res <- deconvo_tme(eset = exp, method = "estimate")
estimate_res
```
### Part VI: IPS
```R
ips_res <- deconvo_tme(eset = exprSet, method = "ips", plot= F)
head(ips_res)
```
> 值得一提的是这个IPS，其是指immunophenotype (免疫表型评分)。 里面一共评估四个主要参数分别是：MHC分子(MHC molecules，MHC) ，免疫调节分子(Immunomodulators，CP) ，效应细胞(Effector cells，EC) ，抑制细胞(Suppressor cells，SC) 。使用者可以根据得到的结果联合生存数据/临床参数分析<br>

### pheatmap热图
#### 整合数据-挑选了三个绝对丰度的结果
```R
data_total <- cbind(ciber_res[,-c(24:26)],epic_res[,-1],quantiseq_res[,-1])
data_total <- as.data.frame(t(data_total))
colnames(data_total) <- data_total[1,]
data_total <- dplyr::slice(data_total,-1)
head(data_total)[1:5,1:5]
```
#### 检查数据类型是否正确转换为数值型
```R
data_total <- data_total %>% mutate_all(as.numeric)
```
#### 数据标准化
* 方法1
```R
stand_fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
                  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
                  if (!is.null(halfwidth)) {
                     outdata[outdata > halfwidth] = halfwidth
                     outdata[outdata < (-halfwidth)] = -halfwidth
  }
  return(outdata)
}
data_total <- stand_fun(data_total,halfwidth = 2)
```
* 方法2
```R
data_total <- scale(data_total)
```
#### 重新对列进行排序
> group_list 是一个向量，表示每个样本的类型(PAH 或 control)<br>
```R
sorted_index <- order(group_list, decreasing = FALSE)
```
#### 使用排序索引对 data_total 重新排序列
```R
data_total <- data_total[, sorted_index]
```
#### 重新对 group_list 进行排序，以匹配列的顺序
```R
group_list <- group_list[sorted_index]
```
#### 创建注释
* 列注释
```
annCol <- data.frame(Type = group_list,
                     row.names = colnames(data_total),
                     stringsAsFactors = FALSE)
```
* 行注释
从行名中提取方法
```R
methods <- sub('.*_', '', rownames(data_total))
```
* 创建行注释数据框

annRow <- data.frame(Methods = factor(methods, levels = unique(methods)),
                     row.names = rownames(data_total),
                     stringsAsFactors = FALSE)


* 调整颜色梯度
```R
breaksList = seq(-5, 5, by = 0.1)
colors <- colorRampPalette(c("#336699", "white", "tomato"))(length(breaksList))
```
#### 绘制热图
```R
pheatmap(data_total,
         annotation_col = annCol,
         annotation_row = annRow,
         color = colors,
         breaks = breaksList,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         gaps_col = cumsum(table(annCol$Type)),  # 使用排序后的列分割点
         gaps_row = cumsum(table(annRow$Methods)), # 行分割
         fontsize_row = 6,
         fontsize_col = 6,
         annotation_names_row = FALSE
         )
```
### complexheatmap
#### 整合数据-挑选了三个绝对丰度的结果
```R
data_total <- cbind(ciber_res[,-c(24:26)],epic_res[,-1],quantiseq_res[,-1])
data_total <- as.data.frame(t(data_total))
colnames(data_total) <- data_total[1,]
data_total <- dplyr::slice(data_total,-1)
head(data_total)[1:5,1:5]
```
#### 检查数据类型是否正确转换为数值型
```R
data_total <- data_total %>% mutate_all(as.numeric)
```
* 标准化-方法3(范围限定在-3至3)
```R 
data_total <- scale(data_total)
data_total[data_total > 3] <- 3
data_total[data_total + 3 < 0] <- -3
```

#### 定义颜色渐变(也可以用circlize创建颜色映射)
```R
library(circlize)
color_fun <- colorRamp2(c(-3, 0, 3), c("#336699", "white", "tomato"))
```
#### pd$type 和 pd$samples 是注释数据怎么上色
```R
type_colors <- c("control" = "green", "PAH" = "red")
samples_colors <- c("all PAH" = "blue", "CTEPH patient" = "orange", "FD" = "purple",
                    "idiopathic PAH patient" = "brown", "normal control" = "yellow",
                    "patient with PAH and CHD" = "pink", "patient with PAH and CTD" = "magenta",
                    "pulmonary arterial hypertension (PAH) patient" = "cyan")
```
```R
columnAnno <- HeatmapAnnotation(type = pd$group,
                                samples = pd$samples,
                                col = list(type = type_colors, samples = samples_colors))
```
#### 绘制热图
```R
ComplexHeatmap::Heatmap(data_total, 
                        na_col = "white",
                        col = color_fun,  # 添加颜色映射函数
                        show_column_names = F,
                        row_names_side = "right",
                        name = "fraction",
                        column_order = c(colnames(data_total)[c(grep("control",pd$group),
                                                                grep("PAH",pd$group))]),
                        column_split = pd$group, 
                        column_title = NULL,
                        cluster_columns = F,
                        top_annotation = columnAnno,
                        heatmap_width = unit(20, "cm"),  # 调整热图宽度
                        row_dend_width = unit(1, "cm"),  # 调整聚类树宽度
                        cluster_rows = T, # 行聚类
                        row_names_gp = gpar(fontsize = 8)  # 调整行名字体大小
                        #cluster_columns = FALSE # 列聚类。                     
                        )
```





