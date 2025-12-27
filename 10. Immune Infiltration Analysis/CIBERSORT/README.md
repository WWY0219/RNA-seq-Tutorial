# CIBERSORT
**CIBERSORT**，The goal of CIBERSORT is to run the CIBERSORT flow.<br>
## DESCRIPTION
### Version
CIBERSORT
### Learning
* Github :<https://github.com/Moonerss/CIBERSORT>
* ZhiHu-1:<https://zhuanlan.zhihu.com/p/555473921>
## Usage
### Part Ⅰ: 加载包并准备环境
```R
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("workspace")
getwd()
library(CIBERSORT)
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
### Part Ⅱ: 加载数据
数据格式为第一列为基因名，不能有空值和重复，重复基因取平均值或最大值处理，第一行为样本名称，数据为data.frame。
```R
data <- read.table(file="tpm.txt",header=TRUE, sep="\t",check.names=F,quote="")
head(data,n=3)
```
#### 读取包自带的LM22文件（免疫细胞特征基因文件）
```R
sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
```
### Part Ⅲ: 运行cibersot
```R
results <- cibersort(sig_matrix, "tpm.txt",perm = 0,QN = T)
```
> perm置换次数=1000，QN分位数归一化=TRUE<br>

#### 得到的结果是第一列为样本，第一行为细胞类群
```R
head(results[,1:4],n=12)
```

### Part Ⅳ: 结果可视化
#### 箱线图
数据整理
```R
res <- data.frame(results[,1:22])%>%
  mutate(group = c(rep('shControl',4),rep('shPHF8',8)))%>%
  rownames_to_column("sample")%>%
  pivot_longer(cols = colnames(.)[2:23],
               names_to = "cell.type",
               values_to = 'value')
head(res,n=6)
```
绘图
```R
ggplot(res,aes(cell.type,value,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 14,face = "italic",colour = 'black'),
        axis.text.y = element_text(face = "italic",size = 14,colour = 'black'))+
  scale_fill_nejm()+
  stat_compare_means(aes(group = group),label = "p.format",size=3,method = "kruskal.test")
```
#### 热图
```R
normalize <- function(x) {
  if((max(x) - min(x)) == 0){
    return(mean(x))
  }else{
    return((x - min(x)) / (max(x) - min(x)))
  }
}
heat_map_res <- apply(results[,1:22], 2,normalize)
heat_map_res <- t(as.data.frame(heat_map_res))
mycol<-colorRampPalette(c("navy", "white", "firebrick3"))(100)
library(pheatmap)
pheatmap(heat_map_res,color = mycol,cluster_rows = F,cluster_cols = F,scale = 'none')
```

















