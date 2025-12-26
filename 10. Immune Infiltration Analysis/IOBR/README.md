# CellChat
IOBR(Immuno-Oncology Biological Research，免疫肿瘤学生物研究)这个R包整合了8种已发表的方法（CIBERSORT、TIMER、xCell、MCPcounter、ESTIMATE、EPIC、IPS和quanTIseq)，用于解析肿瘤微环境的背景信息<br>
## DESCRIPTION
### Version
IOBR V2.0.0
### Learning
* Github :<https://github.com/IOBR/IOBR>
* ZhiHu-1:<https://blog.csdn.net/zfyyzhys/article/details/141048663>
* Zhihu-2:<https://zhuanlan.zhihu.com/p/715185869>
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
set.seed(1234)
# 查看工作路径下的文件
list.files()
dir.create("./CellChat/)
```




