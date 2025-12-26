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
```R
load("~/Desktop/data/combat.Rdata")
head(exp)[1:4,1:4]
```
数据结构belike
> #           GSM3106326 GSM3106327 GSM3106328 GSM3106329<br>
> # LINC01128   8.000049   8.044479   7.311868   7.610713<br>
> # SAMD11      7.499611   7.562370   7.412649   7.239141<br>
> # KLHL17      7.662076   7.672331   7.627330   7.564212<br>
> # PLEKHN1     7.666884   7.889537   7.640606   7.746440<br>
临床数据
```R
head(pd)
```
> #                                                title                  samples   GSE_num group<br>
> # GSM3106326 Pulmonary arterial hypertension patient 1   idiopathic PAH patient GSE113439   PAH<br>
> # GSM3106327 Pulmonary arterial hypertension patient 2 patient with PAH and CHD GSE113439   PAH<br>
> # GSM3106328 Pulmonary arterial hypertension patient 3 patient with PAH and CTD GSE113439   PAH<br>
> # GSM3106329 Pulmonary arterial hypertension patient 4 patient with PAH and CHD GSE113439   PAH<br>
> # GSM3106330 Pulmonary arterial hypertension patient 5   idiopathic PAH patient GSE113439   PAH<br>
> # GSM3106331 Pulmonary arterial hypertension patient 6   idiopathic PAH patient GSE113439   PAH<br>
确认临床信息与表达数据是否符合
```R
identical(rownames(pd),colnames(exp))
```
内置免疫分析方法
```R
tme_deconvolution_methods
```
> mcpcounter，epic，xcell，cibersort，cibersort_abs，ips，ESTIMATE，SVR，lsei，TIMER，quanTIseq<br>

### Part I：cibersort
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


