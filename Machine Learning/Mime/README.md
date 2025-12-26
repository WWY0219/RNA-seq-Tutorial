# Mime
机器学习相关.<br>
## DESCRIPTION
### Version
Mime1
### Learning
* Github :
* ZhiHu-1:[<https://zhuanlan.zhihu.com/p/1894789522887250489>](https://jishuzhan.net/article/1954880570147778561)
* Zhihu-2:
## Usage
### 01. 加载包并准备环境
```R
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("workspace")
getwd()
library(Mime1)
library(patchwork)
library(ggplot2)
library(ComplexHeatmap)
library(svglite)
library(qs)
library(dplyr)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(harmony)
library(patchwork)
library(RColorBrewer)
set.seed(1234)
# 查看工作路径下的文件
list.files()
dir.create("./Mime1/)
```
### 02. 数据整理


### 03. 构建预后预测模型
res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$Dataset1,
                     list_train_vali_Data = list_train_vali_Data,
                     unicox.filter.for.candi = T,
                     unicox_p_cutoff = 0.05,
                     candidate_genes = genelist,
                     mode = 'all',nodesize =5,seed = 5201314 )
> `ML.Dev.Prog.Sig()`提供了三种模式，包括all、single和double。
>    - ALL:标识使用所有十种算法及其组合。
>    - Single:表示只使用十种算法中的一种。
>    - double:表示使用两种算法的组合。
>    - 在大多数情况下，我们通常会使用all模式来分析数据。<br>
> 如果将unicox.filter.for.candi设置为T (默认值),Mime首先会在训练数据集中对提供的基因进行单变量cox回归，以筛选出预后变量，然后用于构建模型。
