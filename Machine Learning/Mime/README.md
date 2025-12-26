# Mime
机器学习相关.<br>
## DESCRIPTION
### Version
Mime1
### Learning
* Github :<https://github.com/l-magnificence/Mime>
* ZhiHu-1:<https://jishuzhan.net/article/1954880570147778561>
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
#### 选择最优模型
```R
res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$Dataset1,
                       list_train_vali_Data = list_train_vali_Data,
                       unicox.filter.for.candi = T,
                       unicox_p_cutoff = 0.05,
                       candidate_genes = genelist,
                       mode = 'all',
                       nodesize =5,
                       seed = 5201314 )
```
> `ML.Dev.Prog.Sig()`提供了三种模式，包括all、single和double。
>    - ALL:标识使用所有十种算法及其组合。
>    - Single:表示只使用十种算法中的一种。
>    - double:表示使用两种算法的组合。
>    - 在大多数情况下，我们通常会使用all模式来分析数据。<br>

> 如果将unicox.filter.for.candi设置为T (默认值),Mime首先会在训练数据集中对提供的基因进行单变量cox回归，以筛选出预后变量，然后用于构建模型。
#### 绘制每个模型的C-index
```R
cindex_dis_all(res,validate_set = names(list_train_vali_Data)[-1],
order =names(list_train_vali_Data),width = 0.35)
```
#### 在不同数据集之间绘制特定模型的C-index
```R
cindex_dis_select(res,
                  model="StepCox[forward] + plsRcox",
                  order= names(list_train_vali_Data))
```
> 如果输入对象res来自ML.Dev.Prog.Sig()中使用的所有模型，则应将model定义为特定的模型名称，同时将model定义为SOD<br>

#### 根据特定模型在不同数据集中计算的风险评分绘制患者生存曲线
```R
survplot <- vector("list",2) 
for (i in c(1:2)) {
  print(survplot[[i]]<-rs_sur(res, model_name = "StepCox[forward] + plsRcox",dataset = names(list_train_vali_Data)[i],
                              #color=c("blue","green"),
                              median.line = "hv",
                              cutoff = 0.5,
                              conf.int = T,
                              xlab="Day",pval.coord=c(1000,0.9)))
}
aplot::plot_list(gglist=survplot,ncol=2)
```
#### 计算每个模型的AUC scores
```R
all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Dataset1"]],
                            inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 1,
                            auc_cal_method="KM")
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Dataset1"]],
                            inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 3,
                            auc_cal_method="KM")
all.auc.5y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Dataset1"]],
                            inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 5,
                            auc_cal_method="KM")
```
> Cal_AUC_ml_res()提供了三种模式，这些模式应该与 ML.Dev.Prog.Sig () 使用的模式一致。<br>
> AUC_time for 1 year, 2 years, 3 years..., We recommend using the shortest survival time among all datasets.<br>

* only plot 1-year AUC predicted by all models:
```R
auc_dis_all(all.auc.1y,

        dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.35,
            year=1)
```
#### 绘制不同数据集中特定模型的ROC图
```R
roc_vis(all.auc.1y,
        model_name = "StepCox[forward] + plsRcox",
        dataset = names(list_train_vali_Data),
        order= names(list_train_vali_Data),
        anno_position=c(0.65,0.55),
        year=1)
```
* Plot 1, 3, and 5-year AUC of specific model among different datasets:
```R
auc_dis_select(list(all.auc.1y,all.auc.3y,all.auc.5y),
               model_name="StepCox[forward] + plsRcox",
               dataset = names(list_train_vali_Data),
               order= names(list_train_vali_Data),
               year=c(1,3,5))
```
#### 特定模型的单变量 Cox 回归元分析
```R
unicox.rs.res <- cal_unicox_ml_res(res.by.ML.Dev.Prog.Sig = res,optimal.model = "StepCox[forward] + plsRcox",type ='categorical')
metamodel <- cal_unicox_meta_ml_res(input = unicox.rs.res)
meta_unicox_vis(metamodel,
                dataset = names(list_train_vali_Data))
```
> type includes categorical and continuous.

#### 与先前发表的模型进行比较
```R
rs.glioma.lgg.gbm <- cal_RS_pre.prog.sig(use_your_own_collected_sig = F,type.sig = c('LGG','GBM','Glioma'),
                                        list_input_data = list_train_vali_Data)
```
> cal_RS_pre.prog.sig()will calculate the risk score based on the signatures from previous papers.

* 特定模型的HR与之前发表的模型进行比较
```R
HR_com(rs.glioma.lgg.gbm,
       res,
       model_name="StepCox[forward] + plsRcox",
       dataset=names(list_train_vali_Data),
       type = "categorical")
```
```R
cc.glioma.lgg.gbm <- cal_cindex_pre.prog.sig(use_your_own_collected_sig = F,type.sig = c('Glioma','LGG','GBM'),
                                            list_input_data = list_train_vali_Data)
```
> cal_cindex_pre.prog.sig()will calculate the C-index based on the signatures from previous papers like the fuction cal_RS_pre.prog.sig()

* 将特定模型的C-index与之前发表的模型进行比较
```R
cindex_comp(cc.glioma.lgg.gbm,
            res,
            model_name="StepCox[forward] + plsRcox",
            dataset=names(list_train_vali_Data))
```
```R
auc.glioma.lgg.gbm.1 <- cal_auc_pre.prog.sig(use_your_own_collected_sig = F,
                                            type.sig = c('Glioma','LGG','GBM'),
                                            list_input_data = list_train_vali_Data,AUC_time = 1,
                                            auc_cal_method = 'KM')
```
> Call_auc_pre.prog.sig()将基于以前论文的签名计算 AUC, 例如 cal_RS_pre.prog.sig()函数。<br>
> AUC_time类似于cal_AUC_ml_res()的要求<br>

#### 与特定模型的AUC与之前发布的模型进行比较
```R
auc_comp(auc.glioma.lgg.gbm.1,
         all.auc.1y,
         model_name="StepCox[forward] + plsRcox",
         dataset=names(list_train_vali_Data))
```

#### 免疫浸润分析
```R
devo <- TME_deconvolution_all(list_train_vali_Data)
```
> TME_deconvolution_all() includes 10 deconvolution methods ("quantiseq", "xcell", "epic", "abis", "mcp_counter", "estimate", "cibersort", "cibersort_abs", "timer", "consensus_tme") from immunedeconv::deconvolution_methods. By default, deconvolution methods are set as ("xcell", "epic", "abis", "estimate", "cibersort", "cibersort_abs").

```R
immuno_heatmap(res,
               devo,
               model_name="StepCox[backward] + plsRcox",
               dataset="Dataset1")
```

### 04. 构建响应的预测模型
```R
load("./Example.ici.Rdata")
load("./genelist.Rdata")
res.ici <- ML.Dev.Pred.Category.Sig(train_data = list_train_vali_Data$training,
                                      list_train_vali_Data = list_train_vali_Data,
                                      candidate_genes = genelist,
                                      methods = c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass'),
                                      seed = 5201314,
                                      cores_for_parallel = 60
)
```
> ML.Dev.Pred.Category.Sig()使用机器学习算法开发二进制变量的预测模型<br>

#### 绘制不同数据集中不同方法的AUC:
```R
auc_vis_category_all(res.ici,dataset = c("training","validation"),
  order= c("training","validation"))
```
#### 绘制特定方法在不同数据集之间的ROC
```R
plot_list<-list()
methods <- c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass')
for (i in methods) {
  plot_list[[i]]<-roc_vis_category(res.ici,model_name = i,dataset = c("training","validation"),
                                   order= c("training","validation"),
                                   anno_position=c(0.4,0.25))
}
aplot::plot_list(gglist=plot_list,ncol=3)
```
#### 将AUC与其他已发表的与免疫治疗反应相关的模型进行比较：
```R
auc.other.pre <- cal_auc_previous_sig(list_train_vali_Data = list_train_vali_Data,seed = 5201314,
                                      train_data = list_train_vali_Data$training,
                                      cores_for_parallel = 32)
```
> cal_auc_previous_sig(): 将根据先前论文的签名计算免疫治疗反应的AUC。<br>
> cores_for_parallel 是指您可以选择用于并行操作的核心。如果多核心条件是错误的，请将cores_for_parallel设置为 1。<br>

#### 绘制特定模型的比较结果图
```R
auc_category_comp(res.ici,
                  auc.other.pre,
                  model_name="svmRadialWeights",
                  dataset=names(list_train_vali_Data))
```

### 05. 核心特征选择
```R
load("./Example.cohort.Rdata")
load("./genelist.Rdata")
res.feature.all <- ML.Corefeature.Prog.Screen(InputMatrix = list_train_vali_Data$Dataset1,
                                            candidate_genes = genelist,
                                            mode = "all",nodesize =5,seed = 5201314 )
```
> ML.Corefeature.Prog.Screen()提供了三种模式，包括all、single和all_without_SVM。
>   - all模式意味着使用所有八种方法 ("RSF"、"Enet"、"Boruta"、"Xgboost"、"SVM-REF"、"Lasso"、"CoxBoost"、"StepCox") 进行选择。
>   - single模式意味着只使用一种方法进行运行。如果使用single模式，则应在八种方法中指定single_ml。
>   - 由于SVM耗时太长，我们定义了其他七种用于选择的方法作为all_without_SVM模式。

#### 通过不同方法筛选的基因的Upset plot
```R
core_feature_select(res.feature.all)
```
#### 绘制不同方法过滤的基因序列
```R
core_feature_rank(res.feature.all, top=20)
```
#### 随机选择前两个基因来分析它们的相关性：
```R
dataset_col<-c("#3182BDFF","#E6550DFF")
corplot <- list()
for (i in c(1:2)) {
  print(corplot[[i]]<-cor_plot(list_train_vali_Data[[i]],
                               dataset=names(list_train_vali_Data)[i],
                               color = dataset_col[i],
                               feature1="PSEN2",
                               feature2="WNT5B",
                               method="pearson"))
}
aplot::plot_list(gglist=corplot,ncol=2)
```
#### 根据不同数据集中特定基因的中位表达水平绘制患者生存曲线：
```R
survplot <- vector("list",2) 
for (i in c(1:2)) {
  print(survplot[[i]]<-core_feature_sur("PSEN2", 
                                        InputMatrix=list_train_vali_Data[[i]],
                                        dataset = names(list_train_vali_Data)[i],
                                        #color=c("blue","green"),
                                        median.line = "hv",
                                        cutoff = 0.5,
                                        conf.int = T,
                                        xlab="Day",pval.coord=c(1000,0.9)))
}
aplot::plot_list(gglist=survplot,ncol=2)
```















































































