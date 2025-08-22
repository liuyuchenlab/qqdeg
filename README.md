# qqdeg！ 
# 一键获得DEG和enrich!
## 安装

```
  install.packages("devtools")

  library(devtools)  
```

# 如果连接失败：  
## 1.尝试修复Hosts配置  
## 2.尝试Win+R，输入inetcpl.cpl 直接打开Internet选项。打开后，在高级中勾选使用TLS 1.0、使用TLS 1.1、使用TLS 1.2、使用TLS 1.3。

```
  devtools::install_github('liuyuchenlab/qqdeg')  


  library(qqdeg)  

```

# 示例代码
###### 处理小鼠数据或者人类数据，想用其他物种，可以改下代码，非常方便 
###### 数据可以用xlsx | csv | txt文件
###### 可以直接承接TECOUNTS的counts矩阵，选择gene或者te即可提取相应矩阵
###### 其他软件的gene定量选gene即可
###### 最重要的是数据格式，基因的列名必须为gene_id，样本名可以为A-1，A_1，A.1三种形式之一

![image](https://github.com/user-attachments/assets/4499d333-b5a1-4bf3-8051-7435f5d0cf97)

###### 第一个是实验组，第二个是对照组，fc是foldchange，默认是1.5
###### GO和KEGG默认p.adj < 0.05，TOP 5，删除重复项

```
result <- qqdeg("mouse.xlsx","gene","male-ko","male-wt",fc_threshold = 1.5,species = "mouse")
```

## result中存有大量对象，可自行保存或者修改

### 输出一张主成分分析图

<img width="936" height="515" alt="image" src="https://github.com/user-attachments/assets/4013145c-9145-4961-ba05-9837fcc0bfa6" />

### 输出一张火山图

<img width="936" height="515" alt="image" src="https://github.com/user-attachments/assets/64ae0365-6315-44b7-97f7-0cd1dd0694e8" />

###### 可使用输出的DEG文件自行进行富集分析

### GO BP 富集

<img width="936" height="515" alt="image" src="https://github.com/user-attachments/assets/54e9ace1-7387-4952-96a2-6f6e9a2402f9" />

### KEGG 富集

<img width="936" height="515" alt="image" src="https://github.com/user-attachments/assets/920d42e5-94db-4375-9785-9e58d12f39a0" />

### HALLMARK GSEA 打分排序

<img width="953" height="712" alt="image" src="https://github.com/user-attachments/assets/c2894c34-f89c-453a-8d1f-471551c09045" />

### 运行函数会自动新建文件夹存放文件

<img width="772" height="565" alt="image" src="https://github.com/user-attachments/assets/872ee939-0638-4ca5-9d1a-9730c5d90c3c" />















