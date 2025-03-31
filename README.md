# qqdeg！ 
###### 一键获得DEG和enrich!
##### 安装

```
  install.packages("devtools")

  library(devtools)  
```

###### 如果连接失败：  
###### 1.尝试修复Hosts配置  
###### 2.尝试Win+R，输入inetcpl.cpl 直接打开Internet选项。打开后，在高级中勾选使用TLS 1.0、使用TLS 1.1、使用TLS 1.2、使用TLS 1.3。

```
  devtools::install_github('liuyuchenlab/qqdeg')  


  library(qqdeg)  

```
###### 数据格式只能用xlsx文件
###### 可以直接承接TECOUNTS的counts矩阵，选择gene或者te即可
###### 其他软件的gene定量选gene即可
###### 最重要的是数据格式，基因的列名为gene_id，样本名为A-1，A-2，B-1，B-2等等

![image](https://github.com/user-attachments/assets/4499d333-b5a1-4bf3-8051-7435f5d0cf97)




###### 示例代码，自己用的，处理小鼠数据，想用其他物种，或者改其他参数可以改下代码，非常方便 

###### 第一个是实验组，第二个是对照组，fc是foldchange，默认是1.5
###### GO和KEGG默认p.adj < 0.05，TOP 5，删除重复项


```
  result <- qqdeg("rlim.xlsx","gene","male-ko","male-wt",fc_threshold = 1.5)
```

###### result中存有大量对象，可自行保存或者修改

###### 输出一张火山图
![image](https://github.com/user-attachments/assets/eeb526cb-933f-4b71-afa7-a81a553e9546)



###### 可使用输出的DEG文件自行进行富集分析，如果没出问题那就最好了

##### GO BP 富集

![image](https://github.com/user-attachments/assets/a661df4f-8841-4aba-b783-b73e72909eb0)


##### KEGG 富集

![image](https://github.com/user-attachments/assets/8e1e5d48-f04b-4e57-99ba-59bf48f0e085)


##### HALLMARK GSEA 富集

![image](https://github.com/user-attachments/assets/0f8659c8-a025-4f58-b003-f4ba096fb48c)


##### 运行函数会自动新建文件夹存放文件

![image](https://github.com/user-attachments/assets/5eef6abb-2f9d-4e12-b90a-91c852324645)














