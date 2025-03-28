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
###### 数据格式
###### 可以直接承接TECOUNTS的counts矩阵，选择gene或者te即可，其他软件的gene定量选gene即可

![image](https://github.com/user-attachments/assets/7f48fb48-c7e4-4c7d-94b6-edd3760a83aa)



###### 示例代码，自己用的，只能用xlsx文件，处理小鼠数据，想用其他物种，或者改其他参数可以改下代码，非常方便 。

###### 第一个是实验组，第二个是对照组，fc是foldchange，默认是1.5，result里存有很多对象
###### GO和KEGG默认p.adj < 0.05，TOP 5,删除重复项。


```
  result <- qqdeg("rlim.xlsx","gene","female-ko","female-wt",fc_threshold = 1.5)
```

###### result中存有大量对象，可自行保存或者修改

###### 输出一张火山图
![image](https://github.com/user-attachments/assets/8443dde9-1d17-47a0-bcff-d22ea42c6b49)


###### 可使用输出的DEG文件自行进行富集分析，如果没出问题那就最好了

###### GO BP 富集
![image](https://github.com/user-attachments/assets/612d569c-f38d-4bd3-902d-a6b0ed72022a)

###### KEGG 富集
![image](https://github.com/user-attachments/assets/246a62b7-0951-4246-88a9-e8083a37bc60)

###### 运行函数会自动新建文件夹存放文件
![image](https://github.com/user-attachments/assets/fa96f77d-620d-4fe6-8693-72fb5d5bb3ca)












