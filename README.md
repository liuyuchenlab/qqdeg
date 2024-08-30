# qqdeg！ 
###### 差异分析重复做太麻烦，可以一键获得DEG和enrich
##### 安装

install.packages("devtools")

library(devtools)  

###### 如果连接失败：  
###### 1.尝试修复Hosts配置  
###### 2.尝试Win+R，输入inetcpl.cpl 直接打开Internet选项。打开后，在高级中勾选使用TLS 1.0、使用TLS 1.1、使用TLS 1.2、使用TLS 1.3。

devtools::install_github('liuyuchenlab/qqdeg')  


library(qqdeg)  

###### 示例代码，自己用的，只能用xlsx文件，处理小鼠数据，可处理gene数据和TE数据，第一个是实验组，第二个是对照组，fc是foldchang，默认是1.5倍，result里存有很多对象。

result <- qqdeg("rlim.xlsx","gene","female-ko","female-wt",fc_threshold = 1.5)

###### 输出一张火山图
![image](https://github.com/user-attachments/assets/8443dde9-1d17-47a0-bcff-d22ea42c6b49)

###### GO BP 富集
![image](https://github.com/user-attachments/assets/612d569c-f38d-4bd3-902d-a6b0ed72022a)

###### KEGG 富集
![image](https://github.com/user-attachments/assets/246a62b7-0951-4246-88a9-e8083a37bc60)

###### 还会新建文件夹存放文件
![image](https://github.com/user-attachments/assets/4e87eed0-7c75-4008-92e2-0123e2cafec6)











