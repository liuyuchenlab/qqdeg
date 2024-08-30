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

###### 示例代码

result <- qqdeg("rlim.xlsx","gene","female-ko","female-wt",fc_threshold = 1.5)






