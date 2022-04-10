# crystal_identification
从点阵数据中分析金属晶体的结构

## 缩写:
BAA: Bond angle analysis
DBSCAN: Density-based spatial clustering of applications with noise
BAAm: Bond angle anaysis(modification)
CSP: Center symmetry parameter

## 备注:  
1.BBAm.m文件里面使用的pdist函数, 需要在matlab的附加功能里面安装Statistics and Machine Learning Toolbox  
  
## .m文件介绍：  
main.m: 主要的脚本程序,所有的函数都在这里执行  
BAA.m: 将BAA方法的封装函数,识别晶体结构,[输入:原子位置,输出:结构类型]  
BAAm.m: 对BAA.m进行改进并注释的版本  
DBSCAN.m: 对DBSCAN方法的封装函数,识别单晶的方向  
Neg.m: 寻找邻居原子的函数  
Orien.m: 计算晶胞的晶向向量的函数  
p_norm.m: 求空间中点阵的两两点之间的距离  

## .mat文件介绍:  
RefDataBall.mat: 球状晶体数据  
RefDataCell.mat: 晶体原胞数据  
RefDataSquare.mat: 方形晶体数据  
  
## RawData是从lammps导出的原始数据  