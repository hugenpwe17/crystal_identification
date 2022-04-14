# crystal_identification  
从点阵数据中分析金属晶体的结构  

## 备注  
BBAm.m文件里面使用的pdist函数, 需要在matlab的附加功能里面安装Statistics and Machine Learning Toolbox  

## 缩写
BAA: Bond angle analysis  
DBSCAN: Density-based spatial clustering of applications with noise  
BAAm: Bond angle anaysis(modification)  
CSP: Center symmetry parameter  
  
### script info
main.m: 主要的脚本程序,所有的函数都在这里执行  
my_script_1.m: 计算在使用rand模拟的热扰动下的键角余弦的分布情况

### function info
BAA.m: 将BAA方法的封装函数,识别晶体结构,[输入:原子位置,输出:结构类型]  
BAAm.m: 对BAA.m进行改进并注释的版本  
BAAm_2.m: 使用向量化参数改进的BBA方法
BAAm_2_orth.m: 对BBAm_2.m中的参数向量进行正交化改进
DBSCAN.m: 对DBSCAN方法的封装函数,识别单晶的方向  
Neg.m: 寻找邻居原子的函数  
Orien.m: 计算晶胞的晶向向量的函数  
p_norm.m: 求空间中点阵的两两点之间的距离  
myorth.m: 对三个向量进行施密特正交化处理

## .mat info
RefDataBall.mat: 球状晶体数据  
RefDataCell.mat: 晶体原胞数据  
RefDataSquare.mat: 方形晶体数据  
20210416: bcc/fcc/hcp的参数向量对/bcc/fcc/hcp结构在不同形变下的识别效果
  
## .fig info
RefBCC.fig: BCC单胞在无扰动情况下的键角余弦分布情况
RefFCC.fig: FCC单胞在无扰动情况下的键角余弦分布情况
RefHCP.fig: HCP单胞在无扰动情况下的键角余弦分布情况
20210606: 理想模型在算法下的边缘检测
20210607: BBA和BBAm在不同形变下的识别效果


## RawData是从lammps导出的原始数据  

./
├── BAA.m
├── BAAm_2.m
├── BAAm_2_orth.m
├── BAAm.m
├── DBSCAN.m
├── figure
│   ├── 20210606
│   │   ├── figure_edge_bcc.eps
│   │   ├── figure_edge_bcc.jpg
│   │   ├── figure_edge_combine.eps
│   │   ├── figure_edge_combine.jpg
│   │   ├── figure_edge_fcc.eps
│   │   ├── figure_edge_fcc.jpg
│   │   ├── figure_edge_hcp.eps
│   │   └── figure_edge_hcp.jpg
│   ├── 20210607
│   │   ├── figure_delta0.03bcc_modify.eps
│   │   ├── figure_delta0.03bcc_modify.jpg
│   │   ├── figure_delta0.03bcc_origin.eps
│   │   ├── figure_delta0.03bcc_origin.jpg
│   │   ├── figure_delta0.06bcc_modify.eps
│   │   ├── figure_delta0.06bcc_modify.jpg
│   │   ├── figure_delta0.06bcc_origin.eps
│   │   ├── figure_delta0.06bcc_origin.jpg
│   │   ├── figure_delta0.09bcc_modify.eps
│   │   ├── figure_delta0.09bcc_modify.jpg
│   │   ├── figure_delta0.09bcc_origin.eps
│   │   ├── figure_delta0.09bcc_origin.jpg
│   │   ├── hist_delta0.03bcc_modify.eps
│   │   ├── hist_delta0.03bcc_modify.jpg
│   │   ├── hist_delta0.03bcc_origin.eps
│   │   ├── hist_delta0.03bcc_origin.jpg
│   │   ├── hist_delta0.06bcc_modify.eps
│   │   ├── hist_delta0.06bcc_modify.jpg
│   │   ├── hist_delta0.06bcc_origin.eps
│   │   ├── hist_delta0.06bcc_origin.jpg
│   │   ├── hist_delta0.09bcc_modify.eps
│   │   ├── hist_delta0.09bcc_modify.jpg
│   │   ├── hist_delta0.09bcc_origin.eps
│   │   ├── hist_delta0.09bcc_origin.png
│   │   └── Untitled.ipynb
│   ├── RefBCC.fig
│   ├── RefFCC.fig
│   └── RefHCP.fig
├── list
├── main.m
├── Mat
│   ├── 20210410
│   │   ├── BAA_delta_bcc_in_bcc.mat
│   │   ├── BAA_delta_bcc_in_cp.mat
│   │   ├── BAA_delta_bcc_in_fcc.mat
│   │   ├── BAA_delta_bcc_in_hcp.mat
│   │   ├── BAA_delta_fcc_in_bcc.mat
│   │   ├── BAA_delta_fcc_in_cp.mat
│   │   ├── BAA_delta_fcc_in_fcc.mat
│   │   ├── BAA_delta_fcc_in_hcp.mat
│   │   ├── BAA_delta_hcp_in_bcc2.mat
│   │   ├── BAA_delta_hcp_in_bcc.mat
│   │   ├── BAA_delta_hcp_in_cp2.mat
│   │   ├── BAA_delta_hcp_in_cp.mat
│   │   ├── BAA_delta_hcp_in_fcc2.mat
│   │   ├── BAA_delta_hcp_in_fcc.mat
│   │   ├── BAA_delta_hcp_in_hcp2.mat
│   │   ├── BAA_delta_hcp_in_hcp.mat
│   │   ├── delta_bcc_in_bcc.mat
│   │   ├── delta_bcc_in_fcc.mat
│   │   ├── delta_bcc_in_hcp.mat
│   │   ├── delta_fcc_in_bcc.mat
│   │   ├── delta_fcc_in_fcc.mat
│   │   ├── delta_fcc_in_hcp2.mat
│   │   ├── delta_fcc_in_hcp.mat
│   │   ├── delta_hcp_in_bcc.mat
│   │   ├── delta_hcp_in_fcc.mat
│   │   ├── delta_hcp_in_hcp2.mat
│   │   ├── delta_hcp_in_hcp.mat
│   │   └── fit_cft_bcc.mat
│   ├── 20210416
│   │   ├── delta_bcc_in_bcc.mat
│   │   ├── delta_bcc_in_fcc.mat
│   │   ├── delta_bcc_in_hcp.mat
│   │   ├── delta_fcc_in_bcc.mat
│   │   ├── delta_fcc_in_fcc.mat
│   │   ├── delta_fcc_in_hcp.mat
│   │   ├── delta_hcp_in_bcc.mat
│   │   ├── delta_hcp_in_fcc.mat
│   │   └── delta_hcp_in_hcp.mat
│   ├── 5nm-Au_data.mat
│   ├── RefDataBall.mat
│   ├── RefDataCell.mat
│   └── RefDataSquare.mat
├── myorth.m
├── my_script_1.m
├── Neg.m
├── Orien.m
├── p_norm.m
├── RawData
│   ├── 5nm-Au_data.txt
│   ├── Au-8.16nm.txt
│   ├── Au-sphere_diameter4.896nm-initial.xyz
│   ├── Au-sphere_diameter4.896nm-nvt-room.xyz
│   ├── Au-sphere_diameter9.792nm-initial.xyz
│   ├── Au-sphere_diameter9.792nm-nvt-room.xyz
│   ├── bak
│   │   ├── Au-sphere_diameter9.792nm-initial - 副本.xyz
│   │   └── Au-sphere_diameter9.792nm-nvt-room - 副本.xyz
│   └── README.md
├── README
└── references
    └── 0.ackland2006.pdf

9 directories, 100 files
