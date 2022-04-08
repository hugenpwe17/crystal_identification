%---------------------------主程序-----------------------------------------
%% 导入数据
Atm_All                 =   load('data/Au-8.16nm.txt');
Atm_Dat                 =   Atm_All(:,3:5);

%% 参数区
parameter_N             =   length(Atm_All);
%粒子数
parameter_latconst      =   4.0783;
%晶格常数
parameter_conum         =   12;
%配位数
parameter_Phase         =   2;
%晶相 (1=BCC;2=FCC;3=HCP)
parameter_corfac        =   1.2;
parameter_Eps           =   (parameter_latconst/sqrt(2))*parameter_corfac;
%搜索半径Eps(&a/sqrt(2)=FCC,BCC;?=HCP)
parameter_MinPts     	=   parameter_conum;
%?密度阈值(4=FCC,5=BCC,?=HCP)

%% 识别晶相(1=BCC,2=FCC,3=HCP,0=Other)
Atm_Dat                 =   BAA(Atm_Dat);
%此程序运行时间较久
Atm_Dat_single          =   Atm_Dat(Atm_Dat(:,4)==parameter_Phase,:);

%% 单相晶体邻居粒子
Atm_Dat_single_neg=Neg(Atm_Dat_single,parameter_conum);

%% 消除晶界处原子

temp                    =   Atm_Dat_single;
Atm_Dat_single_pair     =   zeros(length(temp),1);
for n                   =   1:length(temp)
    Vec                 =   temp(Atm_Dat_single_neg(n,:),:);
    Vec                 =   Vec-Vec(1,:);
    for i=2:1+parameter_conum
        for j=i+1:1+parameter_conum
            cos_nij     =   dot(Vec(i,:),Vec(j,:))/norm(Vec(i,:))/norm(Vec(j,:));
            if ((cos_nij>=-1.0)&&(cos_nij<-0.945))
                %求(-1.000,-0.945)角度余弦分布
                Atm_Dat_single_pair(n)     =   Atm_Dat_single_pair(n)+1;
            end
        end
    end
end
clear temp n i j V cos_nij 

Atm_Dat_grain           =   unique(Atm_Dat_single...
    (Atm_Dat_single_neg(Atm_Dat_single_pair==(parameter_conum/2),:),:),'rows');
%将数据集中化

%% 提取晶粒
temp                    =   Atm_Dat_grain(:,1:3);
dis_mat                 =   p_norm(length(temp),temp,2);
Atm_Dat_grain(:,4)      =   DBSCAN(dis_mat,parameter_Eps,parameter_MinPts);
clear temp dis_mat
%DBSCAN识别晶粒
Atm_Dat_grain           =   Atm_Dat_grain(Atm_Dat_grain(:,4)~=0,:);
%去除噪音点

%% 输出图像
figure(1);
scatter3(Atm_Dat_grain(:,1),Atm_Dat_grain(:,2),Atm_Dat_grain(:,3),10,Atm_Dat_grain(:,4),'filled')
figure(2);
histogram(Atm_Dat_grain(:,4));

%% 求各晶粒的最对称晶胞
Atm_Dat_grain_cell      =   zeros(parameter_conum+1,3,max(Atm_Dat_grain(:,4)));
for i=1:max(Atm_Dat_grain(:,4))
    temp                =   Atm_Dat_grain(Atm_Dat_grain(:,4)==i,:);
    temp_neg            =   Neg(temp,parameter_conum);
    Vec_norm            =   zeros(1,length(temp));
    for j=1:length(temp)
        Vec             =   temp(temp_neg(j,:),1:3);
        Vec             =   Vec-Vec(1,:);
        Vec_norm(j)     =   norm(sum(Vec));
        %中心对称参数
    end
    index               =   Vec_norm==min(Vec_norm);
    Atm_Dat_grain_cell(:,:,i)  ... 
                        =   temp(temp_neg(index,:),1:3);
end
    clear i j index Vec Vec_norm temp temp_neg
%% 晶胞中心移向原点
for i=1:size(Atm_Dat_grain_cell,3)
    Atm_Dat_grain_cell_mod(:,:,i)...
                        =   Atm_Dat_grain_cell(:,:,i)-Atm_Dat_grain_cell(1,:,i);
end
clear i

%%
%分配晶向数组
Atm_Dat_grain_cell_orien...
                        =   zeros(1,3,size(Atm_Dat_grain_cell_mod,3));
for n=1:size(Atm_Dat_grain_cell_mod,3)
Atm_Dat_grain_cell_orien(:,:,n)...
                        =   Orien(n,Atm_Dat_grain_cell_mod);
end
clear n

%% 三角剖分画图
%选择n号晶胞
n                       =   1;
T                       =   Atm_Dat_grain_cell_mod(:,:,n);
Tes                     =   delaunayn(T);
tetramesh(Tes,T)
hold on
quiver3(0,0,0,Atm_Dat_grain_cell_orien(:,1,n),Atm_Dat_grain_cell_orien(:,2,n),Atm_Dat_grain_cell_orien(:,3,n))
clear n T Tes