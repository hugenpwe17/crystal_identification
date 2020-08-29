function origin_index = tracesource(Atm_All,Atm_Dat_grain_cell,k)
%myFun - Description
%   输入:形变坐标,输出:初始坐标的索引
% Sycontection output = myFun(input)
%
% Long description
origin_index=zeros(11,1);    
for i=1:13
[M,~]=find(Atm_All(:,3:5)==Atm_Dat_grain_cell(i,:,k));
temp=max(tabulate(M(:)));
origin_index(i,1)=Atm_All(temp(1),1);%对应的lammps输出的原子ID
%origin_index(i,1)=temp(1);         %对应matlab里面的索引
end