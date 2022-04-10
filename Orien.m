% 计算晶胞的晶向向量
function [Atm_Dat_grain_cell_orien]=Orien(n,Atm_Dat_grain_cell)
count=0;
temp=Atm_Dat_grain_cell;
for i=2:13
    for j=2:13
        angle=dot(temp(i,:,n),temp(j,:,n))/(norm(temp(i,:,n))*norm(temp(j,:,n)));
        if angle>-0.2&&angle<0.2
            count=count+1;
            if count==2
                break
            end
        end
    end
    count=0;
end
clear count angle
Atm_Dat_grain_cell_orien=cross(temp(i,:,n),temp(j,:,n));
end