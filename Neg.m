%���ھ���������,[�ھ�����]=Neb[��������,��λ��]
function [neighbor]=Neg(Dat_position,parameter_conum)
temp            =   Dat_position;
neighbor        =   zeros(length(temp),1+parameter_conum);
for i=1:length(temp)
    
    temp                    =   temp-temp(i,:);
    dis                     =   sum(abs(temp).^2,2).^(1/2);
    disort                  =   sort(dis);
    [~,neighbor(i,:)]       =   ismember(disort(1:parameter_conum+1,:),dis);
    

end
end