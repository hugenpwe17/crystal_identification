%[�������]=p_norm[������,��������,��������]
%��ռ��е����������֮��ķ���
%p��ȡ2,�Ա�ʾ����

function [dis_mat]=p_norm(N,data,p)
    dis_mat=zeros(N,N);
    for i=1:N
    V=data;
    Vn=V-data(i,:);
    dis_mat(:,i)=sum(abs(Vn).^p,2).^(1/p);
    end
end