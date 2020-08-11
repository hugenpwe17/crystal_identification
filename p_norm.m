%[距离矩阵]=p_norm[粒子数,坐标数据,范数阶数]
%求空间中点阵的两两点之间的范数
%p常取2,以表示距离

function [dis_mat]=p_norm(N,data,p)
    dis_mat=zeros(N,N);
    for i=1:N
    V=data;
    Vn=V-data(i,:);
    dis_mat(:,i)=sum(abs(Vn).^p,2).^(1/p);
    end
end