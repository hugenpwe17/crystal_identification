%---------------------------������-----------------------------------------
%% ��������
Atm_All                 =   load('C:\Users\hideo\Documents\��ѧѵ���γ�\new code\Au-8.16nm.txt');
Atm_Dat                 =   Atm_All(:,3:5);

%% ������
parameter_N             =   length(Atm_All);
%������
parameter_latconst      =   4.0783;
%������
parameter_conum         =   12;
%��λ��
parameter_Phase         =   2;
%���� (1=BCC;2=FCC;3=HCP)
parameter_corfac        =   1.2;
parameter_Eps           =   (parameter_latconst/sqrt(2))*parameter_corfac;
%�����뾶Eps(&a/sqrt(2)=FCC,BCC;?=HCP)
parameter_MinPts     	=   parameter_conum;
%?�ܶ���ֵ(4=FCC,5=BCC,?=HCP)

%% ʶ����(1=BCC,2=FCC,3=HCP,0=Other)
Atm_Dat                 =   BAA(Atm_Dat);
%�˳�������ʱ��Ͼ�,��ֻ������һ��,֮���ظ��������ע��
Atm_Dat_single          =   Atm_Dat(Atm_Dat(:,4)==parameter_Phase,:);

%% ���ྦྷ���ھ�����
Atm_Dat_single_neg=Neg(Atm_Dat_single,parameter_conum);

%% �������紦ԭ��,�������

temp                    =   Atm_Dat_single;
Atm_Dat_single_pair     =   zeros(length(temp),1);
for n                   =   1:length(temp)
    Vec                 =   temp(Atm_Dat_single_neg(n,:),:);
    Vec                 =   Vec-Vec(1,:);
    for i=2:1+parameter_conum
        for j=i+1:1+parameter_conum
            cos_nij     =   dot(Vec(i,:),Vec(j,:))/norm(Vec(i,:))/norm(Vec(j,:));
            if ((cos_nij>=-1.0)&&(cos_nij<-0.945))
                %��(-1.000,-0.945)�Ƕ����ҷֲ�
                Atm_Dat_single_pair(n)     =   Atm_Dat_single_pair(n)+1;
            end
        end
    end
end
clear temp n i j V cos_nij 

Atm_Dat_grain           =   unique(Atm_Dat_single...
    (Atm_Dat_single_neg(Atm_Dat_single_pair==(parameter_conum/2),:),:),'rows');
%�����ݼ��л�

%% ��ȡ����
temp                    =   Atm_Dat_grain(:,1:3);
dis_mat                 =   p_norm(length(temp),temp,2);
Atm_Dat_grain(:,4)      =   DBSCAN(dis_mat,parameter_Eps,parameter_MinPts);
clear temp dis_mat
%DBSCANʶ����
Atm_Dat_grain           =   Atm_Dat_grain(Atm_Dat_grain(:,4)~=0,:);
%ȥ��������

%% ���ͼ��
figure(1);
scatter3(Atm_Dat_grain(:,1),Atm_Dat_grain(:,2),Atm_Dat_grain(:,3),50,Atm_Dat_grain(:,4),'filled')
figure(2);
histogram(Atm_Dat_grain(:,4));

%% �����������Գƾ���
Atm_Dat_grain_cell      =   zeros(parameter_conum+1,3,max(Atm_Dat_grain(:,4)));
for i=1:max(Atm_Dat_grain(:,4))
    temp                =   Atm_Dat_grain(Atm_Dat_grain(:,4)==i,:);
    temp_neg            =   Neg(temp,parameter_conum);
    Vec_norm            =   zeros(1,length(temp));
    for j=1:length(temp)
        Vec             =   temp(temp_neg(j,:),1:3);
        Vec             =   Vec-Vec(1,:);
        Vec_norm(j)     =   norm(sum(Vec));
        %���ĶԳƲ���
    end
    index               =   Vec_norm==min(Vec_norm);
    Atm_Dat_grain_cell(:,:,i)  ... 
                        =   temp(temp_neg(index,:),1:3);
end
    clear i j index Vec Vec_norm temp temp_neg
%% ������������ԭ��
for i=1:length(Atm_Dat_grain_cell)
    Atm_Dat_grain_cell(:,:,i)...
                        =   Atm_Dat_grain_cell(:,:,i)-Atm_Dat_grain_cell(1,:,i);
end
clear i

%%
%���侧������
Atm_Dat_grain_cell_orien...
                        =   zeros(1,3,length(Atm_Dat_grain_cell));
for n=1:length(Atm_Dat_grain_cell)
Atm_Dat_grain_cell_orien(:,:,n)...
                        =   Orien(n,Atm_Dat_grain_cell);
end
clear n

%% �����ʷֻ�ͼ
n                       =   8;
%ѡ��n�ž���
T                       =   Atm_Dat_grain_cell(:,:,n);
Tes                     =   delaunayn(T);
tetramesh(Tes,T)
hold on
quiver3(0,0,0,Atm_Dat_grain_cell_orien(:,1,n),Atm_Dat_grain_cell_orien(:,2,n),Atm_Dat_grain_cell_orien(:,3,n))
clear n T Tes

