%[��������]=DBSCAN[�������,�����뾶,�ܶ���ֵ]
%-----------------------------------------------------------------------------
%DBSCAN�����㷨�����������dismat����С�����뾶eps����������Сֵminpts���õ�����ֵ
%�����ܶȵ��㷨,�������������Ĵ��Ϳռ������з��ּ�Ⱥ
%û��ʵ�ֿռ���ʷ���,����ʱ�临�Ӷ�O(N^2)
%-----------------------------------------------------------------------------
%����:      �ռ�������DistMat, �����뾶Eps, �ܶ���ֵ:MinPts
%DistMat:   һ��N*N�������, (i,j)�ϵ�Ԫ�ر�ʾ��i����j�ľ���
%Eps:       �������������뾶
%MinPts:    Eps�����п��Ա������������ٵĵ������
%-----------------------------------------------------------------------------
%���:      ��������Clust
%Clust:  һ��N*1������,��¼ÿ���������ĸ���,0��ʾ������
%-----------------------------------------------------------------------------
%���������ӳ�ʼ��Ϊ-1,��ʾδ����
%-----------------------------------------------------------------------------
function [Clust] = DBSCAN(DistMat,Eps,MinPts)

Clust=zeros(size(DistMat,1),1)-1;
ClusterId=1;

%�������
VisitSequence=randperm(length(Clust));

for i=1:length(Clust)
    %��ÿ�������Ƿ�δ�����
    pt=VisitSequence(i);
    if Clust(pt)==-1
        %ͨ���ܶȿɴ��Ե�����չ��
        [Clust,isnoise]=ExpandCluster(DistMat,pt,ClusterId,Eps,MinPts,Clust);
        if ~isnoise
            ClusterId=ClusterId+1;
        end
    end
end

end
%-----------------------------------------------------------------------------
function [Clust,isnoise]=ExpandCluster(DistMat,pt,ClusterId,Eps,MinPts,Clust)

%�����ѯ
seeds=find(DistMat(:,pt)<=Eps);     %��ȡ����Eps�ĵ������

if length(seeds)<MinPts             %��Ⲣ��������
    Clust(pt)=0;                    %����Ϊ����
    isnoise=true;
    return
else
    Clust(seeds)=ClusterId;         %�������
    seeds=setxor(seeds,pt);         %ȥ�����ĵ�
    while ~isempty(seeds)           %��ʼ����,��������û����������������
        currentP=seeds(1);
        %������
        result=find(DistMat(:,currentP)<=Eps);
        if length(result)>=MinPts
            for i=1:length(result)
                resultP=result(i);
                if Clust(resultP)==-1||Clust(resultP)==0 %�����δ�����������
                    if Clust(resultP)==-1                %ɸȥ����
                        seeds=[seeds(:);resultP];        %������������
                    end
                    Clust(resultP)=ClusterId;
                end
                
            end
        end
        seeds=setxor(seeds,currentP);
    end
    isnoise=false;
    return 
end
end