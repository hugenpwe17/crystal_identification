%[聚类数据]=DBSCAN[距离矩阵,搜索半径,密度阈值]
%-----------------------------------------------------------------------------
%DBSCAN聚类算法，输入坐标的dismat，最小搜索半径eps，和所含最小值minpts，得到聚类值
%基于密度的算法,用于在有噪声的大型空间数据中发现集群
%没有实现空间访问方法,运行时间复杂度O(N^2)
%-----------------------------------------------------------------------------
%输入:      空间距离矩阵DistMat, 搜索半径Eps, 密度阈值:MinPts
%DistMat:   一个N*N距离矩阵, (i,j)上的元素表示点i到点j的距离
%Eps:       粒子邻域搜索半径
%MinPts:    Eps邻域中可以被搜索到的最少的点的数量
%-----------------------------------------------------------------------------
%输出:      聚类数据Clust
%Clust:  一个N*1的向量,记录每个点属于哪个簇,0表示是噪音
%-----------------------------------------------------------------------------
%将所有粒子初始化为-1,表示未分类
%-----------------------------------------------------------------------------
function [Clust] = DBSCAN(DistMat,Eps,MinPts)

Clust=zeros(size(DistMat,1),1)-1;
ClusterId=1;

%随机遍历
VisitSequence=randperm(length(Clust));

for i=1:length(Clust)
    %对每个点检查是否未被检测
    pt=VisitSequence(i);
    if Clust(pt)==-1
        %通过密度可达性迭代扩展簇
        [Clust,isnoise]=ExpandCluster(DistMat,pt,ClusterId,Eps,MinPts,Clust);
        if ~isnoise
            ClusterId=ClusterId+1;
        end
    end
end

end
%-----------------------------------------------------------------------------
function [Clust,isnoise]=ExpandCluster(DistMat,pt,ClusterId,Eps,MinPts,Clust)

%区域查询
seeds=find(DistMat(:,pt)<=Eps);     %获取满足Eps的点的索引

if length(seeds)<MinPts             %检测并处理噪音
    Clust(pt)=0;                    %保留为噪音
    isnoise=true;
    return
else
    Clust(seeds)=ClusterId;         %加入簇中
    seeds=setxor(seeds,pt);         %去除中心点
    while ~isempty(seeds)           %开始迭代,结束条件没有满足条件的粒子
        currentP=seeds(1);
        %区域检测
        result=find(DistMat(:,currentP)<=Eps);
        if length(result)>=MinPts
            for i=1:length(result)
                resultP=result(i);
                if Clust(resultP)==-1||Clust(resultP)==0 %检测是未定义或者噪音
                    if Clust(resultP)==-1                %筛去噪音
                        seeds=[seeds(:);resultP];        %点加入待检测簇内
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