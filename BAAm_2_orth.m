function [Data,chimat,str,delta] = Copy_of_BAAm(data)
% idengtification of structure type with bond-angle analysis
% [structure type] = BAA[point coords]
% measure the number of atom in crystal
N = length(data);
%
str = zeros(N,1);
% define norm function to calculate the row value
norm2 = @(x)sum(abs(x).^2,2).^(1/2);

% Neighbor matrix including n point 
% nNeimat: rowid is atom id,colid is nei atom id with ascengding distance
NeiNum  = 6;
% Reference radius (Six-atom dist average)
R0 = zeros(size(data,1),1);
% First shell particle set (fcc)
N0 = cell(size(data,1),1);
% Second shell particle set (bcc)
N1 = N0;
for i = 1:N
    % distance
    DistNei = norm2(data-data(i,:));
    % ascending order distance
    DistOrd = sort(DistNei);
    % reference distance of a atom
%     R0(i)   = mean(DistOrd(2:NeiNum+1));
    R0(i)   = sqrt(sum(DistOrd(2:NeiNum+1).^2)/NeiNum);
    % N0: {[x,y]|R < (1+sqrt(2))/2 *R0}
    N0(i)   = {(find(DistNei <= (sqrt(1.45) *R0(i))))'};
    % delete itself
    N0{i}(N0{i}==i)= [];
    % N1: {[x,y]|R < (1+sqrt(2))/2 *R0}
    N1(i)   = {(find(DistNei <= (sqrt(1.55) *R0(i))))'};
    % delete itself
    N1{i}(N1{i}==i)= [];
end
% angle distribution
chi = cell(N,1);
% angle interval
Y     = [-0.945 -0.915 -0.755 -0.195 0.195 0.245 0.795 1];
for i = 1:N
    % vectorization
    NeiofAtom = data(N0{i},:)-data(i,:);
    % pairwise angle consine value
    AngleJug  = 1 - pdist(NeiofAtom,'cosine');
    % Interval distribution (>=,<)
    chi(i,:)= {sum(cumsum((Y-AngleJug')>0,2)==1)};

end
% cell to matrix
chimat = cell2mat(chi);

% vec
% % bcc
% bccvec =[7 0 0 36 12 0 36 0]; 
% % fcc
% fccvec =[6 0 0 24 12 0 24 0];
% % hcp
% hcpvec =[3 0 6 21 12 0 24 0];

% bcc
bccvec =[-1.6041 0 -0.0406  1.5228 -5.2487 0  1.5025 0]; 
% fcc
fccvec =[ 1.5553 0 -1.4777 -0.2304  2.6918 0 -0.9693 0];
% hcp
hcpvec =[-2.5794 0  6.0000 -1.3925  0.8598 0  1.6075 0];
delta = zeros(N,3);
for i = 1:N
delta(i,1) = abs(dot(chimat(i,:),bccvec)/norm(bccvec)/norm(bccvec)-1);
delta(i,2) = abs(dot(chimat(i,:),fccvec)/norm(fccvec)/norm(fccvec)-1);
delta(i,3) = abs(dot(chimat(i,:),hcpvec)/norm(hcpvec)/norm(hcpvec)-1);
end

% Judgment module
% define the structure type
% none = 0;bcc = 1;fcc = 2;hcp = 3;cp = 4;ico = 4;

% Preliminary judgment
for i =1:N
    str(i,1)=(delta(i,:)==min(delta(i,:)))*(1:3)';
end
% Removes atoms that are not complete in the local area
for i=1:N
    if length(N1{i}) < 12
        str(i) = 0;
    elseif min(delta(i,:))>0.6
        str(i) = 0;
    end
end
%



% add the structure type column in data matrix
Data = data;
Data(:,4) = str;
return 
end