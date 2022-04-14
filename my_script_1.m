%% Calculate the distribution of bond angle cosines under thermal disturbances (simulated using rand)
% load bcc/fcc/hcp cell unit
load('RefDataCell.mat');
% set test data
data = bcc;
% minimal distance (data(1,:) is the center point)
dist_list = sum(abs(data - data(1,:)).^2,2).^(1/2);
dist_list(dist_list==0) = [];
dist_min = min(dist_list);
% rand displacement (disp form -1 to 1)
disp = -1 + 2 * rand(size(data,1),3);
% normlization disp 
for i = 1:size(data,1)
    disp (i,:) = disp (i,:) / norm(disp(i,:));
end
clear i
% set displacement
disp_factor = 0.1;
datam = data + disp_factor * disp * dist_min;
% total number of particle
N = length(datam);
% define norm function to calculate the row value
norm2 = @(x)sum(abs(x).^2,2).^(1/2);
% Neighbor matrix including n point 
% nNeimat: rowid is atom id,colid is nei atom id with ascengding distance
NeiNum  = 6;
% Reference radius (Six-atom dist average)
R0 = zeros(size(datam,1),1);
% calculate N0
N0 = cell(size(datam,1),1);
for i = 1:N
    % distance
    dist_nei = norm2(datam-datam(i,:));
    % ascending order distance
    dist_order = sort(dist_nei);
    % reference distance of a atom
    R0(i) = mean(dist_order(2:NeiNum+1));
    % N0: {[x,y]|R < (1+sqrt(2))/2 *R0}
    N0(i) = {(find(dist_nei <= ((1+sqrt(2))/2) *R0(i)))'};
    % delete itself
    N0{i}(N0{i}==i)= [];
end
clear i

for i = 1:N
    % vectorization
    NeiofAtom = datam(N0{i},:)-datam(i,:);
    % pairwise angle consine value
    AngleJug = 1 - pdist(NeiofAtom,'cosine');
    % Interval distribution (>=,<)
    A{i} = AngleJug;
end
clear i

A1 = cell2mat(A);
edges = (-1:0.02:1);
h1 = histogram(A1,edges);
% frequency
L = length(h1.BinEdges);
x1 = zeros(size(h1.BinEdges,1)-1,1);
for i = 1:L-1
    x1(i) = (h1.BinEdges(i) + h1.BinEdges(i+1))/2;
end
clear i

%y1 = h1.Values/size(data,1)/h1.BinWidth;
%figure(2)
%plot(x1,y1,'-','linewidth',2,'MarkerSize',10)