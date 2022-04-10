function [Data,chimat] = BAAm(data)
% Idengtification of structure type with bond-angle analysis
% [structure type] = BAA[point coords]

% Measure the total number of atoms in the (crystal) dataset
N = length(data);

% Initialize crystal structure data
str = zeros(N,1);

% define norm function to calculate the row value
norm2 = @(x)sum(abs(x).^2,2).^(1/2);

% The average distance of the nearest six neighboring atoms is used as the reference distance
RefNeiNum  = 6;
% Initialize the reference radius matrix (Six-atom dist average)
R0 = zeros(N,1);

% Initialize the first shell particle set matrix (fcc\hcp:12)
N0 = cell(N,1);
% Initialize the second shell particle set matrix (bcc:14)
N1 = N0;

% Find neighbor atoms of each atom
% row index is atom id,coloum index is neighbor atoms' id (with ascengding distance)
for i = 1:N
    % Calculate the distance of each atom relative to other atoms
    DistNei = norm2(data - data(i,:));
    % ascending order distance
    DistOrd = sort(DistNei);
    
    % Calculate R0
    R0(i)   = sqrt(sum(DistOrd(2:RefNeiNum+1).^2)/RefNeiNum);
    
    % Calculate N0, method 1
    N0(i)   = {(find(DistNei <= ( sqrt(1.45) *R0(i))))'};
    % Calculate N0, method 2 
    % ((1 + sqrt(2))/2 ~= sqrt(1.45))
    % N0(i) = {(find(DistNei <= ((1+sqrt(2))/2 * R0(i))))'};
    
    % Calculate N1, method 1
    N1(i)   = {(find(DistNei <= (sqrt(1.55) * R0(i))))'};
    % Calculate N1, method 2
    % N1(i) = {(find(DistNei <= (((1 + sqrt(2))/2 + 1) * R0(i))))'};
    % Calculate N1, method 3
    % N1(i) = {(find(DistNei <= ((1 + sqrt(2))/2 * (8/14 + 6/14 * 2/sqrt(3)) * R0(i))))'};
    
    % Remove element 0s (the distance of an atom to itself) from the distance data to neighboring atoms
    N0{i}(N0{i} == i) = [];
    N1{i}(N1{i} == i) = [];
end

% Initialize the angle cosine distribution matrix
chi = cell(N,1);

% Cos value matrix of characteristic angles
Y   = [-0.945 -0.915 -0.755 -0.195 0.195 0.245 0.795 1];

% Calculate the angle cosine distribution matrix
for i = 1:N
    % vectorization
    NeiofAtom = data(N0{i},:) - data(i,:);
    % pairwise atoms' angle consine value
    AngleJug  = 1 - pdist(NeiofAtom, 'cosine');
    % Interval distribution (left open, right closed)
    chi(i,:)  = {sum(cumsum((Y - AngleJug')>0, 2) == 1)};
end

% Cell to matrix
chimat = cell2mat(chi);

% Delta
% bcc
delta(:,1)=0.35*chimat(:,5)./(chimat(:,6)+chimat(:,7)-chimat(:,5));        
% fcc
delta(:,2)=0.61*(abs(chimat(:,1)+chimat(:,2)-6)+chimat(:,3))/6;        
% hcp
delta(:,3)=(abs(chimat(:,1)-3)+abs(chimat(:,1)+chimat(:,2)+chimat(:,3)-9))/12;   
% cp
delta(:,4)=abs(1-chimat(:,7)./24);                                       

% Judgment module
% define the structure type
% none = 0;bcc = 1;fcc = 2;hcp = 3;cp = 4;ico = 4;
% Removes atoms that are not complete in the local area
for i=1:N
    if length(N0{i})<6
        str(i)=0;
    end
end
% Preliminary judgment
delta(chimat(:,1)==7, 1) = 0;
delta(chimat(:,1)==6, 2) = 0;
delta(chimat(:,1)<=3, 3) = 0;
str  (chimat(:,1)==7, 1) = 1;
str  (chimat(:,1)==6, 1) = 2;
str  (chimat(:,1)<=3, 1) = 3;
% judge flow
for i=1:N
    if chimat(i,8)>0
         str(i) = 0;
    elseif chimat(i,5)<3
        if (length(N1{i})>13||length(N1{i})<11)
            str(i) = 0;
        else
            str(i) = 0;%(or ico ==4)
        end
    elseif (delta(i,1)<=delta(i,4))
        if (length(N1{i})<11)
            str(i) = 0;
        else
            str(i) = 1;
        end
    elseif (length(N1{i})>12)||(length(N1{i})<11)
            str(i) = 0;
    elseif (delta(i,2)<delta(i,3))
            str(i) = 2;
    else
            str(i) = 3;
    end
end
% add the structure type column in data matrix
Data = data;
Data(:,4) = str;
return 
end