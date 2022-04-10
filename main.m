%---------------------------main script-----------------------------------------
%% Import Data
% Atm_All is raw data
Atm_All = load('RawData/Au-8.16nm.txt');
% Atm_Dat is coords ([x,y,z]) of each atom in crystal
Atm_Dat = Atm_All(:,3:5);

%% Default parameters
% BAA: Identify crystal phases
% The total number of particles contained in Atm_Dat
parameter_N = length(Atm_All);
% Lattice constants defined in the in file of lammps
parameter_lat_const = 4.0783;
% The coordination number of the target structure to be retrieved (BCC=>8(14);FCC=>12;HCP=>12)
parameter_coord_num = 12;
% The crystal phase type of the target structure to be retrieved (BCC=>1;FCC=>2;HCP=>3)
parameter_phase_type = 2;

% DBSCAN: Identify crystal orientation
% search radius correction factor
parameter_correction_factor = 1.2;
% Search radius (clustering method)
parameter_Eps = (parameter_lat_const / sqrt(2) ) * parameter_correction_factor;
% Density threshold
parameter_MinPts = parameter_coord_num;

%% Identify crystal phases
% main function (longer run time)
Atm_Dat = BAAm(Atm_Dat);
% Filter out single-phase structures (1=>BCC,2=>FCC,3=>HCP,0=>Other)
Atm_Dat_single = Atm_Dat(Atm_Dat(:,4)==parameter_phase_type,:);

%% neighbor particles set of single-phase crystal
Atm_Dat_single_neg = Neg(Atm_Dat_single, parameter_coord_num);

%% Delete atoms at grain boundaries
temp                = Atm_Dat_single;
Atm_Dat_single_pair = zeros(length(temp),1);
for n               = 1:length(temp)
    Vec             = temp(Atm_Dat_single_neg(n,:),:);
    Vec             = Vec-Vec(1,:);
    for i = 2 : 1+parameter_coord_num
        for j = i+1 : 1+parameter_coord_num
            cos_nij = dot(Vec(i,:),Vec(j,:))/norm(Vec(i,:))/norm(Vec(j,:));
            % not grain boundary atomic condition
            if ((cos_nij>=-1.0)&&(cos_nij<-0.945))
                Atm_Dat_single_pair(n) = Atm_Dat_single_pair(n)+1;
            end
        end
    end
end
clear temp n i j V cos_nij

% deduplicate data
Atm_Dat_grain = unique(Atm_Dat_single(Atm_Dat_single_neg(Atm_Dat_single_pair==(parameter_coord_num/2),:),:),'rows');

%% Filter out the desired grains
temp               = Atm_Dat_grain(:,1:3);
dis_mat            = p_norm(length(temp),temp,2);
% DBSCAN identifies the grain
Atm_Dat_grain(:,4) = DBSCAN(dis_mat,parameter_Eps,parameter_MinPts);
clear temp dis_mat
% remove noise points
Atm_Dat_grain      = Atm_Dat_grain(Atm_Dat_grain(:,4)~=0,:);

%% Print images
% 3D image of grains
figure(1);
scatter3(Atm_Dat_grain(:,1),Atm_Dat_grain(:,2),Atm_Dat_grain(:,3),10,Atm_Dat_grain(:,4),'filled')
% histograms of atom number within the grains
figure(2);
histogram(Atm_Dat_grain(:,4));

%% Find the most probable symmetric unit cell for each grain
Atm_Dat_grain_cell      =   zeros(parameter_coord_num+1,3,max(Atm_Dat_grain(:,4)));
for i=1:max(Atm_Dat_grain(:,4))
    temp                =   Atm_Dat_grain(Atm_Dat_grain(:,4)==i,:);
    temp_neg            =   Neg(temp,parameter_coord_num);
    Vec_norm            =   zeros(1,length(temp));
    for j=1:length(temp)
        Vec             =   temp(temp_neg(j,:),1:3);
        Vec             =   Vec-Vec(1,:);
        % Center symmetry parameter
        Vec_norm(j)     =   norm(sum(Vec));
    end
    index               =   Vec_norm==min(Vec_norm);
    Atm_Dat_grain_cell(:,:,i)  ...
        =   temp(temp_neg(index,:),1:3);
end
clear i j index Vec Vec_norm temp temp_neg

%% The center of the unit cell moves to the origin of the coordinates
Atm_Dat_grain_cell_mod = Atm_Dat_grain_cell;
for i=1:size(Atm_Dat_grain_cell,3)
    Atm_Dat_grain_cell_mod(:,:,i)...
        =   Atm_Dat_grain_cell(:,:,i)-Atm_Dat_grain_cell(1,:,i);
end
clear i

%% Calculate the orientation vector
Atm_Dat_grain_cell_orien...
    =   zeros(1,3,size(Atm_Dat_grain_cell_mod,3));
for n=1:size(Atm_Dat_grain_cell_mod,3)
    Atm_Dat_grain_cell_orien(:,:,n)...
        =   Orien(n,Atm_Dat_grain_cell_mod);
end
clear n
disp(['可选择的晶胞数量为:',num2str(size(Atm_Dat_grain_cell_orien,3))])

%% Print image
% Select the n(th) unit cell
n   = 18;
T   = Atm_Dat_grain_cell_mod(:,:,n);
Tes = delaunayn(T);
tetramesh(Tes,T)
disp([num2str(n), '号晶胞晶向是:']);
disp(Atm_Dat_grain_cell_orien(:,:,1));
disp(['可选择的晶胞序号为:','1-',num2str(size(Atm_Dat_grain_cell_orien,3))])
hold on
quiver3(0,0,0,Atm_Dat_grain_cell_orien(:,1,n),Atm_Dat_grain_cell_orien(:,2,n),Atm_Dat_grain_cell_orien(:,3,n))
hold off
clear n T Tes