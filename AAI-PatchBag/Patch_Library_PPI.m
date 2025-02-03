%%%% Special Edition for 6 plus 1 Patches
%%%% Patch Preprocessing
% Load Patch Merge
load('Reduced_Sample_Int_Res_Results/Unique_Patch_Merge_PPI.mat');

% Unique Res_Type Patches
Res_Type=[Unique_Patch_Merge.Res_Type];
Reshpaed_Res_Type=reshape(Res_Type,7,[])';
[~, Idx]=unique(Reshpaed_Res_Type, 'rows');
Unique_Res_Type_Patch=Unique_Patch_Merge(Idx);

Reshpaed_Unique_Res_Type=reshape([Unique_Res_Type_Patch.Res_Type],7,[])';
Unique_NonX_Res_Type_Patch=Unique_Res_Type_Patch(find(~sum(Reshpaed_Unique_Res_Type==0,2)));

%%%% % Patch Library
% Library Size
Count=8000;

% Patch Library
Seed=RandStream('mrg32k3a');
Idx=randsample(Seed, size(Unique_NonX_Res_Type_Patch,2), Count);

% Reshpaed_Norm_V=reshape(Norm_V,3,[])';
% plot3(Reshpaed_Norm_V(Idx,1), Reshpaed_Norm_V(Idx,2), Reshpaed_Norm_V(Idx,3), '.')
Patch_Library_8000=Unique_NonX_Res_Type_Patch(Idx);
clearvars -except Count Patch_Library_8000

% Distance Matrix (DISM) of Patch Library
% parpool (8)
P=perms(1:size(Patch_Library_8000(1).Coord,1)-1);
nm_RMSD_DISM=[];
Error_ID=[];

parfor n=1:(Count-1)
    P1=Patch_Library_8000(n).Coord;
    m_RMSD_DISM=zeros(Count-n, 1);

    for m=n+1:Count
        RMSD_Table=zeros(size(P,1), 1);
        Angle_Table=zeros(size(P,1), 1);

        for r=1:size(P,1)
            P2=Patch_Library_8000(m).Coord([P(r,:), 7],:);

            [R,~,eRMSD] = CoordiExam_AC(P1, P2);
            RMSD_Table(r,1)=eRMSD;
            Angle_Table(r,1)=acosd(dot(Patch_Library_8000(n).Norm_V*R, Patch_Library_8000(m).Norm_V));
        end

        if isempty(RMSD_Table(Angle_Table<90))
            Error_ID=[Error_ID; [n,m]];
            m_RMSD_DISM(m-n,1)=20;
        else
            m_RMSD_DISM(m-n,1)=min(RMSD_Table(Angle_Table<90),[],'all');
        end
    end
    nm_RMSD_DISM=[nm_RMSD_DISM; m_RMSD_DISM];
    n
end

Patch_Library_DISM=squareform(nm_RMSD_DISM);
save('Reduced_Sample_Int_Res_Results/Patch_Library_DISM.mat', 'Patch_Library_DISM');

%

%%%% Clustering of Patch Library
% Load PDB ID of Qualified Reduced  Samples
Filename='Qualified_Reduced_Sample_3737(Exclude 4YKN).txt';
ID_List=fopen(Filename);
sp=1;Sample_Info(1).ID=[];
while (1)
    line=fgetl(ID_List);
    if line==-1, break, end
    Sample_Info(sp).ID=sscanf(line,'%c');
    sp=sp+1;
end
fclose(ID_List);

% Average Count No.
% mean(Count_Table)
Count_Table=zeros(size(Sample_Info));
for i=1:size(Sample_Info,1)
    PDB_ID=lower(Sample_Info(i).ID);
    File=table2array(readtable([PDB_ID, '_Patch.txt']));
    Count_Table(i)=size(File,1);
end


%%%% K-medoid Clustering
% Load Patch_Library_DISM.mat
load('Reduced_Sample_Int_Res_Results/Patch_Library_DISM.mat')

% Distance matrix
DISM=Patch_Library_DISM;
N = size(DISM,1);

Max_Cluster_No=300;
[idx, C] = kmedoids((1:N)', Max_Cluster_No, 'Distance', @(row, col) DISM(row, col), 'Replicates', 5);

% Medoid
load('Reduced_Sample_Int_Res_Results/Patch_Library_8000.mat')
Clust_medoids=Patch_Library_8000(C);

save('Reduced_Sample_Int_Res_Results/Clust_medoids_300.mat', 'Clust_medoids');
