%%%% PPI
%%%% Patch to Vector
% Load Qualified Reduced Sample ID
% Filename='Reduced_Sample_Training_Test_Dataset.xlsx';
Filename='target_entry.txt';
ID_List=fopen(Filename);
sp=1;Sample_Info(1).ID=[];
while (1)
    line=fgetl(ID_List);
    if line==-1, break, end
    Sample_Info(sp).ID=sscanf(line,'%c');
    sp=sp+1;
end
fclose(ID_List);

% Load Merged Error ID
% load('Merged_Error_ID.mat')
% Sample_Info=Sample_Info(setdiff([1:2472], Merged_Error_ID));

% Load Clust Info (Representative Patch Library)
Max_Cluster_No=300;
load(['Package/Clust_medoids_', num2str(Max_Cluster_No), '_PPI.mat'])

% Patch Classification
Selected_Sample=Sample_Info;
Selected_Sample(1).Patch_code=[];
Selected_Sample(1).Patch_Res=[];
Selected_Sample(1).Bar_code=[];

Blank_ID=[];
P=perms(1:6);
for i=1:size(Selected_Sample,1)
    try
        % Patch Files
        PDB_ID=lower(Selected_Sample(i).ID);
        Patch_File=table2array(readtable(['Reduced_Sample_Int_Res_Results/Patch_Results_6plus1/',PDB_ID,'_Patch.txt']));
        [Unique_Patch, ~, Uidx]=unique(Patch_File, 'rows', 'stable');

        % Patch Codes
        Patch_code=zeros(size(Unique_Patch,1),1);

        parfor j=1:size(Unique_Patch,1)
            P1=reshape(Unique_Patch(j,1:21),7,3);
            RMSD_DISM=zeros(Max_Cluster_No, 1);

            for m=1:Max_Cluster_No
                RMSD_Table=zeros(size(P,1), 1);
                Angle_Table=zeros(size(P,1), 1);

                for r=1:size(P,1)
                    P2=Clust_medoids(m).Coord([P(r,:), 7],:);

                    [R,~,eRMSD] = CoordiExam_AC(P1, P2);
                    RMSD_Table(r,1)=eRMSD;
                    Angle_Table(r,1)=acosd(dot(Unique_Patch(j,22:24)*R, Clust_medoids(m).Norm_V));
                end

                if isempty(RMSD_Table(Angle_Table<90))
                    RMSD_DISM(m,1)=min(RMSD_Table,[],'all');
                else
                    RMSD_DISM(m,1)=min(RMSD_Table(Angle_Table<90),[],'all');
                end
            end
            [~,c]=min(RMSD_DISM);
            Patch_code(j,1)=c;
        end

        Unique_Count=hist(Uidx, (1:max(Uidx))');
        Selected_Sample(i).Patch_code=[Patch_code, Unique_Count'];
        Selected_Sample(i).Patch_Res=Unique_Patch(:,25:31);

        % Patch Vectorlization
        Bar_code=zeros(Max_Cluster_No,1);
        for u=1:Max_Cluster_No
            Bar_code(u,1)= sum(Unique_Count(Patch_code==u));
        end
        Selected_Sample(i).Bar_code=Bar_code;
    catch
        Blank_ID=[Blank_ID; i];
    end
    i
end

% save('Blank_ID.mat', 'Blank_ID');
save('Reduced_Sample_Int_Res_Results/Selected_Sample.mat', 'Selected_Sample');

%
%
%

load('Package/AA_Index.mat')
load('Reduced_Sample_Int_Res_Results/Selected_Sample.mat')

% Selected_Sample.Bar_code -> Spatial Distance
Count=[];
Merged_Table=[];
Selected_Sample=Selected_Sample';

parfor i=1:size(Selected_Sample,2)
    Patch_Res=Selected_Sample(i).Patch_Res;
    if ~isempty(Patch_Res)
        for j=1:size(Patch_Res,1)
            Idx=Patch_Res(j,:);
            Idx(Idx==0)=8;
            Merged_Table=[Merged_Table, mean(AA_Index(:,Idx),2)];
        end
        Count=[Count; Selected_Sample(i).Patch_code(:,2)];
    end
    i
end

Completed_Table=[];
parfor k=1:size(Count,1)
    Completed_Table=[Completed_Table, repmat(Merged_Table(:,k),1,Count(k))];
    k
end
% histogram(Completed_Table(1,:))
% Ruler=[quantile(Completed_Table,0.2,2), quantile(Completed_Table,0.4,2), quantile(Completed_Table,0.6,2), quantile(Completed_Table,0.8,2)];
load('Package/Ruler_PPI.mat')

% MEAN=mean(Completed_Table,2);
% STD=std(Completed_Table,[],2);
% Ruler=[MEAN-3*STD, MEAN-2.5*STD, MEAN-2*STD, MEAN-1.5*STD, MEAN-1*STD, MEAN-0.5*STD, MEAN, MEAN+0.5*STD, MEAN+1*STD, MEAN+1.5*STD, MEAN+2*STD, MEAN+2.5*STD, MEAN+3*STD];
% Ruler=[MEAN-3*STD, MEAN-2*STD, MEAN-1*STD, MEAN, MEAN+1*STD, MEAN+2*STD, MEAN+3*STD];

Hist=zeros(size(Merged_Table));
for m=1:size(Merged_Table,2)
    Det=Merged_Table(:,m)<Ruler;
    for n=1:size(Merged_Table,1)
        if isempty(find(Det(n,:),1))
            Hist(n,m)=size(Ruler,2)+1;
        else
            Hist(n,m)=find(Det(n,:),1);
        end
    end
    m
end

Selected_Sample(1).Merged_AAIdx=[]; Acc=0;
for p=1:size(Selected_Sample,2)
    if ~isempty(Selected_Sample(p).Patch_code)
        Acc_Start=Acc+1;
        Acc=Acc+size(Selected_Sample(p).Patch_code,1);
        Selected_Sample(p).Merged_AAIdx=[Selected_Sample(p).Patch_code(:,2), Selected_Sample(p).Patch_code(:,1), Hist(:,Acc_Start:Acc)'];
    end
    p
end

Merged_AAIdx=Selected_Sample;
% save('Merged_AAIdx.mat','Merged_AAIdx')

% Patch Vectorlization
Max_Cluster_No=300; Bin=size(Ruler,2)+1;

clear AA_Idx_Patch_code;
AA_Idx_Patch_code(1).ID=[];
AA_Idx_Patch_code(1).Patch_code=[];

for u=1:158
    for v=1:size(Merged_AAIdx,2)
        if ~isempty(Merged_AAIdx(v).Patch_code)
            Data=Merged_AAIdx(v).Merged_AAIdx(:,[1,2,u+2]);
            Patch_code=zeros(Max_Cluster_No, Bin);
            for x=1:size(Data,1)
                Patch_code(Data(x,2), Data(x,3))=Patch_code(Data(x,2), Data(x,3))+Data(x,1);
            end
            AA_Idx_Patch_code(v).ID=Merged_AAIdx(v).ID;
            AA_Idx_Patch_code(v).Patch_code=sum(Patch_code,2);
            AA_Idx_Patch_code(v).(['Patch_code_', num2str(u)])=Patch_code;
        end
    end
    u
end
save('Reduced_Sample_Int_Res_Results/AA_Idx_Patch_code.mat','AA_Idx_Patch_code')

%
%
%

%%%% AA_Idx Reduction
% 1-R2 as Similarity Distance 
load('Package/AA_Index.mat')
R2=corrcoef(AA_Index').^2;

% Hierachical Clustering (R2<0.3 as Cutoff)
Z=linkage(squareform(1-R2),'average');
T=cluster(Z,'cutoff',0.3,'criterion','distance');

% Representative AA_Idx Feature
Representative_Clust=zeros(max(T),1);
for t=1:max(T)
    Ind=find(T==t);
    if size(Ind,1)==1
        Representative_Clust(t,1)=Ind;
    elseif size(Ind,1)==2
        Ind1=sort(R2(:,Ind(1)),'descend');
        Ind2=sort(R2(:,Ind(2)),'descend');
        if Ind1(3)>Ind2(3)
            Representative_Clust(t,1)=Ind(2);
        else
            Representative_Clust(t,1)=Ind(1);
        end
    else
         [~,c]=min(sum(R2(Ind, Ind)));
         Representative_Clust(t,1)=Ind(c);
    end
end

%%%% AA_Idx_Patch_code Table Rearrangement
load ('Reduced_Sample_Int_Res_Results/AA_Idx_Patch_code.mat')

clear AA_Idx_Patch_code_RA;
AA_Idx_Patch_code_RA(1).ID=[];
AA_Idx_Patch_code_RA(1).Patch_code=[];

for u=1:size(AA_Idx_Patch_code,2)
    AA_Idx_Patch_code_RA(u).ID=AA_Idx_Patch_code(u).ID;
    AA_Idx_Patch_code_RA(u).Patch_code=AA_Idx_Patch_code(u).Patch_code;
    for v=1:max(T)
        AA_Idx_Patch_code_RA(u).(['Patch_code_', num2str(v)])=AA_Idx_Patch_code(u).(['Patch_code_', num2str(Representative_Clust(v))]);
    end
    u
end
save('Reduced_Sample_Int_Res_Results/AAI_PatchBag_PPI.mat','AA_Idx_Patch_code_RA')

%
%
%

%%%% PSI

%%%% Patch to Vector
% Load Qualified Reduced Sample ID
% Filename='Reduced_Sample_Training_Test_Dataset.xlsx';
Filename='target_entry.txt';
ID_List=fopen(Filename);
sp=1;Sample_Info(1).ID=[];
while (1)
    line=fgetl(ID_List);
    if line==-1, break, end
    Sample_Info(sp).ID=sscanf(line,'%c');
    sp=sp+1;
end
fclose(ID_List);

% Load Merged Error ID
% load('Merged_Error_ID.mat')
% Sample_Info=Sample_Info(setdiff([1:2472], Merged_Error_ID));

% Load Clust Info (Representative Patch Library)
Max_Cluster_No=300;
load(['Package/Clust_medoids_', num2str(Max_Cluster_No), '_PSI.mat'])

% Patch Classification
Selected_Sample=Sample_Info;
Selected_Sample(1).Patch_code=[];
Selected_Sample(1).Patch_Res=[];
Selected_Sample(1).Bar_code=[];

Blank_ID=[];
P=perms(1:6);
for i=1:size(Selected_Sample,1)
    try
        % Patch Files
        PDB_ID=lower(Selected_Sample(i).ID);
        Patch_File=table2array(readtable(['Reduced_Sample_nInt_Res_Results/Patch_Results_nInt/',PDB_ID,'_Patch.txt']));
        [Unique_Patch, ~, Uidx]=unique(Patch_File, 'rows', 'stable');

        % Patch Codes
        Patch_code=zeros(size(Unique_Patch,1),1);

        parfor j=1:size(Unique_Patch,1)
            P1=reshape(Unique_Patch(j,1:18),6,3);
            RMSD_DISM=zeros(Max_Cluster_No, 1);

            for m=1:Max_Cluster_No
                RMSD_Table=zeros(size(P,1), 1);
                Angle_Table=zeros(size(P,1), 1);

                for r=1:size(P,1)
                    P2=Clust_medoids(m).Coord(P(r,:),:);

                    [R,~,eRMSD] = CoordiExam_AC(P1, P2);
                    RMSD_Table(r,1)=eRMSD;
                    Angle_Table(r,1)=acosd(dot(Unique_Patch(j,19:21)*R, Clust_medoids(m).Norm_V));
                end

                if isempty(RMSD_Table(Angle_Table<90))
                    RMSD_DISM(m,1)=min(RMSD_Table,[],'all');
                else
                    RMSD_DISM(m,1)=min(RMSD_Table(Angle_Table<90),[],'all');
                end
            end
            [~,c]=min(RMSD_DISM);
            Patch_code(j,1)=c;
        end

        Unique_Count=hist(Uidx, (1:max(Uidx))');
        Selected_Sample(i).Patch_code=[Patch_code, Unique_Count'];
        Selected_Sample(i).Patch_Res=Unique_Patch(:,22:27);

        % Patch Vectorlization
        Bar_code=zeros(Max_Cluster_No,1);
        for u=1:Max_Cluster_No
            Bar_code(u,1)= sum(Unique_Count(Patch_code==u));
        end
        Selected_Sample(i).Bar_code=Bar_code;
    catch
        Blank_ID=[Blank_ID; i];
    end
    i
end

% save('Blank_ID.mat', 'Blank_ID');
save('Reduced_Sample_nInt_Res_Results/Selected_Sample.mat', 'Selected_Sample');

%
%
%

load('Package/AA_Index.mat')
load('Reduced_Sample_nInt_Res_Results/Selected_Sample.mat')

% Selected_Sample.Bar_code -> Spatial Distance
Count=[];
Merged_Table=[];
Selected_Sample=Selected_Sample';

parfor i=1:size(Selected_Sample,2)
    Patch_Res=Selected_Sample(i).Patch_Res;
    if ~isempty(Patch_Res)
        for j=1:size(Patch_Res,1)
            Idx=Patch_Res(j,:);
            Idx(Idx==0)=8;
            Merged_Table=[Merged_Table, mean(AA_Index(:,Idx),2)];
        end
        Count=[Count; Selected_Sample(i).Patch_code(:,2)];
    end
    i
end

Completed_Table=[];
parfor k=1:size(Count,1)
    Completed_Table=[Completed_Table, repmat(Merged_Table(:,k),1,Count(k))];
    k
end
% histogram(Completed_Table(1,:))
% Ruler=[quantile(Completed_Table,0.2,2), quantile(Completed_Table,0.4,2), quantile(Completed_Table,0.6,2), quantile(Completed_Table,0.8,2)];
load('Package/Ruler_PSI.mat')

% MEAN=mean(Completed_Table,2);
% STD=std(Completed_Table,[],2);
% Ruler=[MEAN-3*STD, MEAN-2.5*STD, MEAN-2*STD, MEAN-1.5*STD, MEAN-1*STD, MEAN-0.5*STD, MEAN, MEAN+0.5*STD, MEAN+1*STD, MEAN+1.5*STD, MEAN+2*STD, MEAN+2.5*STD, MEAN+3*STD];
% Ruler=[MEAN-3*STD, MEAN-2*STD, MEAN-1*STD, MEAN, MEAN+1*STD, MEAN+2*STD, MEAN+3*STD];

Hist=zeros(size(Merged_Table));
for m=1:size(Merged_Table,2)
    Det=Merged_Table(:,m)<Ruler;
    for n=1:size(Merged_Table,1)
        if isempty(find(Det(n,:),1))
            Hist(n,m)=size(Ruler,2)+1;
        else
            Hist(n,m)=find(Det(n,:),1);
        end
    end
    m
end

Selected_Sample(1).Merged_AAIdx=[]; Acc=0;
for p=1:size(Selected_Sample,2)
    if ~isempty(Selected_Sample(p).Patch_code)
        Acc_Start=Acc+1;
        Acc=Acc+size(Selected_Sample(p).Patch_code,1);
        Selected_Sample(p).Merged_AAIdx=[Selected_Sample(p).Patch_code(:,2), Selected_Sample(p).Patch_code(:,1), Hist(:,Acc_Start:Acc)'];
    end
    p
end

Merged_AAIdx=Selected_Sample;
% save('Merged_AAIdx.mat','Merged_AAIdx')

% Patch Vectorlization
Max_Cluster_No=300; Bin=size(Ruler,2)+1;

clear AA_Idx_Patch_code;
AA_Idx_Patch_code(1).ID=[];
AA_Idx_Patch_code(1).Patch_code=[];

for u=1:158
    for v=1:size(Merged_AAIdx,2)
        if ~isempty(Merged_AAIdx(v).Patch_code)
            Data=Merged_AAIdx(v).Merged_AAIdx(:,[1,2,u+2]);
            Patch_code=zeros(Max_Cluster_No, Bin);
            for x=1:size(Data,1)
                Patch_code(Data(x,2), Data(x,3))=Patch_code(Data(x,2), Data(x,3))+Data(x,1);
            end
            AA_Idx_Patch_code(v).ID=Merged_AAIdx(v).ID;
            AA_Idx_Patch_code(v).Patch_code=sum(Patch_code,2);
            AA_Idx_Patch_code(v).(['Patch_code_', num2str(u)])=Patch_code;
        end
    end
    u
end
save('Reduced_Sample_nInt_Res_Results/AA_Idx_Patch_code.mat','AA_Idx_Patch_code')

%
%
%

%%%% AA_Idx Reduction
% 1-R2 as Similarity Distance 
load('Package/AA_Index.mat')
R2=corrcoef(AA_Index').^2;

% Hierachical Clustering (R2<0.3 as Cutoff)
Z=linkage(squareform(1-R2),'average');
T=cluster(Z,'cutoff',0.3,'criterion','distance');

% Representative AA_Idx Feature
Representative_Clust=zeros(max(T),1);
for t=1:max(T)
    Ind=find(T==t);
    if size(Ind,1)==1
        Representative_Clust(t,1)=Ind;
    elseif size(Ind,1)==2
        Ind1=sort(R2(:,Ind(1)),'descend');
        Ind2=sort(R2(:,Ind(2)),'descend');
        if Ind1(3)>Ind2(3)
            Representative_Clust(t,1)=Ind(2);
        else
            Representative_Clust(t,1)=Ind(1);
        end
    else
         [~,c]=min(sum(R2(Ind, Ind)));
         Representative_Clust(t,1)=Ind(c);
    end
end

%%%% AA_Idx_Patch_code Table Rearrangement
load ('Reduced_Sample_nInt_Res_Results/AA_Idx_Patch_code.mat')

clear AA_Idx_Patch_code_RA;
AA_Idx_Patch_code_RA(1).ID=[];
AA_Idx_Patch_code_RA(1).Patch_code=[];

for u=1:size(AA_Idx_Patch_code,2)
    AA_Idx_Patch_code_RA(u).ID=AA_Idx_Patch_code(u).ID;
    AA_Idx_Patch_code_RA(u).Patch_code=AA_Idx_Patch_code(u).Patch_code;
    for v=1:max(T)
        AA_Idx_Patch_code_RA(u).(['Patch_code_', num2str(v)])=AA_Idx_Patch_code(u).(['Patch_code_', num2str(Representative_Clust(v))]);
    end
    u
end
save('Reduced_Sample_nInt_Res_Results/AAI_PatchBag_PSI.mat','AA_Idx_Patch_code_RA')

%
%
%

clear all

PPI=load('Reduced_Sample_Int_Res_Results/AAI_PatchBag_PPI.mat');
PSI=load('Reduced_Sample_nInt_Res_Results/AAI_PatchBag_PSI.mat');

AAI_PatchBag.PPI=PPI.AA_Idx_Patch_code_RA;
AAI_PatchBag.PSI=PSI.AA_Idx_Patch_code_RA;

save('AAI_PatchBag.mat','AAI_PatchBag')