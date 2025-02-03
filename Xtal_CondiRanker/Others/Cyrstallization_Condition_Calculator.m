%%%% Transform Reagent Condition Table into Component Condition Table
PDB_ID='1lyy';

% Load AAI_PatchBag_Ranking
load(['AAI_PatchBag_Ranking_',PDB_ID,'.mat'])
Input=readtable(['Input_',PDB_ID,'.xlsx'],'VariableNamingRule','preserve');

% Load Database_Sample_Condition_2472
load('Package/Sample_Condition_2472.mat')

Idx=[];
for s=1:size(AAI_PatchBag_Ranking,1)
    Idx=[Idx, find(strcmp(Sample_Condition_2472.ID,AAI_PatchBag_Ranking{s}))];
end
Condition_Table=table2array([Input(1,12:end); Sample_Condition_2472(Idx,12:end)]);
Qualified_Sample_ID=table2struct([Input(1,1); Sample_Condition_2472(Idx,1)]);

% Load Decomposition_Table
Decom_Table=readtable('Package/Decomposition_Table(Modified).csv');
Component_list=readcell('Package/Unique_Component.csv');
Salt_Idx=cell2mat(Component_list(:,2));

% Decomposition Transformation
Result=zeros(size(Condition_Table,1), size(Component_list,1));
for Idx_r=1:size(Condition_Table,1)
    Index=find(Condition_Table(Idx_r,:));
    for i=1:size(Index,2)
        for j=1:4
            Component=table2cell(Decom_Table(Index(i),j+1));
            if ~isempty(Component{:})
                [row, ~]=find(ismember(Component_list(:,1), Component));
                Result(Idx_r,row)=sum([Result(Idx_r,row), Condition_Table(Idx_r,Index(i))*table2array(Decom_Table(Index(i),j+5))]);
            end
        end
    end
end
clearvars -except Input Sample_Condition_2472 Qualified_Sample_ID Salt_Idx Result

%

%%%% Scaling of Component Condition Table
% Load the Component_Scaling_Factor.mat for following Robust Scaling
load('Package/Component_Scaling_Factor.mat')

% Separate Additives from Principal Components
Total_Size=length(Component_Scaling_Factor)+sum(Component_Scaling_Factor(:,1)>0);
Scaled_Result=zeros(size(Result,1), Total_Size);

% Robust Scaling
Median_Divergence_Ratio=[Component_Scaling_Factor(:,2)./Component_Scaling_Factor(:,3); Component_Scaling_Factor(:,4)./Component_Scaling_Factor(:,5)];
Dist=Median_Divergence_Ratio(logical(1-(isinf(Median_Divergence_Ratio)+isnan(Median_Divergence_Ratio))));
[f,xi]=ksdensity(Dist, 'Bandwidth', std(Dist)/2);

Cali_Value=median(Dist);

t=0; 
for k=1:size(Result,2)
    if ~isnan(Component_Scaling_Factor(k,2))
        if Component_Scaling_Factor(k,1)>0
            if Component_Scaling_Factor(k,3)>0
                Scaled_Result(:,k+t)=(Result(:,k)-Component_Scaling_Factor(k,2))/Component_Scaling_Factor(k,3);
            else
                Scaled_Result(:,k+t)=(Result(:,k)-Component_Scaling_Factor(k,2))/(Component_Scaling_Factor(k,2)/Cali_Value);
            end

            t=t+1;
            if Component_Scaling_Factor(k,5)>0
                Scaled_Result(:,k+t)=(Result(:,k)-Component_Scaling_Factor(k,4))/Component_Scaling_Factor(k,5);
            else
                Scaled_Result(:,k+t)=(Result(:,k)-Component_Scaling_Factor(k,4))/(Component_Scaling_Factor(k,4)/Cali_Value);
            end       
        else
            if Component_Scaling_Factor(k,3)>0
                Scaled_Result(:,k+t)=(Result(:,k)-Component_Scaling_Factor(k,2))/Component_Scaling_Factor(k,3);
            else
                Scaled_Result(:,k+t)=(Result(:,k)-Component_Scaling_Factor(k,2))/(Component_Scaling_Factor(k,2)/Cali_Value);
            end
        end
    end
end

%

%%%% Environmental_factor


Environmental_factor=table2array([Input(1,10:11); Sample_Condition_2472(:,10:11)]);
Scaled_Environmental_factor=zeros(size(Environmental_factor,1),2);

% Temprature
AVG_Temp=median(nonzeros(Environmental_factor(:,1)));
Environmental_factor(Environmental_factor(:,1)==0,1)=AVG_Temp;
Environmental_factor(Environmental_factor(:,1)>320)=AVG_Temp;
Environmental_factor(Environmental_factor(:,1)<=273)=273+4;
Scaled_Environmental_factor(Environmental_factor(:,1)>(273.15+10),1)=Cali_Value;

% pH
AVG_pH=median(nonzeros(Environmental_factor(:,2)));
Scaled_pH=iqr(nonzeros(Environmental_factor(:,2)));
Environmental_factor(Environmental_factor(:,2)==0,2)=AVG_pH;
Scaled_Environmental_factor(:,2)=(Environmental_factor(:,2)-AVG_pH)/Scaled_pH;

% Scaled_Table
RS_Table=[Scaled_Environmental_factor, Scaled_Result];

% Remove w/o Unused Components
Used_Idx=sum(RS_Table-mode(RS_Table))~=0;
Reduced_RS_Table=RS_Table(:,Used_Idx);

%
%
%

%%%% Weighted-Euclidean Distance
% Load Component Similarity Table (Completed_AP_DISM)
Completed_AP_DISM=readcell('Package/Completed_AP_DISM(Modified).xlsx');
DISM=cell2mat(Completed_AP_DISM(2:end,2:end));
clear Completed_AP_DISM;

% PC Index
SR_Idx=Component_Scaling_Factor(:,[2,4])';
SR_Idx(SR_Idx~=0)=1:sum(SR_Idx~=0,'all');

t=0;
DISM_Idx=zeros(1, Total_Size);
for k=1:size(SR_Idx,2)
    if SR_Idx(2,k)~=0
        t=t+1;
        DISM_Idx(1,[k+t-1, k+t])=k;
    end
    DISM_Idx(1,k+t)=k;
end
DISM_En_Idx=[0, 0, DISM_Idx];
Reduced_DISM_Idx=DISM_En_Idx(Used_Idx);

Salt_Cali_M=eye(size(DISM_En_Idx,2));
Scaled_Salt_Idx=Salt_Idx(DISM_Idx);
for ssi=1:size(Scaled_Salt_Idx,1)
    % Salt
    if Scaled_Salt_Idx(ssi)==1
        Salt_Cali_M(ssi+2,ssi+2)=0.5;
    % Buffer
    elseif Scaled_Salt_Idx(ssi)==2
        Salt_Cali_M(ssi+2,ssi+2)=0.3;
    % Organic
    else
        Salt_Cali_M(ssi+2,ssi+2)=1;
    end
end
Reduced_Salt_Cali_M=Salt_Cali_M(Used_Idx, Used_Idx);

%

%%%% Weighted-Euclidean Distance
% Preallocate Memory for Condition DISM Storage & Reduced DISM Index
[row, col]=size(Reduced_RS_Table);

Condition_DISM=1:row*(row-1)/2;
Reduced_DISM_Ind=tril(ones(row),-1);
Reduced_DISM_Ind(Reduced_DISM_Ind==1)=Condition_DISM;

% Weighted-Euclidean Distance
[Ridx, Cidx]=find(Reduced_DISM_Ind);

[Count, UQ_DISM_Idx]=hist(Reduced_DISM_Idx, unique(Reduced_DISM_Idx));
Duplicated_Idx=setdiff(UQ_DISM_Idx(Count==2),0);
parfor m=1:length(Condition_DISM)
    V=Reduced_RS_Table(Ridx(m),:)-Reduced_RS_Table(Cidx(m),:);
    for r=Duplicated_Idx
        Duplicated_Value=V(Reduced_DISM_Idx==r);
        if Duplicated_Value(1)~=0 && Duplicated_Value(2)~=0
            [~,maxow]=max(abs(Duplicated_Value));
            Vidx=find(Reduced_DISM_Idx==r);
            V(Vidx(maxow))=0;
        end
    end
    Condition_DISM(m)=sqrt(V*Reduced_Salt_Cali_M*V');
    m
end

save('Condition_DISM_2472.mat','Condition_DISM')
A=squareform(Condition_DISM);
plot(A(1,:),'.')
ylim([0 30])