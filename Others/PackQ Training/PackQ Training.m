%%%% Merge of Dock Quality Assessment Results
% Filename of Dock Quality Assessment Results
Filename_List = {'Merged_Result_P41212_MXRGA.mat','Merged_Result_P41212_SurfOrca.mat','Merged_Result_P41212_Tetramer.mat',...
                 'Merged_Result_P43212_MXRGA.mat','Merged_Result_P43212_SurfOrca.mat','Merged_Result_P43212_Tetramer.mat'};

% Merge of Dock Quality Assessment Results
M_Fnat=[];
M_iRMS=[];
M_RMSD_95=[];

for f=1:6
    % Load Corresponding Files
    load(Filename_List{f})

    % Merge Results
    M_Fnat=[M_Fnat, [Merged_Result.Fnat]];          %#ok<AGROW> 
    M_iRMS=[M_iRMS, [Merged_Result.iRMS]];          %#ok<AGROW> 
    M_RMSD_95=[M_RMSD_95, [Merged_Result.RMSD_95]]; %#ok<AGROW> 
end

% Combination (Remove NaN Values)
Merged_Pmt=[M_Fnat; M_iRMS; M_RMSD_95]';
Reduced_Pmt=Merged_Pmt(~isnan(Merged_Pmt(:,1)) & ~isnan(Merged_Pmt(:,2)) & ~isnan(Merged_Pmt(:,3)),:);

%
%
%

%%%% Pack Quality Claasification (Class Refined by KJ)
for r=1:size(Reduced_Pmt,1)
    
    Fnat=Reduced_Pmt(r,1);
    iRMS=Reduced_Pmt(r,2);
    RMSD_95=Reduced_Pmt(r,3);

    if Fnat < 0.3 || (RMSD_95 > 8.0 && iRMS > 4.0)
        % Incorrect
        Reduced_Pmt(r,4)=0;
    elseif (Fnat >= 0.3 && Fnat < 0.5) && (RMSD_95 <= 8.0 || iRMS <= 4.0) || (Fnat >= 0.5 && RMSD_95 > 5.0 && iRMS > 2.0)
        % Acceptable
        Reduced_Pmt(r,4)=1;
    elseif (Fnat >= 0.5 && Fnat < 0.7) && (RMSD_95 <= 5.0 || iRMS <= 2.0) || (Fnat >= 0.7 && RMSD_95 > 3.0 && iRMS > 1.0)
        % Medium
        Reduced_Pmt(r,4)=2;
    elseif Fnat >= 0.7 && (RMSD_95 <= 3.0 || iRMS <= 1.0)
        % High
        Reduced_Pmt(r,4)=3;
    else
        % Undef
        Reduced_Pmt(r,4)=-1;
    end
end

% find(Reduced_Pmt(:,4)==-1)
% Please Check if there is Undefiened Cases (-1)

%
%
%

s=RandStream('mrg32k3a');
s=RandStream('mt19937ar');
s=RandStream('mlfg6331_64');

%%%% Division of Training/Testing Sets
% Divide into Training/Testing Sets
Training_Idx=randsample(s, 1:size(Reduced_Pmt,1), floor(size(Reduced_Pmt,1)*0.7));
Testing_Idx=setdiff(1:size(Reduced_Pmt,1), Training_Idx,'stable');

% Training Sets
Training_Feature=Reduced_Pmt(Training_Idx,:);

% Testing Sets
Testing_Feature=Reduced_Pmt(Testing_Idx,:);

%%%% Training of PackQ Pmts
% Rough: Determine Cutoff to Seperate Incorrect/Acceptable(C1) & Acceptable/Medium(C2) & Medium/High(C3)
t=0; Pmt=zeros(5000000,12);
for C1=0.1:0.03:0.5
    for C2=0.4:0.03:0.7
        for C3=0.6:0.03:0.95
            C = [C1, C2, C3];
            for d1=0.5:0.5:5
                for d2=0.5:0.5:10

                    % Record of PackQ and Class
                    Table=nan(size(Training_Feature,1),2);
                    for r=1:size(Training_Feature,1)
                        Fnat=Training_Feature(r,1);
                        iRMS=Training_Feature(r,2);
                        RMSD_95=Training_Feature(r,3);
                        Table(r,1) = (Fnat + 1/(1+(iRMS/d1)^2) + 1/(1+(RMSD_95/d2)^2))/3;   % PackQ
                        Table(r,2) = Training_Feature(r,4);                                 % Class
                    end

                    t=t+1;
                    Pmt(t,1)=d1;    % d1
                    Pmt(t,2)=d2;    % d2
                    Pmt(t,3)=C1;    % Incorrect/Acceptable(C1)
                    Pmt(t,4)=C2;    % Acceptable/Medium(C2)
                    Pmt(t,5)=C3;    % Medium/High(C3)

                    for Class_Value=1:3
                        Reduced_Table=Table((Table(:,2)==Class_Value|Table(:,2)==(Class_Value-1)),:);
                        TP= sum(Reduced_Table(:,1)>=C(Class_Value) & Reduced_Table(:,2)==Class_Value);
                        FN= sum(Reduced_Table(:,1) <C(Class_Value) & Reduced_Table(:,2)==Class_Value);
                        TN= sum(Reduced_Table(:,1) <C(Class_Value) & Reduced_Table(:,2)==(Class_Value-1));
                        FP= sum(Reduced_Table(:,1)>=C(Class_Value) & Reduced_Table(:,2)==(Class_Value-1));
                        
                        PPV=TP/(TP+FP);                             % Precision
                        TPR=TP/(TP+FN);                             % Recall
                        Pmt(t,Class_Value+5)=2*PPV*TPR/(PPV+TPR);   % F1 Score
                        Pmt(t,Class_Value+8)=sum(Training_Feature(:,4)==Class_Value|Training_Feature(:,4)==(Class_Value-1));    % Count
                    end

                    Pmt(t,12)=(Pmt(t,6)*Pmt(t,9)+ Pmt(t,7)*Pmt(t,10)+ Pmt(t,8)*Pmt(t,11))/sum(Pmt(t,9:11));                      % Weighted F1 Score
                end
            end
        end
    end
    C1 %#ok<NOPTS> 
end
Pmt(t+1,:)=[];

[~,row]=max(Pmt(:,12));
Pmt(row,[1:8,12]) %#ok<NOPTS> 

% mrg32k3a    : 1.0000    9.5000    0.3400    0.4900    0.6300    0.7344    0.7901    0.8571    0.7383
% mt19937ar   : 1.0000    9.5000    0.3400    0.4900    0.6300    0.7315    0.8152    0.7619    0.7363
% mlfg6331_64 : 1.0000    9.5000    0.3400    0.4900    0.6300    0.7342    0.8092    0.8148    0.7389               
% Overall     : 1.0000    9.5000    0.3400    0.4900    0.6300    0.7360    0.8000    0.8571    0.7977

%

% Detail: Determine Cutoff to Seperate Incorrect/Acceptable(C1) & Acceptable/Medium(C2) & Medium/High(C3)
C1_int=Pmt(row,3);
C2_int=Pmt(row,4);
C3_int=Pmt(row,5);

t=0; Pmt=zeros(5000000,12);
for C1=C1_int-0.05:0.01:C1_int+0.05
    for C2=C2_int-0.05:0.01:C2_int+0.05
        for C3=C3_int-0.05:0.01:C3_int+0.05
            C = [C1, C2, C3];
            for d1=0.5:0.5:5
                for d2=0.5:0.5:10

                    % Record of PackQ and Class
                    Table=nan(size(Training_Feature,1),2);
                    for r=1:size(Training_Feature,1)
                        Fnat=Training_Feature(r,1);
                        iRMS=Training_Feature(r,2);
                        RMSD_95=Training_Feature(r,3);
                        Table(r,1) = (Fnat + 1/(1+(iRMS/d1)^2) + 1/(1+(RMSD_95/d2)^2))/3;   % PackQ
                        Table(r,2) = Training_Feature(r,4);                                 % Class
                    end

                    t=t+1;
                    Pmt(t,1)=d1;    % d1
                    Pmt(t,2)=d2;    % d2
                    Pmt(t,3)=C1;    % Incorrect/Acceptable(C1)
                    Pmt(t,4)=C2;    % Acceptable/Medium(C2)
                    Pmt(t,5)=C3;    % Medium/High(C3)

                    for Class_Value=1:3
                        Reduced_Table=Table((Table(:,2)==Class_Value|Table(:,2)==(Class_Value-1)),:);
                        TP= sum(Reduced_Table(:,1)>=C(Class_Value) & Reduced_Table(:,2)==Class_Value);
                        FN= sum(Reduced_Table(:,1) <C(Class_Value) & Reduced_Table(:,2)==Class_Value);
                        TN= sum(Reduced_Table(:,1) <C(Class_Value) & Reduced_Table(:,2)==(Class_Value-1));
                        FP= sum(Reduced_Table(:,1)>=C(Class_Value) & Reduced_Table(:,2)==(Class_Value-1));
                        
                        PPV=TP/(TP+FP);                             % Precision
                        TPR=TP/(TP+FN);                             % Recall
                        Pmt(t,Class_Value+5)=2*PPV*TPR/(PPV+TPR);   % F1 Score
                        Pmt(t,Class_Value+8)=sum(Training_Feature(:,4)==Class_Value|Training_Feature(:,4)==(Class_Value-1));    % Count
                    end

                    Pmt(t,12)=(Pmt(t,6)*Pmt(t,9)+ Pmt(t,7)*Pmt(t,10)+ Pmt(t,8)*Pmt(t,11))/sum(Pmt(t,9:11));                      % Weighted F1 Score
                end
            end
        end
    end
    C1 %#ok<NOPTS> 
end
Pmt(t+1,:)=[];

[~,row]=max(Pmt(:,12));
Pmt(row,[1:8,12]) %#ok<NOPTS> 

% mrg32k3a    : 1.0000    9.5000    0.3400    0.5000    0.6500    0.7344    0.7922    0.8889    0.7387
% mt19937ar   : 1.0000    9.5000    0.3400    0.4900    0.6300    0.7315    0.8152    0.7619    0.7363
% mlfg6331_64 : 1.5000   10.0000    0.3600    0.5400    0.6700    0.7370    0.8050    0.8276    0.7415                
% Overall     : 1.0000    9.0000    0.3300    0.4800    0.6200    0.7317    0.7984    0.8649    0.7983

%
%
%

%%%% Performance Assessment (Testing_Set)
best_d1=1.0; best_d2=9.5; C=[0.34, 0.50, 0.65]; % mrg32k3a 
best_d1=1.0; best_d2=9.5; C=[0.34, 0.49, 0.63]; % mt19937ar
best_d1=1.5; best_d2=10.0;C=[0.36, 0.54, 0.67]; % mlfg6331_64

% Record of PackQ and Class
Table=nan(size(Testing_Feature,1),2);
for r=1:size(Testing_Feature,1)
    Fnat=Testing_Feature(r,1);
    iRMS=Testing_Feature(r,2);
    RMSD_95=Testing_Feature(r,3);
    Table(r,1) = (Fnat + 1/(1+(iRMS/best_d1)^2) + 1/(1+(RMSD_95/best_d2)^2))/3;     % PackQ
    Table(r,2) = Testing_Feature(r,4);                                              % Class
end

F1_Score=zeros(1,3);
Count=zeros(1,3);
for Class_Value=1:3
    Reduced_Table=Table((Table(:,2)==Class_Value|Table(:,2)==(Class_Value-1)),:);
    TP= sum(Reduced_Table(:,1)>=C(Class_Value) & Reduced_Table(:,2)==Class_Value);
    FN= sum(Reduced_Table(:,1) <C(Class_Value) & Reduced_Table(:,2)==Class_Value);
    TN= sum(Reduced_Table(:,1) <C(Class_Value) & Reduced_Table(:,2)==(Class_Value-1));
    FP= sum(Reduced_Table(:,1)>=C(Class_Value) & Reduced_Table(:,2)==(Class_Value-1));

    PPV=TP/(TP+FP);                                % Precision
    TPR=TP/(TP+FN);                                % Recall
    F1_Score(1,Class_Value)=2*PPV*TPR/(PPV+TPR);   % F1 Score
    Count(1,Class_Value)=sum(Testing_Feature(:,4)==Class_Value|Testing_Feature(:,4)==(Class_Value-1));
end

Average_F1_Score=F1_Score*Count'/sum(Count);

% mrg32k3a   : 0.7395   :   0.7398    0.7606    0.6000
% mt19937ar  : 0.7484   :   0.7463    0.7500    1.0000
% mlfg6331_64: 0.7277   :   0.7261    0.7302    0.8889

%
%
%

%%%% Report
% Choose mlfg6331_64
% best_d1=1.5; best_d2=10.0; C=[0.36, 0.54, 0.67];
% mrg32k3a   : 0.7367 (0.7325    0.7974    0.8333) / 0.7390 (0.7368    0.7536    0.8571)
% mt19937ar  : 0.7335 (0.7298    0.7953    0.7500) / 0.7451 (0.7431    0.7451    1.0000)
% mlfg6331_64: 0.7415 (0.7370    0.8050    0.8276) / 0.7277 (0.7261    0.7302    0.8889)
% OVerall    : 0.7373 (0.7338    0.7838    0.8421)

% Pmts Setting
best_d1=1.5;
best_d2=10.0; 
C=[0.36, 0.54, 0.67];

% Load Reduced Pmt
load('Reduced_Pmt.mat')
s=RandStream('mlfg6331_64');

% Divide into Training/Testing Sets
Training_Idx=randsample(s, 1:size(Reduced_Pmt,1), floor(size(Reduced_Pmt,1)*0.7));
Testing_Idx=setdiff(1:size(Reduced_Pmt,1), Training_Idx,'stable');

Training_Feature=Reduced_Pmt(Training_Idx,:);
Testing_Feature=Reduced_Pmt(Testing_Idx,:);

%%% Boxplot
PackQ = (Training_Feature(:,1) + 1./(1+(Training_Feature(:,2)/best_d1).^2) + 1./(1+(Training_Feature(:,3)/best_d2).^2))/3;
Class = Training_Feature(:,4);

P0=PackQ(Class==0); g0=ones(size(P0));
P1=PackQ(Class==1); g1=2*ones(size(P1));
P2=PackQ(Class==2); g2=3*ones(size(P2));
P3=PackQ(Class==3); g3=4*ones(size(P3));

boxplot([P0;P1;P2;P3],[g0;g1;g2;g3]);

%

PackQ = (Testing_Feature(:,1) + 1./(1+(Testing_Feature(:,2)/best_d1).^2) + 1./(1+(Testing_Feature(:,3)/best_d2).^2))/3;
Class = Testing_Feature(:,4);

P0=PackQ(Class==0); G0=ones(size(P0));
P1=PackQ(Class==1); G1=2*ones(size(P1));
P2=PackQ(Class==2); G2=3*ones(size(P2));
P3=PackQ(Class==3); G3=4*ones(size(P3));

Data=[P0;P1;P2;P3];
% boxplot(Data,[G0;G1;G2;G3]);

%%% Violin Plot
% Data & Group
Group=[];
for g=0:3
    Group=[Group; repmat({['G', num2str(g)]},sum(Class==g),1)];
end

% Violin Plot
GroupList={'G0','G1','G2','G3'};
Colorlist = [0         0.4470    0.7410;
             0.8500    0.3250    0.0980;
             0.9290    0.6940    0.1250;
             0.4940    0.1840    0.5560];

for Ind=1:length(GroupList)
    vs = Violin({Data(strcmp(Group, GroupList{Ind}))},...
        Ind+0.5,...
        'EdgeColor', [0.65 0.65 0.65],...
        'HalfViolin','full',...           % right,left, full
        'QuartileStyle','boxplot',...     % shadow, boxplot, none
        'DataStyle', 'none',...           % histogram, scatter, none
        'ShowData', true,...
        'ShowBox',  true,...
        'ShowMean', false,...
        'ShowMedian',  true,...
        'ShowNotches', false,...
        'ViolinColor', {Colorlist(Ind,:)});;
end

xlim([1,5])
xticks([1.5:1:5])
xticklabels({'Incorrect', 'Acceptable', 'Medium', 'High'})

ylim([0,1])
yticks([0:0.1:1])

% C = [0.36, 0.54, 0.67]
% Median: 0.1829  0.4068  0.5799  0.7303
% Mean:   0.1939  0.4179  0.5865  0.7499
% p-value_Incorrect_vs_Acceptable: 0
% p-value_Acceptable_vs_Medium: 1.8276e-83
% p-value_Medium_vs_High: 3.3870e-15