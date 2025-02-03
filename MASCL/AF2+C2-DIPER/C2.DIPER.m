%%%% Parallel Computation
% Activate the Parallel Cores (default: 8 cores)
parpool(8)

%%%% Default Settings
% Load PDB ID
Filename='target_entry.txt';
List=fopen(Filename);
sp=1;Sample(1).ID=[];
while (1)
    line=fgetl(List);
    if line==-1, break, end
    Sample(sp).ID=sscanf(line,'%c');
    sp=sp+1;
end
fclose(List);

% Load PIPER C2-Symmetry/Rotation Matrix
Rots_list=table2array(readtable('PIPER_Package/C2_rots.txt'));

%%%% DIPER Analysis
% Combination List for Translation Table
itr=0; l=20; List=zeros((l*2+1)^3,3);
for p=-l:l
    for q=-l:l
        for r=-l:l
            itr=itr+1;
            List(itr,1:3)=[p,q,r];
        end
    end
end

% Boundary Condition
% (Default: Limit(1800); Depth(600))
Limit=1800; Depth=600;

% DIPER Analysis
Sample(1).Lowest_RMSD_95=[];
Sample(1).Result_Idx=[]; Sample(1).RMSD_95=[];
Sample(1).Output_Idx=[]; Sample(1).Output_RMSD_95=[];

parfor r=1:size(Sample,2)

    %%%% Cull Crystal Homodimers for RMSD Examination
    % Load PDB File
    PDB=pdbread(['PDB/',lower(Sample(r).ID),'.pdb']);
    File=fopen(['PDB/',lower(Sample(r).ID),'.pdb']);

    % Symmetry Operator
    Operator=[];
    while(1)
        line=fgetl(File);
        if isequal(sscanf(line(1:10),'%c'),'REMARK 300'), break, end

        if isequal(sscanf(line(1:10),'%c'),'REMARK 290')
            Operator=[Operator; sscanf(line(25:79),'%f %f %f')'];
        end
    end

    % Unit Cell Parameter
    Unit_Cell=[];
    while(1)
        line=fgetl(File);
        if isequal(sscanf(line(1:3),'%c'),'CRY')
            Unit_Cell=sscanf(line(8:55),'%f %f %f');
            break,
        end
    end
    fclose(File);

    % Extract Coordinates of Alpha Carbon (CA)
    PDB_Model=[PDB.Model(1).Atom.X; PDB.Model(1).Atom.Y; PDB.Model(1).Atom.Z]';

    PDB_resSeq_Idx=[PDB.Model(1).Atom.resSeq]';
    pseudo_PDB_CA_Idx=strcmp({PDB.Model(1).Atom.AtomName}, 'CA')';
    PDB_resSeq_Unique_Idx=unique(PDB_resSeq_Idx(pseudo_PDB_CA_Idx));

    PDB_CA_Idx=zeros(size(PDB_resSeq_Unique_Idx,1),1);
    for RSUidx=1:size(PDB_resSeq_Unique_Idx,1)
        PDB_CA_Idx(RSUidx,1)=find((PDB_resSeq_Idx==PDB_resSeq_Unique_Idx(RSUidx))&pseudo_PDB_CA_Idx, 1, 'first');
    end
    Asy_Unit=PDB_Model(PDB_CA_Idx,:);

    % Sequence
    PDB_Seq=aminolookup(strrep(strrep([PDB.Model(1).Atom(PDB_CA_Idx).resName], 'SEC', 'CYS'), 'UNK', 'GLY'));

    %%% Construct Reduced Unit Box
    % Basis: (a, b, c) to (x, y, z)
    basis=eye(3);
    basis(1,2)=cosd(Unit_Cell(6));
    basis(2,2)=sind(Unit_Cell(6));
    basis(1,3)=cosd(Unit_Cell(5));
    basis(2,3)=(cosd(Unit_Cell(4))-basis(1,2)*basis(1,3))/basis(2,2);
    basis(3,3)=sqrt(1-basis(1,3)^2 -basis(2,3)^2);
    Unit_Cell_extent=basis.*repmat(Unit_Cell(1:3)',3,1);

    % Set of Symmertic Mates
    Unit=[];
    Sym_No=size(Operator,1)/3;
    for j=1:Sym_No
        SMTRY_Unit=Asy_Unit*Operator(3*j-2:3*j,1:3)'+ Operator(3*j-2:3*j,4)';

        % Shift
        Shift=List*Unit_Cell_extent';
        Translation_Table=mean(Asy_Unit)-(mean(SMTRY_Unit)+Shift);

        % Store Shifted SMTRY_Unit
        [~,idx]=min(sqrt(sum(Translation_Table.^2,2)));
        SMTRY_Unit=SMTRY_Unit+Shift(idx,:);
        Unit=[Unit; SMTRY_Unit];
    end

    % Unit Box
    Unit_Box=Unit;
    for x=-2:2
        for y=-2:2
            for z=-2:2
                if ~(x==0 && y==0 && z==0)
                    NB_Unit=Unit+[x,y,z]*Unit_Cell_extent';
                    Unit_Box=[Unit_Box; NB_Unit];
                end
            end
        end
    end

    % Identify Closest Packing Neighbors (PDB_NBs)
    NB_Cutoff=12;
    while(1)
        i=0; PDB_NB=struct('Coord', []);
        for j=2:size(Unit_Box,1)/size(Asy_Unit,1)
            DISM=pdist2(Asy_Unit, Unit_Box(size(Asy_Unit,1)*(j-1)+1:size(Asy_Unit,1)*j,:));
            if min(min(DISM))<NB_Cutoff
                i=i+1;
                PDB_NB(i).Coord=Unit_Box(size(Asy_Unit,1)*(j-1)+1:size(Asy_Unit,1)*j,:);
            end
        end

        NeiB_Idx=[];
        RMSD_Table=zeros(size(PDB_NB,2));
        for k=1:size(PDB_NB,2)
            for l=1:size(PDB_NB,2)
                [~,~, eRMSD_1]=CoordiExam([Asy_Unit; PDB_NB(k).Coord], [Asy_Unit; PDB_NB(l).Coord]);
                [~,~, eRMSD_2]=CoordiExam([Asy_Unit; PDB_NB(k).Coord], [PDB_NB(l).Coord; Asy_Unit]);
                RMSD_Table(k,l)=min(eRMSD_1, eRMSD_2);
            end
            NeiB_Idx=[NeiB_Idx, find(RMSD_Table(k,:)<1, 1)];
        end

        % Generate Closest/"Unique" Packing Neighbors (PDB_UNBs)
        UNB=unique(NeiB_Idx);
        m=0; PDB_UNB=struct('Coord',[]);
        for n=1:length(UNB)
            R=CoordiExam(Asy_Unit, PDB_NB(UNB(n)).Coord);
            Eigen=sortrows([real(eig(R)), imag(eig(R))])*[1;sqrt(-1)];

            % Examine C2 Symmetry (Angle)
            if norm(Eigen-[-1; -1; 1])<6*pi/180
                [V,D]=eig(R);
                [eig_v,eig_idx]=max(real(diag(D)'));

                % Examine C2 Symmetry (Axis: loose criteria (<5) for tolerance of CA mass center shift)
                if abs(acosd((mean(Asy_Unit)-mean(PDB_NB(UNB(n)).Coord))/norm(mean(Asy_Unit)-mean(PDB_NB(UNB(n)).Coord))*V(:,eig_idx))-90)<5
                    m=m+1;
                    PDB_UNB(m).Coord=[Asy_Unit; PDB_NB(UNB(n)).Coord];
                end
            end
        end

        % Termination
        if m>0
            break,
        else
            % Expand Searching Region
            NB_Cutoff=NB_Cutoff+6;
        end
    end

    %

    %%%% Construct C2-DIPER Homodimers
    % Extract CA Coordinate of Truncated Alphafold2 Monomer
    AF=pdbread([lower(Sample(r).ID),'_trun_AF.pdb']);
    Rcpt_Model=[AF.Model.Atom.X; AF.Model.Atom.Y; AF.Model.Atom.Z]';

    % Sequence
    AF_CA_Idx=strcmp({AF.Model.Atom.AtomName}, 'CA')';
    AF_Seq=aminolookup([AF.Model.Atom(AF_CA_Idx).resName]);

    %

    %%%% Sequence Alignment for Further RMSD 95 Calculation
    % Extract 1-letter Sequences
    seq1=[PDB_Seq, PDB_Seq]; seq2=[AF_Seq, AF_Seq];

    % Global Sequence Alignment to Extract Superimposed CA for RMSD Calculation
    [~, Alignment] = nwalign(seq1,seq2);
    alignidex=find((Alignment(1,:)~='-')&(Alignment(3,:)~='-'));
    aseq1=find(Alignment(1,:)=='-');
    aseq2=find(Alignment(3,:)=='-');

    aseq1_idx=[]; aseq2_idx=[];
    for o=1:length(alignidex)
        aseq1_idx=[aseq1_idx, alignidex(o)-sum(aseq1<alignidex(o))];
        aseq2_idx=[aseq2_idx, alignidex(o)-sum(aseq2<alignidex(o))];
    end

    %

    %%%% C2.DIPER Top Rank Output Analysis
    % Load PIPER Docking Results
    Result=fopen(['C2.DIPER_Result/',lower(Sample(r).ID),'/ft.000.00']);
    PIPER_result=zeros(Limit, 10);
    for p=1:Limit
        line=fgetl(Result);
        PIPER_result(p,1:10)=sscanf(line,'%f')';
    end
    fclose(Result);

    % PIPER Top Rank Output Analysis
    Table=zeros(Limit,10);
    Count=0; Result_Idx=0;
    while(1)
        % Construct Lgnd
        Result_Idx=Result_Idx+1;
        Trans_M=PIPER_result(Result_Idx,2:4);
        Rot_Info=Rots_list(PIPER_result(Result_Idx,1)+1,:);
        Rot_M=[Rot_Info(2:4); Rot_Info(5:7); Rot_Info(8:10)];
        Lgnd_Model=(Rcpt_Model-mean(Rcpt_Model))*Rot_M'+(Trans_M+mean(Rcpt_Model));

        % Examine the C2 symmertry
        Eigen=sortrows([real(eig(Rot_M)), imag(eig(Rot_M))])*[1;sqrt(-1)];
        if norm(Eigen-[-1; -1; 1])<6*pi/180

            % Examine the C2 axis
            [V,D]=eig(Rot_M);
            [eig_v,eig_idx]=max(real(diag(D)'));
            if abs(acosd((mean(Rcpt_Model)-mean(Lgnd_Model))/norm(mean(Rcpt_Model)-mean(Lgnd_Model))*V(:,eig_idx))-90)<3

                % Record
                Count=Count+1;
                Table(Count,1)=Result_Idx;
                Table(Count,2:4)=rotm2eul(Rot_M);                   % Rotation Matrix
                Table(Count,5:7)=mean(Lgnd_Model);                  % Lgnd Mass Center
                Table(Count,8)=PIPER_result(Result_Idx,5);          % PIPER Energy

                % RMSD 95 calculation (ATTTENTION: ORDER OF MONOMERS)
                Rcpt=Rcpt_Model(AF_CA_Idx,:); Lgnd=Lgnd_Model(AF_CA_Idx,:);
                PDB_1=[Rcpt; Lgnd]; PDB_CA_1=PDB_1(aseq2_idx,:);
                PDB_2=[Lgnd; Rcpt]; PDB_CA_2=PDB_2(aseq2_idx,:);
                len=size(PDB_CA_1,1);

                RMSD_95_Table=[];
                for q=1:size(PDB_UNB,2)
                    [~,~,~,~, Sort_Sum_Square_1]=CoordiExam(PDB_UNB(q).Coord(aseq1_idx,:), PDB_CA_1);
                    [~,~,~,~, Sort_Sum_Square_2]=CoordiExam(PDB_UNB(q).Coord(aseq1_idx,:), PDB_CA_2);
                    RMSD_95_1=sqrt(sum(Sort_Sum_Square_1(1:ceil(len*0.95)))/ceil(len*0.95));
                    RMSD_95_2=sqrt(sum(Sort_Sum_Square_2(1:ceil(len*0.95)))/ceil(len*0.95));
                    RMSD_95_Table(q)=min(RMSD_95_1, RMSD_95_2);
                end
                Table(Count,9)=min(RMSD_95_Table,[],'all');
            end
        end

        % Termination
        if Count==Depth||Result_Idx==Limit, break; end            % Enough Counts (600)||Depth Limit (1800)
    end
    Table(Count+1:end,:)=[];
    Sample(r).Lowest_RMSD_95=min(Table(:,9));

    %

    %%%% Agglomerative Hierarchical Clustering
    % oRMSD DISM (for Representative Pose Selection)
    oRMSD_DISM=zeros(1,Count*(Count-1)/2);
    t=0;
    for u=1:(Count-1)
        u_Trans_M=PIPER_result(Table(u,1), 2:4);
        u_Rot_Info=Rots_list(PIPER_result(Table(u,1),1)+1,:);
        u_Rot_M=[u_Rot_Info(2:4); u_Rot_Info(5:7); u_Rot_Info(8:10)];
        fromXYZ_Model=(Rcpt_Model-mean(Rcpt_Model))*u_Rot_M'+(u_Trans_M+mean(Rcpt_Model));
        fromXYZ=fromXYZ_Model(AF_CA_Idx,:);

        for v=u+1:Count
            t=t+1;
            v_Trans_M=PIPER_result(Table(v,1), 2:4);
            v_Rot_Info=Rots_list(PIPER_result(Table(v,1),1)+1,:);
            v_Rot_M=[v_Rot_Info(2:4); v_Rot_Info(5:7); v_Rot_Info(8:10)];
            toXYZ_Model=(Rcpt_Model-mean(Rcpt_Model))*v_Rot_M'+(v_Trans_M+mean(Rcpt_Model));
            toXYZ=toXYZ_Model(AF_CA_Idx,:);

            % oRMSD
            [~,~,oRMSD]=CoordiExam([Rcpt; fromXYZ], [Rcpt; toXYZ]);
            oRMSD_DISM(1,t)=oRMSD;
        end
    end
    Squared_oRMSD_DISM=squareform(oRMSD_DISM);

    % Agglomerative Hierarchical Clustering
    oRMSD_Cutoff=5;
    while (1)
        % Linkage Clustering (Farest Within-Cluster Distance)
        oRMSD_Link=linkage(oRMSD_DISM,'complete');
        Table(:,10)=cluster(oRMSD_Link,'cutoff',oRMSD_Cutoff,'Criterion','distance');

        % No of Unique Clustering
        Clust_Idx=unique(Table(:,10));

        % Termination
        if size(Clust_Idx,1)>=15
            break,
        else
            oRMSD_Cutoff=oRMSD_Cutoff-0.5;
        end
        if oRMSD_Cutoff==0, break; end
    end

    %

    %%%% Rank Clusters in Consideration of Energy & Cluster Info
    % Record Info of Clusters
    Clust_Info=zeros(size(Clust_Idx,1),6);
    for g=1:size(Clust_Idx,1)
        Clust_Info(g,1)=Clust_Idx(g);
        Clust_Info(g,2)=sum(Table(:,10)==Clust_Idx(g));                   % Cluster Size
        Clust_Info(g,3)=median(Table(Table(:,10)==Clust_Idx(g),8));       % Average PIPER Energy of Cluster

        % Select Representative Pose in Consideration of PIPER Energy & oRMSD
        Pose_Idx=find(Table(:,10)==Clust_Idx(g));
        Pose_Info=zeros(size(Pose_Idx,1),3);
        
        for h=1:size(Pose_Idx,1)
            Pose_Info(h,1)=Table(Pose_Idx(h),8);
            Pose_Info(h,2)=sum(Squared_oRMSD_DISM(Pose_Idx(h), Pose_Idx))/size(Pose_Idx,1);
        end

        [~,Energy_Pose_Rank]=sort(Pose_Info(:,1));
        [~,Energy_Pose_Order]=sort(Energy_Pose_Rank);
        [~,oRMSD_Pose_Rank]=sort(Pose_Info(:,2));
        [~,oRMSD_Pose_Order]=sort(oRMSD_Pose_Rank);

        % *** Tradeoff of PIPER Energy & oRMSD
        Pose_Info(:,3)=oRMSD_Pose_Order*0.6+Energy_Pose_Order*0.4;
        [~,col]=min(Pose_Info(:,3));

        Clust_Info(g,4)=Table(Pose_Idx(col),1);                           % Result_Idx
        Clust_Info(g,5)=Table(Pose_Idx(col),8);                           % PIPER Energy
        Clust_Info(g,6)=Table(Pose_Idx(col),9);                           % RMSD_95
    end

    % Ranked by PIPER Energy (PIPER_Rank)
    [~,PIPER_Idx]=sort(Clust_Info(:,5), 'ascend');
    PIPER_Rank=Clust_Info(PIPER_Idx,:);

    % Ranking by PIPER Energy & Number/Average PIPER Energy of Cluster Size
    [~,NeiP_No_PIdx]=sort(PIPER_Rank(:,2), 'descend');                     % Number of Cluster Size
    [~,NeiP_En_PIdx]=sort(PIPER_Rank(:,3), 'ascend');                      % Average PIPER Energy of Cluster
    NeiP_RankP=[NeiP_No_PIdx, NeiP_En_PIdx];

    NeiP_PIdx=[];
    for PR=1:size(PIPER_Rank,1)
        [row,~]=find(NeiP_RankP==PR);
        NeiP_PIdx(PR,1)=min(row)+mean(row)/(Limit*10);
    end

    [~, NeiP_Idx_RankP]=sort(NeiP_PIdx);
    [~, Nei_OrderP]=sort(NeiP_Idx_RankP);

    % *** Tradeoff of PIPER Energy & Cluster
    PC_Rank=(1:size(PIPER_Rank,1))'*0.6+Nei_OrderP*0.4;
    [~,PC_Idx]=sort(PC_Rank);

    Sample(r).Result_Idx=PIPER_Rank(PC_Idx,4);
    Sample(r).RMSD_95=PIPER_Rank(PC_Idx,6);

    %

    %%%% Output
    PC_Result=PIPER_Rank(PC_Idx,:);

    if size(PC_Result,1)>15
        % Once More If Too Much Clusterings
        OM_Result_Idx=zeros(size(PC_Result,1),1);
        for CR=1:size(PC_Result,1)
            OM_Result_Idx(CR,1)=find(Table(:,1)==PC_Result(CR,4));
        end

        OM_Result_oRMSD=Squared_oRMSD_DISM(OM_Result_Idx, OM_Result_Idx);
        warning('off','all');
        Z=linkage(squareform(OM_Result_oRMSD),'centroid');

        Default_Cutoff=5;
        while (1)
            T=cluster(Z,'cutoff',Default_Cutoff,'Criterion','distance');
            if max(T)>=15
                break,
            else
                Default_Cutoff=Default_Cutoff-1;
            end
            if Default_Cutoff==0, break; end
        end

        % Output
        Output=zeros(max(T),7);
        for f=1:max(T)
            Output_List=find(T==f);
            if sum(T==f)<=2
                Output(f,1)=min(Output_List);
            else
                [~,Output_Idx]=min(sum(OM_Result_oRMSD(T==f,T==f)));
                Output(f,1)=Output_List(Output_Idx);
            end
            Output(f,2:7)=PC_Result(Output(f,1),1:6);
        end

        [~, Output_Order]=sort(Output(:,1));
        Output_Rank=Output(Output_Order,:);

        Sample(r).Output_Idx=Output_Rank(:,5);
        Sample(r).Output_RMSD_95=Output_Rank(:,7);
    else
        Sample(r).Output_Idx=Sample(r).Result_Idx;
        Sample(r).Output_RMSD_95=Sample(r).RMSD_95;
    end
    r
end

%
%
%

%%%% Extract C2_TopRank_Result
C2_TopRank_Result(1).ID=[];
C2_TopRank_Result(1).Lowest_RMSD_95=[];
C2_TopRank_Result(1).Output_Idx=[];
C2_TopRank_Result(1).Output_RMSD_95=[];
C2_TopRank_Result(1).Best_Idx=[];
C2_TopRank_Result(1).Best_RMSD_95=[];

for r=1:size(Sample,2)
    % Extract Top 25 Results
    Pose_No=size(Sample(r).Output_Idx,1);
    if Pose_No>25, Pose_No=25; end

    % Record ID, Idx, RMSD95 & Best Result
    C2_TopRank_Result(r).ID=Sample(r).ID;
    C2_TopRank_Result(r).Lowest_RMSD_95=Sample(r).Lowest_RMSD_95;
    C2_TopRank_Result(r).Output_Idx=Sample(r).Output_Idx(1:Pose_No);
    C2_TopRank_Result(r).Output_RMSD_95=Sample(r).Output_RMSD_95(1:Pose_No);

    % Best Result
    [val, row]=min(C2_TopRank_Result(r).Output_RMSD_95);
    C2_TopRank_Result(r).Best_Idx=C2_TopRank_Result(r).Output_Idx(row);
    C2_TopRank_Result(r).Best_RMSD_95=val;
end

save('C2.DIPER_Result/C2_TopRank_Result.mat', 'C2_TopRank_Result');

%
%
%

%%%% Construct C2_TopRank_Result Homodimers
% Load PIPER C2-Symmetry/Rotation Matrix
Rots_list=table2array(readtable('PIPER_Package/C2_rots.txt'));

% Default Parameters
prm=0; set=['ft.00', num2str(prm), '.00']; Limit=1800;

% Combination List for Translation Table
itr=0; l=20; List=zeros((l*2+1)^3,3);
for a=-l:l
    for b=-l:l
        for c=-l:l
            itr=itr+1;
            List(itr,1:3)=[a,b,c];
        end
    end
end

% Establish Deposited File of C2.DIPER Top 25 Outputs
mkdir('C2.DIPER_Result/C2.Top25_Pose');

for r=1:size(C2_TopRank_Result,2)

    %%%% Cull Crystal Homodimers for RMSD Examination
    % Load PDB File
    PDB=pdbread(['PDB/',lower(C2_TopRank_Result(r).ID),'.pdb']);
    File=fopen(['PDB/',lower(C2_TopRank_Result(r).ID),'.pdb']);
    Ans.Model.Atom=repmat(PDB.Model.Atom,1,2);

    % Symmetry Operator
    Operator=[];
    while(1)
        line=fgetl(File);
        if isequal(sscanf(line(1:10),'%c'),'REMARK 300'), break, end

        if isequal(sscanf(line(1:10),'%c'),'REMARK 290')
            Operator=[Operator; sscanf(line(25:79),'%f %f %f')'];
        end
    end

    % Unit Cell Parameter
    Unit_Cell=[];
    while(1)
        line=fgetl(File);
        if isequal(sscanf(line(1:3),'%c'),'CRY')
            Unit_Cell=sscanf(line(8:55),'%f %f %f');
            break,
        end
    end
    fclose(File);

    % Extract Coordinates of Alpha Carbon (CA)
    PDB_Model=[PDB.Model(1).Atom.X; PDB.Model(1).Atom.Y; PDB.Model(1).Atom.Z]';

    PDB_resSeq_Idx=[PDB.Model(1).Atom.resSeq]';
    pseudo_PDB_CA_Idx=strcmp({PDB.Model(1).Atom.AtomName}, 'CA')';
    PDB_resSeq_Unique_Idx=unique(PDB_resSeq_Idx(pseudo_PDB_CA_Idx));

    PDB_CA_Idx=zeros(size(PDB_resSeq_Unique_Idx,1),1);
    for RSUidx=1:size(PDB_resSeq_Unique_Idx,1)
        PDB_CA_Idx(RSUidx,1)=find((PDB_resSeq_Idx==PDB_resSeq_Unique_Idx(RSUidx))&pseudo_PDB_CA_Idx, 1, 'first');
    end
    Asy_Unit=PDB_Model(PDB_CA_Idx,:);

    %%% Construct Reduced Unit Box
    % Basis: (a, b, c) to (x, y, z)
    basis=eye(3);
    basis(1,2)=cosd(Unit_Cell(6));
    basis(2,2)=sind(Unit_Cell(6));
    basis(1,3)=cosd(Unit_Cell(5));
    basis(2,3)=(cosd(Unit_Cell(4))-basis(1,2)*basis(1,3))/basis(2,2);
    basis(3,3)=sqrt(1-basis(1,3)^2 -basis(2,3)^2);
    Unit_Cell_extent=basis.*repmat(Unit_Cell(1:3)',3,1);

    % Set of Symmertic Mates
    Unit=[];
    Sym_No=size(Operator,1)/3;
    for j=1:Sym_No
        SMTRY_Unit=Asy_Unit*Operator(3*j-2:3*j,1:3)'+ Operator(3*j-2:3*j,4)';

        % Shift
        Shift=List*Unit_Cell_extent';
        Translation_Table=mean(Asy_Unit)-(mean(SMTRY_Unit)+Shift);

        % Store Shifted SMTRY_Unit
        [~,idx]=min(sqrt(sum(Translation_Table.^2,2)));
        SMTRY_Unit=SMTRY_Unit+Shift(idx,:);
        Unit=[Unit; SMTRY_Unit];
    end

    % Unit Box
    Unit_Box=Unit;
    for x=-2:2
        for y=-2:2
            for z=-2:2
                if ~(x==0 && y==0 && z==0)
                    NB_Unit=Unit+[x,y,z]*Unit_Cell_extent';
                    Unit_Box=[Unit_Box; NB_Unit];
                end
            end
        end
    end

    % Identify Closest Packing Neighbors (PDB_NBs)
    NB_Cutoff=12;
    while(1)
        i=0; PDB_NB=struct('Coord', []);
        for j=2:size(Unit_Box,1)/size(Asy_Unit,1)
            DISM=pdist2(Asy_Unit, Unit_Box(size(Asy_Unit,1)*(j-1)+1:size(Asy_Unit,1)*j,:));
            if min(min(DISM))<NB_Cutoff
                i=i+1;
                PDB_NB(i).Coord=Unit_Box(size(Asy_Unit,1)*(j-1)+1:size(Asy_Unit,1)*j,:);
            end
        end

        NeiB_Idx=[];
        RMSD_Table=zeros(size(PDB_NB,2));
        for k=1:size(PDB_NB,2)
            for l=1:size(PDB_NB,2)
                [~,~, eRMSD_1]=CoordiExam([Asy_Unit; PDB_NB(k).Coord], [Asy_Unit; PDB_NB(l).Coord]);
                [~,~, eRMSD_2]=CoordiExam([Asy_Unit; PDB_NB(k).Coord], [PDB_NB(l).Coord; Asy_Unit]);
                RMSD_Table(k,l)=min(eRMSD_1, eRMSD_2);
            end
            NeiB_Idx=[NeiB_Idx, find(RMSD_Table(k,:)<1, 1)];
        end

        % Generate Closest/"Unique" Packing Neighbors (PDB_UNBs)
        UNB=unique(NeiB_Idx);
        m=0; PDB_UNB=struct('Coord',[]);
        for n=1:length(UNB)
            R=CoordiExam(Asy_Unit, PDB_NB(UNB(n)).Coord);
            Eigen=sortrows([real(eig(R)), imag(eig(R))])*[1;sqrt(-1)];

            % Examine C2 Symmetry (Angle)
            if norm(Eigen-[-1; -1; 1])<6*pi/180
                [V,D]=eig(R);
                [eig_v,eig_idx]=max(real(diag(D)'));

                % Examine C2 Symmetry (Axis: loose criteria (<5) for tolerance of CA mass center shift)
                if abs(acosd((mean(Asy_Unit)-mean(PDB_NB(UNB(n)).Coord))/norm(mean(Asy_Unit)-mean(PDB_NB(UNB(n)).Coord))*V(:,eig_idx))-90)<5
                    m=m+1;
                    PDB_UNB(m).Coord=[Asy_Unit; PDB_NB(UNB(n)).Coord];
                end
            end
        end

        % Termination
        if m>0
            break,
        else
            % Expand Searching Region
            NB_Cutoff=NB_Cutoff+6;
        end
    end

    %
    
    %%%% Construct C2-DIPER Homodimers
    % Establish Deposited File
    mkdir(['C2.DIPER_Result/C2.Top25_Pose/', lower(C2_TopRank_Result(r).ID)]);

    % Load PIPER Docking Result
    DIPER=fopen(['C2.DIPER_Result/',lower(C2_TopRank_Result(r).ID),'/',set]);
    PIPER_result=zeros(Limit, 10);
    for p=1:Limit
        line=fgetl(DIPER);
        PIPER_result(p,1:10)=sscanf(line,'%f')';
    end
    fclose(DIPER);

    % Extract CA Coordinate of Truncated Alphafold2 Monomer
    AF=pdbread([lower(C2_TopRank_Result(r).ID),'_trun_AF.pdb']);
    Model=[AF.Model.Atom.X; AF.Model.Atom.Y; AF.Model.Atom.Z]';
    AF_CA_Idx=strcmp({AF.Model.Atom.AtomName}, 'CA')';
    Rcpt=Model(AF_CA_Idx,:);

    % Construct Representative Pose (RePo)
    AF.Model.Atom=repmat(AF.Model.Atom,1,2);
    Pose_No=size(C2_TopRank_Result(r).Output_Idx,1);

    % Extract 1-letter Sequences
    PDB_Seq=aminolookup(strrep(strrep([PDB.Model(1).Atom(PDB_CA_Idx).resName], 'SEC', 'CYS'), 'UNK', 'GLY'));
    AF_Seq=aminolookup([AF.Model.Atom(AF_CA_Idx).resName]);

    seq1=[PDB_Seq, PDB_Seq]; seq2=[AF_Seq, AF_Seq];

    % Global Sequence Alignment to Extract Superimposed CA for RMSD Calculation
    [~, Alignment] = nwalign(seq1,seq2);
    alignidex=find((Alignment(1,:)~='-')&(Alignment(3,:)~='-'));
    aseq1=find(Alignment(1,:)=='-');
    aseq2=find(Alignment(3,:)=='-');

    aseq1_idx=[]; aseq2_idx=[];
    for o=1:length(alignidex)
        aseq1_idx=[aseq1_idx, alignidex(o)-sum(aseq1<alignidex(o))];
        aseq2_idx=[aseq2_idx, alignidex(o)-sum(aseq2<alignidex(o))];
    end

    %

    % Construct C2-DIPER Homodimers
    for q=1:Pose_No
        Result_Idx=C2_TopRank_Result(r).Output_Idx(q);
        Trans_M=PIPER_result(Result_Idx,2:4);
        Rot_Info=Rots_list(PIPER_result(Result_Idx,1)+1,:);
        Rot_M=[Rot_Info(2:4); Rot_Info(5:7); Rot_Info(8:10)];
        Lgnd=(Rcpt-mean(Rcpt))*Rot_M'+repmat(Trans_M+mean(Rcpt),size(Rcpt,1),1);

        % Output C2-Symmetric PDB File
        [R_AF,T_AF]=CoordiExam(Rcpt, Lgnd);
        Model_Lgnd=Model*R_AF+repmat(T_AF,size(Model,1),1);

        len_AF=size(Model_Lgnd,1);
        for s=1:len_AF
            AF.Model.Atom(s).chainID='A';
            AF.Model.Atom(len_AF+s).chainID='B';
            AF.Model.Atom(s).AtomSerNo=s;
            AF.Model.Atom(len_AF+s).AtomSerNo=len_AF+s;
            AF.Model.Atom(len_AF+s).X=Model_Lgnd(s,1);
            AF.Model.Atom(len_AF+s).Y=Model_Lgnd(s,2);
            AF.Model.Atom(len_AF+s).Z=Model_Lgnd(s,3);
        end

        warning('off','all');
        pdbwrite(['C2.DIPER_Result/C2.Top25_Pose/', lower(C2_TopRank_Result(r).ID), '/', lower(C2_TopRank_Result(r).ID), set(3:7), num2str(Result_Idx), '.pdb'], AF);

        %

        % RMSD 95 calculation (ATTTENTION: ORDER OF MONOMERS)
        PDB_1=[Rcpt; Lgnd]; PDB_CA_1=PDB_1(aseq2_idx,:);
        PDB_2=[Lgnd; Rcpt]; PDB_CA_2=PDB_2(aseq2_idx,:);
        len_CA=size(PDB_CA_1,1);

        RMSD_95_Table=[];
        for w=1:size(PDB_UNB,2)
            [~,~,~,~, Sort_Sum_Square_1]=CoordiExam(PDB_UNB(w).Coord(aseq1_idx,:), PDB_CA_1);
            [~,~,~,~, Sort_Sum_Square_2]=CoordiExam(PDB_UNB(w).Coord(aseq1_idx,:), PDB_CA_2);
            RMSD_95_1=sqrt(sum(Sort_Sum_Square_1(1:ceil(len_CA*0.95)))/ceil(len_CA*0.95));
            RMSD_95_2=sqrt(sum(Sort_Sum_Square_2(1:ceil(len_CA*0.95)))/ceil(len_CA*0.95));
            RMSD_95_Table(w)=min(RMSD_95_1, RMSD_95_2);
        end
        [~, row]=min(RMSD_95_Table,[],'all');

        % Output Ans PDB File
        [R_PDB,T_PDB]=CoordiExam(Asy_Unit, PDB_UNB(row).Coord(length(Asy_Unit)+1:end,:));
        PDB_Model_Pair=PDB_Model*R_PDB+repmat(T_PDB,size(PDB_Model,1),1);

        len_PDB=size(PDB_Model_Pair,1);
        for t=1:len_PDB
            Ans.Model.Atom(t).chainID='A';
            Ans.Model.Atom(len_PDB+t).chainID='B';
            Ans.Model.Atom(t).AtomSerNo=t;
            Ans.Model.Atom(len_PDB+t).AtomSerNo=len_PDB+t;
            Ans.Model.Atom(len_PDB+t).X=PDB_Model_Pair(t,1);
            Ans.Model.Atom(len_PDB+t).Y=PDB_Model_Pair(t,2);
            Ans.Model.Atom(len_PDB+t).Z=PDB_Model_Pair(t,3);
        end

        warning('off','all');
        pdbwrite(['C2.DIPER_Result/C2.Top25_Pose/', lower(C2_TopRank_Result(r).ID), '/', lower(C2_TopRank_Result(r).ID), set(3:7), num2str(Result_Idx), '_Ans.pdb'], Ans);
    end
    r
end

% Shut Down Idle Cores
delete(gcp('nocreate'));