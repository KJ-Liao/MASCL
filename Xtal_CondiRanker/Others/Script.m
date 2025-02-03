%%%% Default Settings
% Parallel Computation
parpool(6)

%%%% Input
PDB_ID='1f0w';

%
%
%

% Establish Output File
mkdir(PDB_ID);

% Residue Table
Res={'GLY';'ALA';'VAL';'ILE';'LEU';
    'SER';'THR';'ASP';'ASN';'GLU';
    'GLN';'LYS';'ARG';'CYS';'MET';
    'PHE';'TYR';'TRP';'HIS';'PRO'};

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

%

%%%%%%%% Interaction txt %%%%%%%%

% Establish Output File
mkdir([PDB_ID,'/PatchBag_Result_',PDB_ID]);

% Distance Cutoff (A) of Interacting Residues
Threshold=10;

%%%% Int Part

try
    %%% Load Corresponding PDB File
    File=fopen([PDB_ID,'.pdb']);

    % Symmetry Operator
    Operator=[];
    while(1)
        line=fgetl(File);
        if isequal(sscanf(line(1:10),'%c'),'REMARK 300'), break, end

        if isequal(sscanf(line(1:10),'%c'),'REMARK 290')
            SMTRY=sscanf(line(25:79),'%f %f %f');
            Operator=[Operator; SMTRY'];
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
    PDB=pdbread([PDB_ID,'.pdb']);
    PDB_Model=[PDB.Model(1).Atom.X; PDB.Model(1).Atom.Y; PDB.Model(1).Atom.Z]';

    PDB_resSeq_Idx=[PDB.Model(1).Atom.resSeq]';
    pseudo_PDB_CA_Idx=strcmp({PDB.Model(1).Atom.AtomName}, 'CA')';
    PDB_resSeq_Unique_Idx=unique(PDB_resSeq_Idx(pseudo_PDB_CA_Idx));

    PDB_CA_Idx=zeros(size(PDB_resSeq_Unique_Idx,1),1);
    for RSUidx=1:size(PDB_resSeq_Unique_Idx,1)
        PDB_CA_Idx(RSUidx,1)=find((PDB_resSeq_Idx==PDB_resSeq_Unique_Idx(RSUidx))&pseudo_PDB_CA_Idx, 1, 'first');
    end
    Asy_Unit=PDB_Model(PDB_CA_Idx,:);

    %

    if PDB.Sequence(1).NumOfResidues<20
        disp('Please check the protein sequence.');
    elseif size(Asy_Unit,1)<15
        disp('Please check the protein sequence.');
    else

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

        % Reduced Unit Box
        Up_B=max(Asy_Unit,[],1)+Threshold;
        Low_B=min(Asy_Unit,[],1)-Threshold;
        X=Unit_Box(:,1)>Low_B(1)&Unit_Box(:,1)<Up_B(1);
        Y=Unit_Box(:,2)>Low_B(2)&Unit_Box(:,2)<Up_B(2);
        Z=Unit_Box(:,3)>Low_B(3)&Unit_Box(:,3)<Up_B(3);

        InBox=find(X+Y+Z==3);
        Reduced_Box=Unit_Box(InBox,:);
        Reduced_Neighbor_Box=Reduced_Box(size(Asy_Unit,1)+1:end,:);

        %

        %%% Identify Interacting Residues
        % Index of Interacting Residues
        DISM=pdist2(Asy_Unit, Reduced_Neighbor_Box);
        Asy_Unit_Idx=find(sum(DISM<Threshold,2));       % Row-wise
        Neighbor_Idx=find(sum(DISM<Threshold,1));       % Column-wise

        % Residue Type List of Interacting Residues
        PDB_Res={PDB.Model(1).Atom.resName}';
        Res_type=PDB_Res(PDB_CA_Idx,:);
        Unit_Box_Res_type=repmat(Res_type,5*5*5*size(Operator,1)/3,1);

        Reduced_Box_Res_type=Unit_Box_Res_type(InBox);
        Reduced_Neighbor_Box_Res_type=Reduced_Box_Res_type(size(Asy_Unit,1)+1:end);

        %

        %%% Output
        % Coordinate of Interacting Residues
        Asy_Coord=Asy_Unit(Asy_Unit_Idx,:);
        NB_Coord=Reduced_Neighbor_Box(Neighbor_Idx,:);

        % Residue Type of Interacting Residues
        Asy_Rtype=Res_type(Asy_Unit_Idx);
        NB_Rtype=Reduced_Neighbor_Box_Res_type(Neighbor_Idx);

        % Res_Type_Seperated Structure
        Coord=[Asy_Coord; NB_Coord];
        Mo_Idx=[ones(size(Asy_Coord,1),1); zeros(size(NB_Coord,1),1)];
        Rtype=[Asy_Rtype; NB_Rtype];
        Rtype_Idx=zeros(size(Rtype));

        for r=1:20
            Rtype_Idx=Rtype_Idx+strcmp(Rtype,Res{r})*r;
        end

        % Output
        Output=[];
        Output.X=Coord(:,1);
        Output.Y=Coord(:,2);
        Output.Z=Coord(:,3);
        Output.Res=Rtype_Idx;
        Output.Mo=Mo_Idx;
        writetable(struct2table(Output),[PDB_ID,'/PatchBag_Result_',PDB_ID,'/',PDB_ID,'_Int.txt']);
        
    end
catch
    disp('Please check the protein pdb file.');
end

clearvars -except PDB_ID Res List Threshold

%%%% nInt Part

try
    % Extract Coordinates of Alpha Carbon (CA)
    PDB=pdbread([PDB_ID,'.pdb']);
    PDB_Model=[PDB.Model(1).Atom.X; PDB.Model(1).Atom.Y; PDB.Model(1).Atom.Z]';

    PDB_CA_Idx=strcmp({PDB.Model(1).Atom.AtomName}, 'CA')';
    PDB_resSeq_Idx=[PDB.Model(1).Atom.resSeq]';

    Res_ID_Idx=zeros(size(PDB_resSeq_Idx,1),1);
    Sub_Idx=sort(unique(PDB_resSeq_Idx));
    for Sidx=1:size(Sub_Idx,1)
        Res_ID_Idx(PDB_resSeq_Idx==Sub_Idx(Sidx))=Sidx;
    end

    % Make Ply file
    try
        [status, ~] = system(['python3 generate_3dzd.py ', PDB_ID, '.pdb ',PDB_ID, '/']);
    catch
        disp('Please check generate_3dzd.py.');
    end
    
    % Load Protein Molecular Surface of Corresponding PDB ID
    Ply_file = fopen([PDB_ID,'/',PDB_ID, '.ply']);
    while(1)
        line=fgetl(Ply_file);
        if length(line)>14 && strcmp(line(1:14),'element vertex')
            Vertex_No=str2double(line(16:end));
        elseif length(line)>12 && strcmp(line(1:12),'element face')
            Face_No=str2double(line(14:end));
        elseif length(line)==10 && strcmp(line(1:10),'end_header')
            RawText = fread(Ply_file,inf,'*char');
            splitLines = cell2mat(textscan(RawText, '%f %f %f %f %*[^\n]'));
            break;
        else
            continue;
        end
    end
    fclose(Ply_file);
    Vertex_Table=splitLines(1:Vertex_No,1:3);

    % Extract Surface Residues (Dedault Threshold: 1.8A)
    DISM=pdist2(PDB_Model, Vertex_Table);
    Surf_Atom=find(min(DISM, [], 2)<1.8);
    Surf_Res=unique(Res_ID_Idx(Surf_Atom));
    ref_Surf_Res=intersect(Surf_Res, unique(Res_ID_Idx(PDB_CA_Idx)));

    Surf_Res_Sym_Idx=zeros(size(ref_Surf_Res,1),1);
    for RSRidx=1:size(ref_Surf_Res,1)
        Surf_Res_Sym_Idx(RSRidx,1)=find((Res_ID_Idx==ref_Surf_Res(RSRidx))&PDB_CA_Idx, 1, 'first');
    end
    PDB_Res={PDB.Model(1).Atom.resName}';
    Surf_Res_Sym=PDB_Res(Surf_Res_Sym_Idx,:);

    Surf_Res_Idx=zeros(size(ref_Surf_Res,1),1);
    for r=1:20
        Surf_Res_Idx=Surf_Res_Idx+strcmp(Surf_Res_Sym,Res{r})*r;
    end
    Surf_Res_Info=[PDB_Model(Surf_Res_Sym_Idx,:), Surf_Res_Idx];

    % non_Int_Surf_Res Coordiinates
    Int_File=table2array(readtable([PDB_ID,'/PatchBag_Result_',PDB_ID,'/',PDB_ID,'_Int.txt']));
    Int_Res_Info=Int_File(Int_File(:,5)==1,1:4);

    % Output n_Int Surf Res
    n_Int_Res_Info=setdiff([Surf_Res_Info; Int_Res_Info], Int_Res_Info, 'rows');

    n_Int_Output=[];
    n_Int_Output.X=n_Int_Res_Info(:,1);
    n_Int_Output.Y=n_Int_Res_Info(:,2);
    n_Int_Output.Z=n_Int_Res_Info(:,3);
    n_Int_Output.Res=n_Int_Res_Info(:,4);
    writetable(struct2table(n_Int_Output),[PDB_ID,'/PatchBag_Result_',PDB_ID,'/',PDB_ID,'_nInt.txt']);

catch
    disp('Please check the protein pdb file.');
end

clearvars -except PDB_ID 

%
%
%

%%%%%%%% Patch Extraction %%%%%%%%

% Log Postion of Alpha Carbon and Corresponding Index
PDB=pdbread([PDB_ID,'/',PDB_ID,'.pdb']);
PDB_resSeq_Idx=[PDB.Model(1).Atom.resSeq]';
PDB_Model=[PDB.Model(1).Atom.X; PDB.Model(1).Atom.Y; PDB.Model(1).Atom.Z]';

pseudo_PDB_CA_Idx=strcmp({PDB.Model(1).Atom.AtomName}, 'CA')';
PDB_resSeq_Unique_Idx=unique(PDB_resSeq_Idx(pseudo_PDB_CA_Idx));

PDB_CA_Idx=zeros(size(PDB_resSeq_Unique_Idx,1),1);
for RSUidx=1:size(PDB_resSeq_Unique_Idx,1)
    PDB_CA_Idx(RSUidx,1)=find((PDB_resSeq_Idx==PDB_resSeq_Unique_Idx(RSUidx))&pseudo_PDB_CA_Idx, 1, 'first');
end

Pool=[PDB_Model, PDB_resSeq_Idx];
CA_Coord=Pool(PDB_CA_Idx,:);

% Patch Size (No of Alpha Carbons)
Patch_Size=6;

% Pre-allocate RAM for Patch Struct
fields = {'Coord','Norm_V','Res_Type'};

%%%% Int Part

try
    % Load Info of Interacting Residues
    Int_File=table2array(readtable([PDB_ID,'/PatchBag_Result_',PDB_ID,'/',PDB_ID,'_Int.txt']));
    Int_Asy_Unit=Int_File(Int_File(:,5)==1,1:4);

    if size(Int_Asy_Unit,1)<6
        Patch = cell2struct(cell(length(fields), []), fields);
        writetable(struct2table(Patch),[PDB_ID,'/PatchBag_Result_',PDB_ID,'/',PDB_ID,'_Int_Patch.txt']);
    else
        % Read Molecule Surface
        Ply_file = fopen([PDB_ID,'/',PDB_ID, '.ply']);
        while(1)
            line=fgetl(Ply_file);
            if length(line)>14 && strcmp(line(1:14),'element vertex')
                Vertex_No=str2double(line(16:end));
            elseif length(line)>12 && strcmp(line(1:12),'element face')
                Face_No=str2double(line(14:end));
            elseif length(line)==10 && strcmp(line(1:10),'end_header')
                RawText = fread(Ply_file,inf,'*char');
                splitLines = cell2mat(textscan(RawText, '%f %f %f %f %*[^\n]'));
                break;
            else
                continue;
            end
        end
        fclose(Ply_file);

        % Normal Vector Calculation
        Vertex_Table=splitLines(1:Vertex_No,1:3);
        CA_Idx_Table=pdist2(Int_Asy_Unit(:,1:3), CA_Coord(:,1:3));

        for k=1:size(CA_Idx_Table,1)
            Res_Idx=CA_Coord(CA_Idx_Table(k,:)==0,4);
            CA_Position=CA_Coord(CA_Coord(:,4)==Res_Idx,1:3);
            Residue_Surface_DISM=pdist2(Pool(Pool(:,4)==Res_Idx,1:3), Vertex_Table);

            Cutoff=2;
            while(1)
                if Cutoff>8
                    disp('Please check the Cutoff.');
                    break;
                end

                if sum(sum(Residue_Surface_DISM<Cutoff)>0)>10
                    break;
                end
                Cutoff=Cutoff+1;
            end

            V=mean(Vertex_Table(sum(Residue_Surface_DISM<Cutoff)>0,:))-CA_Position;
            Int_Asy_Unit(k,5:7)=V/norm(V);
        end

        %

        % Patch Extraction
        DISM=squareform(pdist(Int_Asy_Unit(:,1:3)));
        Patch_Ridx=zeros(size(Int_Asy_Unit,1),Patch_Size);
        for p=1:size(Int_Asy_Unit,1)
            [~, Pidx]=mink(DISM(p,:), Patch_Size);
            Patch_Ridx(p,:)=sort(Pidx);
        end

        %%%% NEED TO BE IMPROVED %%%%
        Patch = cell2struct(cell(length(fields), []), fields);
        n_Int_File=Int_File(Int_File(:,5)==0,1:4);
        for j=1:size(Patch_Ridx,1)
            Nei_Dist=pdist2(Int_Asy_Unit(j,1:3), n_Int_File(:,1:3));
            [~, Nei_ID]=min(Nei_Dist);
            Patch(j).Coord=[Int_Asy_Unit(Patch_Ridx(j,:),1:3); n_Int_File(Nei_ID, 1:3)];

            B=sum(Int_Asy_Unit(Patch_Ridx(j,:),5:7));
            Patch(j).Norm_V=B/norm(B);
            Patch(j).Res_Type=[Int_Asy_Unit(Patch_Ridx(j,:),4); n_Int_File(Nei_ID, 4)];
        end

        % Output
        writetable(struct2table(Patch),[PDB_ID,'/PatchBag_Result_',PDB_ID,'/',PDB_ID,'_Int_Patch.txt']);
    end
catch
    disp('Please check the Int.txt.');
end

clearvars -except PDB_ID CA_Coord Pool Patch_Size fields

%%%% nInt Part

try
    % Load Info of Interacting Residues
    nInt_File=table2array(readtable([PDB_ID,'/PatchBag_Result_',PDB_ID,'/',PDB_ID,'_nInt.txt']));
    nInt_Asy_Unit=nInt_File(:,1:4);

    if size(nInt_Asy_Unit,1)<6
        Patch = cell2struct(cell(length(fields), []), fields);
        writetable(struct2table(Patch),[PDB_ID,'/PatchBag_Result_',PDB_ID,'/',PDB_ID,'_nInt_Patch.txt']);
    else
        % Read Molecule Surface
        Ply_file = fopen([PDB_ID,'/',PDB_ID, '.ply']);
        while(1)
            line=fgetl(Ply_file);
            if length(line)>14 && strcmp(line(1:14),'element vertex')
                Vertex_No=str2double(line(16:end));
            elseif length(line)>12 && strcmp(line(1:12),'element face')
                Face_No=str2double(line(14:end));
            elseif length(line)==10 && strcmp(line(1:10),'end_header')
                RawText = fread(Ply_file,inf,'*char');
                splitLines = cell2mat(textscan(RawText, '%f %f %f %f %*[^\n]'));
                break;
            else
                continue;
            end
        end
        fclose(Ply_file);

        % Normal Vector Calculation
        Vertex_Table=splitLines(1:Vertex_No,1:3);
        CA_Idx_Table=pdist2(nInt_Asy_Unit(:,1:3), CA_Coord(:,1:3));

        for k=1:size(CA_Idx_Table,1)
            Res_Idx=CA_Coord(CA_Idx_Table(k,:)==0,4);
            CA_Position=CA_Coord(CA_Coord(:,4)==Res_Idx,1:3);
            Residue_Surface_DISM=pdist2(Pool(Pool(:,4)==Res_Idx,1:3), Vertex_Table);

            Cutoff=2;
            while(1)
                if Cutoff>8
                    disp('Please check the Cutoff.');
                    break;
                end

                if sum(sum(Residue_Surface_DISM<Cutoff)>0)>10
                    break;
                end
                Cutoff=Cutoff+1;
            end

            V=mean(Vertex_Table(sum(Residue_Surface_DISM<Cutoff)>0,:))-CA_Position;
            nInt_Asy_Unit(k,5:7)=V/norm(V);
        end

        %

        % Patch Extraction
        DISM=squareform(pdist(nInt_Asy_Unit(:,1:3)));
        Patch_Ridx=zeros(size(nInt_Asy_Unit,1),Patch_Size);
        for p=1:size(nInt_Asy_Unit,1)
            [~, Pidx]=mink(DISM(p,:), Patch_Size);
            Patch_Ridx(p,:)=sort(Pidx);
        end

        %%%% NEED TO BE IMPROVED %%%%
        Patch = cell2struct(cell(length(fields), []), fields);
        for j=1:size(Patch_Ridx,1)
            Patch(j).Coord=nInt_Asy_Unit(Patch_Ridx(j,:),1:3);
            B=sum(nInt_Asy_Unit(Patch_Ridx(j,:),5:7));
            Patch(j).Norm_V=B/norm(B);
            Patch(j).Res_Type=nInt_Asy_Unit(Patch_Ridx(j,:),4);
        end

        % Output
        writetable(struct2table(Patch),[PDB_ID,'/PatchBag_Result_',PDB_ID,'/',PDB_ID,'_nInt_Patch.txt']);
    end
catch
    disp('Please check the nInt.txt.');
end

clearvars -except PDB_ID

%
%
%

%%%%%%%% AAI-PatchBag %%%%%%%%

%%% AA_Idx Reduction
% 1-R2 as Similarity Distance 
load('Database/AA Index/AA_Index.mat')
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

%%% Default Setting
Max_Cluster_No=300;
P=perms(1:6);

%%%% Int Part

% Load Clust Info (Representative Patch Library)
load(['Database/Int_Patch/Clust_medoids_', num2str(Max_Cluster_No), '.mat'])

% Patch Classification
Selected_Sample(1).ID=PDB_ID;
Selected_Sample(1).Patch_code=[];
Selected_Sample(1).Patch_Res=[];
Selected_Sample(1).Bar_code=[];

try
    % Patch Files
    Patch_File=table2array(readtable([PDB_ID,'/PatchBag_Result_',PDB_ID,'/',PDB_ID,'_Int_Patch.txt']));
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
    Selected_Sample(1).Patch_code=[Patch_code, Unique_Count'];
    Selected_Sample(1).Patch_Res=Unique_Patch(:,25:31);

    % Patch Vectorlization
    Bar_code=zeros(Max_Cluster_No,1);
    for u=1:Max_Cluster_No
        Bar_code(u,1)= sum(Unique_Count(Patch_code==u));
    end
    Selected_Sample(1).Bar_code=Bar_code;
catch
    disp('Please check the Int_Patch.txt.');
end

% 

% Patch Vectorlization
load('Database/Int_Patch/Ruler.mat');

Count=[];
Merged_Table=[];
Patch_Res=Selected_Sample(1).Patch_Res;
if ~isempty(Patch_Res)
    for j=1:size(Patch_Res,1)
        Idx=Patch_Res(j,:);
        Idx(Idx==0)=8;
        Merged_Table=[Merged_Table, mean(AA_Index(:,Idx),2)];
    end
    Count=[Count; Selected_Sample(1).Patch_code(:,2)];
end

Completed_Table=[];
parfor k=1:size(Count,1)
    Completed_Table=[Completed_Table, repmat(Merged_Table(:,k),1,Count(k))];
end

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
end

Selected_Sample(1).Merged_AAIdx=[]; Acc=0;
if ~isempty(Selected_Sample(1).Patch_code)
    Acc_Start=Acc+1;
    Acc=Acc+size(Selected_Sample(1).Patch_code,1);
    Selected_Sample(1).Merged_AAIdx=[Selected_Sample(1).Patch_code(:,2), Selected_Sample(1).Patch_code(:,1), Hist(:,Acc_Start:Acc)'];
end

Merged_AAIdx=Selected_Sample;
Bin=size(Ruler,2)+1;

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
end

clear AA_Idx_Patch_code_RA;
AA_Idx_Patch_code_RA(1).ID=[];
AA_Idx_Patch_code_RA(1).Patch_code=[];

AA_Idx_Patch_code_RA(1).ID=AA_Idx_Patch_code(1).ID;
AA_Idx_Patch_code_RA(1).Patch_code=AA_Idx_Patch_code(1).Patch_code;
for v=1:size(Representative_Clust,1)
    AA_Idx_Patch_code_RA(1).(['Patch_code_', num2str(v)])=AA_Idx_Patch_code(1).(['Patch_code_', num2str(Representative_Clust(v))]);
end

save([PDB_ID,'/','PatchBag_Result_',PDB_ID,'/AA_Idx_Patch_code_',PDB_ID,'_Int_RA.mat'],'AA_Idx_Patch_code_RA')

clearvars -except PDB_ID AA_Index Max_Cluster_No P Representative_Clust

%

%%%% nInt Part

% Load Clust Info (Representative Patch Library)
load(['Database/nInt_Patch/Clust_medoids_', num2str(Max_Cluster_No), '.mat'])

% Patch Classification
Selected_Sample(1).ID=PDB_ID;
Selected_Sample(1).Patch_code=[];
Selected_Sample(1).Patch_Res=[];
Selected_Sample(1).Bar_code=[];

try
    % Patch Files
    Patch_File=table2array(readtable([PDB_ID,'/PatchBag_Result_',PDB_ID,'/',PDB_ID,'_nInt_Patch.txt']));
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
    Selected_Sample(1).Patch_code=[Patch_code, Unique_Count'];
    Selected_Sample(1).Patch_Res=Unique_Patch(:,22:27);

    % Patch Vectorlization
    Bar_code=zeros(Max_Cluster_No,1);
    for u=1:Max_Cluster_No
        Bar_code(u,1)= sum(Unique_Count(Patch_code==u));
    end
    Selected_Sample(1).Bar_code=Bar_code;
catch
    disp('Please check the nInt_Patch.txt.');
end

% 

% Patch Vectorlization
load('Database/nInt_Patch/Ruler.mat');

Count=[];
Merged_Table=[];
Patch_Res=Selected_Sample(1).Patch_Res;
if ~isempty(Patch_Res)
    for j=1:size(Patch_Res,1)
        Idx=Patch_Res(j,:);
        Idx(Idx==0)=8;
        Merged_Table=[Merged_Table, mean(AA_Index(:,Idx),2)];
    end
    Count=[Count; Selected_Sample(1).Patch_code(:,2)];
end

Completed_Table=[];
parfor k=1:size(Count,1)
    Completed_Table=[Completed_Table, repmat(Merged_Table(:,k),1,Count(k))];
end

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
end

Selected_Sample(1).Merged_AAIdx=[]; Acc=0;
if ~isempty(Selected_Sample(1).Patch_code)
    Acc_Start=Acc+1;
    Acc=Acc+size(Selected_Sample(1).Patch_code,1);
    Selected_Sample(1).Merged_AAIdx=[Selected_Sample(1).Patch_code(:,2), Selected_Sample(1).Patch_code(:,1), Hist(:,Acc_Start:Acc)'];
end

Merged_AAIdx=Selected_Sample;
Bin=size(Ruler,2)+1;

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
end

clear AA_Idx_Patch_code_RA;
AA_Idx_Patch_code_RA(1).ID=[];
AA_Idx_Patch_code_RA(1).Patch_code=[];

AA_Idx_Patch_code_RA(1).ID=AA_Idx_Patch_code(1).ID;
AA_Idx_Patch_code_RA(1).Patch_code=AA_Idx_Patch_code(1).Patch_code;
for v=1:size(Representative_Clust,1)
    AA_Idx_Patch_code_RA(1).(['Patch_code_', num2str(v)])=AA_Idx_Patch_code(1).(['Patch_code_', num2str(Representative_Clust(v))]);
end

save([PDB_ID,'/','PatchBag_Result_',PDB_ID,'/AA_Idx_Patch_code_',PDB_ID,'_nInt_RA.mat'],'AA_Idx_Patch_code_RA')

clearvars -except PDB_ID

%
%
%

%%%%%%%% Patch Distance Calculation %%%%%%%%

%%%% Int Part
% Load Input PDB ID
load([PDB_ID,'/','PatchBag_Result_',PDB_ID,'/AA_Idx_Patch_code_',PDB_ID,'_Int_RA.mat'])
Input=AA_Idx_Patch_code_RA;
clear AA_Idx_Patch_code_RA

% Load Database of AA_Idx_Patch_code.mat
load('Database/Int_Patch/AA_Idx_Patch_code_Int_5bin_RA.mat')
Dataset=AA_Idx_Patch_code_RA;
clear AA_Idx_Patch_code_RA

Count=size(Dataset,2);
Int_Result=zeros(length(fieldnames(Dataset))-2,Count);

% Int_Result
for r=1:(length(fieldnames(Dataset))-2)
    P1=Input(1).(['Patch_code_', num2str(r)]);
    if isempty(P1)
        P1=zeros(300,5);
    end
    t=0;
    localResult=zeros(1, Count);
    for m=1:Count
        P2=Dataset(m).(['Patch_code_', num2str(r)]);
        if isempty(P2)
            P2=zeros(300,5);
        end

        Diff=P2-P1;

        t=t+1;
        localResult(t) = norm(Diff);
    end
    Int_Result(r, :) = localResult;
end

%

%%%% nInt Part
% Load Input PDB ID
load([PDB_ID,'/','PatchBag_Result_',PDB_ID,'/AA_Idx_Patch_code_',PDB_ID,'_nInt_RA.mat'])
Input=AA_Idx_Patch_code_RA;
clear AA_Idx_Patch_code_RA

% Load Database of AA_Idx_Patch_code.mat
load('Database/nInt_Patch/AA_Idx_Patch_code_nInt_5bin_RA.mat')
Dataset=AA_Idx_Patch_code_RA;
clear AA_Idx_Patch_code_RA

Count=size(Dataset,2);
nInt_Result=zeros(length(fieldnames(Dataset))-2,Count);

% nInt_Result
for r=1:(length(fieldnames(Dataset))-2)
    P1=Input(1).(['Patch_code_', num2str(r)]);
    if isempty(P1)
        P1=zeros(300,5);
    end
    t=0;
    localResult=zeros(1, Count);
    for m=1:Count
        P2=Dataset(m).(['Patch_code_', num2str(r)]);
        if isempty(P2)
            P2=zeros(300,5);
        end

        Diff=P2-P1;

        t=t+1;
        localResult(t) = norm(Diff);
    end
    nInt_Result(r, :) = localResult;
end

clearvars -except PDB_ID Int_Result nInt_Result

%
%
%

%%%%%%%% Feature Collection %%%%%%%%
%%% 2472 Database
load('Database/Sample_Info_2472.mat')

%

%%% Sequence DISM Calculation
% Input Sequence
PDB=pdbread([PDB_ID, '/', PDB_ID,'.pdb']);
Input_Seq=strrep(strrep(PDB.Sequence(1).Sequence,'U','C'),'X','G');

% Seq DISM Calculation
Seq_Linear_DISM=zeros(1,size(Sample_Info_2472,2));
for n=1:size(Sample_Info_2472,2)
    Seq_Linear_DISM(1,n)=seqpdist({Input_Seq; Sample_Info_2472(n).Seq},'Method', 'p-distance');
end

%

%%% RMSD DISM Calculation
% Input Sequence Alpha Carbon Coordinate
% Sequence and Alpha Carbon Position Extraction
File=[PDB_ID, '/', PDB_ID,'.pdb'];
Input(1).CA_Coord=PDB_CA_Extracting(File);

% RMSD DISM Calculation
RMSD_Linear_DISM=zeros(1,size(Sample_Info_2472,2));

Seq_1=strrep(strrep(Input(1).CA_Coord.resName{1},'U','C'),'X','G');
Coord_1=[Input(1).CA_Coord.X, Input(1).CA_Coord.Y, Input(1).CA_Coord.Z];

for q=1:size(Sample_Info_2472,2)

    Sample_Info_2472(1).CA_Coord;

    Seq_2=Sample_Info_2472(q).CA_Coord.resName{1};
    Coord_2=[Sample_Info_2472(q).CA_Coord.X, Sample_Info_2472(q).CA_Coord.Y, Sample_Info_2472(q).CA_Coord.Z];

    % Global Sequence Alignment to Extract Superimposed CA for RMSD Calculation
    [~, Alignment] = nwalign(Seq_1,Seq_2);
    alignidex=find((Alignment(1,:)~='-')&(Alignment(3,:)~='-'));
    aseq1=find(Alignment(1,:)=='-');
    aseq2=find(Alignment(3,:)=='-');

    aseq1_idx=[]; aseq2_idx=[];
    for o=1:length(alignidex)
        aseq1_idx=[aseq1_idx, alignidex(o)-sum(aseq1<alignidex(o))];
        aseq2_idx=[aseq2_idx, alignidex(o)-sum(aseq2<alignidex(o))];
    end

    % RMSD
    [~,~,RMSD]=CoordiExam_AC(Coord_1(aseq1_idx,:), Coord_2(aseq2_idx,:));
    RMSD_Linear_DISM(1,q)=RMSD;
end

%

%%% Jerhoud DISM Calculation
% Calculate Jerhoud 3DZD of Input Protein
%%1 Ply2off
Ply_file = fopen([PDB_ID, '/', PDB_ID, '.ply']);
while(1)
    line=fgetl(Ply_file);
    if length(line)>14 && strcmp(line(1:14),'element vertex')
        Vertex_No=str2double(line(16:end));
    elseif length(line)>12 && strcmp(line(1:12),'element face')
        Face_No=str2double(line(14:end));
    elseif length(line)==10 && strcmp(line(1:10),'end_header')
        RawText = fread(Ply_file,inf,'*char');
        splitLines = cell2mat(textscan(RawText, '%f %f %f %f %*[^\n]'));
        break;
    else
        continue;
    end
end
fclose(Ply_file);

Vertex_Table=splitLines(1:Vertex_No,1:3);
Face_Table=splitLines(Vertex_No+1:end,1:4);

% Output Off File
Off_file = fopen([PDB_ID, '/', PDB_ID, '.off'],'w');
fprintf(Off_file, 'OFF\n');
fprintf(Off_file, [num2str(Vertex_No), ' ', num2str(Face_No), ' 0\n']);
Vertex = sprintf('%.4f %.4f %.4f\n', Vertex_Table');
Face = sprintf('%d %d %d %d\n', Face_Table');
fprintf(Off_file, [Vertex Face]);
fclose(Off_file);

%%2 MakeShpae
MS_CMD=['./MakeShape -l ', PDB_ID, '/', PDB_ID, '.off -cr 0.8 -o ', PDB_ID, '/', PDB_ID, '_0.8.off'];
try 
   [status, ~] = system(MS_CMD);
catch
    disp('Please install the Jernoud 3DZD program and check the glibcxx version.');
end

%%3 Shape2Zernike
S2Z_CMD=['./Shape2Zernike -o ', PDB_ID, '/', PDB_ID, '_0.8_50_a8.zm -t 8 -a 8 50 ', PDB_ID, '/', PDB_ID, '_0.8.off -v'];
try 
    [status, ~]=system(S2Z_CMD);
catch
    disp('Please install the Jernoud 3DZD program and check the glibcxx version.');
end

%%4 Jerhoud ZM Restoration
% Zernike Moment Table
scale=0.8;
Count = 0; Order=50;
Table = zeros((Order+1)*(Order+2)*(Order+3)/6,5);

for n=0:Order
    l_start=mod(n,2);
    for l=l_start:2:n
        for m=-l:l
            Count = Count + 1;
            Table(Count,1) = n;
            Table(Count,2) = l;
            Table(Count,3) = m;
        end
    end
end

% Jerhoud Zernike Moment
Jerhoud_ZM=fopen([PDB_ID, '/', PDB_ID, '_', num2str(scale), '_', num2str(Order), '_a8.zm'],'r');
while(1)
    line=fgetl(Jerhoud_ZM);
    if length(line)>5 && strcmp(line(1:5),'ORTHO')
        RawText = fread(Jerhoud_ZM,inf,'*char');
        splitLines = cell2mat(textscan(RawText, '%f %f %f %f %f %*[^\n]'));
        break;
    else
        continue;
    end
end
splitLines(isnan(splitLines))=0;
fclose(Jerhoud_ZM);

% Zernike Moment Restoration
for r=1:size(Table,1)
    Idx_raw=find(splitLines(:,1)==Table(r,1) & splitLines(:,2)==Table(r,2) & splitLines(:,3)==Table(r,3));
    Idx_abs=find(splitLines(:,1)==Table(r,1) & splitLines(:,2)==Table(r,2) & splitLines(:,3)==abs(Table(r,3)));

    if ~isempty(Idx_raw)
        Table(r,4:5)=splitLines(Idx_raw,4:5);
    elseif isempty(Idx_raw) && ~isempty(Idx_abs)
        Table(r,4:5)=(-1)^(splitLines(Idx_abs,3))*[splitLines(Idx_abs,4), -splitLines(Idx_abs,5)];
    else
        Table(r,4:5)=[0, 0];
    end
end

% Output
Jerhoud_ZM_file = fopen([PDB_ID, '/', PDB_ID, '_Jerhoud_All.zm'],'w');
fprintf(Jerhoud_ZM_file, 'ZM\n');
fprintf(Jerhoud_ZM_file, ['ORTHO ', num2str(Order), ' COMPLEX\n']);
ZM = sprintf('%d %d %d %.8e %.8e\n', Table');
fprintf(Jerhoud_ZM_file, ZM);
fclose(Jerhoud_ZM_file);

%%5 Moments2Descriptors
Zernike_Moments=Table;

% Re-rder Zernike Moments
Ridx=Zernike_Moments(:,3)<=0;
Ordered_ZM=Zernike_Moments(Ridx,:);
Ordered_ZM(:,3)=Ordered_ZM(:,2)+Ordered_ZM(:,3);

% Compute Descriptors from Zernile Moments
Order=max(Zernike_Moments(:,1));
nInvariants=floor(((Order+2)/2)^2);
Input_Descriptor = zeros(nInvariants,1);

Count=0;
for n=0:Order
    for l=mod(n,2):2:n

        % Sum_Moments
        sum_mom=0;

        for m=-l:l

            % The ZC_nlm moment
            Idx=find(Ordered_ZM(:,1)==n & Ordered_ZM(:,2)==l & Ordered_ZM(:,3)==abs(m));
            mom = complex(Ordered_ZM(Idx,4), Ordered_ZM(Idx,5));

            % Conjugate if m is negative
            if m<0
                mom = conj(mom);
                if mod(abs(m),2)
                    mom = -1*mom;
                end
            end
            sum_mom = sum_mom + norm(mom)^2;
        end

        % Norm of Sum_Moments
        Count=Count+1;
        Input_Descriptor(Count,1)=sqrt(sum_mom);
    end
end

clearvars -except PDB_ID Int_Result nInt_Result Sample_Info_2472 Seq_Linear_DISM RMSD_Linear_DISM Input_Descriptor

% Jerhoud 3DZD DISM Calculation
Jerhoud_Linear_DISM=zeros(1,size(Sample_Info_2472,2));

ZCD_1=Input_Descriptor;
for n=1:size(Sample_Info_2472,2)

    ZCD_2=Sample_Info_2472(n).Jerhoud_InvD;

    %%%% ZCD Distance
    % the Euclidean distance between two shape descriptors of given structure
    Jerhoud_Linear_DISM(1,n) = norm(ZCD_1-ZCD_2, 2);
end

clearvars -except PDB_ID Int_Result nInt_Result Sample_Info_2472 Seq_Linear_DISM RMSD_Linear_DISM Descriptors Jerhoud_Linear_DISM

%

%%% Kihara DISM Calculation
% Calculate Kihara 3DZD of Input Protein
Inv_Descriptors=readtable([PDB_ID, '/', PDB_ID, '.obj.grid.inv'], 'FileType', 'text');
Input_Kihara_Inv=Inv_Descriptors(2:end,1);

% Kihara 3DZD DISM Calculation
Kihara_Linear_DISM=zeros(1,size(Sample_Info_2472,2));

ZCD_1=table2array(Input_Kihara_Inv);
for n=1:size(Sample_Info_2472,2)

    ZCD_2=Sample_Info_2472(n).Kihara_InvD;

    %%%% ZCD Distance
    % the Euclidean distance between two shape descriptors of given structure
    Kihara_Linear_DISM(1,n) = norm(ZCD_1-ZCD_2, 2);
end

clearvars -except PDB_ID Int_Result nInt_Result Sample_Info_2472 Seq_Linear_DISM RMSD_Linear_DISM Descriptors Jerhoud_Linear_DISM Kihara_Linear_DISM

%

%%% Molecular Weight / Patch Num / Interaction Matrix Calculation
%%1 Molecular Weight
% Input Molecular Weight
PDB=pdbread([PDB_ID, '/', PDB_ID,'.pdb']);
Input_Seq=PDB.Sequence.Sequence;
Input_MW=molweight(strrep(strrep(Input_Seq,'U','C'),'X','G'));

%%2 Patch Num
load([PDB_ID,'/','PatchBag_Result_',PDB_ID,'/AA_Idx_Patch_code_',PDB_ID,'_Int_RA.mat'])
Input_Int_Patch_No=sum(AA_Idx_Patch_code_RA(1).Patch_code);

load([PDB_ID,'/','PatchBag_Result_',PDB_ID,'/AA_Idx_Patch_code_',PDB_ID,'_nInt_RA.mat'])
Input_nInt_Patch_No=sum(AA_Idx_Patch_code_RA(1).Patch_code);

Input_MWP=[Input_MW; Input_Int_Patch_No; Input_nInt_Patch_No];

%

Count=size(Sample_Info_2472,2);
MWP_Result=zeros(4,Count);

for r=1:size(Input_MWP,1)
    P1=Input_MWP(r);
    for n=1:Count
        P2=Sample_Info_2472(n).MW_Patch_Num(r);
        Distance=abs(P1-P2);
        MWP_Result(r,n)=Distance;
    end
end

%%3 Interaction Matrix
% Input_ResM
Int_File=table2array(readtable([PDB_ID,'/PatchBag_Result_',PDB_ID,'/',PDB_ID,'_Int.txt']));
Int_File_Res=Int_File(:,4); Int_File_Res(Int_File_Res==0)=1;

if ~isempty(Int_File)
    Patch_File=table2array(readtable([PDB_ID,'/PatchBag_Result_',PDB_ID,'/',PDB_ID,'_Int_Patch.txt']));

    if ~isempty(Patch_File)
        Res_Info=Patch_File(:,25:31); Res_Info(Res_Info==0)=1;

        Res_Matrix=zeros(20);
        for j=1:size(Patch_File,1)
            Res_Matrix(Int_File_Res(j,1),Res_Info(j,7))=Res_Matrix(Int_File_Res(j,1),Res_Info(j,7))+1;
        end
        Input_ResM=Res_Matrix;
    else
        disp('Please check the Int.txt and Int_Patch.txt.');
    end
else
    disp('Please check the Int.txt and Int_Patch.txt.');
end

% Output
P1=Input_ResM;
if isempty(P1)
    P1=zeros(20,20);
end
for n=1:size(Sample_Info_2472,2)
    P2=Sample_Info_2472(n).Res_Matrix;
    if isempty(P2)
        P2=zeros(20,20);
    end

    Distance=norm(P1-P2);
    MWP_Result(4,n)=Distance;
end

clearvars -except PDB_ID Int_Result nInt_Result Seq_Linear_DISM RMSD_Linear_DISM Descriptors Jerhoud_Linear_DISM Kihara_Linear_DISM MWP_Result

%

%%%% Feature Merge
Feature=[ones(1,size(Int_Result,2));
        Seq_Linear_DISM; RMSD_Linear_DISM; Jerhoud_Linear_DISM; Kihara_Linear_DISM; 
        MWP_Result; 
        Int_Result; nInt_Result];

clearvars -except PDB_ID Feature

%
%
%

%%%% Ranking
% AAI-PatchBag-Based Dissimilarity (Score)
load('Database/Rank_Coeff.mat')
load('Database/Significant_Coeff.mat')

Reduced_Feature=Feature(p_value,:);
AAI_PatchBag_Score=-Coeff'*Reduced_Feature;
[score,row]=sort(AAI_PatchBag_Score,'ascend');

% Database of 2472 Qualified Proteins
load('Database/Sample_Info_2472.mat')

Datacase_ID={Sample_Info_2472.ID}';
AAI_PatchBag_Ranking=Datacase_ID(row);

save('AAI_PatchBag_Ranking.mat','AAI_PatchBag_Ranking');
