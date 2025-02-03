% The Only Exception
% I 41/a Space Group

%%%% Default Settings
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

% Distance Cutoff (A) of Interacting Residues
Threshold=10;

% Establish Deposited Folder
mkdir('Reduced_Sample_Int_Res_Results');

%

%%%% Cull Coordinates of Interacting Residues
% List of Sample ID
% Filename='Qualified_Reduced_Sample_3737(Exclude 4YKN).txt';
Filename='target_entry.txt';
ID_List=fopen(Filename);
sp=1;Sample(1).ID=[];
while (1)
    line=fgetl(ID_List);
    if line==-1, break, end
    Sample(sp).ID=sscanf(line,'%c');
    sp=sp+1;
end
fclose(ID_List);

Scale_Table=Sample;
Scale_Table(1).Scale_Factor=[];
Error_ID=[];

parfor t=1:size(Sample,1)

    try
        %%% Load Corresponding PDB File
        PDB_ID=lower(Sample(t).ID);
        File=fopen(['PDB/', PDB_ID,'.pdb']);

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
        PDB=pdbread(['PDB/', PDB_ID,'.pdb']);
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

        if PDB.Sequence.NumOfResidues<20
            Error_ID=[Error_ID; t];
        elseif size(Asy_Unit,1)<15
            Error_ID=[Error_ID; t];
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
            writetable(struct2table(Output),['Reduced_Sample_Int_Res_Results/',lower(Sample(t).ID),'.txt']);

            % Scale Record
            Asy_Unit_Scale=pdist2(mean(Asy_Unit,1),Asy_Unit);
            Asy_Coord_Scale=pdist2(mean(Asy_Coord,1),Asy_Coord);
            NB_Coord_Scale=pdist2(mean(NB_Coord,1),NB_Coord);
            Scale_Table(t).Scale_Factor=[max(Asy_Unit_Scale), max(Asy_Coord_Scale), max(NB_Coord_Scale)];
        end
    catch
        Error_ID=[Error_ID; t];
    end
    t
end

writetable(struct2table(Scale_Table),'Reduced_Sample_Int_Res_Results/Scale_Table.xlsx');