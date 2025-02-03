%%%% Packing Simulation
%%% Default Settings
% Establish Deposited File of Packing Results
mkdir('AF3_MPCP_Result/P43212_Packing_Pose');

% Combination List for Translation Table
itr=0; lmt=20;
List=zeros((lmt*2+1)^3,3);
for o=-lmt:lmt
    for p=-lmt:lmt
        for q=-lmt:lmt
            itr=itr+1;
            List(itr,1:3)=[o,p,q];
        end
    end
end

% Distance Cutoff (A) of Interacting Residues
Threshold=10;

% Chain ID List
ChainID_List = {'A';'B';'C';'D';'E';'F';'G';'I';'J';'K';
                'L';'M';'N';'O';'P';'Q';'R';'S';'T';'U';
                'V';'W';'X';'Y';'Z';'#';'$';'%';'/';'\'};

%%% Packing Simulation
% Load MXRGA Result
load('AF3_SurfOrca_Result/AF3_P43212_MXRGA_Result.mat')

% Erroneous Index
Error_ID=[];

% Packing Simulation
for r=1:size(MXRGA,2)

    if ~isempty(MXRGA(r).SurfOrca_Cell_Pmt)

        % Establish Deposited File
        if exist(['AF3_MPCP_Result/P43212_Packing_Pose/',MXRGA(r).PDB_ID],'dir')==0
            mkdir(['AF3_MPCP_Result/P43212_Packing_Pose/',MXRGA(r).PDB_ID]);
        end

        %%% Load Corresponding PDB File
        File=fopen(['PDB/', MXRGA(r).PDB_ID,'.pdb']);

        % Symmetry Operator
        Operator=[];
        while(1)
            line=fgetl(File);
            if isequal(sscanf(line(1:10),'%c'),'REMARK 290')
                SMTRY=sscanf(line(25:79),'%f %f %f');
                Operator=[Operator; SMTRY'];
            elseif isequal(sscanf(line(1:10),'%c'),'REMARK 300')
                break;
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
        PDB=pdbread(['PDB/', MXRGA(r).PDB_ID,'.pdb']);
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
        for i=1:Sym_No
            SMTRY_Unit=Asy_Unit*Operator(3*i-2:3*i,1:3)'+ Operator(3*i-2:3*i,4)';

            % Shift
            Shift=List*Unit_Cell_extent';
            Translation_Table=mean(Asy_Unit)-(mean(SMTRY_Unit)+Shift);

            % Store Shifted SMTRY_Unit
            [~,idx]=min(sqrt(sum(Translation_Table.^2,2)));
            SMTRY_Unit=SMTRY_Unit+Shift(idx,:);
            Unit=[Unit; SMTRY_Unit];
        end

        % Unit Box
        s=2;
        Unit_Size=size(Unit,1);
        Unit_Box=zeros(Unit_Size*(s*2+1)^3,4);
        Unit_Box(1:Unit_Size,1:3)=Unit;

        t=1;
        for x=-s:s
            for y=-s:s
                for z=-s:s
                    if ~(x==0 && y==0 && z==0)
                        NB_Unit=Unit+[x,y,z]*Unit_Cell_extent';
                        Unit_Box(Unit_Size*t+1:Unit_Size*(t+1),1:3)=NB_Unit;
                        Unit_Box(Unit_Size*t+1:Unit_Size*(t+1),4)=t;
                        t=t+1;
                    end
                end
            end
        end


        % Reduced Unit Box
        ASU_Up_B=max(Asy_Unit,[],1)+Threshold;
        ASU_Low_B=min(Asy_Unit,[],1)-Threshold;
        X=Unit_Box(:,1)>ASU_Low_B(1)&Unit_Box(:,1)<ASU_Up_B(1);
        Y=Unit_Box(:,2)>ASU_Low_B(2)&Unit_Box(:,2)<ASU_Up_B(2);
        Z=Unit_Box(:,3)>ASU_Low_B(3)&Unit_Box(:,3)<ASU_Up_B(3);

        % Identify Interacting Neighbors
        AU_Size=size(Asy_Unit,1);
        Mo_Idx=reshape(repmat(1:Sym_No*(s*2+1)^3,AU_Size,1),AU_Size*Sym_No*(s*2+1)^3,1);

        Reduced_Unit_Box=Unit_Box(X&Y&Z,1:3);
        Reduced_Neighbor_Box=Reduced_Unit_Box(AU_Size+1:end,:);

        Reduced_Mo_Idx=Mo_Idx(X&Y&Z);
        Reduced_Neighbor_Mo_Idx=Reduced_Mo_Idx(AU_Size+1:end,:);

        % Unique Neighbor ID (Column-wise)
        DISM=pdist2(Asy_Unit, Reduced_Neighbor_Box);

        Neighbor_Idx=sum(DISM<Threshold,1);
        UNI=unique(Reduced_Neighbor_Mo_Idx(Neighbor_Idx~=0));

        Xtal_Packing_Set=zeros(AU_Size*(size(UNI,1)+1),4);
        Xtal_Packing_Set(1:AU_Size,1:3)=Asy_Unit;

        Cent_Pos=zeros(size(UNI,1)+1,5);
        Cent_Pos(1,1:4)=[mean(Asy_Unit), 1];

        for m=1:size(UNI,1)
            Xtal_Packing_Set(AU_Size*m+1:AU_Size*(m+1),1:3)=Unit_Box(Mo_Idx==UNI(m),1:3);
            Xtal_Packing_Set(AU_Size*m+1:AU_Size*(m+1),4)=UNI(m);

            Cent_Pos(m+1,1:3)=mean(Xtal_Packing_Set(AU_Size*m+1:AU_Size*(m+1),1:3));
            Cent_Pos(m+1,4)=UNI(m);
        end
        Cent_Pos(:,5)=pdist2(Cent_Pos(1,1:3), Cent_Pos(:,1:3))';

        %
        %
        %

        % Load TopRank PDB File
        Sim_PDB=pdbread(['AF3_SurfOrca_Result/P43212_Refined_TopRank_300_Pose/',MXRGA(r).PDB_ID,'/',MXRGA(r).SiXt_ID,'_Refined.pdb']);

        % Extract Coordinates from Output TopRank PDB File
        CA_Idx=strcmp({Sim_PDB.Model.Atom.AtomName}, 'CA');
        Model=[Sim_PDB.Model.Atom.X; Sim_PDB.Model.Atom.Y; Sim_PDB.Model.Atom.Z]';

        % CA Coordinate of Opt_Results ([toXYZ; toXYZ_pair; fromXYZ; fromXYZ_pair])
        Sim_PDB_CA_Idx=strcmp({Sim_PDB.Model(1).Atom.AtomName}, 'CA')';
        CA_Model=Model(Sim_PDB_CA_Idx,:);
        len=size(CA_Model,1)/4;

        toXYZ=CA_Model(len*0+1:len*1,:);
        toXYZ_pair=CA_Model(len*1+1:len*2,:);

        fromXYZ=CA_Model(len*2+1:len*3,:);
        fromXYZ_pair=CA_Model(len*3+1:len*4,:);

        % CA Coordinate of Packing Symmetric Unit
        [R1,T1]=CoordiExam(fromXYZ, toXYZ_pair);
        Opp_Pairs_1=[toXYZ; toXYZ_pair]*R1+T1;

        [R2,T2]=CoordiExam(toXYZ, fromXYZ_pair);
        Opp_Pairs_2=[fromXYZ; fromXYZ_pair]*R2+T2;

        Packing_Sym=[toXYZ; toXYZ_pair; fromXYZ; fromXYZ_pair; Opp_Pairs_1; Opp_Pairs_2];

        % Unit Box
        Sim_Size=size(Packing_Sym,1);
        Sim_Unit_Box=zeros(Sim_Size*3^3,4);
        Sim_Unit_Box(1:Sim_Size,1:3)=Packing_Sym;

        u=1;
        for h=-1:1
            for k=-1:1
                for l=-1:1
                    if ~(h==0 && k==0 && l==0)
                        NB_Unit=Packing_Sym+[h,k,l].*MXRGA(r).SurfOrca_Cell_Pmt;
                        Sim_Unit_Box(Sim_Size*u+1:Sim_Size*(u+1),1:3)=NB_Unit;
                        Sim_Unit_Box(Sim_Size*u+1:Sim_Size*(u+1),4)=u;
                        u=u+1;
                    end
                end
            end
        end

        % Reduced Sim Unit Box
        Sim_Up_B=max(toXYZ,[],1)+Threshold*1.5;
        Sim_Low_B=min(toXYZ,[],1)-Threshold*1.5;
        Reduced_X=Sim_Unit_Box(:,1)>Sim_Low_B(1)&Sim_Unit_Box(:,1)<Sim_Up_B(1);
        Reduced_Y=Sim_Unit_Box(:,2)>Sim_Low_B(2)&Sim_Unit_Box(:,2)<Sim_Up_B(2);
        Reduced_Z=Sim_Unit_Box(:,3)>Sim_Low_B(3)&Sim_Unit_Box(:,3)<Sim_Up_B(3);

        % Identify Interacting Neighbors
        Reduced_Sim_Unit_Box=Sim_Unit_Box(Reduced_X&Reduced_Y&Reduced_Z,1:3);
        Sim_Mo_Idx=reshape(repmat(1:Sym_No*3^3,len,1),Sim_Size*3^3,1);

        Reduced_Sim_Mo_Idx=Sim_Mo_Idx(Reduced_X&Reduced_Y&Reduced_Z);
        Sim_DISM=pdist2(toXYZ, Reduced_Sim_Unit_Box);

        URMI=unique(Reduced_Sim_Mo_Idx);
        for n=1:size(URMI,1)
            URMI(n,2)=min(Sim_DISM(:,Reduced_Sim_Mo_Idx==URMI(n)),[],'all');
        end

        %

        %%%% Sequence Alignment for Further RMSD Calculation
        % Sequence Extraction
        PDB_Seq=aminolookup(strrep(strrep([PDB.Model(1).Atom(PDB_CA_Idx).resName], 'SEC', 'CYS'), 'UNK', 'GLY'));
        Sim_PDB_Seq=aminolookup([Sim_PDB.Model.Atom(Sim_PDB_CA_Idx).resName]);
        Mono_PDB_Seq=Sim_PDB_Seq(1:len);

        % Global Sequence Alignment to Extract Superimposed CA for RMSD Calculation
        [~, Alignment] = nwalign(PDB_Seq, Mono_PDB_Seq);
        alignidex=find((Alignment(1,:)~='-')&(Alignment(3,:)~='-'));
        aseq1=find(Alignment(1,:)=='-');
        aseq2=find(Alignment(3,:)=='-');

        aseq1_idx=[]; aseq2_idx=[];
        for o=1:length(alignidex)
            aseq1_idx=[aseq1_idx, alignidex(o)-sum(aseq1<alignidex(o))];
            aseq2_idx=[aseq2_idx, alignidex(o)-sum(aseq2<alignidex(o))];
        end

        %

        [~,Ord]=sort(Cent_Pos(:,5));
        Rank=Cent_Pos(Ord,4);

        Count=1;
        Template=Asy_Unit(aseq1_idx,:);
        Target=toXYZ(aseq2_idx,:);

        Potentail_Neighbor=URMI(logical([0; URMI(2:end,2)<25]),:);

        try
            while(1)
                Count=Count+1;
                if Count == (m+2), break; end

                Addon_Template=Xtal_Packing_Set(Xtal_Packing_Set(:,4)==Rank(Count),1:3);
                Template=[Template; Addon_Template(aseq1_idx,:)]; %#ok<AGROW>

                for w=1:size(Potentail_Neighbor,1)
                    Addon_Target=Sim_Unit_Box(Sim_Mo_Idx==Potentail_Neighbor(w,1),1:3);
                    [~,~,RMSD]=CoordiExam(Template, [Target; Addon_Target(aseq2_idx,:)]);
                    Potentail_Neighbor(w,2)=RMSD;
                end

                [val,pos]=min(Potentail_Neighbor(:,2));
                Elong_Target=Sim_Unit_Box(Sim_Mo_Idx==Potentail_Neighbor(pos,1),1:3);
                Target=[Target; Elong_Target(aseq2_idx,:)]; %#ok<AGROW>
                Potentail_Neighbor(pos,:)=[];
            end

            %
            %
            %

            %%%% Output
            %%% Crystal Packing
            % Benefit for Following PDB writing
            Xtal_Output.Model.Atom=PDB.Model(1).Atom;

            % Remove Non-Alignment Residues
            Xtal_Output_Res_No=[Xtal_Output.Model.Atom.resSeq];
            Xtal_Res_No_List=unique(Xtal_Output_Res_No);

            RM_ID_1=setdiff(1:size(PDB_CA_Idx,1),aseq1_idx);
            if ~isempty(RM_ID_1)
                Logic_ID1=zeros(length(RM_ID_1),size(Xtal_Output_Res_No,2));
                for idx_1=1:length(RM_ID_1)
                    Logic_ID1(idx_1,:)=Xtal_Output_Res_No==Xtal_Res_No_List(RM_ID_1(idx_1));
                end
                Xtal_Output.Model.Atom(sum(Logic_ID1)~=0)=[];
            end

            Xtal_Packing_Output.Model.Atom=repmat(Xtal_Output.Model.Atom(:),m+1,1)';

            % Re-write AtomSerNo/chainID/Coordinate
            Xtal_Model=[Xtal_Output.Model(1).Atom.X; Xtal_Output.Model(1).Atom.Y; Xtal_Output.Model(1).Atom.Z]';
            len_aseq1=length(aseq1_idx);
            len_xtal=size(Xtal_Model,1);
            for TA=0:m
                % Atom Serial Number
                SerNo=num2cell(len_xtal*TA+1:len_xtal*(TA+1));
                [Xtal_Packing_Output.Model.Atom(len_xtal*TA+1:len_xtal*(TA+1)).AtomSerNo]=SerNo{:};

                % Chain ID
                ChainID=repmat(ChainID_List(TA+1),len_xtal,1);
                [Xtal_Packing_Output.Model.Atom(len_xtal*TA+1:len_xtal*(TA+1)).chainID]=ChainID{:};

                % Atom Coordinate
                [R,T]=CoordiExam(Asy_Unit(aseq1_idx,:), Template(len_aseq1*TA+1:len_aseq1*(TA+1),:));

                Coordinate=num2cell(Xtal_Model*R+T);
                [Xtal_Packing_Output.Model.Atom(len_xtal*TA+1:len_xtal*(TA+1)).X]=Coordinate{:,1};
                [Xtal_Packing_Output.Model.Atom(len_xtal*TA+1:len_xtal*(TA+1)).Y]=Coordinate{:,2};
                [Xtal_Packing_Output.Model.Atom(len_xtal*TA+1:len_xtal*(TA+1)).Z]=Coordinate{:,3};
            end

            % New Seq_No
            Seq_No=[Xtal_Packing_Output.Model(1).Atom.resSeq];
            UniSeq_No=unique(Seq_No);

            New_Seq_No=zeros(size(Seq_No));
            for UN=1:size(UniSeq_No,2)
                New_Seq_No(Seq_No==UniSeq_No(UN))=UN;
            end

            New_Seq_No_cell=num2cell(New_Seq_No);
            [Xtal_Packing_Output.Model.Atom.resSeq]=New_Seq_No_cell{:};

            warning('off','all');
            pdbwrite(['AF3_MPCP_Result/P43212_Packing_Pose/',MXRGA(r).PDB_ID,'/',MXRGA(r).SiXt_ID, '_Xtal_Packing.pdb'], Xtal_Packing_Output);

            %
            %
            %

            %%% Simulated Packing
            % Benefit for Following PDB writing
            Mono_Size=size(Sim_PDB.Model.Atom,2)/4;
            Mono_Output.Model.Atom=Sim_PDB.Model.Atom(1:Mono_Size);

            % Remove Non-Alignment Residues
            Mono_Output_Res_No=[Mono_Output.Model.Atom.resSeq];
            Res_No_List=unique(Mono_Output_Res_No);

            RM_ID2=setdiff(1:len,aseq2_idx);
            if ~isempty(RM_ID2)
                Logic_ID2=zeros(length(RM_ID2),size(Mono_Output_Res_No,2));
                for idx2=1:length(RM_ID2)
                    Logic_ID2(idx2,:)=Mono_Output_Res_No==Res_No_List(RM_ID2(idx2));
                end
                Mono_Output.Model.Atom(sum(Logic_ID2)~=0)=[];
            end

            Packing_Output.Model.Atom=repmat(Mono_Output.Model.Atom(:),m+1,1)';

            % Re-write AtomSerNo/chainID/Coordinate
            Mono_Model=[Mono_Output.Model(1).Atom.X; Mono_Output.Model(1).Atom.Y; Mono_Output.Model(1).Atom.Z]';
            len_aseq2=length(aseq2_idx);
            len_momo=size(Mono_Model,1);
            for TA=0:m
                % Atom Serial Number
                SerNo=num2cell(len_momo*TA+1:len_momo*(TA+1));
                [Packing_Output.Model.Atom(len_momo*TA+1:len_momo*(TA+1)).AtomSerNo]=SerNo{:};

                % Chain ID
                ChainID=repmat(ChainID_List(TA+1),len_momo,1);
                [Packing_Output.Model.Atom(len_momo*TA+1:len_momo*(TA+1)).chainID]=ChainID{:};

                % Atom Coordinate
                [R,T]=CoordiExam(toXYZ(aseq2_idx,:), Target(len_aseq2*TA+1:len_aseq2*(TA+1),:));

                Coordinate=num2cell(Mono_Model*R+T);
                [Packing_Output.Model.Atom(len_momo*TA+1:len_momo*(TA+1)).X]=Coordinate{:,1};
                [Packing_Output.Model.Atom(len_momo*TA+1:len_momo*(TA+1)).Y]=Coordinate{:,2};
                [Packing_Output.Model.Atom(len_momo*TA+1:len_momo*(TA+1)).Z]=Coordinate{:,3};
            end

            % New Seq_No
            Seq_No=[Packing_Output.Model(1).Atom.resSeq];
            UniSeq_No=unique(Seq_No);

            New_Seq_No=zeros(size(Seq_No));
            for UN=1:size(UniSeq_No,2)
                New_Seq_No(Seq_No==UniSeq_No(UN))=UN;
            end

            New_Seq_No_cell=num2cell(New_Seq_No);
            [Packing_Output.Model.Atom.resSeq]=New_Seq_No_cell{:};

            warning('off','all');
            pdbwrite(['AF3_MPCP_Result/P43212_Packing_Pose/',MXRGA(r).PDB_ID,'/',MXRGA(r).SiXt_ID, '_Sim_Packing.pdb'], Packing_Output);
        catch
            Error_ID=[Error_ID; r];
        end
    else
        Error_ID=[Error_ID; r];
    end
    r
end