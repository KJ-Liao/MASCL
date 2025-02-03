% Elapsed time is 1945.294211 seconds.
% The Only Exception
% I 41/a Space Group

%%%% Default Settings
% Residue Table
Res={'GLY';'ALA';'VAL';'ILE';'LEU';
    'SER';'THR';'ASP';'ASN';'GLU';
    'GLN';'LYS';'ARG';'CYS';'MET';
    'PHE';'TYR';'TRP';'HIS';'PRO'};

% Establish Deposited Folder
mkdir('Reduced_Sample_nInt_Res_Results');

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

% Memory Issue for 4YKN
Error_ID=[];

% Parallel Computation Failure due to Memory Restriction
for t=1:size(Sample,1) % Example: t=1

    try
        % Extract Coordinates of Alpha Carbon (CA)
        PDB_ID=lower(Sample(t).ID);
        PDB=pdbread(['PDB/', PDB_ID,'.pdb']);
        PDB_Model=[PDB.Model(1).Atom.X; PDB.Model(1).Atom.Y; PDB.Model(1).Atom.Z]';

        PDB_CA_Idx=strcmp({PDB.Model(1).Atom.AtomName}, 'CA')';
        PDB_resSeq_Idx=[PDB.Model(1).Atom.resSeq]';

        Res_ID_Idx=zeros(size(PDB_resSeq_Idx,1),1);
        Sub_Idx=sort(unique(PDB_resSeq_Idx));
        for Sidx=1:size(Sub_Idx,1)
            Res_ID_Idx(PDB_resSeq_Idx==Sub_Idx(Sidx))=Sidx;
        end

        % Load Protein Molecular Surface of Corresponding PDB ID
        Ply_file = fopen(['Surf/',PDB_ID, '.ply']);
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
        Int_File=table2array(readtable(['Reduced_Sample_Int_Res_Results/', PDB_ID,'.txt']));
        Int_Res_Info=Int_File(Int_File(:,5)==1,1:4);

        % Output All Surf Res
        % Union_Surf_Res_Info=union(Surf_Res_Info, Int_Res_Info, 'rows');     
        % Surf_Output=[];
        % Surf_Output.X=Union_Surf_Res_Info(:,1);
        % Surf_Output.Y=Union_Surf_Res_Info(:,2);
        % Surf_Output.Z=Union_Surf_Res_Info(:,3);
        % Surf_Output.Res=Union_Surf_Res_Info(:,4);
        % writetable(struct2table(Surf_Output),['Reduced_Sample_Surf_Res_Results/', PDB_ID,'.txt']);

        % Output n_Int Surf Res
        n_Int_Res_Info=setdiff([Surf_Res_Info; Int_Res_Info], Int_Res_Info, 'rows');

        n_Int_Output=[];
        n_Int_Output.X=n_Int_Res_Info(:,1);
        n_Int_Output.Y=n_Int_Res_Info(:,2);
        n_Int_Output.Z=n_Int_Res_Info(:,3);
        n_Int_Output.Res=n_Int_Res_Info(:,4);
        writetable(struct2table(n_Int_Output),['Reduced_Sample_nInt_Res_Results/', PDB_ID,'.txt']);

    catch
        Error_ID=[Error_ID; t];
    end
    t
end