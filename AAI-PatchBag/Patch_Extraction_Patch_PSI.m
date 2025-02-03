%%%% Default Setting
% Parallel Computation
parpool(6)

% Establish Deposited Folder
mkdir('Reduced_Sample_nInt_Res_Results/Patch_Results_nInt/');

% Patch Size (No of Alpha Carbons)
Patch_Size=6;

% Pre-allocate RAM for Error ID & Patch Struct
Error_ID=[]; fields = {'Coord','Norm_V','Res_Type'};

% Load PDB ID of Qualified Reduced  Samples
% Filename='Qualified_Reduced_Sample_3737(Exclude 4YKN).txt';
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

% Patch Extraction
for i=1:size(Sample_Info,1) % Example: i=1
    try
        % Load Info of Interacting Residues
        PDB_ID=lower(Sample_Info(i).ID);
        nInt_File=table2array(readtable(['Reduced_Sample_nInt_Res_Results/', PDB_ID,'.txt']));
        nInt_Asy_Unit=nInt_File(:,1:4);

        if size(nInt_Asy_Unit,1)<6
            Patch = cell2struct(cell(length(fields), []), fields);
            writetable(struct2table(Patch),['Reduced_Sample_nInt_Res_Results/Patch_Results_nInt/',PDB_ID, '_Patch.txt']);
        else
            % Read Molecule Surface
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

            % Log Postion of Alpha Carbon and Corresponding Index
            PDB=pdbread(['PDB/', PDB_ID,'.pdb']);
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
                        Error_ID=[Error_ID; i];
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
            writetable(struct2table(Patch),['Reduced_Sample_nInt_Res_Results/Patch_Results_nInt/',PDB_ID, '_Patch.txt']);
        end
    catch
        Error_ID=[Error_ID; i];
    end
    i
end