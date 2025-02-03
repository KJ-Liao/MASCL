%%%% Patch Merge
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

% Pre-allocate RAM for Error ID & Patch Struct
fields = {'Coord','Norm_V','Res_Type'};
Patch_Merge = cell2struct(cell(length(fields), []), fields);

% Merge
pm=0; Blank_ID=[];
for si=1:size(Sample_Info,1)
    PDB_ID=lower(Sample_Info(si).ID);
    Patch_File=table2array(readtable(['Reduced_Sample_nInt_Res_Results/Patch_Results_nInt/',PDB_ID,'_Patch.txt']));
    if ~isempty(Patch_File)
        for sj=1:size(Patch_File,1)
            pm=pm+1;
            Patch_Merge(pm).Coord    = reshape(Patch_File(sj,1:18),6,3);
            Patch_Merge(pm).Norm_V   = Patch_File(sj,19:21);
            Patch_Merge(pm).Res_Type = Patch_File(sj,22:27);
        end
        if isnan(Patch_Merge(pm).Norm_V(1))
            Error_ID=[Error_ID; si];
        end
    else
        % si=1272: '1xs7'
        Blank_ID=[Blank_ID; si];
    end
    si
end

% Output
Coord=[Patch_Merge.Coord];
Reshpaed_Coord=reshape(Coord,18,[])';
[~, Idx]=unique(Reshpaed_Coord, 'rows');

Unique_Patch_Merge=Patch_Merge(Idx);
save('Reduced_Sample_nInt_Res_Results/Unique_Patch_Merge_PSI.mat', 'Unique_Patch_Merge')