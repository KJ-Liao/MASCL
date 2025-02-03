%%%% PDBID2AF_Structure
% Build up List of PDB Sample ID
Filename='target_entry.txt';
File=fopen(Filename);
Sample_Info(1).ID=[];i=1;
while (1)
    line=fgetl(File);
    if line==-1, break, end 
    Sample_Info(i).ID=sscanf(line,'%c');
    i=i+1;         
end
fclose(File);

% Extract the Uniprot ID and Start/End Site of Seq
o=1;
Error(1).ID=[];
Sample_Info(1).UNP=[];
Sample_Info(1).Res_Sta=[];
Sample_Info(1).Res_End=[];
for j=1:size(Sample_Info,2)

    % Read Corresponding PDB
    PDB_Info=pdbread(['PDB/',lower(Sample_Info(j).ID),'.pdb']);

    % Extract the Uniprot ID and Start/End Site of Seq
    if isfield(PDB_Info,'DBReferences') & ~isempty(PDB_Info.DBReferences)
        DBREF_UNP=find(strcmp(deblank({PDB_Info.DBReferences.database}),'UNP'));
        if ~isempty(DBREF_UNP) & (size(DBREF_UNP,2)==1)
            Sample_Info(j).UNP=deblank(PDB_Info.DBReferences(DBREF_UNP).dbAccession);
            Sample_Info(j).Res_Sta=PDB_Info.DBReferences(DBREF_UNP).dbseqBegin;
            Sample_Info(j).Res_End=PDB_Info.DBReferences(DBREF_UNP).dbseqEnd;
        else
            Error(o).ID=upper(Sample_Info(j).ID);
            Sample_Info(j).UNP=upper(Sample_Info(j).ID);
            Sample_Info(j).Res_Sta=1;
            Sample_Info(j).Res_End=PDB_Info.Sequence.NumOfResidues;        
            o=o+1;
        end
    end
    Reamined_No=size(Sample_Info,2)-j;
end

%
%
%

%%%% Download Alphafold2_V3 Structure
% Download Structure (Version 3)
for p=1:size(Sample_Info,2) 
    URL=strcat('https://alphafold.ebi.ac.uk/files/AF-', Sample_Info(p).UNP,'-F1-model_v3.pdb');
    websave([lower(Sample_Info(p).ID), '_AF'],URL);
    Reamined_No=size(Sample_Info,2)-p;
end

%
%
%

%%%% Remove pLDDT<70 Residues
% Residue Index with pLDDT>70
Sample_Info(1).Sta=[];
Sample_Info(1).End=[];

% Calculate pLDDT>70 Residue Index
for q=1:size(Sample_Info,2)

    % Read Corresponding Alphafold2_V3 PDB
    AF_PDB=pdbread([lower(Sample_Info(q).ID),'_AF.pdb']);

    % Record pLDDT Info
    Data=[AF_PDB.Model.Atom.resSeq; AF_PDB.Model.Atom.tempFactor]';
    [~,idx,~]=unique([AF_PDB.Model.Atom.resSeq],'stable');
    idx(length(idx)+1)=length([AF_PDB.Model.Atom.resSeq])+1;
    pLDDT=Data(idx(1:end-1),:);
    
    % Extract the pLDDT>70 Residues (Remove Uncertain Structure from Sta/End)
    % Start
    i=0;
    while(1)
        if min([pLDDT(Sample_Info(q).Res_Sta+i,2), pLDDT(Sample_Info(q).Res_Sta+1+i,2), pLDDT(Sample_Info(q).Res_Sta+2+i,2), pLDDT(Sample_Info(q).Res_Sta+3+i,2), pLDDT(Sample_Info(q).Res_Sta+4+i,2)])>70, break, end
        i=i+1;
    end
    
    % End
    j=0;
    while(1)
        if min([pLDDT(Sample_Info(q).Res_End-j,2), pLDDT(Sample_Info(q).Res_End-1-j,2), pLDDT(Sample_Info(q).Res_End-2-j,2), pLDDT(Sample_Info(q).Res_End-3-j,2), pLDDT(Sample_Info(q).Res_End-4-j,2)])>70, break, end
        j=j+1;
    end

    Sample_Info(q).Sta=Sample_Info(q).Res_Sta+i;
    Sample_Info(q).End=Sample_Info(q).Res_End-j;  
    Seq=[Sample_Info(q).Sta:Sample_Info(q).End];

    % Extract the pLDDT>70 Residues (Remove Uncertain Structure of Intra-loop)
    for p=Sample_Info(q).Sta+2:Sample_Info(q).End-2
        if mean(pLDDT(p-2:p+2,2))<70,
            Seq=setdiff(Seq,p,'stable');
        end
    end
    
    % Output Seq_Index with pLDDT>70
    Seq_Index=[];
    for l=1:length(Seq)
        Seq_Index=[Seq_Index; find(Data(:,1)==Seq(l))];
    end
    
    %%%% Generate Truncated (pLDDT>70) Alphafold2_V3 PDB (_trun_AF,pdb)
    % Loading Simplified Template PDB
    PDBStr=pdbread('Package/Simplified temp.pdb');

    % Replace Default Template with Model Coordinates
    PDBStr.Model.Atom=AF_PDB.Model.Atom(Seq_Index);
    AtomSerNo=num2cell(1:length(Seq_Index));
    [PDBStr.Model.Atom.AtomSerNo]=AtomSerNo{:};

    % Output PDB file
    warning('off','all')
    pdbwrite([lower(Sample_Info(q).ID),'_trun_AF.pdb'], PDBStr);
    
    Reamined_No=size(Sample_Info,2)-q;
end