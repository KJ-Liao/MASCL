%%%% DLoad PDB ID
File=fopen('target_entry.txt');
sp=1;Sample_Info(1).ID=[];
while (1)
    line=fgetl(File);
    if line==-1, break, end
    Sample_Info(sp).ID=sscanf(line,'%c');
    sp=sp+1;
end
fclose(File);

%%%% Remove pLDDT<50 Residues
% Residue Index with pLDDT>50
Sample_Info(1).Res_Sta=[];
Sample_Info(1).Res_End=[];
Sample_Info(1).Sta=[];
Sample_Info(1).End=[];


% Calculate pLDDT>50 Residue Index
for q=1:size(Sample_Info,2)

    % Establish Deposited File
    mkdir(['AF3_trun_Dimer/', lower(Sample_Info(q).ID)]);

    for n=1:5
        % Read Corresponding Alphafold2_V3 PDB
        AF_PDB=pdbread(['fold_', lower(Sample_Info(q).ID), '/model_', num2str(n-1), '.pdb']);

        % Record pLDDT Info
        Data=[AF_PDB.Model.Atom.resSeq; AF_PDB.Model.Atom.tempFactor]';
        [~,idx,~]=unique([AF_PDB.Model.Atom.resSeq],'stable');
        idx(length(idx)+1)=length([AF_PDB.Model.Atom.resSeq])+1;
        pLDDT=Data(idx(1:end-1),:);

        % Extract the pLDDT>50 Residues (Remove Uncertain Structure from Sta/End)
        Sample_Info(q).Res_Sta=1;
        Sample_Info(q).Res_End=size(pLDDT,1);

        % Start
        i=0;
        while(1)
            if min([pLDDT(Sample_Info(q).Res_Sta+i,2), pLDDT(Sample_Info(q).Res_Sta+1+i,2), pLDDT(Sample_Info(q).Res_Sta+2+i,2), pLDDT(Sample_Info(q).Res_Sta+3+i,2), pLDDT(Sample_Info(q).Res_Sta+4+i,2)])>50, break, end
            i=i+1;
        end

        % End
        j=0;
        while(1)
            if min([pLDDT(Sample_Info(q).Res_End-j,2), pLDDT(Sample_Info(q).Res_End-1-j,2), pLDDT(Sample_Info(q).Res_End-2-j,2), pLDDT(Sample_Info(q).Res_End-3-j,2), pLDDT(Sample_Info(q).Res_End-4-j,2)])>50, break, end
            j=j+1;
        end

        Sample_Info(q).Sta=Sample_Info(q).Res_Sta+i;
        Sample_Info(q).End=Sample_Info(q).Res_End-j;
        Seq=[Sample_Info(q).Sta:Sample_Info(q).End];

        % Extract the pLDDT>50 Residues (Remove Uncertain Structure of Intra-loop)
        for p=Sample_Info(q).Sta+2:Sample_Info(q).End-2
            if mean(pLDDT(p-2:p+2,2))<50,
                Seq=setdiff(Seq,p,'stable');
            end
        end

        % Output Seq_Index with pLDDT>50
        Seq_Index=[];
        for l=1:length(Seq)
            Seq_Index=[Seq_Index; find(Data(:,1)==Seq(l))];
        end
        Seq_Index=sort(Seq_Index);

        %%%% Generate Truncated (pLDDT>50) Alphafold2_V3 PDB (_trun_AF,pdb)
        % Loading Simplified Template PDB
        PDBStr=pdbread('Package/Simplified temp.pdb');

        % Replace Default Template with Model Coordinates
        PDBStr.Model.Atom=AF_PDB.Model.Atom(Seq_Index);
        AtomSerNo=num2cell(1:length(Seq_Index));
        [PDBStr.Model.Atom.AtomSerNo]=AtomSerNo{:};

        % Output PDB file
        warning('off','all')
        pdbwrite(['AF3_trun_Dimer/', lower(Sample_Info(q).ID), '/', lower(Sample_Info(q).ID), '_trun_AF_model_', num2str(n-1), '.pdb'], PDBStr);
    end
end