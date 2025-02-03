%%%% Load PDB ID
File=fopen('target_entry.txt');
sp=1;Sample_Info(1).ID=[];
while (1)
    line=fgetl(File);
    if line==-1, break, end
    Sample_Info(sp).ID=sscanf(line,'%c');
    sp=sp+1;
end
fclose(File);
                
%

%%%% Sample Preparation: Sample_prep.txt
load('AF3_Dimer_Result/AF3_Dimer_Result.mat')

Script_Prep = fopen([pwd,'/Sample_prep.txt'],'w');
fprintf(Script_Prep,'#!/bin/bash');
fprintf(Script_Prep,'\n');
for r=1:size(Sample_Info,2)
    for s=1:5
        if Table(r,s) < 6
            fprintf(Script_Prep,'\n');
            fprintf(Script_Prep,['./protein_prep/prepare.py ', pwd,'/AF3_trun_Dimer/',lower(Sample_Info(r).ID),'/',lower(Sample_Info(r).ID),'_trun_AF_model_',num2str(s-1),'.pdb']);   
        end 
    end
end
fclose(Script_Prep);

%

%%%% Script fot PIPER: Sample_PIPER.txt

% Please correct the path of you pdb for PIPER by yourself!!!!!
% Please correct the path of you pdb for PIPER by yourself!!!!!
% Please correct the path of you pdb for PIPER by yourself!!!!!

% Write Sample_PIPER.txt
Script_PIPER = fopen([pwd,'/Sample_PIPER.txt'],'w');
fprintf(Script_PIPER,'#!/bin/bash');
fprintf(Script_PIPER,'\n');
for r=1:size(Sample_Info,2)
    for s=1:5
        if Table(r,s) < 6
            fprintf(Script_PIPER,'\n');
            fprintf(Script_PIPER,['./piper -vv -c1.0 -k4 --msur_k=1.0 --maskr=1.0 -T FFTW_EXHAUSTIVE -R 3395 -t 1 -p ./prms/atoms.prm -f ./prms/Coeff_DIPER.prm -r ./prms/C4_rots.prm',' ./','AF_Dimer4Docking/',lower(Sample_Info(r).ID),'_trun_AF_model_',num2str(s-1),'_pnon.pdb',' ./','AF_Dimer4Docking/',lower(Sample_Info(r).ID),'_trun_AF_model_',num2str(s-1),'_pnon.pdb',' ','--o ./AF_Dimer4Docking/',lower(Sample_Info(r).ID),'_trun_AF_model_',num2str(s-1)]);
        end  
    end
end
fclose(Script_PIPER);