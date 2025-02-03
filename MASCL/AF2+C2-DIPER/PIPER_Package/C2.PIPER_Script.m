%%%% Create the Directory Pathes
% Load Input_Sample.txt
% Filename='P41212_Reduced_Sample.txt';
Filename='target_entry.txt';
File=fopen(Filename);

% Create the Deposited File
s=1;
Sample_Info(1).ID=[];
while (1)
    line=fgetl(File);
    if line==-1, break, end 
    Sample_Info(s).ID=sscanf(line,'%c');
    mkdir(lower(Sample_Info(s).ID(1:4)));
    s=s+1;         
end
fclose(File);

%
%
%

%%%% Sample Preparation: Sample_prep.txt
Script_Prep = fopen([pwd,'/Sample_prep.txt'],'w');
fprintf(Script_Prep,'#!/bin/bash');
fprintf(Script_Prep,'\n');
for r=1:size(Sample_Info,2)
    fprintf(Script_Prep,'\n');
    fprintf(Script_Prep,['../protein_prep/prepare.py',' ',lower(Sample_Info(r).ID(1:4)),'_trun_AF.pdb']);   
end  
fclose(Script_Prep);

%
%
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
    fprintf(Script_PIPER,'\n');
    fprintf(Script_PIPER,['../piper -vv -c1.0 -k4 --msur_k=1.0 --maskr=1.0 -T FFTW_EXHAUSTIVE -R 3281 -t 1 -p ../prms/atoms.prm -f ../prms/Coeff_DIPER.prm -r ../prms/C2_rots.prm',' ',lower(Sample_Info(r).ID(1:4)),'_trun_AF_pnon.pdb',' ',lower(Sample_Info(r).ID(1:4)),'_trun_AF_pnon.pdb',' ','--o ./C2.DIPER_Result/',lower(Sample_Info(r).ID(1:4)),'/']);
end  
fclose(Script_PIPER);