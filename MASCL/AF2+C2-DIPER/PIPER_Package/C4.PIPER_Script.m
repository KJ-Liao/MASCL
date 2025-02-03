%%%% Load Acceptable Pose
load ('C2.DIPER_Result/C2_TopRank_Result.mat');
                       
%

%%%% Sample Preparation: Sample_prep.txt
Set_Idx='.000.'; 
Script_Prep = fopen([pwd,'/C4_Sample_prep.txt'],'w');
fprintf(Script_Prep,'#!/bin/bash');
fprintf(Script_Prep,'\n');
for r=1:size(C2_TopRank_Result,2)
    for s=1:size(C2_TopRank_Result(r).Output_Idx,1)
        fprintf(Script_Prep,'\n');
        fprintf(Script_Prep,['../protein_prep/prepare.py ',lower(C2_TopRank_Result(r).ID),Set_Idx,num2str(C2_TopRank_Result(r).Output_Idx(s)),'.pdb']);   
    end  
end
fclose(Script_Prep);

%

%%%% Script fot PIPER: Sample_PIPER.txt

% Please correct the path of you pdb for PIPER by yourself!!!!!
% Please correct the path of you pdb for PIPER by yourself!!!!!
% Please correct the path of you pdb for PIPER by yourself!!!!!

% Write Sample_PIPER.txt
Script_PIPER = fopen([pwd,'/C4_Sample_PIPER.txt'],'w');
fprintf(Script_PIPER,'#!/bin/bash');
fprintf(Script_PIPER,'\n');
for r=1:size(C2_TopRank_Result,2)
    for s=1:size(C2_TopRank_Result(r).Output_Idx,1)
        fprintf(Script_PIPER,'\n');
        fprintf(Script_PIPER,['../piper -vv -c1.0 -k4 --msur_k=1.0 --maskr=1.0 -T FFTW_EXHAUSTIVE -R 3395 -t 1 -p ../prms/atoms.prm -f ../prms/Coeff_DIPER.prm -r ../prms/C4_rots.prm',' ',lower(C2_TopRank_Result(r).ID),Set_Idx,num2str(C2_TopRank_Result(r).Output_Idx(s)),'_pnon.pdb',' ',lower(C2_TopRank_Result(r).ID),Set_Idx,num2str(C2_TopRank_Result(r).Output_Idx(s)),'_pnon.pdb',' ','--o ./C4.DIPER_Result/',lower(C2_TopRank_Result(r).ID),Set_Idx,num2str(C2_TopRank_Result(r).Output_Idx(s)),'/']);
    end  
end
fclose(Script_PIPER);