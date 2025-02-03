%%%% Surf-Orca: Symmetry Correction & Local Docking
% Establish Deposited File of Refined TopRank Poses
mkdir('SurfOrca_Result/P43212_Refined_TopRank_300_Pose');

% Load TopRank Result
load('C4.DIPER_Result/P43212_Sum_Table.mat')
load('C4.DIPER_Result/P43212_TopRank_Result.mat');

% Docking Grid Cloud
Radius=8;
Ind=0; Grid_Cloud=zeros(Radius*2*2+1,3);
for x=-Radius:Radius
    for y=-Radius:Radius
        for z=-Radius:Radius
            if sqrt(x^2+y^2+z^2)<Radius
                Ind=Ind+1;
                Grid_Cloud(Ind,:)=[x,y,z];
            end
        end
    end
end
Grid_Cloud(Ind+1:end,:,:)=[];

%

%%%% MXRGA_Results
% Experimental Data
MXRGA(1).PDB_ID=[];
MXRGA(1).Xtal_Cell_Pmt=[];

% DIPER-Simulated Result
MXRGA(1).SiXt_ID=[];
MXRGA(1).Original_Cell_Pmt=[];
MXRGA(1).Original_RMSD95=[];

% SurfOrca-Refined Result
MXRGA(1).SurfOrca_Cell_Pmt=[];
MXRGA(1).SurfOrca_RMSD95=[];

% Surf-Orca
for r=1:size(Sum_Table,2) % Example: r=7

    % Establish Deposited File
    if exist(['SurfOrca_Result/P43212_Refined_TopRank_300_Pose/',Sum_Table(r).ID(1:4)],'dir')==0
        mkdir(['SurfOrca_Result/P43212_Refined_TopRank_300_Pose/',Sum_Table(r).ID(1:4)]);
    end

    %%%% Record PDB ID & Cell Parameters of Corresponding Crystal
    % PDB_ID
    MXRGA(r).PDB_ID=Sum_Table(r).ID(1:4);

    % Cell Parameter
    File=fopen(['PDB/',Sum_Table(r).ID(1:4),'.pdb']);
    Unit_Cell=[];
    while(1)
        line=fgetl(File);
        if isequal(sscanf(line(1:3),'%c'),'CRY')
            Unit_Cell=sscanf(line(8:55),'%f %f %f');
            break,
        end
    end
    fclose(File);
    MXRGA(r).Xtal_Cell_Pmt=Unit_Cell(1:3);

    %

    %%%% Load TopRank PDB File & Extract Coordinates
    % MXRGA.SiXt_ID
    MXRGA(r).SiXt_ID=Sum_Table(r).ID;

    % Load TopRank PDB File
    Filename=['C4.DIPER_Result/P43212_TopRank_300_Pose/',Sum_Table(r).ID(1:4),'/',Sum_Table(r).ID,'.pdb'];
    AF=pdbread(Filename);

    % Benefit for Following PDB writing
    Output_AF=AF;

    % Extract Coordinates from Output TopRank PDB File
    CA_Idx=strcmp({AF.Model.Atom.AtomName}, 'CA');
    Model=[AF.Model.Atom.X; AF.Model.Atom.Y; AF.Model.Atom.Z]';

    % SiXt/Xtal in Order of [toXYZ; toXYZ_pair; fromXYZ; fromXYZ_pair]
    len=sum(CA_Idx)/4;
    TopRank_DIPER=Model(CA_Idx,:);

    % Divide C4_DIPER into two Monomers: Rcpt, Lgnd
    Rcpt=TopRank_DIPER(1:len*2,:);
    Lgnd=TopRank_DIPER(len*2+1:end,:);

    % Divide Rcpt/Lgnd into Rcpt_A, Rcpt_B, Lgnd_A, Lgnd_B
    toXYZ=Rcpt(1:len,:);
    toXYZ_pair=Rcpt(len+1:end,:);
    fromXYZ=Lgnd(1:len,:);
    fromXYZ_pair=Lgnd(len+1:end,:);

    % Calculate Closest Neighbor Protein Distance (Dist_thres) for Further Pose Selection
    Dist_thres=min(pdist([mean(toXYZ); mean(toXYZ_pair); mean(fromXYZ); mean(fromXYZ_pair)]));

    % Sequence
    AF_Seq=aminolookup([AF.Model.Atom(CA_Idx).resName]);

    %

    %%%% Jugle (Juicy Angle, the Angle for Matched Packing) Screening
    % P 41 21 2 Setting
    Jable=zeros(360,6);
    theo_eul=rotm2eul([1, 0, 0; 0, -1, 0; 0, 0, -1]); ;
    Rz=[cos(pi/4), -sin(pi/4), 0; sin(pi/4), cos(pi/4), 0; 0, 0, 1];

    % Jugle Screening
    for theta=1:360

        % Spin Along with y=x Axis, z=0
        Rot_factor=Rz*[1, 0, 0; 0, cosd(theta), -sind(theta); 0, sind(theta), cosd(theta)]*Rz';
        f_Mo=fromXYZ*Rot_factor;
        t_Mo=toXYZ*Rot_factor;
        t_Mo_pair=toXYZ_pair*Rot_factor;

        % Calculate Origin_oRMSD of f_Mo (O_f_Mo) Relative to t_Mo (O_t_Mo)of DIPER Homodimer
        O_f_Mo=f_Mo-repmat(mean(f_Mo), len, 1);
        O_t_Mo=t_Mo-repmat(mean(t_Mo), len, 1);
        oRMSD=sqrt(sum(sum((O_t_Mo*[1, 0, 0; 0, -1, 0; 0, 0, -1]-O_f_Mo).^2))/len);

        % Find Eular Angle (Rz->Ry->Rx) of f_Mo Relative to t_Mo of DIPERX Homodimer
        f_R=CoordiExam(f_Mo, t_Mo);
        Angle_Det=abs(rotm2eul(f_R)-theo_eul);
        Det_Idx=find(Angle_Det>pi,1);
        Angle_Det(Det_Idx)=2*pi-Angle_Det(Det_Idx);
        MRE=sqrt(mean(Angle_Det.^2))*180/pi;

        % Calculate Cell Parameter
        Cell_pmt_1=abs((mean(f_Mo)-mean(t_Mo))*[2; 0; 0]);
        Cell_pmt_2=(mean(f_Mo)-mean(t_Mo_pair))*[0; 0; 4];
        Cell_pmt=[Cell_pmt_1, Cell_pmt_1, Cell_pmt_2];

        % Jable Record
        Jable(theta,1)=theta;
        Jable(theta,2)=oRMSD;                           % Origin_oRMSD
        Jable(theta,3)=MRE;                             % Mean Rotation Error
        Jable(theta,4:6)=Cell_pmt;                      % Cell Parameter
    end

    % Pose Examination via the Origin_oRMSD (Default<10), Mean Rotation Error (MRE, Default<20) & Dist_thres
    Qualified_Pose=(Jable(:,2)<10)&(Jable(:,3)<20)&(sum(Jable(:,4:6)>Dist_thres,2)==3);

    % Pose Examination via the Lowest Origin_oRMSD
    Qualified_Table=Jable(Qualified_Pose,:);
    Jugle=Qualified_Table(Qualified_Table(:,2)==min(Qualified_Table(:,2)),1);

    % MXRGA_Original_Cell_Pmt
    MXRGA(r).Original_Cell_Pmt=Jable(Jugle,4:6);
    MXRGA(r).Original_RMSD95=Sum_Table(r).Overall_RMSD95;

    %

    if ~isempty(Jugle)
        %%%% Surf_Orca (Grid Dock)
        % Determine Closest Docking Partners
        tf_Dist=mean(mink(min(pdist2(toXYZ, fromXYZ)),round(size(toXYZ,1)*0.1)));
        tf_pair_Dist=mean(mink(min(pdist2(toXYZ_pair, fromXYZ_pair)),round(size(toXYZ_pair,1)*0.1)));

        % Assign the Shell Boundary
        CC_min_lengh=2.8;

        % Assign the Angle Resolution
        Angle_Res=1;
        Jugle_list_Up_B=Jugle+15;
        Jugle_list_Low_B=Jugle-15;

        % Find Refined Jugle
        i=0; Grid_Ocean=[];

        parfor Local_theta=Jugle_list_Low_B:Angle_Res:Jugle_list_Up_B

            % Spin Along with y=x Axis, z=0
            Rot_factor=Rz*[1, 0, 0; 0, cosd(Local_theta), -sind(Local_theta); 0, sind(Local_theta), cosd(Local_theta)]*Rz';

            t_Mo=toXYZ*Rot_factor;
            t_Mo_pair=toXYZ_pair*Rot_factor;
            t_Mo_set=[t_Mo; t_Mo_pair];

            f_Mo=fromXYZ*Rot_factor;
            f_Mo_pair=fromXYZ_pair*Rot_factor;

            f_Mo_Pose=t_Mo*[1, 0, 0; 0, -1, 0; 0, 0, -1];
            f_Mo_pair_Pose=t_Mo_pair*[1, 0, 0; 0, -1, 0; 0, 0, -1];

            if tf_Dist < tf_pair_Dist
                Theo_f_Mo=f_Mo_Pose-mean(f_Mo_Pose)+mean(f_Mo);
                Theo_f_Mo_pair=f_Mo_pair_Pose-mean(f_Mo_Pose)+mean(f_Mo);
            else
                Theo_f_Mo=f_Mo_Pose-mean(f_Mo_pair_Pose)+mean(f_Mo_pair);
                Theo_f_Mo_pair=f_Mo_pair_Pose-mean(f_Mo_pair_Pose)+mean(f_Mo_pair);
            end

            % Dock Grid Cloud Screening
            for m=1:size(Grid_Cloud,1)

                % f_Mo_set
                R_f_Mo=Theo_f_Mo+Grid_Cloud(m,1:3);
                R_f_Mo_pair=Theo_f_Mo_pair+Grid_Cloud(m,1:3);
                R_f_Mo_set=[R_f_Mo; R_f_Mo_pair];

                % Calculate Cell Parameter
                Cell_pmt_1=abs((mean(R_f_Mo)-mean(t_Mo))*[2; 0; 0]);
                Cell_pmt_2=(mean(R_f_Mo)-mean(t_Mo_pair))*[0; 0; 4];
                Cell_pmt=[Cell_pmt_1, Cell_pmt_1, Cell_pmt_2];

                % Cell Parameters & Local/Global RMSD Cutoff
                if all(Cell_pmt>Dist_thres)
                    ASU=[t_Mo_set; R_f_Mo_set];
                    [~,~,RMSD_Global]=CoordiExam([t_Mo_set; f_Mo; f_Mo_pair], ASU);

                    if RMSD_Global < 5
                        if tf_Dist < tf_pair_Dist
                            [~,~,RMSD_Local]=CoordiExam([t_Mo; f_Mo], [t_Mo; R_f_Mo]);
                        else
                            [~,~,RMSD_Local]=CoordiExam([t_Mo_pair; f_Mo_pair], [t_Mo_pair; R_f_Mo_pair]);
                        end

                        % Interaction & Collision Properties of the Dimer-Dimer (within 12 A)
                        if RMSD_Local < 5
                            Up_B=max(t_Mo_set)+12;
                            Low_B=min(t_Mo_set)-12;
                            X=R_f_Mo_set(:,1)>Low_B(1)&R_f_Mo_set(:,1)<Up_B(1);
                            Y=R_f_Mo_set(:,2)>Low_B(2)&R_f_Mo_set(:,2)<Up_B(2);
                            Z=R_f_Mo_set(:,3)>Low_B(3)&R_f_Mo_set(:,3)<Up_B(3);
                            Inbox_f_Mo_set=R_f_Mo_set(X&Y&Z,:);

                            Up_B=max(R_f_Mo_set)+12;
                            Low_B=min(R_f_Mo_set)-12;
                            X=t_Mo_set(:,1)>Low_B(1)&t_Mo_set(:,1)<Up_B(1);
                            Y=t_Mo_set(:,2)>Low_B(2)&t_Mo_set(:,2)<Up_B(2);
                            Z=t_Mo_set(:,3)>Low_B(3)&t_Mo_set(:,3)<Up_B(3);
                            Inbox_t_Mo_set=t_Mo_set(X&Y&Z,:);

                            DISM=pdist2(Inbox_t_Mo_set, Inbox_f_Mo_set);

                            % Exclude No-Interaction Dimer-Dimer (within 12 A)
                            if min(DISM,[],'all') < 12
                                Collided_Res_no=sum(DISM<=CC_min_lengh,'all');

                                % Collision Properties of the (Dimer-Dimer)-(Dimer-Dimer)
                                if Collided_Res_no==0
                                    Pack=[t_Mo_set; R_f_Mo_set; R_f_Mo_set*[0, 1, 0; 1, 0, 0; 0, 0, -1]];
                                    
                                    Pack_Up_B=max(Pack)+CC_min_lengh;
                                    Pack_Low_B=min(Pack)-CC_min_lengh;

                                    Nerb_1=[Pack(:,1), Pack(:,2), Pack(:,3)+Cell_pmt(3)];
                                    X=Nerb_1(:,1)>Pack_Low_B(1)&Nerb_1(:,1)<Pack_Up_B(1);
                                    Y=Nerb_1(:,2)>Pack_Low_B(2)&Nerb_1(:,2)<Pack_Up_B(2);
                                    Z=Nerb_1(:,3)>Pack_Low_B(3)&Nerb_1(:,3)<Pack_Up_B(3);
                                    Inbox_Nerb1=Nerb_1(X&Y&Z,:);

                                    Up_B=max(Nerb_1)+CC_min_lengh;
                                    Low_B=min(Nerb_1)-CC_min_lengh;
                                    X=Pack(:,1)>Low_B(1)&Pack(:,1)<Up_B(1);
                                    Y=Pack(:,2)>Low_B(2)&Pack(:,2)<Up_B(2);
                                    Z=Pack(:,3)>Low_B(3)&Pack(:,3)<Up_B(3);
                                    Inbox_Pack=Pack(X&Y&Z,:);

                                    if ~any(pdist2(Inbox_Pack, Inbox_Nerb1)<=CC_min_lengh,'all')
                                        Nerb_2=[Pack(:,1), Pack(:,2)+Cell_pmt(2), Pack(:,3)];
                                        X=Nerb_2(:,1)>Pack_Low_B(1)&Nerb_2(:,1)<Pack_Up_B(1);
                                        Y=Nerb_2(:,2)>Pack_Low_B(2)&Nerb_2(:,2)<Pack_Up_B(2);
                                        Z=Nerb_2(:,3)>Pack_Low_B(3)&Nerb_2(:,3)<Pack_Up_B(3);
                                        Inbox_Nerb2=Nerb_2(X&Y&Z,:);

                                        Up_B=max(Nerb_2)+CC_min_lengh;
                                        Low_B=min(Nerb_2)-CC_min_lengh;
                                        X=Pack(:,1)>Low_B(1)&Pack(:,1)<Up_B(1);
                                        Y=Pack(:,2)>Low_B(2)&Pack(:,2)<Up_B(2);
                                        Z=Pack(:,3)>Low_B(3)&Pack(:,3)<Up_B(3);
                                        Inbox_Pack=Pack(X&Y&Z,:);

                                        if ~any(pdist2(Inbox_Pack, Inbox_Nerb2)<=CC_min_lengh,'all')
                                            Nerb_3=[Pack(:,1)+Cell_pmt(1), Pack(:,2), Pack(:,3)];
                                            X=Nerb_3(:,1)>Pack_Low_B(1)&Nerb_3(:,1)<Pack_Up_B(1);
                                            Y=Nerb_3(:,2)>Pack_Low_B(2)&Nerb_3(:,2)<Pack_Up_B(2);
                                            Z=Nerb_3(:,3)>Pack_Low_B(3)&Nerb_3(:,3)<Pack_Up_B(3);
                                            Inbox_Nerb3=Nerb_3(X&Y&Z,:);

                                            Up_B=max(Nerb_3)+CC_min_lengh;
                                            Low_B=min(Nerb_3)-CC_min_lengh;
                                            X=Pack(:,1)>Low_B(1)&Pack(:,1)<Up_B(1);
                                            Y=Pack(:,2)>Low_B(2)&Pack(:,2)<Up_B(2);
                                            Z=Pack(:,3)>Low_B(3)&Pack(:,3)<Up_B(3);
                                            Inbox_Pack=Pack(X&Y&Z,:);

                                            if ~any(pdist2(Inbox_Pack, Inbox_Nerb3)<=CC_min_lengh,'all')
                                                Grid_Ocean=[Grid_Ocean; [Local_theta, m, Grid_Cloud(m,1:3), Cell_pmt, RMSD_Local, RMSD_Global]];
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        %
        %
        %

        %%%% Scoring Complexes of Docking Results
        if ~isempty(Grid_Ocean)

            % Selected via Popular Angles & Smallest Global RMSD
            [~, pos]=min(Grid_Ocean(:,10));

            % Refined Final Pose
            Opt_Pose = Grid_Ocean(pos,:);

            %

            %%%% Output (Write the PDB Files of Opt_Pose)
            % Spin Along with y=x Axis, z=0
            Refined_theta=Opt_Pose(1,1);
            Rot_factor=Rz*[1, 0, 0; 0, cosd(Refined_theta), -sind(Refined_theta); 0, sind(Refined_theta), cosd(Refined_theta)]*Rz';

            t_Mo=toXYZ*Rot_factor;
            t_Mo_pair=toXYZ_pair*Rot_factor;

            f_Mo=fromXYZ*Rot_factor;
            f_Mo_pair=fromXYZ_pair*Rot_factor;

            f_Mo_Pose=t_Mo*[1, 0, 0; 0, -1, 0; 0, 0, -1];
            f_Mo_pair_Pose=t_Mo_pair*[1, 0, 0; 0, -1, 0; 0, 0, -1];

            % Refined/Normalized Theo_f_Mo, Theo_f_Mo_pair
            if tf_Dist < tf_pair_Dist
                Theo_f_Mo=f_Mo_Pose-mean(f_Mo_Pose)+mean(f_Mo);
                Theo_f_Mo_pair=f_Mo_pair_Pose-mean(f_Mo_Pose)+mean(f_Mo);
            else
                Theo_f_Mo=f_Mo_Pose-mean(f_Mo_pair_Pose)+mean(f_Mo_pair);
                Theo_f_Mo_pair=f_Mo_pair_Pose-mean(f_Mo_pair_Pose)+mean(f_Mo_pair);
            end

            R_f_Mo=Theo_f_Mo+Opt_Pose(3:5);
            R_f_Mo_pair=Theo_f_Mo_pair+Opt_Pose(3:5);

            % CA Coordinate of Opt_Results ([toXYZ; toXYZ_pair; fromXYZ; fromXYZ_pair])
            Opt_Results_Coord=[t_Mo; t_Mo_pair; R_f_Mo; R_f_Mo_pair];

            % Construct Coordinate of Opt_Results (Output_Results_Coord)
            [R1,T1]=CoordiExam(toXYZ, t_Mo);
            [R2,T2]=CoordiExam(toXYZ, t_Mo_pair);
            [R3,T3]=CoordiExam(toXYZ, R_f_Mo);
            [R4,T4]=CoordiExam(toXYZ, R_f_Mo_pair);

            Opt_Results_Model=Model(1:size(Model,1)/4,:);
            Output_Results_Coord=[Opt_Results_Model*R1+T1; Opt_Results_Model*R2+T2; Opt_Results_Model*R3+T3; Opt_Results_Model*R4+T4];

            X=num2cell(Output_Results_Coord(:,1)); [Output_AF.Model.Atom.X]=X{:};
            Y=num2cell(Output_Results_Coord(:,2)); [Output_AF.Model.Atom.Y]=Y{:};
            Z=num2cell(Output_Results_Coord(:,3)); [Output_AF.Model.Atom.Z]=Z{:};

            warning('off','all');
            pdbwrite(['SurfOrca_Result/P43212_Refined_TopRank_300_Pose/', MXRGA(r).SiXt_ID(1:4), '/', MXRGA(r).SiXt_ID, '_Refined.pdb'], Output_AF);

            %

            %%%% RMSD 95 Calculation
            % Load Corresponding Ans PDB (Xtal_PDB) for RMSD Calculation
            Filename=['C4.DIPER_Result/P43212_TopRank_300_Pose/', MXRGA(r).SiXt_ID(1:4), '/', MXRGA(r).SiXt_ID, '_Ans.pdb'];
            PDB=pdbread(Filename);

            % Extract Coordinates of Alpha Carbon (CA)
            PDB_Model=[PDB.Model.Atom.X; PDB.Model.Atom.Y; PDB.Model.Atom.Z]';
            CA_Idx=strcmp({PDB.Model.Atom.AtomName}, 'CA')';
            Xtal_PDB=PDB_Model(CA_Idx,:);

            % Sequence
            PDB_Seq=aminolookup([PDB.Model.Atom(CA_Idx).resName]);

            % Extract 1-letter Sequences
            seq1=AF_Seq; seq2=PDB_Seq;

            % Global Sequence Alignment to Extract Superimposed CA for RMSD Calculation
            [~, Alignment] = nwalign(seq1,seq2);
            alignidex=find((Alignment(1,:)~='-')&(Alignment(3,:)~='-'));
            aseq1=find(Alignment(1,:)=='-');
            aseq2=find(Alignment(3,:)=='-');

            aseq1_idx=[]; aseq2_idx=[];
            for a=1:length(alignidex)
                aseq1_idx=[aseq1_idx, alignidex(a)-sum(aseq1<alignidex(a))];
                aseq2_idx=[aseq2_idx, alignidex(a)-sum(aseq2<alignidex(a))];
            end

            % RMSD 95 Calculation
            Opt_fromXYZ=Opt_Results_Coord(aseq1_idx,:);
            Opt_toXYZ=Xtal_PDB(aseq2_idx,:);
            [~,~,~,Sum_Square]=CoordiExam(Opt_fromXYZ, Opt_toXYZ);
            RMSD_95=sqrt(sum(Sum_Square(1:ceil(size(Opt_toXYZ,1)*0.95)))/ceil(size(Opt_toXYZ,1)*0.95));

            MXRGA(r).SurfOrca_RMSD95=RMSD_95;
            MXRGA(r).SurfOrca_Cell_Pmt=Opt_Pose(6:8);
        end
    end
    r
end

save('SurfOrca_Result/P43212_MXRGA_Result.mat','MXRGA');
% writetable(struct2table(MXRGA),'MXRGA_Result.xlsx');