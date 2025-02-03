% Loading PDB ID of Accessible_SCH_Sample (16123)
Sample_Info=table2struct(readtable('SCHP_Xtal_Cond_Parser_(XCP)_Dataset.xlsx'));

% Load Condition Dictionary
M=readtable('Package/Extraction_Dictionary.csv');
E=readtable('Package/Exception.csv');

% Build Condition MetaTable for Record
Condition_MetaTable=zeros(size(M,1)+1, size(Sample_Info,1));

% Record Other Reagents & Missing Info
i=1; Other_Reagent(1).Reagent=[];
j=1; Other_Reagent(1).Protein=[];
k=1; Missing_Info(1).Content=[];

% Extract Crystallization Conditions
for si=1:size(Sample_Info,1)

    %%%% Record Crystallization Condition
    File=fopen([lower(Sample_Info(si).ID),'.pdb']);
    try
        while(1)
            line=fgetl(File);
            if isequal(sscanf(line(1:37),'%c'),'REMARK 280 CRYSTALLIZATION CONDITIONS')
                Condition=line(39:end);
                line=fgetl(File);
                while isequal(sscanf(line(1:10),'%c'),'REMARK 280')
                    Condition=strcat(Condition,line(12:end));
                    line=fgetl(File);
                    if ~isequal(sscanf(line(1:10),'%c'),'REMARK 280'), break, end
                end
            end
            if isequal(sscanf(line(1:10),'%c'),'REMARK 290')||isequal(sscanf(line(1:3),'%c'),'END'), break, end
        end
    catch
        % Record ID with Problematic PDB Format
        Condition_MetaTable(1,si)=1;
    end
    fclose(File);

    %

    %%%% Condition Formatting
    if Condition_MetaTable(1,si)==0

        %%% Format Preprocessing
        % Remove Redundant Blank & '.'/','/';' in Both End
        Condition=lower(strtrim(Condition));

        % Empty (Null) Examination
        if ~isempty(Condition)

            if  any(ismember({':', ';', '.', ','}, Condition(1)))
                Condition=strtrim(Condition(2:end));
            end

            if  any(ismember({':', ';', '.', ','}, Condition(end)))
                Condition=strtrim(Condition(1:end-1));
            end

            % Null Examination
            if strcmp(Condition, 'null')

                % Null Condition
                Condition_MetaTable(2,si)=1;
            else

                %%% Condition Formatting: Label Anchor Point
                % Label Anchor Point
                Condition=strtrim(Condition);
                Condition=[', ',Condition,','];
                Condition=strrep(Condition,'  ',' ');

                %

                %%% Condition Formatting: Erase Unncessary Words
                % Remove Common Words
                Common_Words={' was ', ' were ', ' the ', ' grown ', ' grew ', ' grows ', ' grow ', 'fresh', 'obtained', 'obtains', 'obtain', 'harvested', 'harvests', 'harvest', 'crystallized', 'crystals',  'crystal', 'crystallization', 'reservoirs', 'reservoir', 'buffered', 'buffers', 'buffer'}; 
                for CW=1:size(Common_Words,2)
                    Condition=strrep(Condition, Common_Words{CW}, ' ');
                    Condition=strrep(Condition,'  ',' ');
                    Condition=strrep(Condition,'  ',' ');
                end
                Condition=strrep(Condition,'  ',' ');

                % Remove Preposition
                Preposition_List={' at ', ' by ', ' in ', ' on ', ' of ', ' for ', ' with ', ' from ', 'using', 'containing'};
                for PL=1:size(Preposition_List,2)
                    Condition=strrep(Condition, Preposition_List{PL}, ', ');
                    Condition=strrep(Condition,'  ',' ');
                end

                % Concentration Formatting (%)
                Condition=strrep(Condition,'mg/ml',' mg/ml ');  % 'mg/ml' to ' mg/ml '
                Condition=strrep(Condition,'mg/l',' mg/l ');    % 'mg/l'  to ' mg/l '
                if ~contains(Condition,'mg/ml')
                    Condition=strrep(Condition,'g/ml',' g/ml ');    % 'm/ml'  to ' g/ml '
                end

                Condition=strrep(Condition,'  ',' ');           % '  ' to ' '
                Condition=strrep(Condition,'%.','%');           % '%.' to '%'
                Condition=strrep(Condition,'%;','%');           % '%;' to '%'
                Condition=strrep(Condition,'%,','%');           % '%,' to '%'
                Condition=strrep(Condition,'%-','-');           % '%-' to '-'

                Condition=strrep(Condition,'/ ','/');           % '/ ' to '/'
                Condition=strrep(Condition,' /','/');           % ' /' to '/'
                Condition=strrep(Condition,'( ','(');           % '( ' to '('
                Condition=strrep(Condition,' )',')');           % ' )' to ')'

                Conc_List={'weight', 'volume', 'wt.', 'wt', 'w', 'vol.', 'vol', 'v', 'mt.', 'mt', 'm'};
                for CL_1=1:size(Conc_List,2)
                    for CL_2=1:size(Conc_List,2)
                        if ~(CL_1==size(Conc_List,2)&&CL_2==size(Conc_List,2))
                            Condition=erase(Condition, ['(', Conc_List{CL_1}, '/', Conc_List{CL_2}, ')']);
                            Condition=erase(Condition, ['(', Conc_List{CL_1}, ' ', Conc_List{CL_2}, ')']);
                            Condition=erase(Condition, [Conc_List{CL_1}, '/', Conc_List{CL_2}]);
                            Condition=strrep(Condition,'  ',' ');
                        end
                    end
                end

                Condition=strrep(Condition, ' saturated', ', saturated% ');      % ' saturated ' to ' saturated%'
                Condition=strrep(Condition,'%','% ');                            % '%,' to '%'
                Condition=strrep(Condition,'  ',' ');                            % '  ' to ' '

                %

                %%% Condition Formatting: Frequent Typo/ Number Format/Protection/Others
                % Frequent Typo
                Condition=strrep(Condition,'o.','0.');
                Condition=strrep(Condition,'.o','.0');
                Condition=strrep(Condition,'?',', ');
                Condition=strrep(Condition,'  ',' ');

                % Number Format
                Condition=strrep(Condition,'1,500','1500');                 % PEG 1,500     
                Condition=strrep(Condition,'3,350','3350');                 % PEG 3,350
                Condition=strrep(Condition,',000','000');
                Condition=strrep(Condition,'  ',' ');

                % Protection
                for p=1:size(E,1)
                    Condition=strrep(Condition, E.Origin{p}, E.Protect{p});
                end

                % Others
                Condition=strrep(Condition,'lipidic cubic phase',', lipidic cubic phase');
                Condition=strrep(Condition,'vapor diffusion',', vapor diffusion');
                Condition=strrep(Condition,'sitting',', sitting');
                Condition=strrep(Condition,'hanging',', hanging');
                

                %

                %%% Condition Formatting: Format Separating/(Abbreviated pH) Layout
                % Format Separating/Space Layout
                Condition=strrep(Condition,':',' ');            % ':' to ' '
                Condition=strrep(Condition,'=',' ');            % '=' to ' '
                Condition=strrep(Condition,'  ',' ');           % '  ' to ' '

                Condition=strrep(Condition,'(+/-)','($/-)');
                Condition=strrep(Condition,'nadp+','nadp$');
                Condition=strrep(Condition,'nad+','nad$');

                Condition=strrep(Condition,'+',', ');           % '+' to ','
                Condition=strrep(Condition,'&',', ');           % '&' to ','
                Condition=strrep(Condition,' and ',', ');       % ' and ' to ','

                Condition=strrep(Condition,'$','+');

                Condition=strrep(Condition,' to ', ' - ');      % ' to ' to '-'
                Condition=strrep(Condition,'--', ', ');         % '--' to ', '
                Condition=strrep(Condition,'~','-');            % '~'  to '-'
                Condition=strrep(Condition,' -','-');           % ' -' to '-'
                Condition=strrep(Condition,'- ','-');           % '- ' to '-'
                Condition=strrep(Condition,'  ',' ');           % '  ' to ' '

                Condition=strrep(Condition,';',', ');           % ';'   to ', '
                Condition=strrep(Condition,',.',', ');          % ',.'  to ', '
                Condition=strrep(Condition,'. ',', ');          % '. '  to ', '
                Condition=strrep(Condition,'.,',', ');          % '.,'  to ', '
                Condition=strrep(Condition,';,',', ');          % ';,'  to ', '
                Condition=strrep(Condition,',,',', ');          % ',,'  to ', '
                Condition=strrep(Condition,'. ,',', ');         % '. ,' to ', '
                Condition=strrep(Condition,'; ,',', ');         % '; ,' to ', '
                Condition=strrep(Condition,', ,',', ');         % ', ,' to ', '
                Condition=strrep(Condition,'  ',' ');           % '  ' to ' '

                Condition=strrep(Condition,' , ',', ');         % ' , ' to ', '
                Condition=strrep(Condition,' ,',', ');          % ' ,' to ', '
                Condition=strrep(Condition,'  ',' ');           % '  ' to ' '

                % (Abbreviated pH)
                Num_List={'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'};
                for NL_1=1:13
                    Condition=strrep(Condition, ['(', Num_List{NL_1}, ')'], [', ph ', Num_List{NL_1}, ',']);     % e.g. (7) to ph 7

                    for NL_2=1:13
                        pH_1=[Num_List{NL_1}, '.', Num_List{NL_2}];
                        Condition=strrep(Condition, ['(', pH_1, ')'], [', ph ', pH_1, ',']);                     % e.g. (7.0) to ph 7.0

                        for NL_3=1:13
                            for NL_4=1:10
                                pH_2=[Num_List{NL_3}, '.', Num_List{NL_4}];
                                Condition=strrep(Condition, ['(', pH_1, '-', pH_2, ')'], [', ph ', pH_1, '-', pH_2, ',']);     % e.g. (6.8-7.5) to ph 6.8-7.5
                            end
                        end
                    end
                end
                Condition=strrep(Condition,'  ',' ');

                Condition=strrep(Condition,',',', ');           % ',' to ', '
                Condition=strrep(Condition,'  ',' ');           % '  ' to ' '
                Condition=strrep(Condition,',,',', ');          % ',,'  to ', '
                Condition=strrep(Condition,', ,',', ');         % ', ,' to ', '
                Condition=strrep(Condition,' , ',', ');         % ' , ' to ', '
                Condition=strrep(Condition,'  ',' ');           % '  '  to ' '

                %

                %%% Correct Molar(m)/Micro-Molar(mm)/Temperature/pH Layout
                % Correct Molar(m)/Micro-Molar(mm) Layout
                m_Table={'m.', 'm,', 'm;', 'mm.', 'mm,', 'mm;'};
                for mT=1:size(m_Table,2)
                    m_Pos=strfind(Condition, m_Table{mT});
                    for mP=1:size(m_Pos,2)
                        if (~isnan(str2double(Condition(m_Pos(mP)-1)))&isreal(str2double(Condition(m_Pos(mP)-1))))||...
                                (~isnan(str2double(Condition(m_Pos(mP)-2)))&isreal(str2double(Condition(m_Pos(mP)-2))))
                            if mT>3
                                Condition(m_Pos(mP)+2)=' ';
                            else
                                Condition(m_Pos(mP)+1)=' ';
                            end
                        end
                    end
                end
                Condition=strrep(Condition,'  ',' ');

                % Correct Temperature Layout
                Condition=strrep(Condition,'temp.','temp ');
                Condition=strrep(Condition,'temp,','temp ');
                Condition=strrep(Condition,'temperature.','temperature ');
                Condition=strrep(Condition,'temperature,','temperature ');
                Condition=strrep(Condition,'  ',' ');

                if contains(Condition,'room temperature')
                    Condition=strrep(Condition,'room temperature',', room temperature, ');
                else
                    Condition=strrep(Condition,'temp',', temp');
                end
                Condition=strrep(Condition,'  ',' ');

                % Correct pH Layout
                Condition=strrep(Condition,'( ph',', ph');
                Condition=strrep(Condition,'(ph',', ph');
                Condition=strrep(Condition,'ph.','ph ');
                Condition=strrep(Condition,'ph,','ph ');
                Condition=strrep(Condition,'ph:','ph ');
                Condition=strrep(Condition,'ph/','ph ');
                Condition=strrep(Condition,'  ',' ');

                pH_Table={',ph', '-ph', ' ph'};
                for pHT=1:size(pH_Table,2)
                    pH_Pos=strfind(Condition, pH_Table{pHT});
                    for pHP=1:size(pH_Pos,2)
                        if (~isnan(str2double(Condition(pH_Pos(pHP)+3)))&isreal(str2double(Condition(pH_Pos(pHP)+3))))||...
                                (~isnan(str2double(Condition(pH_Pos(pHP)+4)))&isreal(str2double(Condition(pH_Pos(pHP)+4))))
                            Condition(pH_Pos(pHP))='$';
                        end
                    end
                end
                Condition=strrep(Condition,'$ph',', ph');

                Condition=strrep(Condition,',',', ');           % ',' to ', '
                Condition=strrep(Condition,'  ',' ');           % '  ' to ' '
                Condition=strrep(Condition,',,',', ');          % ',,'  to ', '
                Condition=strrep(Condition,', ,',', ');         % ', ,' to ', '
                Condition=strrep(Condition,' , ',', ');         % ' , ' to ', '
                Condition=strrep(Condition,'  ',' ');           % '  '  to ' '

                %

                %%% Content Fragmentation
                % All Gap
                Condition=strtrim(Condition);
                Gap_Anc=strfind(Condition, ' ');

                % Anchor Exception
                % Comma
                Comma_Anc=strfind(Condition,', ')+1;

                % Percentage
                Percent_Anc=strfind(Condition,'% ')+1;

                % Reagent
                PEG_Anc=[strfind(Condition,'peg ')+3, strfind(Condition,'(peg) ')+5, strfind(Condition,'polyethylene glycol ')+19, strfind(Condition,'poly ethylene glycol ')+20, strfind(Condition,'poly-ethylene glycol ')+20];
                PPG_Anc=[strfind(Condition,'ppg ')+3, strfind(Condition,'propylene glycol ')+16, strfind(Condition,'propylene glycol p ')+18];
                MPEG_Anc=[strfind(Condition,'mpeg ')+4, strfind(Condition,'pegmme ')+6, strfind(Condition,'peg-mme ')+7, strfind(Condition,'peg mme ')+7, strfind(Condition,'monomethyl ether ')+16];
                PPGBA_Anc=strfind(Condition,'ppgba ')+5;

                Jeffamine_Anc=[strfind(Condition,'jeffamine m ')+11, strfind(Condition,'jeffamine-m ')+11, strfind(Condition,'jeffamine ed ')+12, strfind(Condition,'jeffamine-ed ')+12];
                PVPK_Anc=[strfind(Condition,'polyvinylpyrrolidone k ')+22, strfind(Condition,'polyvinylpyrrolidone-k ')+22];
                Sokalan_Anc=[strfind(Condition,'cp ')+2, strfind(Condition,'hp ')+2, strfind(Condition,'pa ')+2];
                PAASalt_Anc=strfind(Condition,'salt) ')+5;
                NDSB_Anc=strfind(Condition,'ndsb ')+4;

                % Temp/pH
                Temp_Anc=[strfind(Condition,'temp ')+4, strfind(Condition,'temperature ')+11];
                pH_Anc=strfind(Condition,'ph ')+2;

                % Segregation Anchor Point (Break String Number with ,)
                Exception_Anc=unique([Comma_Anc, Percent_Anc, PEG_Anc, PPG_Anc, MPEG_Anc, PPGBA_Anc, Jeffamine_Anc, PVPK_Anc, Sokalan_Anc, PAASalt_Anc, NDSB_Anc, Temp_Anc, pH_Anc]);
                Seg_Anc=setdiff(Gap_Anc, Exception_Anc);

                % isdouble/isreal/not cmp with ,-/not oxylate/ not pga
                for Anc=1:size(Seg_Anc,2)
                    if (~isnan(str2double(Condition(Seg_Anc(Anc)+1))))&&isreal(str2double(Condition(Seg_Anc(Anc)+1)))&&...
                            ((Condition(Seg_Anc(Anc)+2))~=',')&&((Condition(Seg_Anc(Anc)+2))~='-')&&...
                            ((Condition(Seg_Anc(Anc)+2))~='/')&&((Condition(Seg_Anc(Anc)+3))~='/')&&...
                            ((Condition(Seg_Anc(Anc)+2))~='p')
                        Condition(Seg_Anc(Anc))='$';
                    end
                end
                Condition=strrep(Condition,'$',', ');

                Condition=strrep(Condition,',',', ');           % ',' to ', '
                Condition=strrep(Condition,'  ',' ');           % '  ' to ' '
                Condition=strrep(Condition,',,',', ');          % ',,'  to ', '
                Condition=strrep(Condition,', ,',', ');         % ', ,' to ', '
                Condition=strrep(Condition,' , ',', ');         % ' , ' to ', '
                Condition=strrep(Condition,'  ',' ');           % '  '  to ' '

                 % De-preotection
                Condition=strrep(Condition,'#',',');
                Condition=strrep(Condition,'^',' ');
                Condition=strrep(Condition,'  ',' '); 

                %
                %
                %

                %%%% Fill in the Condition_MetaTable
                % Condition Fragmentation Labelled with ', '
                Frag_L=strfind(Condition,', ');

                % Fill in the Condition_MetaTable
                for k=1:(size(Frag_L,2)-1)

                    % Reagent Content
                    Reagent=strtrim(Condition(Frag_L(k):Frag_L(k+1)-1));
                    if ~contains(Reagent,'(')||~contains(Reagent,')')
                        Reagent=strrep(Reagent,'(','');
                        Reagent=strrep(Reagent,')','');
                    end
                    
                    % Check Protein Status  (58: In-list Proteins)
                    protein_num=0;
                    for q=2:58
                        protein_num=protein_num+contains(Reagent, char(M{6,q}));
                    end 
                    
                    %%%% Protein
                    if protein_num~=0
                        Condition_MetaTable(6,si)=1;
                    elseif contains(Reagent,'mutant')||contains(Reagent,'mutation')||contains(Reagent,'purified')||contains(Reagent,'recombinant')||contains(Reagent,'protein')||contains(Reagent,'proteins')||contains(Reagent,'protien')||contains(Reagent,'protiens')
                        Condition_MetaTable(6,si)=1;

                    %%%% Missing Info
                    elseif isempty(strtrim(Reagent(2:end)))
                        Condition_MetaTable(3,si)=1;

                    %%%% Missing Conc. Value
                    elseif any(ismember(M{:,1:size(M,2)}, strtrim(Reagent(2:end))), 'all')
                        [row, ~]=find(ismember(M{:,1:size(M,2)}, strtrim(Reagent(2:end))));
                        Condition_MetaTable(row,si)=777;
                        Condition_MetaTable(5,si)=1;

                    %%%% Percentage (%)
                    elseif contains(Reagent,'% ')

                        % Reagent Content
                        Anchor=strfind(Reagent,'% ');
                        KEYWORD=strtrim(Reagent(Anchor+1:end));

                        % Missing Info
                        if isempty(KEYWORD)
                            Condition_MetaTable(3,si)=1;
                        else
                            % Reagent Examination
                            Index=ismember(M{:,1:size(M,2)},KEYWORD);
                            [row, ~]=find(Index);

                            if ~any(Index, 'all')                    
                                % Not-In_List Reagent
                                Condition_MetaTable(4,si)=1;
                                Other_Reagent(i).Reagent=KEYWORD;
                                i=i+1;
                            else
                                % Conc. Value
                                Num=strtrim(Reagent(2:Anchor-1));

                                if isempty(Num)
                                    % Missing Conc. Value
                                    Condition_MetaTable(row,si)=777;
                                    Condition_MetaTable(5,si)=1;
                                elseif strcmp(Num, 'saturate')||strcmp(Num, 'saturated')
                                    Condition_MetaTable(row,si)=666;                                                              % Saturated
                                else
                                    Num=strrep(Num, 'o', '0');
                                    Num=strrep(Num, 'l', '1');
                                    
                                    if ~contains(Num,'-')
                                        Condition_MetaTable(row,si)=str2double(Num);                                              % Specific Value
                                    else
                                        Sub=strfind(Num,'-');
                                        Condition_MetaTable(row,si)=mean([str2double(Num(1:Sub-1)), str2double(Num(Sub+1:end))]); % Range Average
                                    end
                                end
                            end
                        end

                    %%%% Micro-Molar (mm)
                    elseif contains(Reagent,'mm ')&&...
                            ((~isnan(str2double(Reagent(min(strfind(Reagent,'mm '))-1)))&...
                            isreal(str2double(Reagent(min(strfind(Reagent,'mm '))-1))))||...
                            strcmp(Reagent(min(strfind(Reagent,'mm '))-1),' '))                 
                        
                        % Reagent Content
                        Anchor=min(strfind(Reagent,'mm '));
                        KEYWORD=strtrim(Reagent(Anchor+2:end));

                        % Missing Info
                        if isempty(KEYWORD)
                            Condition_MetaTable(3,si)=1;
                        else
                            % Reagent Examination
                            Index=ismember(M{:,1:size(M,2)},KEYWORD);
                            [row, ~]=find(Index);

                            if ~any(Index, 'all')                    
                                % Not-In_List Reagent
                                Condition_MetaTable(4,si)=1;
                                Other_Reagent(i).Reagent=KEYWORD;
                                i=i+1;
                            else
                                % Conc. Value
                                Num=strtrim(Reagent(2:Anchor-1));

                                if isempty(Num)
                                    % Missing Conc. Value
                                    Condition_MetaTable(row,si)=777;
                                    Condition_MetaTable(5,si)=1;
                                else
                                    Num=strrep(Num, 'o', '0');
                                    Num=strrep(Num, 'l', '1');

                                    if ~contains(Num,'-')
                                        Condition_MetaTable(row,si)=str2double(Num);                                              % Specific Value
                                    else
                                        Sub=strfind(Num,'-');
                                        Condition_MetaTable(row,si)=mean([str2double(Num(1:Sub-1)), str2double(Num(Sub+1:end))]); % Range Average
                                    end
                                end
                            end
                        end

                    %%%% Molar (m)
                    elseif contains(Reagent,'m ')&&...
                            ((~isnan(str2double(Reagent(min(strfind(Reagent,'m '))-1)))&...
                            isreal(str2double(Reagent(min(strfind(Reagent,'m '))-1))))||...
                            strcmp(Reagent(min(strfind(Reagent,'m '))-1),' '))     
                        
                        % Reagent Content
                        Anchor=min(strfind(Reagent,'m '));
                        KEYWORD=strtrim(Reagent(Anchor+1:end));

                        % Missing Info
                        if isempty(KEYWORD)
                            Condition_MetaTable(3,si)=1;
                        else
                            % Reagent Examination
                            Index=ismember(M{:,1:size(M,2)},KEYWORD);
                            [row, ~]=find(Index);

                            if ~any(Index, 'all')                    
                                % Not-In_List Reagent
                                Condition_MetaTable(4,si)=1;
                                Other_Reagent(i).Reagent=KEYWORD;
                                i=i+1;
                            else
                                % Conc. Value
                                Num=strtrim(Reagent(2:Anchor-1));

                                if isempty(Num)
                                    % Missing Conc. Value
                                    Condition_MetaTable(row,si)=777;
                                    Condition_MetaTable(5,si)=1;
                                else
                                    Num=strrep(Num, 'o', '0');
                                    Num=strrep(Num, 'l', '1');

                                    if ~contains(Num,'-')
                                        Condition_MetaTable(row,si)=str2double(Num)*1000;                                               % Specific Value
                                    else
                                        Sub=strfind(Num,'-');
                                        Condition_MetaTable(row,si)=mean([str2double(Num(1:Sub-1)), str2double(Num(Sub+1:end))])*1000;  % Range Average
                                    end
                                end
                            end
                        end

                    %%%% w/v (mg/ml, mg/l, g/ml)
                    elseif contains(Reagent,'mg/ml')||contains(Reagent,'mg/l')||contains(Reagent,'g/ml')
   
                        % Reagent Content
                        Anchor=min(strfind(Reagent,'/'));
                        KEYWORD=strtrim(Reagent(Anchor+3:end));

                        % Missing Info
                        if isempty(KEYWORD)
                            Condition_MetaTable(3,si)=1;
                        else
                            % Reagent Examination
                            Index=ismember(M{:,1:size(M,2)},KEYWORD);
                            [row, ~]=find(Index);

                            if ~any(Index, 'all')                    
                                % Not-In_List Reagent/Protein
                                Condition_MetaTable(4,si)=1;
                                Other_Reagent(i).Reagent=KEYWORD;
                                i=i+1;
                                Other_Reagent(j).Protein=KEYWORD;
                                j=j+1;
                            else
                                % Conc. Value
                                Num=strtrim(Reagent(2:Anchor-3));

                                if isempty(Num)
                                    % Missing Conc. Value
                                    Condition_MetaTable(row,si)=777;
                                    Condition_MetaTable(5,si)=1;
                                else
                                    Num=strrep(Num, 'o', '0');
                                    Num=strrep(Num, 'l', '1');
                                    Condition_MetaTable(end,si)=1;

                                    if ~contains(Num,'-')
                                        Condition_MetaTable(row,si)=str2double(Num);                                              % Specific Value
                                    else
                                        Sub=strfind(Num,'-');
                                        Condition_MetaTable(row,si)=mean([str2double(Num(1:Sub-1)), str2double(Num(Sub+1:end))]); % Range Average
                                    end
                                end
                            end
                        end

                    %%%% Method
                    elseif contains(Reagent,'vapor')||contains(Reagent,'diffusion')||contains(Reagent,'lipidic cubic phase')
                        Condition_MetaTable(7,si)=1;
                    elseif contains(Reagent,'sit')
                        Condition_MetaTable(8,si)=1;
                    elseif contains(Reagent,'hang')
                        Condition_MetaTable(8,si)=2;
                    elseif contains(Reagent,'batch')                            
                        Condition_MetaTable(8,si)=3;

                    %%%% Temperature (temp: room temperature, deg c, k, no unit)
                    elseif contains(Reagent,'temp')

                        % Temp_P
                        if contains(Reagent,'temperature')
                            Temp_P=strfind(Reagent,'temperature')+11;
                        else
                            Temp_P=strfind(Reagent,'temp')+4;
                        end

                        % Temp_Unit_P
                        if contains(Reagent,'room')||contains(Reagent,'rt')     % room temperature
                            Condition_MetaTable(9,si)=298;
                        else
                            if contains(Reagent,'deg')                          % deg C
                                Temp_Unit_P=strfind(Reagent,'deg')-1;
                            elseif contains(Reagent,'oc')                           % deg C
                                Temp_Unit_P=strfind(Reagent,'oc')-1;
                            elseif contains(Reagent,'c')                            % deg C
                                Temp_Unit_P=strfind(Reagent,'c')-1;
                            elseif contains(Reagent,'k')                            % k
                                Temp_Unit_P=strfind(Reagent,'k')-1;
                            else                                                    % no unit
                                Temp_Unit_P=length(Reagent);
                            end

                            % Temp Value
                            Num=strtrim(Reagent(Temp_P:Temp_Unit_P));

                            if isempty(Num)
                                Condition_MetaTable(9,si)=777;
                            else
                                Num=strrep(Num, 'o', '0');
                                Num=strrep(Num, 'l', '1');
                                if ~contains(Num,'-')
                                    if str2double(Num)>200
                                        Condition_MetaTable(9,si)=str2double(Num);                                              % Specific Value
                                    else
                                        Condition_MetaTable(9,si)=str2double(Num)+273;                                          % Specific Value
                                    end
                                else
                                    Sub=strfind(Num,'-');
                                    mean_Num=mean([str2double(Num(1:Sub-1)), str2double(Num(Sub+1:end))]);                      % Range Average
                                    if mean_Num>200
                                        Condition_MetaTable(9,si)=mean_Num;
                                    else
                                        Condition_MetaTable(9,si)=mean_Num+273;
                                    end
                                end
                            end
                        end

                    %%%% Temperature (deg c)                        
                    elseif contains(Reagent,'deg')

                        % Temp Value
                        Num=strtrim(Reagent(3:strfind(Reagent,'deg')-1));

                        if isempty(Num)
                            Condition_MetaTable(9,si)=777;
                        else
                            Num=strrep(Num, 'o', '0');
                            Num=strrep(Num, 'l', '1');
                            if ~contains(Num,'-')
                                if str2double(Num)>200
                                    Condition_MetaTable(9,si)=str2double(Num);                                              % Specific Value
                                else
                                    Condition_MetaTable(9,si)=str2double(Num)+273;                                          % Specific Value
                                end
                            else
                                Sub=strfind(Num,'-');
                                mean_Num=mean([str2double(Num(1:Sub-1)), str2double(Num(Sub+1:end))]);                      % Range Average
                                if mean_Num>200
                                    Condition_MetaTable(9,si)=mean_Num;
                                else
                                    Condition_MetaTable(9,si)=mean_Num+273;
                                end
                            end
                        end

                    %%%% Temperature (oc/c/k)
                    elseif (~isempty(str2num(strrep(strrep(strtrim(erase(Reagent(2:end),'oc')),'o','0'),'l','1')))&isreal(str2num(strrep(strrep(strtrim(erase(Reagent(2:end),'oc')),'o','0'),'l','1'))))||...
                            (~isempty(str2num(strrep(strrep(strtrim(erase(Reagent(2:end),'c')),'o','0'),'l','1')))&isreal(str2num(strrep(strrep(strtrim(erase(Reagent(2:end),'c')),'o','0'),'l','1'))))||...
                            (~isempty(str2num(strrep(strrep(strtrim(erase(Reagent(2:end),'k')),'o','0'),'l','1')))&isreal(str2num(strrep(strrep(strtrim(erase(Reagent(2:end),'k')),'o','0'),'l','1'))))
                        
                        % Temp Value
                        Num=strrep(strrep(strtrim(erase(Reagent(2:end),'oc')),'o','0'),'l','1');
                        Num=erase(Num,'c');
                        Num=erase(Num,'k');
                        Num=strtrim(Num);

                        if isempty(Num)
                            Condition_MetaTable(9,si)=777;
                        else
                            if ~contains(Num,'-')
                                if str2double(Num)>200
                                    Condition_MetaTable(9,si)=str2double(Num);                                              % Specific Value
                                else
                                    Condition_MetaTable(9,si)=str2double(Num)+273;                                          % Specific Value
                                end
                            else
                                Sub=strfind(Num,'-');
                                mean_Num=mean([str2double(Num(1:Sub-1)), str2double(Num(Sub+1:end))]);                      % Range Average
                                if mean_Num>200
                                    Condition_MetaTable(9,si)=mean_Num;
                                else
                                    Condition_MetaTable(9,si)=mean_Num+273;
                                end
                            end
                        end

                    %%%% pH
                    elseif contains(Reagent,', ph')

                        % pH Value
                        Num=strtrim(Reagent(5:end));
                        Num=strrep(Num, 'o', '0');
                        Num=strrep(Num, 'l', '1');

                        if isempty(Num)
                            Condition_MetaTable(10,si)=777;
                        elseif strcmp(Reagent(5),' ')              
                            if ~isreal(str2num(Num))
                                Condition_MetaTable(10,si)=777;
                            else
                                if ~contains(Num,'-')
                                    Condition_MetaTable(10,si)=str2double(Num);
                                else
                                    Sub=strfind(Num,'-');
                                    Condition_MetaTable(10,si)=mean([str2double(Num(1:Sub-1)), str2double(Num(Sub+1:end))]);
                                end
                            end
                        else
                            if ~isreal(str2num(Num))||isempty(str2num(Num))
                                Condition_MetaTable(3,si)=1;
                            else
                                if ~contains(Num,'-')
                                    Condition_MetaTable(10,si)=str2double(Num);
                                else
                                    Sub=strfind(Num,'-');
                                    Condition_MetaTable(10,si)=mean([str2double(Num(1:Sub-1)), str2double(Num(Sub+1:end))]);
                                end
                            end
                        end

                    %

                    %%%% Others
                    % Temperature (rt)
                    elseif strcmp(strtrim(Reagent(2:end)),'rt')||contains(Reagent,' rt ')||contains(Reagent,' rt,')
                        Condition_MetaTable(9,si)=298;

                    % Only Num
                    elseif ~isempty(str2num(strtrim(strrep(strrep(Reagent(2:end), 'o', '0'), 'l', '1'))))&&isreal(str2num(strtrim(strrep(strrep(Reagent(2:end), 'o', '0'), 'l', '1'))))
                        
                        % Value
                        Num=strtrim(strrep(strrep(Reagent(2:end), 'o', '0'), 'l', '1'));

                        if ~contains(Num,'-')
                            mean_Num=str2double(Num);
                        else
                            Sub=strfind(Num,'-');
                            mean_Num=mean([str2double(Num(1:Sub-1)), str2double(Num(Sub+1:end))]);
                        end

                        if mean_Num<12
                            Condition_MetaTable(10,si)=mean_Num;
                        end    

                    % Unexpected Missing_Info
                    else
                        Condition_MetaTable(3,si)=1;
                        Missing_Info(k).Content=Reagent;
                        k=k+1;
                    end
                end
            end
        else
            % Empty (Null) Condition
            Condition_MetaTable(2,si)=1;
        end
    end
    si
end

% sum(sum(Condition_MetaTable(1:5,:))==0)
% writetable(struct2table(Missing_Info), 'Missing_Info.xlsx')
% writetable(struct2table(Other_Reagent), 'Other_Reagent.xlsx')

Condition_MetaTable=real(Condition_MetaTable);
csvwrite('Parsing_Result.csv',Condition_MetaTable)

