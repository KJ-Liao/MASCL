function PDBdata = PDB_CA_Extracting(FileName)

    %%%% [PDBdata] = readPDBfile(FileName);
    % outfile    (the name of the PDB, this is the only input on the command line)
    % recordName (the class or type of atom, such as ATOM, HETATM, SOL, etc)
    % atomName   (elemental identification of the atom)
    % X          (X position of atom)
    % Y          (Y position of atom)
    % Z          (Z position of atom)

    % Initialize file
    FileID = fopen(FileName);
    rawText = fread(FileID,inf,'*char');

    % Parse lines by end-of-lines
    splitLines = strread(rawText, '%s', 'delimiter', '\n');

    % Initialize variables
    numLines = length(splitLines);
    recordName = cell(1,numLines);
    atomName   = cell(1,numLines);
    resName    = cell(1,numLines);
    X          = cell(1,numLines);
    Y          = cell(1,numLines);
    Z          = cell(1,numLines);

    % Read each line
    m = 1;
    for n = 1:numLines
        thisLine = cell2mat(splitLines(n));
        if length(thisLine) > 53 && sum(isstrprop(thisLine(23:53), 'alpha')) == 0
            recordName(m) = {thisLine(1:6)};
            atomName(m)   = {thisLine(13:16)};
            resName(m)    = {thisLine(18:20)};
            X(m)          = {thisLine(31:38)};
            Y(m)          = {thisLine(39:46)};
            Z(m)          = {thisLine(47:54)};
            m = m + 1;
        end
    end

    % Trim exess
    keepData = logical(strcmp(recordName,'ATOM  ') + strcmp(recordName,'HETATM'));
    recordName = recordName(keepData);
    atomName   = atomName(keepData);
    resName    = resName(keepData);
    X          = X(keepData);
    Y          = Y(keepData);
    Z          = Z(keepData);

    % Reformat data for convenience
    
    Alldata.recordName = (strtrim(recordName))';
    Alldata.atomName   = (strtrim(atomName))';
    Alldata.resName    = (strtrim(resName))';
    Alldata.X          = (str2double(X))';
    Alldata.Y          = (str2double(Y))';
    Alldata.Z          = (str2double(Z))';

    fclose(FileID);

    %%%% PDBdata = parsePDBstruct(PDBdata)
    % Filter: without HETATOM (i.e. only ATOM)
    keepListTF = strcmp(Alldata.recordName, 'ATOM');

    % Filter: only preserve CA ATOM
    keepCAatomsTF = strcmp(Alldata.atomName, 'CA');

    % Filter without HETATOM and H ATOM
    keepListTF = (keepListTF + keepCAatomsTF) == 2;

    PDBdata.resName    = aminolookup(strrep(strrep([Alldata.resName(keepListTF)], 'SEC', 'CYS'), 'UNK', 'GLY'));
    PDBdata.X          = Alldata.X(keepListTF);
    PDBdata.Y          = Alldata.Y(keepListTF);
    PDBdata.Z          = Alldata.Z(keepListTF);

end