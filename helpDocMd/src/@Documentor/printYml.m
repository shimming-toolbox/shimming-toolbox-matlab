function [] = printYml( Dr, filePath)
%PRINTYML print yml configuration file according to docFiles
%
% ### Syntax
%
%   [] = PRINTYML( Self, filePath) 
%
% for the Yml configuration file to work, the documentation must be in a
% folder called docs in the same folder

% theme and homepage
theme = struct('name', 'material');
home = struct('Home', 'index.md');

% navigation

% seperate each folder in cells to be easier to work with
parts = {};
nFiles = numel(Dr.docFiles);
for iFile = 1:nFiles
    parts{iFile} = strsplit(Dr.docFiles(iFile), filesep);
    
    % remove first letter when its a class ('@') TODO remove whole cell
    % when the save recursive doesnt save the .md file in the class folder
    fieldName = parts{1,iFile}(end-1);
    if (fieldName{1}(1) == '@')
        fieldName = strip(fieldName,'@');
        parts{1,iFile}(end-1) = fieldName;
    end
end

% add each file one by one
nav = {};
for iFile = 1:nFiles
    nav = addFileLayer(nav, parts, iFile, Dr.extOut, Dr.dirOutTop, Dr.docFiles(iFile));
end

% assemble
yml = struct('site_name', 'Shimming_Toolbox', 'theme', theme, 'nav', {{}});
yml.nav = nav;

% write
YAML.write(filePath, yml)

% Quickfix TODO : write it properly
fixWrite(filePath);
disp('done')
end


function [outNav] = addFileLayer(nav, parts, iFile, ext, dirOutTop, fullpath)
    nParts = size(parts{1,iFile},2);
    
    % find the first different element of parts
    isDifferent = 0;
    partDifferent = 0;
    for iPart = 1:nParts
        if ~isDifferent
            % start at first different name
            oldName = parts{1,1}(iPart);
            for iFile1 = 1:numel(parts)
                name = parts{1,iFile1}(iPart);

                if oldName ~= name
                    isDifferent = 1;
                    partDifferent = iPart;
                    break;
                end
            end
            
            if isDifferent
                break;
            end
        end
    end
    if ~isDifferent
        partDifferent = nParts - 1;
    end
    
    % build "nav" depending on tre structure
    expression = 'nav';

    for iPart1 = partDifferent:nParts

        expressionLastGet = expression;

        % find if field is already present in nav
        found = 0;
        foundWhere = 0;
        firstIsEmpty = 0;
        %find size of last input struct
        expressionSize = strcat('size(', expressionLastGet,',1)');
        nNavX = eval(expressionSize);
        
        if nNavX == 0
            firstIsEmpty = 1;    
        elseif iPart1 < nParts
            for iNavX = 1:nNavX

                expressionIsField = strcat('isfield(', expressionLastGet, '{', string(iNavX), ',1}', ', parts{1,iFile}(', string(iPart1), '))');

                % isfield()
                if eval(expressionIsField)
                    found = 1;
                    foundWhere = iNavX;
                    break;
                end

            end
        end
        
        % if it is not, create it
        if ~found
            %if its the last element aka the file
            if iPart1 >= nParts
                %find path
                relativePath = erase(fullpath,strcat(dirOutTop,filesep));
                structAssign = struct( erase(parts{1,iFile}(iPart1), ext),char(relativePath));
                expressionAssign = strcat(expressionLastGet, '{nNavX+1,1}', ' = structAssign;');
                foundWhere = nNavX+1;
            else
                if firstIsEmpty
                    expressionAssign = strcat(expressionLastGet, '{1,1}.(parts{1,iFile}(', string(iPart1), ')) = {};');
                    %nav{1,1}.(parts{1,iFile}(iPart1)) = {};
                    foundWhere = 1;
                else 
                    expressionAssign = strcat(expressionLastGet, '{', string(nNavX+1), ',1}.(parts{1,iFile}(', string(iPart1), ')) = {};');
                    %nav{nNavX+1,1}.(parts{1,iFile}(iPart1)) = {};
                    foundWhere = nNavX+1;
                end
            end
            %eval is needed because the structure becomes very large so
            %remembering how it was created is necessary
            eval(expressionAssign);
        end
        
        % make sure the fiel has been created or has been found
        assert(foundWhere ~= 0,'Field was not found and was not created')
        
        % prepare next iteration of for loop
        expression = strcat(expressionLastGet, '{', string(foundWhere), ',1}.(parts{1,iFile}(',string(iPart1),'))');
    end
    %
    outNav = nav;

end

function [] = fixWrite(filePath)
    % open file
    [fid, errMsg] = fopen( filePath, 'r' ) ;
    assert( fid~=-1, ['Read failed: ' errMsg], '%s' ) ;
    
    text = fread(fid,'*char');
    
    % Find if the following line of character exists "- -". Deletes "- " if
    % found
    for itext  = 1:size(text,1)
        if text(itext) == '-'
            if text(itext+1) == ' '
                if text(itext+2) == '-'
                    text(itext) = ' ';
                end
            end
        end
    end
    
    % Close file
    fclose(fid);
    
    % Open the same file but now empty
    [fid, errMsg] = fopen( filePath, 'w+' ) ;
    assert( fid~=-1, ['write failed: ' errMsg], '%s' ) ;
    
    % Write the fixed characters
    fwrite(fid,text, '*char');
    
    % close file
    fclose(fid);
end