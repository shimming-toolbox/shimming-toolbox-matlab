function [] = printYml( Dr, filePath, Params)
%PRINTYML print yml configuration file according to docFiles
%
% __Syntax__
%
%   [] = PRINTYML( Self, filePath) 
%   [] = PRINTYML( Self, filePath, Params) 
%
% __OPTIONS__
%
% `filePath` is the path of the file where the .yml file will be generated
%
% `Params.theme` is a character vector containing the name of the them to
% use
%
% `Params.projectName` is a character vector of the name of the project
% (will be displayed as the site name)
%
% % `Params.home` is a character vector of the homepage (default value is
% 'index.md')
%
% `Params.repoURL` is a character vector that specifies the link to the
% remote directory.
%
% for the Yml configuration file to work, the documentation must be in a
% folder called docs in the same folder

if nargin == 2
    
    if isstruct(filePath)
        return ;
    end
    
    Params.theme = 'material' ;
    Params.projectName = 'Your_Project' ;
    
elseif nargin == 3
    
    if ~isfield(Params, 'theme')
        Params.theme = 'material' ;
    else
        mustBeStringOrChar( Params.theme ) ;
        Params.theme = char(Params.theme) ;
    end
    
    if ~isfield(Params, 'projectName')
        Params.projectName = 'Your_Project' ;
    else
        mustBeStringOrChar( Params.projectName ) ;
        Params.projectName = char(Params.projectName) ;
    end
    
    if ~isfield(Params, 'home')
        Params.home = 'index.md' ;
    else
        mustBeStringOrChar( Params.home ) ;
        Params.home = char(Params.home) ;
    end
    
    if ~isfield(Params, 'repoURL')
        Params.repoURL = '' ;
    else
        mustBeStringOrChar( Params.repoURL ) ;
        Params.repoURL = char(Params.repoURL) ;        
    end
else
    return ;
end

% theme and homepage
theme = struct('name', Params.theme) ;
homepage = struct('Home', Params.home) ;
%% navigation

% seperate each folder in cells to be easier to work with
parts = {};
nFiles = numel(Dr.dFiles);
for iFile = 1:nFiles
    parts{iFile} = strsplit(Dr.dFiles(iFile), filesep);
    
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
nav{1} = homepage;
for iFile = 1:nFiles
    nav = addFileLayer(nav, parts, iFile, Dr.extOut, Dr.dirOutTop, Dr.dFiles(iFile));
end

%% assemble
if isempty(Params.repoURL)
    yml = struct('site_name', Params.projectName, 'theme' ,theme, 'nav', {{}});
else
    yml = struct('site_name', Params.projectName, 'theme', theme ,'repo_url', Params.repoURL , 'nav', {{}});
end
    
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