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
% seperate each folder in cells and find longest tree structure
parts = {};
% nLevels = 0;
nFiles = numel(Dr.docFiles);
for iFile = 1:nFiles
    parts{iFile} = strsplit(Dr.docFiles(iFile), filesep);

    % remove first letter when its a class ('@')
    fieldName = parts{1,iFile}(end-1);
    if (fieldName{1}(1) == '@')
        fieldName = strip(fieldName,'@');
        parts{1,iFile}(end-1) = fieldName;
    end
end

% add each file one by one
nav = {{}};
for iFile = 1:nFiles
    nav = addFileLayer(nav, parts, iFile, Dr.extOut);
end

% assemble
yml = struct('site_name', 'Shimming_Toolbox', 'theme', theme, 'nav', {{}});
yml.nav = nav;

% write
YAML.write(filePath, yml)

end


function [outNav] = addFileLayer(nav, parts, iFile, ext)
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
        for iNavX = 1:nNavX

            expressionIsEmpty = strcat(expressionLastGet, '{', string(iNavX), ',1}' );

            if ~isempty(eval(expressionIsEmpty))
                
                expressionIsField = strcat('isfield(', expressionLastGet, '{', string(iNavX), ',1}', ', parts{1,iFile}(', string(iPart1), '))');
                
                % isField()
                if eval(expressionIsField)
                    found = 1;
                    foundWhere = iNavX;
                    break;
                end
                
            else
                firstIsEmpty = 1;
            end
            
        end

        % if it is not, create it
        if ~found
            %if its the last element aka the file
            if iPart1 >= nParts
                structAssign = struct( erase(parts{1,iFile}(iPart1), ext),char(parts{1,iFile}(iPart1)));
                expressionAssign = strcat(expressionLastGet, ' = structAssign');
                foundWhere = nNavX+1;
            else
                if firstIsEmpty
                    expressionAssign = strcat(expressionLastGet, '{1,1}.(parts{1,iFile}(', string(iPart1), ')) = {{}}');
                    %nav{1,1}.(parts{1,iFile}(iPart1)) = {{}};
                    foundWhere = 1;
                else 
                    expressionAssign = strcat(expressionLastGet, '{', string(nNavX+1), ',1}.(parts{1,iFile}(', string(iPart1), ')) = {{}}');
                    %nav{nNavX+1,1}.(parts{1,iFile}(iPart1)) = {{}};
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

