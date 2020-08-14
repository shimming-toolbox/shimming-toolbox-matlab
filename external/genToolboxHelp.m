function FileList = genToolboxHelp(ToolBoxPath, ToolBoxName, IgnoredFolders, evalCode, createDB)
%% GENTOOLBOXHELP Generate a Matlab help for a custom toolbox based on HTML
% This script generates a Matlab F1 help for a custom toolbox. It will run
% through all subfolders of the selected ToolBoxPath directory and create a
% Contents.m file in each subfolder containing .m files. It will then
% publish all .m files to .html files and place them in a hierarchical
% folder structure in a 'ToolBoxPath\Help\ToolBoxName' folder. If the
% respective .m file has already been published and was not changed since
% then, the html file will not be updated. The programm will proceed and
% create a hierarchical 'ToolBoxPath\Help\helptoc.xml' file listing the new
% help structure as well as 'ToolBoxPath\info.xml' file required by the
% Matlab help (see also
% http://de.mathworks.com/help/matlab/matlab_prog/display-custom-documentation.html).
% The help files can be accessed through Matlabs internal help via F1 -
% Supplementary Software - ToolBoxName. A GettingStarted.mlx file is created
% in the 'ToolBoxPath' folder which can be filled with nicely formatted
% content according to
% http://de.mathworks.com/help/matlab/matlab_prog/marking-up-matlab-comments-for-publishing.html
% that serves as a homepage for the F1 help.
% 
% Syntax: 
%
%   FileList = genToolboxHelp(ToolBoxPath, ToolBoxName, IgnoredFolders, evalCode, createDB)
%
% Inputs:
%
%   Directory of the Toolbox, can include subfolders--> ToolBoxPath
%   Desired Name of the Toolbox                   	--> ToolBoxName
%   List of folder names that should be ignored     --> IgnoredFolders
%   Evaluate code for publishing (true/false)       --> evalCode
%   Create a searchable database for the toolbox    --> createDB
%
% Outputs:
% 
%   List of all processed .m files (similar to dir) --> FileList
% 
%   A complete, hierarchical, MATLAB-accessible help of the custom toolbox
%   accessible by pressing F1 within Matlab when in the folder of the
%   custom toolbox. Creates a Help folder with all html help files
%   generated with Matlabs "publish" command.
%
% Example: 
%
%   FileList = genToolboxHelp('C:\Matlab\MyToolbox', 'MyToolbox', {TestScripts, OldScripts}, false, true)
%
% Other m-files required: none
%
% Subfunctions: autoTOC, cleanUpHelp, updateHtmlFile
%
% MAT-files required: none
%
% See also: 
%   http://de.mathworks.com/help/matlab/matlab_prog/display-custom-documentation.html
%   http://de.mathworks.com/help/matlab/ref/publish.html
% Author: Johannes Milvich
% email address: Johannes.Milvich@de.bosch.com
% August 2015; Last revision: 31-Aug-2015
% Copyright 2014-2015 Johannes Milvich
%------------- BEGIN CODE --------------% 
%% 1) Initialisation
% Check inputs - if no ignored folders are specified, only the default ones
% (.git and the Help folder itself) are ignored. By default, the code is not
% evaluated for the published .m files and a searchable database is
% created. Also check if the specified ToolBoxpath exists.
global MAKECONTENTSFILES
MAKECONTENTSFILES = false;  % Create Contents.m file within each folder.
global FOLDERHTML
FOLDERHTML = 'html';  % Folder name for output html files.

if nargin<3
    IgnoredFolders = {};
    evalCode = false;
    createDB = true;
elseif nargin<4
    evalCode = false;
    createDB = true;
elseif nargin<5
    createDB = true;
end
if ~exist(ToolBoxPath, 'dir')
    disp('ToolBoxPath does not exist')
    return
end
fprintf('**********Toolbox Help generation: **********\n');
% Naming
HelpFolderName  = 'doc';                                                   % If you want, you can change the default help folder name
IgnoredFolders  = [IgnoredFolders, {'.git', HelpFolderName}];             	% The files in these folders are not published to html files
% Set Help path and Help HTML path
pathHelpTOC     = [ToolBoxPath filesep HelpFolderName];                       	% Path of the helptoc.xml file
pathHelpFiles   = [pathHelpTOC, filesep, FOLDERHTML];                        	% Root path of the html help files (one level below pathHelpTOC)
% If no help has been created before, initialise the corresponding folder
% structure
if ~exist(pathHelpTOC,'dir')                                              
    mkdir(pathHelpTOC);
end                 
if ~exist(pathHelpFiles,'dir')
    mkdir(pathHelpFiles);
end

%% 2.1) Create a helptoc.xml file in the 'ToolBoxPath\Help' subfolder
% This helptoc.xml file hierarchically orders and links the .html files for
% accessing them correctly in the Matlab F1 help
% See also: http://de.mathworks.com/help/matlab/matlab_prog/display-custom-documentation.html
fid = fopen(fullfile(pathHelpTOC,'helptoc.xml'),'w');                       % Create the file/overwrite existing
fprintf(fid,'<?xml version=''1.0'' encoding="utf-8"?>\n\n<toc version="2.0">\n'); % Print initial lines of code

%% 2.2) Create a GettingStarted.mlx under doc/ if not present
% This GettingStarted.mlx file should be used to introduce people to your
% toolbox. It will be published by the matlab publish function and can have
% some examples and/or graphs etc. Code in here WILL be evaluated to show
% examples or nice graphs.
% See also: http://de.mathworks.com/help/matlab/matlab_prog/marking-up-matlab-comments-for-publishing.html
if ~exist(fullfile(ToolBoxPath, 'doc', 'GettingStarted.mlx'), 'file')
    edit(fullfile(ToolBoxPath, 'doc', 'GettingStarted.mlx'))
    msgbox('Please use GettingStarted.mlx to introduce to your toolbox. The file is getting published and serves as a homepage for the help directory. The code in here will be evaluated, so feel free to add nice examples and graphs that express your toolbox.');
end
% Create the home HTML file of your toolbox, GettingStarted.html. It will
% only be updated if the .m file has been changed.
% JULIEN: modified below for Matlab 2019+
% GSstatus = updateHtmlFile(fullfile(ToolBoxPath, 'doc', 'GettingStarted.mlx'), pathHelpFiles, ToolBoxPath);
matlab.internal.liveeditor.openAndConvert(fullfile(ToolBoxPath, 'doc', 'GettingStarted.mlx'), fullfile(pathHelpTOC, 'GettingStarted.html'))
fprintf('Getting Started html file has been updated\n\n')
% Add a top level table of contents item for GettingStarted.html
fprintf(fid,'\t<tocitem target="GettingStarted.html"> GettingStarted\n');

%% 2.3) Recursive file publishing and hierarchical help generation
% This is the workhorse of the script. It runs through the 'ToolBoxPath'
% and all non-ignored subfolders, publishes the. m files to .html files in
% the 'ToolBoxPath\Help\ToolBoxName' folder (or subfolders) and links the
% output .html files correspondingly in the helptoc.xml file.
fprintf('**********Processing folders: **********\n');
D = autoTOC(ToolBoxPath,[ToolBoxPath filesep '**' filesep '*.m'],0,fid,ToolBoxName,IgnoredFolders,evalCode, HelpFolderName);
fprintf('\n')

%% 3.4) Finish and close the helptoc.xml file
% Close the top-level item of the table of contents (GettingStarted.mlx
% entry), the table of contents as such as well as the helptoc.xml file
fprintf(fid,'\t</tocitem>\n');                                              
fprintf(fid,'</toc>');
fclose(fid);
% JULIEN: I commented the section below because we do not need to udpate
% info.xml at every doc generation. 
% %% 4) Create the necessary info.xml file
% % Includes information about the ToolBoxName and where the helptox.xml file
% % is located ('ToolBoxPath\Help\'). This file is located in the main
% % 'ToolBoxPath' folder. Information aabout the Matlab release and the help
% % type are also included.
% fid2 = fopen(fullfile(ToolBoxPath,'info.xml'),'w');
% fprintf(fid2,'<productinfo xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"\n');
% fprintf(fid2,'\txsi:noNamespaceSchemaLocation="optional">\n');
% fprintf(fid2,'\t<?xml-stylesheet type="text/xsl"href="optional"?>\n\t\n');
% fprintf(fid2,'\t<matlabrelease>%s</matlabrelease>\n',version('-release'));
% fprintf(fid2,'\t<name>%s</name>\n',ToolBoxName);
% fprintf(fid2,'\t<type>toolbox</type>\n');
% fprintf(fid2,'\t<icon></icon>\n');
% fprintf(fid2,'\t<help_location>%s</help_location>\n\n',pathHelpTOC);
% fprintf(fid2,'</productinfo>');
% fclose(fid2);

%% 5) Display feedback
% Let the user know which files have been created
fprintf('**********XML file generation: **********\n')
fprintf('According to the procedure in\nhttp://de.mathworks.com/help/matlab/matlab_prog/display-custom-documentation.html\ntwo .xml files were created:\n')
fprintf('1) info.xml \t\tin \t%s\n', ToolBoxPath) 
fprintf('2) helptoc.xml \t\tin \t%s\n', pathHelpTOC) 
fprintf('The corresponding published HTML files for each .m file can be found in:\n%s\n\n', pathHelpFiles)

%% 6) Create a searchable toolbox help database
% Can be used to search for keywords in your published .m-files
fprintf('**********Searchable database generation: **********\n')
if createDB
    builddocsearchdb(pathHelpTOC);
end

%% 7) Give back list of all processed files
% Similar to 'dir' command
fprintf('\n**********File list generation: **********\n')
fprintf('Display the output argument of this function\nfor a ''dir''-like file list\n')
FileList = '';
pp = {'' 'k' 'M' 'G' 'T'};
    for ii=1:length(D),
        sz = D(ii).bytes;
        if sz<=0,
            FileList = [FileList, sprintf(' %31s %-64s\n','',D(ii).name)];
        else
            ss = min(4,floor(log2(sz)/10));
            FileList = [FileList, sprintf('%4.0f %1sb   %20s   %-64s \n',sz/1024^ss,pp{ss+1},D(ii).date,D(ii).name)];
        end
    end
    
%% 8) Clean up old html files and folders of deleted scripts
% If you re-create the help of your toolbox and deleted various scripts or
% folders, this function will delete the corresponding help html files as
% well.
% removed_files = struct();
% fprintf('\n**********Removed html files (.m files do no longer exist): **********\n');
% JULIEN: I commented out the line below, because it is buggy:
%   Error using rmdir
%   /Users/julien/code/shimming-toolbox/doc/html/Shim_Greg is not a directory.
% TODO: implement html file cleaning BEFORE running this pipeline (not at 
% the end)
% cleanUpHelp(ToolBoxPath,pathHelpFiles,[pathHelpFiles filesep '**' filesep '*.html'],0,removed_files);
end_text = sprintf('*****Press F1 --> Supplementary Software --> %s *****', ToolBoxName);
num_chars = numel(end_text);
fprintf('\n%s',repmat('*',num_chars,1));
fprintf('\n*****Press F1 --> Supplementary Software --> %s *****\n', ToolBoxName);
fprintf('%s\n',repmat('*',num_chars,1));
end

%% RECURSIVE SUBFUNCTION TO GENERATE HELP FILES AND helptoc.xml
function [varargout] = autoTOC(rootDir, currDir, lvl, fid, ToolBoxName, IgnoredFolders, evalCode, HelpFolderName)
%% AUTOTOC Generate helptoc.xml file for a custom Matlab toolbox
% Generates .html help files for all of the scripts/functions in the folder
% currDIR (and its subfolders) and places them in a hierarchical structure
% in a folder called "Help" which will be located in the rootdir. The
% initial call of this function should start with lvl = 0 to indicate the
% top folder. Later, recursive calls to this function will have a higher
% level count. A hierarchical .xml file will be created that generates an
% html table. the .xml file has to be created in a higher script and the
% file identifier fid has to be given. The string in ToolBoxName will label
% the Help section. Any folders whose name is listed in IgnoredFolders will
% not be taken into account.
% Pathnames and wildcards may be used. Wild cards can exist in the pathname
% too. A special case is the double * that will match multiple directory
% levels, e.g. c:\**\*.m. Otherwise a single * will only match one
% directory level. e.g. C:\Program Files\Windows *\
% This function is based on the rdir function that can be found in the
% Matlab file exchange
%
% Syntax:
%   autoTOC(rootDir, currDir, lvl, fid, ToolBoxName, IgnoredFolders, evalCode)
%
% Inputs:
%
%   Root directory of the toolbox                   --> rootDdir
%   Current directory (recursive)                   --> currDir
%   Hierarchical level (0 for initial call)         --> lvl
%   File identifier for the helptoc.xml file		--> fid
%   Name of the toolbox                             --> ToolBoxName
%   List of folders that should not be documented   --> IgnoredFolders
%   Evaluate code for publishing                    --> evalCode
%   HelpFolderName                                  --> String. Output folder for html files.
%
% Outputs:
%
%   Creates a helptoc.xml file
%   if nargout=D: Give back currenty investigated directory (similar to the
%   dir command, with the exception that the name field will include
%   the relative path as well as the name to the file that was found)
%   if nargout=0: Output (dir) is sent to the screen
%
% Example:
%
%   autoTOC('C:\Matlab\MyToolbox','C:\Matlab\MyToolbox\**\*.m',0,fid,'MyToolbox',{TestScripts, OldScripts}, false)
%
% Other m-files required: none
%
% Subfunctions: updateHtmlFile
%
% MAT-files required: none
%
% See also: http://de.mathworks.com/help/matlab/matlab_prog/display-custom-documentation.html
%
% Author: Johannes Milvich
% email address: Johannes.Milvich@de.bosch.com
% August 2015; Last revision: 17-Aug-2015
% Copyright 2014-2015 Johannes Milvich

global MAKECONTENTSFILES 
global FOLDERHTML

%% Get path and wildcard information
% split the file path around the wild card specifiers
prepath = '';       % the path before the wild card
wildpath = '';      % the path wild card
postpath = currDir; % the path after the wild card
% 1. stage: find last \
I = find(currDir==filesep,1,'last');                                        % Find the last file separator (\ on windows) in the current path
if ~isempty(I),
    prepath = currDir(1:I);                                                   % Everything before the last \ is the prepath
    postpath = currDir(I+1:end);                                              % Everything after the last \ is the postpath
    
    % 2. stage: find first wildcard *
    I = find(prepath=='*',1,'first');                                         % Find the first wildcard parameter in the prepath
    if ~isempty(I),
        postpath = [prepath(I:end) postpath];
        prepath = prepath(1:I-1);
        
        % 3. stage: find next \ before wildcard, generate wildpath
        I = find(prepath==filesep,1,'last');
        if ~isempty(I),
            wildpath = prepath(I+1:end);
            prepath = prepath(1:I);
        end
        
        % 4. stage: Add part of postpath before \ to wildpath
        I = find(postpath==filesep,1,'first');
        if ~isempty(I),
            wildpath = [wildpath postpath(1:I-1)];
            postpath = postpath(I:end);
        end
    end
end
% Display current path
disp([' "' prepath '" ~ "' wildpath '" ~ "' postpath '" ']);
hlevel = lvl;
fileid = fid;

%% Generate help files and helptoc.xml of current folder (no wildpath)
if isempty(wildpath),
    % If no directory wildcards then just get file/folder list
    D = dir([prepath postpath]);
    if ~isempty(D) && MAKECONTENTSFILES
        makecontentsfile(prepath,'force');
        D = dir([prepath postpath]);
    end
    D([D.isdir]==1) = [];
    
    % For all file/folder names in the current path
    for ii = 1:length(D),
        % If the current name is a FILE
        if (~D(ii).isdir),
            D(ii).name = [prepath D(ii).name];
            % Make a nice, short, capital Name out of the file name
            shortname = strrep(D(ii).name,prepath,'');shortname = strrep(shortname,'.m','');
            shortname = regexprep(shortname,'(\<[a-z])','${upper($1)}');
            
            hleveltmp = hlevel;
            
            % Generate a folder for the help html files where the current file is
            help_folder_name = [rootDir filesep HelpFolderName filesep FOLDERHTML strrep(prepath,rootDir,'')];
            if ~exist(help_folder_name,'dir')
                mkdir(help_folder_name)
            end
            
            % Generate a html help file with the MATLAB PUBLISH script
            
            updateHtmlFile(D(ii).name, help_folder_name, rootDir);
            % Generate entries in the hierarchical .xml html help file
            if ~strcmp(shortname,'Contents') && ~strcmp(shortname,'GettingStarted')
                
                % If we are not at the top level or talking about a
                % contents.m file, print tabs in the .xml hierarchical file
                % according to the current hierarchical level.
                while hleveltmp>0
                    fprintf(fileid,'\t');
                    hleveltmp = hleveltmp -1;
                end
                
                % print an .xml entry (tocitem) with a link to the recently
                % created html file and the correct path and shortname
                fprintf(fileid,'\t\t<tocitem target="%s"> %s\n', [FOLDERHTML strrep(prepath,rootDir,'') shortname '.html'],shortname);
            end
            
            % Close the entries with </tocitem> according to the current
            % hierarchical level
            if ~strcmp(shortname,'Contents') && ~strcmp(shortname,'GettingStarted')
                hleveltmp = hlevel;
                while hleveltmp>0
                    fprintf(fileid,'\t');
                    hleveltmp = hleveltmp -1;
                end
                
                fprintf(fileid,'\t\t</tocitem>\n');
            end
        end
    end
    
    %% Generate help files and helptoc.xml of all sub-directories in current folder (with wildpath)
elseif strcmp(wildpath,'**'), % a double wild directory means recurs down into sub directories
    
    % first look for files in the current directory (remove extra filesep)
    D = autoTOC(rootDir,[prepath postpath(2:end)],hlevel,fileid, ToolBoxName, IgnoredFolders,evalCode, HelpFolderName);
    
    % then look for sub directories
    Dt = dir('');
    tmp = dir([prepath '*']);
    
    % Do not take into account the following subdirectories
    toDel = [];
    for i=1:numel(tmp);
        for j=1:numel(IgnoredFolders)
            if ~isempty(strfind(tmp(i).name,IgnoredFolders{j}))
                toDel = [toDel, i];
            end
        end
    end
    toDel = unique(toDel);
    
    for i=1:numel(toDel)
        tmp(toDel(end+1-i)) = [];
    end
    
    %% Recursive sub-directory execution
    % process each directory
    for ii = 1:length(tmp),
        if (tmp(ii).isdir && ~strcmpi(tmp(ii).name,'.') && ~strcmpi(tmp(ii).name,'..') ),
            
            hleveltmp = hlevel;
            while hleveltmp>0
                fprintf(fileid,'\t');
                hleveltmp = hleveltmp -1;
            end
            
            dirname = regexprep(tmp(ii).name,'(\<[a-z])','${upper($1)}');
            
            fprintf(fileid,'\t\t<tocitem target="%s"> %s\n', [FOLDERHTML strrep(prepath,rootDir,'') dirname filesep 'Contents.html'],dirname);
            
            Dt = [Dt; autoTOC(rootDir,[prepath tmp(ii).name filesep wildpath postpath],hlevel+1,fileid, ToolBoxName, IgnoredFolders, evalCode, HelpFolderName)];
            
            hleveltmp = hlevel;
            while hleveltmp>0
                fprintf(fileid,'\t');
                hleveltmp = hleveltmp -1;
            end
            
            
            fprintf(fileid,'\t\t</tocitem>\n');
            
            
        end
    end
    D = [D; Dt];
    
    %% Only go in specific subdirectories
else
    % Process directory wild card looking for sub directories that match
    tmp = dir([prepath wildpath]);
    D = dir('');
    % process each directory found
    for ii = 1:length(tmp),
        if (tmp(ii).isdir && ~strcmpi(tmp(ii).name,'.') && ~strcmpi(tmp(ii).name,'..') ),
            D = [D; autoTOC(rootDir,[prepath tmp(ii).name postpath],hlevel+1,fileid, ToolBoxName, IgnoredFolders, evalCode, HelpFolderName)];
        end
    end
end
%% Display listing if no output variables are specified
if nargout==0,
    pp = {'' 'k' 'M' 'G' 'T'};
    for ii=1:length(D),
        sz = D(ii).bytes;
        if sz<=0,
            disp(sprintf(' %31s %-64s','',D(ii).name));
        else
            ss = min(4,floor(log2(sz)/10));
            disp(sprintf('%4.0f %1sb   %20s   %-64s ',sz/1024^ss,pp{ss+1},D(ii).date,D(ii).name));
        end
    end
else
    % send list out
    varargout{1} = D;
end
end
%% RECURSIVE SUBFUNCTION TO CLEAN UP HELP FILES AND FOLDERS
function removed_files = cleanUpHelp(rootDir,pathHelpFiles,currDir,lvl,removed_files)
%% CLEANUP Removes old html files and folders
% Deletes old .html files and folders of .m files that have been removed
% from the toolbox.
%
% Syntax:
%   removed_files = cleanUpHelp(rootDir,pathHelpFiles,currDir,lvl,removed_files)
%
% Inputs:
%
%   Root directory of the toolbox                   --> rootDdir
%   Directory of the help html files                --> pathHelpFiles
%   Files to search for / to compare                --> currDir
%   Hierarchical level (0 for initial call)         --> lvl
%   List of files that have already been removed    --> removed_files
%
% Outputs:
%
%   List of removed html files                      --> removed_files
%
% Example:
%
%   removed_files = cleanUpHelp(ToolBoxPath,pathHelpFiles,[pathHelpFiles '\**\*.html'],0, removed_files);
%
% Other m-files required: none
%
% Subfunctions: none
%
% MAT-files required: none
%
% See also: http://de.mathworks.com/help/matlab/matlab_prog/display-custom-documentation.html
%
% Author: Johannes Milvich
% email address: Johannes.Milvich@de.bosch.com
% August 2015; Last revision: 17-Aug-2015
% Copyright 2014-2015 Johannes Milvich
%% Get path and wildcard information
% split the file path around the wild card specifiers
prepath = '';       % the path before the wild card
wildpath = '';      % the path wild card
postpath = currDir; % the path after the wild card
% 1. stage: find last \
I = find(currDir==filesep,1,'last');                                        % Find the last file separator (\ on windows) in the current path
if ~isempty(I),
    prepath = currDir(1:I);                                                   % Everything before the last \ is the prepath
    postpath = currDir(I+1:end);                                              % Everything after the last \ is the postpath
    
    % 2. stage: find first wildcard *
    I = find(prepath=='*',1,'first');                                         % Find the first wildcard parameter in the prepath
    if ~isempty(I),
        postpath = [prepath(I:end) postpath];
        prepath = prepath(1:I-1);
        
        % 3. stage: find next \ before wildcard, generate wildpath
        I = find(prepath==filesep,1,'last');
        if ~isempty(I),
            wildpath = prepath(I+1:end);
            prepath = prepath(1:I);
        end
        
        % 4. stage: Add part of postpath before \ to wildpath
        I = find(postpath==filesep,1,'first');
        if ~isempty(I),
            wildpath = [wildpath postpath(1:I-1)];
            postpath = postpath(I:end);
        end
    end
end
hlevel = lvl;
%% Generate help files and helptoc.xml of current folder (no wildpath)
if isempty(wildpath),
    % If no directory wildcards then just get file/folder list
    D = dir([prepath postpath]);
    D([D.isdir]==1) = [];
    
    % For all file/folder names in the current path
    for ii = 1:length(D),
        % If the current name is a FILE
        if (~D(ii).isdir),
            m_name = fullfile(rootDir,strrep(prepath,pathHelpFiles,''),strrep(D(ii).name,'.html','.m'));
            D(ii).name = [prepath D(ii).name];
            
            % Delete the help file if no corresponding m file exists
            if ~exist(m_name,'file')
                if isempty(fieldnames(removed_files))
                    removed_files = D(ii);
                else
                    removed_files = [removed_files; D(ii)];
                end
                
                delete(D(ii).name);
            end
        end
    end
    
    %% Generate help files and helptoc.xml of all sub-directories in current folder (with wildpath)
elseif strcmp(wildpath,'**'), % a double wild directory means recurs down into sub directories
    % First, look for files in the current directory (remove extra filesep)
    removed_files = cleanUpHelp(rootDir,pathHelpFiles,[prepath postpath(2:end)],hlevel,removed_files);
    
    % then look for sub directories
    %Dt = dir('');
    tmp = dir([prepath '*']);
    %% Recursive sub-directory execution
    % process each directory
    for ii = 1:length(tmp),
        if (tmp(ii).isdir && ~strcmpi(tmp(ii).name,'.') && ~strcmpi(tmp(ii).name,'..') ),
            dirname = regexprep(tmp(ii).name,'(\<[a-z])','${upper($1)}');
            removed_files = cleanUpHelp(rootDir,pathHelpFiles,[prepath tmp(ii).name filesep wildpath postpath],hlevel+1,removed_files);
            % If the directiory of these help files exists, but no
            % corresponding directory in the m-files section, delete the
            % helpfile folder
            if ~exist(fullfile(rootDir,dirname),'dir')
                rmdir(fullfile(pathHelpFiles,dirname),'s');
            end           
        end
    end
    
    %% Only go in specific subdirectories
else
    % Process directory wild card looking for sub directories that match
    tmp = dir([prepath wildpath]);
    D = dir('');
    % process each directory found
    for ii = 1:length(tmp),
        if (tmp(ii).isdir && ~strcmpi(tmp(ii).name,'.') && ~strcmpi(tmp(ii).name,'..') ),
                removed_files = cleanUpHelp(rootDir,pathHelpFiles,[prepath tmp(ii).name postpath],hlevel+1,removed_files);
        end
    end
end
%% Display listing if no output variables are specified
if nargout==0,
    pp = {'' 'k' 'M' 'G' 'T'};
    if isempty(fieldnames(removed_files))
        fprintf('No files have been removed\n');
    else
    for ii=1:length(removed_files),
        sz =  removed_files(ii).bytes;
        if sz<=0,
            disp(sprintf(' %31s %-64s','', removed_files(ii).name));
        else
            ss = min(4,floor(log2(sz)/10));
            disp(sprintf('%4.0f %1sb   %20s   %-64s ',sz/1024^ss,pp{ss+1}, removed_files(ii).date, removed_files(ii).name));
        end
    end
    end
else
    % send list out
    varargout{1} = removed_files;
end
end
function status = updateHtmlFile(m_file_name, help_folder_name, rootDir)
%% UPDATEHTMLFILE Updates the html help file if the .m file has been changed
% Checks if the .m file has been changed since the last help generation and
% updates the corresponding help file.
%
% Syntax:
%   status = updateHtmlFile(m_file_name, help_folder_name, rootDir)
%
% Inputs:
%
%   File name of the .m file                        --> m_file_name
%   Folder of the html help                         --> help_folder_name
%   ToolBox root path                               --> rootDir
%
% Outputs:
%
%   Flag if the html file has been updated (1 = yes)--> status
%
% Example:
%
%   was_updated = updateHtmlFile('test.m', '\Toolbox\Help\ToolBoxName', '\Toolbox')
%
% Other m-files required: none
%
% Subfunctions: none
%
% MAT-files required: none
%
% See also: http://de.mathworks.com/help/matlab/matlab_prog/display-custom-documentation.html
%
% Author: Johannes Milvich
% email address: Johannes.Milvich@de.bosch.com
% August 2015; Last revision: 31-Aug-2015
% Copyright 2014-2015 Johannes Milvich
help_file_name = strrep(strrep(m_file_name,rootDir,''),'.m','.html');
help_file_name_fp = fullfile(help_folder_name,help_file_name);
if ~exist(help_file_name_fp,'file')
    % If a helpfile of that name does not exist yet, directly
    % publish it to the right folder
    publish(m_file_name, 'format','html','outputDir',help_folder_name,'evalCode',false);
    status = 1;
else
    % If a helpfile of that name already exists, publish to a
    % temporary file, see if it is different from the previous
    % one and if, copy the temporary file to the right folder.
    % otherwise, just delete the temporary file and folder and
    % keep the old one
    tmp_folder = fullfile(help_folder_name,'tmp');
    mkdir(tmp_folder);
    tmp_help_file_name_fp = publish(m_file_name, 'format','html','outputDir',tmp_folder,'evalCode',false);
    
    [status,~] = system(['fc ' tmp_help_file_name_fp ' ' help_file_name_fp]);
    if status
        copyfile(tmp_help_file_name_fp,help_folder_name);
    end
    rmdir(tmp_folder,'s');
end
end
%------------- END CODE --------------% 
