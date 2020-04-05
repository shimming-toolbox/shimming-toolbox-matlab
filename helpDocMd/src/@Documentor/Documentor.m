classdef Documentor < handle
% DOCUMENTOR Custom Matlab documentation into markup/down text files
%
% The Documentor prints custom Matlab documentation to text files 
% (principally, as [Markdown](https://daringfireball.net/projects/markdown/).
% 
% Refer to the README document for basic usage.
% 
% __CONSTRUCTOR SYNTAX__
%
%     Dr = Documentor(  ) ;
%     Dr = Documentor( src ) ;
%     Dr = Documentor( src, params ) ;
% 
% Typically a Documentor instance `Dr` will be created by calling the
% constructor with at least the first argument `src`: a string vector of file
% paths to .m source files and/or directories in which to search for them.
% Upon initialization, documentation can be published using the default
% settings by calling `Dr.printdoc();`
%
% `params.isSearchRecursive` is a scalar logical specifying whether any subdirectories
% of those listed in `src` are to be included in a file search. By default, it
% will be `true` when `src` contains just a single directory, and `false` when
% multiple directories are present. These defaults can be bypassed in either
% case by passing the corresponding boolean as the second argument. The value
% retained as a public property in the returned object (`Dr.isSearchRecursive`)
% for future reference.
%
% `params.outputDir` is a path specifying where the output documentation will be
% generated.
%
% `src` is used to initialize the public property `mFiles`, which contains the
% list of source files to be included in any published output, and which can be
% reconfigured after initialization.
% 
% Note that any invalid or *incompatible* paths suggested by `src` will be
% automatically filtered out from assignments to `mFiles`. That is, by setting
% `Dr.mFiles = src`, the actual assignment will always be 
% `Dr.mFiles = Documentor.findfiles( src, Dr.isSearchRecursive )`. 
% For more info, refer to the documentation entries for Documentor.findfiles
% and Documentor.mFiles.
%
% __ETC__
% 
% To test how a markdown sample will display when reformatted to HTML:
% - <https://daringfireball.net/projects/markdown/dingus>




% __OPTIONS__
% 
% ...TODO
%
% ### Basic usage
% 
% To construct a Documentor instance (e.g. one called `Dr`), call:
%    
%    Dr = Documentor( src ) ;
%
% Input `src` is a string vector of file system paths: these can refer to 
% specific .m files of interest, and/or as source code directories. (The latter are
% searched for compatible .m files).
% 
% To print all of the documentation to file, call:
%    
%    Dr.printdoc( ) ;
%
% #### Demo: Documentor.m 
% 
% To print documentation for the Documentor class: 
%
% 1. Ensure the class is accessible on the MATLAB path by typing 
%   
%   which('Documentor') 
%
% into the command prompt. If this displays 'Documentor not found', then find
% the file Documentor.m and move to its parent folder via `cd`.
% 
% 2. If `which('Documentor')` displays the correct path to Documentor.m, the
% documentation should now be printable, using the default configuration, by calling
%     
%    Dr = Documentor( which('Documentor') ) ;
%
%    Dr.printdoc( ) ;
%
% If successful, the documentation file path is displayed in the command window.
%
% ### General usage
% 
% Read the section on **Basic Usage** first!
%
%
% #### Documenting .m files




properties( Constant )

    % Default file extensions
    Ext = struct( 'matlab', ".m", 'markdown', ".md" ) ;

end

properties( Access=private, AbortSet )
    % Informer object-array: Provides info on each .m file in `mFiles`
    Info Informer ; 
    % Info Informer = Informer( which( "Documentor.m" ) ) ;
end 

properties( AbortSet )

    % Index of next .m file in list of files to document (i.e. `mFiles(iM)`)
    iM(1,1) uint64 {mustBePositive, mustBeInteger} = uint64(1) ;

    % List of .m files to document (string vector of full file paths)
    %
    % To ensure that the list consists solely of *documentable* files, property
    % reassignments (or, *assignments*, in the case of object construction) are
    % mediated (filtered) by an implicit function call: i.e. whenever the
    % property is set (as in `Dr.mFiles = src`) it is, in effect, as a return value
    % (namely, `Dr.mFiles = Documentor.findfiles( src, Dr.isSearchRecursive ) ;`).
    % 
    % __ETC__ 
    %
    % For details regarding the implementation and what constitutes a
    % 'documentable' file, 
    %
    % See also 
    % Documentor.findfiles 
    mFiles {mustBeFileOrFolder} = string( [ mfilename('fullpath') '.m' ] ) ;

    % List of output documentation file paths (string vector of full file paths)
    dFiles {mustBeStringOrCharOrCellstr} = "" ; 

    % Toggle whethers to overwrite existing documentation files (logical column vector with length == numel(mFiles))
    isOverwriting(:,1) {mustBeBoolean} = false ;
    
    % Toggle whether subdirectories are included .m file search 
    % 
    % (i.e when calling `Documentor.mFiles = src` and `src` contains directory paths)
    % 
    % See also
    % -Documentor.findfiles
    % -Documentor.mFiles
    isSearchRecursive(1,1) {mustBeBoolean} = true ;
   
    % Toggle to recreate original directory tree in `dirOutTop` (multiple mFiles case only) 
    %
    % __EXAMPLE__
    % 
    % Given Documentor object instance `Dr`:
    %
    % `if Dr.isSaveRecursive == false` 
    % ...all the output documentation files will be written to the same folder (i.e. `Dr.dirOutTop`)
    % 
    % `else`
    % ...an attempt is made to mirror the original directory tree of the top .m file source folder
    %
    % __ETC__
    %
    % See also 
    % HelpDocMd.isSearchRecursive
    isSaveRecursive(1,1) {mustBeBoolean} = true ;
    
    % Output parent directory for the doc files
    %
    % **Example**
    % 
    % If the top shared directory of the input .m Files is ` "/Users/Buddy/Projects/helpDocMd/src/" `;
    % then, by default, ` dirOutTop = "/Users/Buddy/Projects/helpDocMd/doc/" ;`
    % 
    % See also Documentor.isSaveRecursive
    dirOutTop {mustBeStringOrChar} = "" ;     
    
    % Documentation string for `mFiles(iM)` to be printed to `dFiles(iM)`
    mdDoc string {mustBeStringOrChar} = "" ;
  
    % String specifier for output syntax: "mkd" (for Mkdocs markdown), "mat" (for MATLAB markup)
    %
    % The sole difference between "mkd" and "mat" (for now) is that "mkd" will
    % reformat the style of any embedded links in the comments. 
    syntax(1,1) string { mustBeMember( syntax, ["mat" "mkd"] ) } = "mkd" ;

    % Toggles between basic/user (=false) and detailed/developer documentation (=true) [default: true]
    % 
    % When false, classes and class members with private, protected, or hidden
    % attributes are excluded from the output documentation. [default = true]
    %
    % NOTE/TODO: only partially implemented! (also, a selection of *degrees* of details rather than 0/1 might be better) 
    isDetailed(1,1) {mustBeBoolean} = true ;

    % Output file extension (default = ".md")
    extOut(1,1) {mustBeStringOrChar} = Documentor.Ext.markdown ;
    
end

properties( Access=private, Dependent )

    % Parent folder of next file to document (i.e. `mFiles(iM)`)
    mDir {mustBeFolder} = string( fileparts( mfilename('fullpath') ) ) ;
    
    % Top directory of src mFiles
    dirInTop {mustBeStringOrChar} = "" ;

end

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Dr = Documentor( pathIn, params )

    if nargin == 0
        return ;
    elseif nargin < 2
        params.isSearchRecursive = true ;
    elseif nargin == 2
        
        % assign output directory if not empty
        if isfield( params, 'outputDir' )
            Dr.dirOutTop = params.outputDir ;
        end
        
        % assign default if empty
        if ~isfield( params, 'isSearchRecursive' )
            params.isSearchRecursive = true ;
        end
        
    end
    
    Dr.mFiles = Documentor.findfiles( pathIn, params.isSearchRecursive ) ;
    Dr.Info   = Informer( Dr.mFiles(1) ) ;
    Dr.dFiles = "_DEFAULTS_" ;
    
end
% =========================================================================    
function [dirInTop] = get.dirInTop( Dr )
 
    % find folder with fewest parent directories (i.e. slashes in path): 
    parts = {};
    nFiles = numel(Dr.mFiles);
    for iFile = 1:nFiles
        parts{iFile} = strsplit(Dr.mFiles(iFile), filesep);
    end
    nParts = size(parts{1,iFile},2);
    
    % find the first different element of parts
    isDifferent = 0;
    
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
        partDifferent = nParts;
    end
    
    path = "";
    for iPart = 1:(partDifferent-1)
        if iPart ~= (partDifferent-1)
            path = strcat(path, parts{1,1}(iPart), filesep);
        else
            path = strcat(path, parts{1,1}(iPart));
        end
    end
    dirInTop = path;
end
% =========================================================================    
function [dirOutTop] = get.dirOutTop( Dr )

    if strcmp( Dr.dirOutTop, "" )
        try  % output docs folder in parent directory of ./mFiles
            Dr.dirOutTop = strcat( Dr.dirInTop, filesep, "docs" ) ;
        catch ME
            warning( [ 'Failed to create default dirOutTop directory.\n' ... 
                      'Assign the doc output directory manually.' ], '%s' ) ;
           rethrow(ME) ;
        end 
    end
    
    dirOutTop = Dr.dirOutTop ;

end
% =========================================================================    
function [] = set.dirOutTop( Dr, dirOutTop )

    [isValidPath, Values, msgId] = fileattrib( dirOutTop ) ;
    
    if isValidPath
        assert( Values.directory, 'Invalid assignment: Input is not a directory') ;
        assert( Values.UserWrite, [ 'Invalid permissions: Cannot write to assigned directory.\n', ...
           'Change permissions to the given directory with \n fileattrib( ' dirOutTop ',  +w ) \n' ...
            'or choose a different directory' ], '%s' ) ;
    else
        [isMade, msg, msgId] = mkdir( dirOutTop ) ;
        assert( isMade, msgId, ['Directory creation failed: ' msg ] )
    end
    
    Dr.dirOutTop = string( dirOutTop ) ;
    % Dr.dFiles is set again because if the output directory is changed,
    % the output path names have to be regenerated. This triggers the
    % set.dFiles method.
    Dr.dFiles = "_DEFAULTS_" ; 

end
% =========================================================================    
function [] = set.dFiles( Dr, dFiles )  
        
    dFiles = string( dFiles ) ;

    if isequal( dFiles, "_DEFAULTS_" )

        dFiles = strings( size( Dr.mFiles ) ) ;
        
        for iM = 1 : length( Dr.mFiles )
            
            [folder, docName] = fileparts( Dr.mFiles(iM) ) ;
            
            if Dr.isSaveRecursive % try to recreate subdir structure
                % remove parent dir component
                docFolder = [ Dr.dirOutTop + erase( folder, Dr.dirInTop ) ] ;
                
                % class files are need to bbe documented one level higher
                % in tree because they are one folder deeper
                parts = strsplit(Dr.mFiles(iM), filesep);
                fieldName = parts(end-1);
                if (fieldName{1}(1) == '@')
                    % goes to the parent directory if its a class
                    docFolder = fileparts(docFolder);
                end
                
                [isMade, msg, msgId] = mkdir( docFolder ) ;
                
                assert( isMade, msgId, ['Directory creation failed: ' msg ], '%s' ) ; 
            else
                docFolder = Dr.dirOutTop ;
            end
                
            dFiles(iM) = [ docFolder + filesep + docName + Dr.extOut ] ;

        end
    end
    
    assert( length( unique(dFiles) ) == length( Dr.mFiles ), ...
        "The list of output file paths must possess unique entries for each of the input .m files"  ) ;
    
    Dr.dFiles = dFiles ;

end
% =========================================================================    
function [] = set.iM( Dr, iM )  

    % update Informer with new file:
    Dr.Info.mFile = Dr.mFiles( iM ) ; 
    Dr.iM         = iM ;

end
% =========================================================================    
function [] = set.mFiles( Dr, mFiles )  

    [mFiles] = Documentor.findfiles( mFiles, Dr.isSearchRecursive ) ;
    
    if isempty( mFiles ) || isequal( mFiles, "" )
        warning( 'Nothing assigned: No Documentor-compatible .m files were found using the given parameters.') ;
        return ;
    end

    Dr.mFiles = mFiles ;

    if strcmp( Dr.dirOutTop, "" )
       Dr.dirOutTop = Dr.dirInTop ;
    end

end
% =========================================================================    
function [mdDoc] = get.mdDoc( Dr )  
    
    % TODO: Choose desired formating!
    switch Dr.Info.Attributes.mType{:}
        case 'script'
            mdDoc = Dr.documentbasic( ) ; 
        case 'function'
            mdDoc = Dr.documentfunction( ) ;
        case 'classdef'
            mdDoc = Dr.documentclassdef( ) ;
    end

    if Dr.syntax == "mkd"
        mdDoc = Dr.markuptodown( mdDoc ) ;
    end

end
% =========================================================================    

end
% =========================================================================    
% =========================================================================    
methods
    %.....
    [ docFile ] = printdoc( Dr, iM )
end
% =========================================================================    
% =========================================================================    
methods( Access=private )
    %..... 
    [docStr] = documentbasic( Dr, Att, headingLevel )
    %.....
    [docStr] = documentclassdef( Dr )
    %.....
    [docStr] = documentclassmethods( Dr )
    %.....
    [docStr] = documentclassproperties( Dr )
    %.....
    [ tableStr ] = tableattributes( Dr, Attributes )
end
% =========================================================================
% =========================================================================    
methods( Static )
    %.....
    [mFiles]     = findfiles( pathIn, isSearchRecursive )
    %.....
    [ mdDocStr ] = markuptodown( muDocStr )
end
% =========================================================================    
% =========================================================================    

end
