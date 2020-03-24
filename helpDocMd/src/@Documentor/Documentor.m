classdef Documentor < handle
%DOCUMENTOR Custom Matlab documentation into markup/down text files
% 
% The DOCUMENTOR class serves to print custom Matlab documentation to file 
% (e.g. as <https://daringfireball.net/projects/markdown/ Markdown> text)
% while avoiding external dependencies (e.g. <https://github.com/sphinx-contrib/matlabdomain sphinx>)
% and tagging syntaxes at odds with Matlab's own markup 
% <https://www.mathworks.com/help/matlab/matlab_prog/marking-up-matlab-comments-for-publishing.html style>.
%
% Viz., given a list of 
% <https://www.mathworks.com/help/matlab/matlab_prog/add-help-for-your-program.html properly> 
% commented .m files, a DOCUMENTOR instance outputs simple, readable text files
% which are readily hosted online. 
% **See 
% <https://github.com/neuropoly/realtime_shimming/blob/helpDocMd/helpDocMd/doc/Documentor.md Example>.**
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
% #### Construction syntax
%
%     Dr = Documentor(  ) ;
%     Dr = Documentor( src ) ;
%     Dr = Documentor( src, Options ) ;
% 
% When called without arguments, `Documentor()` constructs a default object,
% assigning the respective defaults to each of its properties; public
% properties can still be reconfigured following construction.
% 
% As described previously in the **Basic usage** section, `src` paths can be
% provided as an argument to the constructor to automatically assign and/or
% search for *documentable* .m files. If `src` is the only input argument,
% default values will be assigned to most properties of the returned object.
%
% **Note:** Any invalid or *incompatible* paths suggested by `src` will be
% automatically omitted (refer to the **Documenting .m files** section for
% details).
% 
% When called with a second argument, `Options` --- 
% 
% By default, if `src` paths include directories that contain subfolders, these too
% are recursively searched for .m files. To bypass the default behaviour and
% restrict the search depth to directories explicitly listed in `src`, the
% constructor should be called with the second argument: A parameters struct `Options` a field Options.isSearchRecursive in the second argument position
% without with a `0` in the second argument position, i.e.
%    
%    isSearchRecursive = false ; 
%    Dr = Documentor( src, isSearchRecursive ) ;
%
% #### Documenting .m files
%
% The list of .m files to be documented is assigned to the value of Documentor property `mFiles`.
%
% To ensure that the list consists exclusively of *documentable* files, assignments to `mFiles` are
% implicitly filtered: i.e. when setting `Dr.mFiles = src`, the actual assignment will be 
% `Dr.mFiles = Documentor.findfiles( src, Dr.isSearchRecursive )`.
%
% (For more info, see the documentation entries for Documentor.findfiles and Documentor.mFiles)
%
% #### Configuring options 
% 
% ...TODO
%
% ### References
% 
% To test how a markdown sample will display when reformatted to HTML:
% - <https://daringfireball.net/projects/markdown/dingus>
%
% Re: hosting documentation online:
% - <https://www.mkdocs.org/ MkDocs>, 
% - <https://pages.github.com/ Github>,
% - <https://docs.readthedocs.io/en/stable/ ReadTheDocs>

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
    % #### References
    %
    % For details re: implementation and what constitutes a 'documentable' file,
    % See also 
    % Documentor.findfiles 
    mFiles {mustBeFileOrFolder} = string( [ mfilename('fullpath') '.m' ] ) ;

    % List of output documentation file paths (string vector of full file paths)
    dFiles {mustBeStringOrCharOrCellstr} = "_DEFAULTS_" ; 

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
    % #### Example
    % 
    % Given Documentor object instance `Dr`:
    %
    % `if Dr.isSaveRecursive == false` 
    % ...all the output documentation files will be written to the same folder (i.e. `Dr.dirOutTop`)
    % 
    % `else`
    % ...an attempt is made to mirror the original directory tree of the top .m file source folder
    %
    % #### References
    %
    % See also 
    % - HelpDocMd.isSearchRecursive
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
function Dr = Documentor( pathIn )

    if nargin == 0
        return ;
    elseif nargin == 2
        if ~isempty( Params.dirOutTop )
            Dr.dirOutTop = Params.dirOutTop ;
        end
    end
    
    Dr.mFiles = Documentor.findfiles( pathIn ) ;
    Dr.Info   = Informer( Dr.mFiles(1) ) ;
    % Dr.Info = {} % or [] or zeros(numel(Dr.mFiles))

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
    Dr.dFiles = "_DEFAULTS_" ; % why set again?

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
