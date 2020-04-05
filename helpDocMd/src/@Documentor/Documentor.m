classdef Documentor < handle
% DOCUMENTOR _.m_bedded `doc`s to .md text files
%
% Documentor (re)publishes embedded source code documentation to text files
% (i.e. as [Markdown](https://daringfireball.net/projects/markdown/)).
% 
% **Refer to the README document for basic usage.**
% 
% __CONSTRUCTOR SYNTAX__
%
%     Dr = Documentor(  ) ;
%     Dr = Documentor( src ) ;
%     Dr = Documentor( src, Params ) ;
% 
% Typically a Documentor instance `Dr` will be created by calling the
% constructor with at least the first argument `src`: a string vector of file
% paths to .m source files and/or directories in which to search for them.
% Upon initialization, documentation can be published using the default
% settings by calling `Dr.printdoc();`
%
% __OPTIONS__
% `Params.isSearchRecursive` is a scalar logical specifying whether any subdirectories
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
% -[test](https://daringfireball.net/projects/markdown/dingus) how a
% sample of Markdown text will display once reformatted to HTML


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
    % For details regarding the implementation and what constitutes a
    % 'documentable' file, 
    %
    % See also 
    % Documentor.findfiles 
    mFiles 
    % mFiles {mustBeFileOrFolder} = string( [ mfilename('fullpath') '.m' ] ) ;

    % List of output documentation file paths (string vector of full file paths)
    dFiles {mustBeStringOrCharOrCellstr} = "" ; 

    % Toggle whethers to overwrite existing documentation files (logical column vector with length == numel(mFiles))
    isOverwriting(:,1) {mustBeBoolean} = false ;
    
    % Toggle whether subdirectories are included when searching for .m files
    % 
    % Applies when calling `Documentor.mFiles = src` and `src` contains directory paths
    % 
    % See also
    % -Documentor.findfiles
    % -Documentor.mFiles
    isSearchRecursive(1,1) {mustBeBoolean} = true ;
   
    % Toggle to recreate original directory tree in `dirOutTop` (multiple mFiles case only) 
    %
    % By default, `isSaveRecursive == true` and an attempt is made to mirror
    % the original directory structure of `.mFiles` within `dirOutTop`.
    % Otherwise, `if isSaveRecursive == false`, all the output documentation files
    % are written to `dirOutTop`.
    %
    % __ETC__  
    % See also 
    % HelpDocMd.isSearchRecursive
    isSaveRecursive(1,1) {mustBeBoolean} = true ;
    
    % Output parent directory for the doc files
    % 
    % If ` "/Users/Buddy/Projects/helpDocMd/src/" ` is the top shared directory
    % of the input `.mFiles`, then, by default
    % `dirOutTop = "/Users/Buddy/Projects/helpDocMd/docs/" ;`
    % 
    % See also Documentor.isSaveRecursive
    dirOutTop {mustBeStringOrChar} = "" ;     
    
    % Documentation string for `mFiles(iM)` to be printed to `dFiles(iM)`
    mdDoc string {mustBeStringOrChar} = "" ;
  
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

properties(  Dependent )
% Access=private,
    
    % Top directory of src mFiles
    dirInTop {mustBeStringOrChar} = "" ;

end
    
properties( Hidden )

    % String specifier for output syntax: "mkd" (for Mkdocs markdown), "mat" (for MATLAB markup)
    %
    % The sole difference between "mkd" and "mat" (for now) is that "mkd" will
    % reformat the style of any embedded links in the comments. 
    %
    % LIKELY NOT NEEDED 
    syntax(1,1) string { mustBeMember( syntax, ["mat" "mkd"] ) } = "mkd" ;
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Dr = Documentor( pathIn, Params )

    if nargin == 0
        return ;
    elseif nargin < 2
        Dr.isSearchRecursive = true ;
    elseif nargin == 2
        
        % assign output directory if not empty
        if isfield( Params, 'outputDir' )
            Dr.dirOutTop = Params.outputDir ;
        end
        
        % assign default if empty
        if ~isfield( Params, 'isSearchRecursive' )
            Dr.isSearchRecursive = Params.isSearchRecursive ;
        end
    end
  
    Dr.mFiles = pathIn ;
    Dr.Info   = Informer( Dr.mFiles(1) ) ;
    Dr.dFiles = "_DEFAULTS_" ;
    
end
% =========================================================================    
function [dirInTop] = get.dirInTop( Dr )

    [~, dirInTop] = parent( Dr.mFiles ) ;

end
% =========================================================================    
function [] = set.dirOutTop( Dr, dirOutTop )

    if isfolder( dirOutTop ) 
        [isValidPath, Values, msgId] = fileattrib( dirOutTop ) ;
    
        assert( Values.directory, 'Invalid assignment: Input is not a directory') ;
        assert( Values.UserWrite, [ 'Invalid permissions: Cannot write to assigned directory.\n', ...
           'Change permissions to the given directory with \n fileattrib( ' dirOutTop ',  +w ) \n' ...
            'or choose a different directory' ], '%s' ) ;
    else
        [isMade, msg, msgId] = mkdir( dirOutTop ) ;
        assert( isMade, msgId, ['Directory creation failed: ' msg ] )
    end
    
    Dr.dirOutTop = string( dirOutTop ) ;
    % reset Dr.dFiles to default (recreate subdirectory tree within Dr.dirOutTop)
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
        warning( 'Nothing assigned: No compatible files were found using the given parameters.') ;
    else
        Dr.mFiles = mFiles ;
    end
    
    % output docs folder in parent directory of ../mFiles
    Dr.dirOutTop = strcat( Dr.dirInTop, filesep, "docs" ) ;

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

    % if Dr.syntax == "mkd" % LIKELY NOT NEEDED...
    %     mdDoc = Dr.markuptodown( mdDoc ) ;
    % end

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
    [ mdDocStr ] = markuptodown( muDocStr ) % LIKELY NOT NEEDED...
end
% =========================================================================    
% =========================================================================    

end
