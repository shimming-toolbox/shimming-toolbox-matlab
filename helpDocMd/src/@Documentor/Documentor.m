classdef Documentor < handle
%DOCUMENTOR Custom MATLAB documentation into markup/down text files
%
% Writes *thorough* Matlab documentation as simple, readable 
% <https://daringfireball.net/projects/markdown/ Markdown> text which is
% 
% 1. Readily hosted online (e.g. <https://www.mkdocs.org/ MkDocs>,
% <https://pages.github.com/>, <https://docs.readthedocs.io/en/stable/>)
% 
% 2. Does not require additional dependencies or different syntax/tagging 
% from Matlab's own markup style (e.g. sphinx)
%
% ### Basic Usage ###
%
% 1. User creates a Documentor instance with the list of .m file paths to
% document:
%      
%     Dr = Documentor( mFiles ) ;
%
% 2. To create the .md documentation, the user calls: 
%      
%     Dr.write ; 
%
% ### Example Output ###
% 
% (TODO: add output to github page or readthedocs)
%
% To see the final output online, see <https://ADD_URL.com this> % 
% 
% ### References ###
% 
% To test how the .md output will appear once reformatted to HTML:
% <https://daringfireball.net/projects/markdown/dingus>

properties( Constant )

    % Default file extensions
    Ext = struct( 'matlab', ".m", 'markdown', ".md" ) ;

end

properties( AbortSet )

    % Informer instance: Provides the info-content to document a given .m file
    Info Informer ; % = Informer( which( "Documentor.m" ) ) ;

    % List of .m files to document (string scalar or vector of full file paths)
    mFiles {mustBeFile} = string( [ mfilename('fullpath') '.m' ] ) ;

    % Index of next .m file in mFiles list to document
    iM(1,1) uint64 {mustBePositive, mustBeInteger} = uint64(1) ;

    % Toggle whether to overwrite existing documentation files
    isOverwriting(1,1) {mustBeNumericOrLogical} = true ;
    
    % Toggle whether subdirectories are included in file search (multiple input case only)
    isSearchRecursive(1,1) {mustBeNumericOrLogical} = true ;
   
    % Recreates original directory tree in dirOut (multiple doc output case only) 
    %
    % See also HelpDocMd.isSeachRecursive
    % TODO use mapdirectorytree.m or something to figure out the subdirectory structure to use and change the default to TRUE.
    % (for now, just dumping all documentation into single folder - dirOutTop
    isSaveRecursive(1,1) {mustBeNumericOrLogical} = false ;
    
    % Output parent directory for the doc files
    % 
    % See also HelpDocMd.isSaveRecursive
    dirOutTop {mustBeStringOrChar} = "" ;     
    
    % Reformated documentation
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
    isDetailed(1,1) {mustBeBoolean} = true ;

    % Input/Matlab file extension 
    extIn {mustBeStringOrChar} = Documentor.Ext.matlab ;

    % Input names (without directory path or file extension)
    nameIn {mustBeStringOrChar} = string( mfilename ) ;
    
    % Output file extension (default = ".md")
    extOut(1,1) {mustBeStringOrChar} = Documentor.Ext.markdown ;
    
end

properties( Access=private, Dependent )

    % parent folder of mFiles(iM)  
    mDir {mustBeFolder} = string( fileparts( mfilename('fullpath') ) ) ;
    
    % top directory of src mFiles
    dirInTop {mustBeStringOrChar} = "" ;

end

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Dr = Documentor( pathIn, Params )

    if nargin == 0
        return ;
    end

    Dr.mFiles = Dr.findfilestodocument( pathIn ) ;
    Dr.Info   = Informer( Dr.mFiles(1) ) ;

end
% =========================================================================    
function [dirInTop] = get.dirInTop( Dr )
 
    if numel( Dr.mFiles ) == 1
       dirInTop = fileparts( Dr.mFiles ) ;

    else % find folder with fewest parent directories (i.e. slashes in path): 
        mDirs = arrayfun( @fileparts, Dr.mFiles )

        [~,iTopDir] = min( count( mDirs, filesep ) ) ;
        dirInTop    = fileparts( mDirs( iTopDir ) ) ;
    end

end
% =========================================================================    
function [dirOutTop] = get.dirOutTop( Dr )

    if strcmp( Dr.dirOutTop, "" )
        try  % output doc folder in parent directory of ./mFiles
            Dr.dirOutTop = strcat( fileparts( Dr.dirInTop ), filesep, "doc" ) ;
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
end
% =========================================================================    
function [] = set.iM( Dr, iM )  

    % update Informer with new file:
    Dr.Info.mPath = Dr.mFiles( iM ) ; 
    Dr.iM         = iM ;

end
% =========================================================================    
function [] = set.mFiles( Dr, mFiles )  

    mType = Dr.Info.mfiletype( mFiles ) ; 

    if any( mType == "NA" )
        badPaths = mFiles( mType == "NA" ) ;
        warning('Excluding the following invalid (possibly unimplemented) .m files from documentation: ')
        display( badPaths ) ;
        mFiles( mType == "NA" ) = [] ;
    end
    
    if isempty( mFiles )
        error( "List of .m files suitable for documentation is empty") ;
    end

    Dr.mFiles = mFiles ;
    
    [ Dr.dirIn, Dr.nameIn, Dr.extIn ] = arrayfun( @fileparts, mFiles ) ;

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
    [mFiles]    = findfilestodocument( Dr, pathIn ) ;
    %.....
    [ pathOut ] = write( Dr )
end
% =========================================================================    
% =========================================================================    
methods( Access=private )
    %..... 
    [docStr] = documentbasic( Dr )
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
    [ mdDocStr ] = markuptodown( muDocStr )
end
% =========================================================================    
% =========================================================================    

end
