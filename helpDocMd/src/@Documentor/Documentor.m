classdef Documentor < handle
%DOCUMENTOR Custom Matlab documentation into markup/down text files
% 
% The DOCUMENTOR class serves to print custom Matlab documentation to file 
% (e.g. as <https://daringfireball.net/projects/markdown/ Markdown>)
% while avoiding external dependencies (e.g. <https://github.com/sphinx-contrib/matlabdomain sphinx>)
% and tagging syntaxes at odds with Matlab's own style of
% <https://www.mathworks.com/help/matlab/matlab_prog/marking-up-matlab-comments-for-publishing.html markup>
%
% Viz., given a list of 
% <https://www.mathworks.com/help/matlab/matlab_prog/add-help-for-your-program.html properly> 
% commented .m files, a DOCUMENTOR instance outputs simple, readable text files
% which are readily hosted online.
%
% ### <https://github.com/neuropoly/realtime_shimming/blob/helpDocMd/helpDocMd/doc/Documentor.md Example>
% 
% ### Usage
%
% Create a Documentor instance with the list of .m file paths of interest, then
% call `printdoc` to write to file: 
%    
%    Dr = Documentor( mFiles ) ;  
% 
%    Dr.printdoc ; 
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
    [mFiles]    = findfilestodocument( Dr, pathIn )
    %.....
    [ pathOut ] = printdoc( Dr )
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
