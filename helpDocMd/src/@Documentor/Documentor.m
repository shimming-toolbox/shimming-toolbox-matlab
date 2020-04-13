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
%     Dr = Documentor( src ) ;      
%     Dr = Documentor( src, Options ) ; 
%
% `src` is a string vector of file paths to .m source files and/or directories
% in which to search for them. 
%
% The constructor calls the methods `draftcontent` and `autoassigndocfiles`, 
% which respectively set the properties:  
% -`docContent`: the documentation text content, and  
% -`docFile`: the output filepath for the published documentation.
%
% Both of which can be reconfigured prior to publishing to file via `printdoc()`.
%
% __OPTIONS__
%
% The constructor accepts a struct of `Options` parameters, with the following
% supported fields. 
%
% -`Options.detailLevel`: Specifies the degree of detail to include when
% generating documentation via `draftcontent()`.
%
% -`Options.outputDir`: specifies the *primary* folder for the printed
% documentation files, used in the call to `autoassigndocfiles()`.
% 
% For more info, refer to the documentation for the corresponding methods.
%
% __TIPS__
% 1. Verify that the assignments to `docFile` are indeed suitable prior
% to calling `printdoc`.
%
% 2. If `src` will contain paths from all over the file system, it may be
% better to construct independent Documentor arrays as the default assignments
% to `docFile` are more likely to work if the original `src` paths stem from a
% common user directory.
%
% 3. To get a better idea of what files will actually be documented, you may
% want to first initialize `src` as `src = Documentor.findfiles( src )`, as
% is done during object construction. The reason being that `src` is used to
% assign the `mFile` property (which is fixed upon object construction) and any
% invalid or *incompatible* paths suggested by `src` will be automatically
% filtered out. For more info, refer to the documentation entries for
% `Documentor.findfiles` and `Documentor.mustBeDocumentable`.
%
% __ETC__  
% -[test](https://daringfireball.net/projects/markdown/dingus) how a
% sample of Markdown text will display once reformatted to HTML

properties( Constant, Hidden )

    % Default file extensions
    Ext = struct( 'matlab', ".m", 'markdown', ".md" ) ;

    % Supported .m file types (script, function, classdef)
    mTypesSupported = [ "script" "function" "classdef" ] ; 

end

properties( SetAccess=immutable ) % `immutable` might simplify object handling; could change if there is reason to.

    % String vector specifying the full file path of the .m file to document
    %
    % For details re: file 'documentability',
    % See also 
    % Documentor.mustBeDocumentable 
    % Documentor.findfiles
    mFile {Documentor.mustBeDocumentable} = string( [ mfilename('fullpath') '.m' ] ) ;

end

properties( AbortSet )
    
    % String vector specifying the full file path for the printed documentation
    docFile string {mustBeStringScalarOrCharVector} = "" ;

    % Documentation text content as a string vector 
    docContent(:,1) string = "" ;
    
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Dr = Documentor( src, Params )
    
    if nargin == 0
        return ;
    end
    
    DEFAULTS.detailLevel = 1 ; % for doc content autofill 
    DEFAULTS.outputDir   = [] ; % for docFile assignment (print output path)
    
    if nargin == 1
        Params = DEFAULTS ; 
    else
        assert( isstruct(Params), 'Optional second argument should be a struct of parameters')
        Params = assignifempty( Params, DEFAULTS ) ;
    end
    
    % Find/filter/assign .m files to document
    mFiles = Documentor.findfiles( src, true ) ;
    
    % Construct object(s)
    for iM = 1 : numel( mFiles )
        Dr(iM).mFile = mFiles(iM) ; 
    end

    Dr.draftdoccontent( Params.detailLevel ) ;
    Dr.autoassigndocfiles( Params.outputDir ) ;

end
% =========================================================================    
end
% =========================================================================    
% =========================================================================    
methods 
    %..... 
    [] = autoassigndocfiles( Dr, outputDir )
    %..... 
    [] = draftdoc( Dr )
    %..... 
    [ docFiles, errMsg ] = printdoc( Dr, Options )
end
% =========================================================================    
% =========================================================================    
methods( Static, Hidden ) 
    %..... 
    [docStr] = documentbasic( Info, headingLevel )
    % ....
    [docStr] = documentclassdef( Info )
    %.....
    [docStr] = documentclassmethods( Info )
    %.....
    [docStr] = documentclassproperties( Info )
    %.....
    [docStr] = documentfunction( Info )   
    %.....
    [ tableStr ] = tableattributes( Attributes )
    %.....
    [ mdDocStr ] = markuptodown( muDocStr ) % LIKELY DEPRECATED...
end
% =========================================================================
% =========================================================================    
methods( Static )
    %.....
    [mFiles] = findfiles( pathIn, isSearchRecursive )
    %.....
    [] = mustBeDocumentable( mFile )
end
% =========================================================================    
% =========================================================================    

end
