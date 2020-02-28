classdef Documentor < handle
%% Documentor Custom MATLAB documentation into markup/down text files
%
% Writes *thorough* Matlab documentation as simple, readable Markdown text
% <https://daringfireball.net/projects/markdown/> which is
% 
% 1. readily hosted online (e.g. <https://www.mkdocs.org>,
% <https://pages.github.com/>, <https://docs.readthedocs.io/en/stable/>)
% 
% 2. does not require additional dependencies or different syntax/tagging 
% from Matlab's own markup style (e.g. sphinx)
%
% ### Syntax ###
% 
% TODO
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

    % Informer object instance: Provides the information content to document a
    % given .m file
    Info Informer = Informer( which( "Documentor.m" ) ) ;

    % Input directory  
    dirIn {mustBeFolder} = string( fileparts( mfilename('fullpath') ) ) ;
    
    % Input/Matlab file extension 
    extIn {mustBeStringOrChar} = Documentor.Ext.matlab ;
    
    % Input names (without directory path or file extension)
    nameIn {mustBeStringOrChar} = string( mfilename ) ;

    % Reformated documentation
    mdDoc string {mustBeStringOrChar} = "" ;
   
    % String specifier for output syntax: "mkd" (for Mkdocs markdown), "mat" (for MATLAB markup)
    syntax(1,1) string { mustBeMember( syntax, ["mat" "mkd"] ) } = "mkd" ;

    % Output parent directory for the doc files
    % 
    % See also HelpDocMd.isSaveRecursive
    dirOutTop {mustBeStringOrChar} = "" ;     
    
    % Output file extension (default = ".md")
    extOut(1,1) {mustBeStringOrChar} = Documentor.Ext.markdown ;
    
    % Toggle whether to overwrite existing documentation files
    isOverwriting(1,1) {mustBeNumericOrLogical} = false ;

    % Recreates original directory tree in dirOut (multiple doc output case only) 
    %
    % See also HelpDocMd.isSeachRecursive
    isSaveRecursive(1,1) {mustBeNumericOrLogical} = true ;

    % Toggle whether subdirectories are included in file search (multiple input case only)
    isSearchRecursive(1,1) {mustBeNumericOrLogical} = true ;
    
end

properties( Dependent )

    % Full input file path(s)
    pathIn ;    

end

properties( Access=private, Dependent )
    
    % top directory
    dirInTop {mustBeStringOrChar} = "" ;

    % Help returned from m file
    mHelp ;

end

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function DrSlf = Documentor( pathIn, Params )

if nargin == 0
    return ;
end

%% check input path
if all( ischar( pathIn ) ) || iscellstr( pathIn )
    pathIn = string( pathIn ) ;

elseif ~isstring( pathIn ) || ndims( pathIn > 2 )
   error( ['First input argument must be a file or directory path string ' ...
            ' or a 1-D array of file path strings']  ) ;
end

% convert to fullpath in case the relative path was provided
pathIn = abspath( strtrim( string( pathIn ) ) ) ;  

if size( pathIn, 2 ) > size( pathIn, 1 )
    pathIn = pathIn' ;
end

pathType = unique( Documentor.quiddity( pathIn ) ) ;

assert( numel( pathType ) == 1, [ 'Multiple paths should all be of the same general type,' ...
      ' such that each returns a value of 2 or 7 when input to exist()'] ) ;

switch pathType 
    case 2 % path to matlab files
        DrSlf.pathIn = pathIn ;
    case 7 % directory path
        DrSlf.dirIn = pathIn ;
    otherwise
        error( ['Invalid input path: \n '...
                'Path string must return value of 2 or 7 when input to exist().' ...
                'Check file or directory exists and update the current search path']  ) ;
end
    
end
% =========================================================================    
function [dirInTop] = get.dirInTop( DrSlf )
    
    [~,iTopDir] = min( count( DrSlf.dirIn, filesep ) ) ;
    dirInTop    = fileparts( DrSlf.dirIn( iTopDir ) ) ;

end
% =========================================================================    
function [dirOutTop] = get.dirOutTop( DrSlf )

    if strcmp( DrSlf.dirOutTop, "" )
        try 
            DrSlf.dirOutTop = strcat( DrSlf.dirInTop, filesep, "doc" ) ;
        catch ME
            warning( [ 'Failed to create default dirOutTop directory.\n' ... 
                      'Assign the doc output directory manually.' ], '%s' ) ;
           rethrow(ME) ;
        end 
    end
    
    dirOutTop = DrSlf.dirOutTop ;
           
end
% =========================================================================    
function [] = set.dirOutTop( DrSlf, dirOutTop )

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
    
    DrSlf.dirOutTop = string( dirOutTop ) ;
 
end
% =========================================================================    
function [] = set.pathIn( DrSlf, pathIn )  

    pathType = unique( Documentor.quiddity( pathIn ) ) ;

    assert( all( pathType == 2 ), ['Invalid input path: \n '...
                'Each path string must return value of 2 or 7 when input to exist().\n' ...
                'Check file or directory exists and update the MATLAB search path'], '%s' ) ;
    
    [ DrSlf.dirIn, DrSlf.nameIn, DrSlf.extIn ] = arrayfun( @fileparts, pathIn ) ;
    
    % if unassigned, set .dirOutTop
    if strcmp( DrSlf.dirOutTop, "" )
       DrSlf.dirOutTop = DrSlf.dirInTop ;
    end 

end
% =========================================================================    
function [pathIn] = get.pathIn( DrSlf )  
  
    pathIn = strcat( DrSlf.dirIn, filesep, DrSlf.nameIn, DrSlf.extIn ) ;

end
% =========================================================================    
function [mdDoc] = get.mdDoc( DrSlf )  
    
    % change dir in case list of files to document includes naming conflicts
    dir0 = pwd ; 
    cd( DrSlf.dirIn ) ;
        
    mdDoc = DrSlf.mHelp ;
    
    % names = strcat( "_", join(names, ", "), "_" ) ; % italicize
    % info  = [ info ; strcat( "- Parent classes: ", names ) ] ;
    %
    % if isempty( Mc.InferiorClasses ) 
    %     names = "(None)" ;
    % else 
    %     % NOTE: 'Mc.InferiorClasses' is a cell array whereas Mc.SuperclassList
    %     % is an object array, hence the for-loop: 
    %     names = string( Mc.InferiorClasses{1}.Name ) ;
    %     for iClass = 2 : numel(Mc.InferiorClasses)
    %         names = [ names string( Mc.InferiorClasses{ iClass }.Name ) ] ;
    %     end
    % end
    %
    % names = strcat( "_", join(names, ", "), "_" ) ; % italicize
    % info  = [ info ; strcat( "- Child classes: ", names ) ] ;
    %
    % if isempty( Mc.ContainingPackage ) 
    %     pkg = "N/A" ;
    % else
    %     pkg = string( Mc.ContainingPackage.Name ) ;
    % end
    %
    % info = [ info ; strcat( "- Containing Package: ", "_", pkg, "_" ) ] ;

    % if strcmp( mfiletype( DrSlf.pathIn ), 'classdef' )
    %
    % mdDoc = DrSlf.mHelp ;
    %     %% Class members: Properties    
    %     mdDoc = [ mdDoc; "" ; "### Members ###" ; ""]; 
    %
    %     Props = Mc.PropertyList
    %
    %     mdDoc = [ mdDoc; "#### Properties ####" ] ;
    %
    %     if isempty( Props )
    %
    %         mdDoc = [ mdDoc; "" ; "_(None)_" ; "" ] ;
    %     else
    %
    %         for iProp = 1 : length( Props )
    %
    %             Prop   = Props( iProp )
    %             mdDoc  = [ mdDoc ; strcat( "#####", Prop.Name, "#####" ) ; "" ] ;
    %             fields = fieldnames( Prop ) ;
    %
    %             for iField = 1 : length( fields )   
    %                 field = fields{ iField } ;
    %
    %                 switch field
    %                     case { 'GetAccess', 'SetAccess', 'Dependent', 'Constant', 'Abstract', 'Transient', 'Hidden',
    %                            'GetObservable', 'SetObservable', 'AbortSet', 'NonCopyable', 'GetMethod', 'SetMethod', 'HasDefault' } 
    %                         if isempty( Prop.(field) )
    %                             entry = "_(None)_" ;
    %                         else
    %                             entry = join( string( Prop.( field ) ) ) ;
    %                         end
    %
    %                         mdDoc = [ mdDoc ; strcat( "- ", fields{iField}, ": ", entry ) ] ;
    %
    %                     otherwise
    %                         % do nothing
    %                 end
    %             end
    %     
    %             % TODO : Printing defaults: what is suitable to print (e.g. if
    %             % default is rand(100,100,100) clearly it would be preferable
    %             % to print the function call itself rather than the value).
    %             % might need to parse the code text itself. 
    %            % if Prop.HasDefault
    %            %        mdDoc = [ mdDoc ; "- Default: " ] ;
    %            %        isPrintingDefault=false ;
    %            %     try
    %            %        defaultString = splitlines( string( Prop.Default ) ) ;
    %            %        mdDoc = [mdDoc ; defaultString ] ; 
    %            %      catch Me
    %            %          defaultString = "_(Not printable)_" ;
    %            %    end
    %            %
    %            %  else
    %            %      mdDoc = [ mdDoc ; "- Default: _(None)_" ] ;
    %            % end
    %
    %
    %             end
    %
    %         end 
    %     end
    %
    %
    %
    %
    % end   
    %

    
    if ~strcmp( DrSlf.syntax, "mkd" ) 
        cd( dir0 ) ;
        return ;
    end

    %% Mkdocs Markdown
    %
    % Reformat links and link-text to Markdown syntax:
    % i.e. MATLAB markup uses: <https://www.thisSite.com text to display>
    % whereas Markdown uses: [text to display](https://www.thisSite.com)
    %
    % NOTE the following works for weblinks, but will need to be elaborated
    % for local links to custom functions & classes: either tags or relative
    % paths could be used for Mkdocs build.
    
    links = extractBetween( mdDoc, "<",">" ) ;

    for iLink = 1 : numel(links)
    
        substrings = split( links(iLink) ) ;

        linkUrl = substrings{1} ;
        
        if numel(substrings) == 1 %URL-only
            linkText=substrings{1} ;
        else
            linkText = strcat( substrings{2:end} ) ;
        end
    
        mdDoc = replace( mdDoc, strcat("<", links{iLink}, ">"), ...
                    strcat( "[", linkText, "](", linkUrl,  ")" ) ) ;
    end
    
    mdDoc = strip( splitlines( mdDoc ) ) ;
    
    % trim any empty terminating lines 
     while( strcmp( mdDoc(end), "" ) )
         mdDoc(end) = [] ;
     end

    % % TODO
    % % Replace links (and link-text) properly.
    % % ie. replace HTML in the 'help'/mHelp output with markup (i.e. Matlab or Markdown)
    % % e.g. if <https://www.neuro.polymtl.ca/home/ NeuroPoly> is placed in the header comment,
    % % MATLAB help with then convert it to a weird quasi-HTML when displayed in the prompt:
    % %
    % % <<a href="matlab:web https://www.neuro.polymtl.ca/home/">https://www.neuro.polymtl.ca/home/</a> NeuroPoly> ;
    % %
    % % which, to get into Markdown, needs to be replaced with:
    % %
    % % [NeuroPoly](https://www.neuro.polymtl.ca/home/)
    %
    % % i _think_ (hope) this pattern works to extract the text link (i.e.
    % % "NeuroPoly" in the above e.g.)
    % textIncluded = '(?<=</a>)(\s|\w)*(?=>)' ;
    % % textDefault  =
    %
    % % component-wise explanation:
    % % (?<=</a>) : find strings following (i.e. not including) the pattern "</a>" 
    % % (\s|\w)* : consisting of any number (*) of words (\w) or (|) spaces (\s)
    % % (?=>) : preceding (but not including) ">" 
    % linkText = regexp( mdDoc, textPattern, 'match' ) ; 
        
end
% =========================================================================    
function [mHelp] = get.mHelp( DrSlf )  
   
    % change dir in case list of files to document includes naming conflicts
    dir0 = pwd ; 
    cd( DrSlf.dirIn ) ;
    
    mHelp = help( DrSlf.nameIn ) ;   
    mHelp = string( mHelp ) ;
    
    cd( dir0 ) ;    

end
% =========================================================================    
function [membersDoc] = docclassmembers( DrSlf )
%DOCCLASSMEMBERS

membersDoc = [] ;
    % Mc = metaclass( MrdiIo ) ;
    %
    % for iProp = 1 : length( Mc.PropertyList ) 
    %     if ~Mc.PropertyList(iProp).Constant && Mc.PropertyList(iProp).HasDefault 
    %         DEFAULTS.( Mc.PropertyList(iProp).Name ) = Mc.PropertyList(iProp).DefaultValue ;
    %     end
    % end
end
% =========================================================================    
function [pathOut] = write( DrSlf )  
%WRITE Write doc contents to file
% 
% [pathOut] = write( Self ) 

    pathOut = [DrSlf.dirOutTop+filesep+DrSlf.nameIn+DrSlf.extOut] ;
    
    assert( DrSlf.isOverwriting || ~exist( pathOut ), ...
        ['Doc file already exists. Assign a different file path for the output,' ...
         'or set ' inputname(1) '.isOverwriting =true to force overwrite'], '%s' );
    
    [fid, errMsg] = fopen( pathOut, 'w+') ;
    assert( fid~=-1, ['Write failed:' errMsg] ) ;

    fprintf( strcat("Writing doc: ", pathOut, "\n") ) ;
    fprintf( fid, '%s\n', DrSlf.mdDoc);

    fclose(fid);

end
% =========================================================================    
function [List] = findmfiles( DrSlf )
%FINDMFILES Returns list of Matlab files from dir search
% 
% Wraps to findfiles and to return a 1-D struct array of Matlab files
% 
% Usage 
%
% List = FINDMFILES( Self )
%
% See also dir, findfiles

%% Check inputs
narginchk(0,1) ;

List = findfiles( DrSlf.dirIn, "*.m", DrSlf.isSearchRecursive ) ;

end
% =========================================================================    

end
% =========================================================================    
% =========================================================================    

% =========================================================================
% =========================================================================    
methods( Static )
% =========================================================================    
function [info] = documentclassattributes( Mc )
%DOCUMENTCLASSATTRIBUTES        
    arguments
        Mc(1,1) meta.class ;
    end

    info = [ "" ; "#### Attributes ####" ; "" ] ; 
    
    fields = fieldnames( Mc ) ;

    %% Basic attributes:
    for iField = 1 : length( fields )
        if islogical( Mc.(fields{iField}) ) 
            % Applies to: "Hidden", "Sealed", "Abstract", "Enumeration",
            % "ConstructOnLoad", "HandleCompatible", "RestrictsSubclassing"
            value = join( string( Mc.( fields{iField} ) ), ", " ) ;
            info  = [ info ; strcat( "- ", fields{iField}, ": ", value ) ] ;
        end
    end

    %% Inheritances and package info:
    %
    %TODO when possible: add links to packages, parents, subclasses

    if isempty( Mc.SuperclassList ) 
        names = "(None)" ;
    else
        names = string( {Mc.SuperclassList.Name} ) ;
    end

    names = strcat( "_", join(names, ", "), "_" ) ; % italicize
    info  = [ info ; strcat( "- Parent classes: ", names ) ] ;

    if isempty( Mc.InferiorClasses ) 
        names = "(None)" ;
    else 
        % NOTE: 'Mc.InferiorClasses' is a cell array whereas Mc.SuperclassList
        % is an object array, hence the for-loop: 
        names = string( Mc.InferiorClasses{1}.Name ) ;
        for iClass = 2 : numel(Mc.InferiorClasses)
            names = [ names string( Mc.InferiorClasses{ iClass }.Name ) ] ;
        end
    end
    
    names = strcat( "_", join(names, ", "), "_" ) ; % italicize
    info  = [ info ; strcat( "- Child classes: ", names ) ] ;

    if isempty( Mc.ContainingPackage ) 
        pkg = "N/A" ;
    else
        pkg = string( Mc.ContainingPackage.Name ) ;
    end

    info = [ info ; strcat( "- Containing Package: ", "_", pkg, "_" ) ] ;

end
% =========================================================================    
function [thingTypes] = quiddity( paths )  
%QUIDDITY Check existence of multiple input paths (wraps to exist for each element)
%
% Given an input array *paths* (consisting solely of strings or chars),
% QUIDDITY returns a double array *thingTypes* wherein each element is the
% corresponding return value from exist()
% 
% ### Usage ### 
%
% thingTypes = Documentor.QUIDDITY( paths ) 
    
% XXX stupid method. delete it.
narginchk(1,1) ;
    thingTypes = arrayfun( @exist, string( paths ) ) ;

end
% =========================================================================    

end
% =========================================================================    
% =========================================================================    
end
    %     %% Class members: Properties    
    %     mdDoc = [ mdDoc; "" ; "### Members ###" ; ""]; 
    %
    %     Props = Mc.PropertyList
    %
    %     mdDoc = [ mdDoc; "#### Properties ####" ] ;
    %
    %     if isempty( Props )
    %
    %         mdDoc = [ mdDoc; "" ; "_(None)_" ; "" ] ;
    %     else
    %
    %         for iProp = 1 : length( Props )
    %
    %             Prop   = Props( iProp )
    %             mdDoc  = [ mdDoc ; strcat( "#####", Prop.Name, "#####" ) ; "" ] ;
    %             fields = fieldnames( Prop ) ;
    %
    %             for iField = 1 : length( fields )   
    %                 field = fields{ iField } ;
    %
    %                 switch field
    %                     case { 'GetAccess', 'SetAccess', 'Dependent', 'Constant', 'Abstract', 'Transient', 'Hidden',
    %                            'GetObservable', 'SetObservable', 'AbortSet', 'NonCopyable', 'GetMethod', 'SetMethod', 'HasDefault' } 
    %                         if isempty( Prop.(field) )
    %                             entry = "_(None)_" ;
    %                         else
    %                             entry = join( string( Prop.( field ) ) ) ;
    %                         end
    %
    %                         mdDoc = [ mdDoc ; strcat( "- ", fields{iField}, ": ", entry ) ] ;
    %
    %                     otherwise
    %                         % do nothing
    %                 end
    %             end
    %     
    %             % TODO : Printing defaults: what is suitable to print (e.g. if
    %             % default is rand(100,100,100) clearly it would be preferable
    %             % to print the function call itself rather than the value).
    %             % might need to parse the code text itself. 
    %            % if Prop.HasDefault
    %            %        mdDoc = [ mdDoc ; "- Default: " ] ;
    %            %        isPrintingDefault=false ;
    %            %     try
    %            %        defaultString = splitlines( string( Prop.Default ) ) ;
    %            %        mdDoc = [mdDoc ; defaultString ] ; 
    %            %      catch Me
    %            %          defaultString = "_(Not printable)_" ;
    %            %    end
    %            %
    %            %  else
    %            %      mdDoc = [ mdDoc ; "- Default: _(None)_" ] ;
    %            % end
    %
    %
    %             end
    %
    %         end 
    %     end
