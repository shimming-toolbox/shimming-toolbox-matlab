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
% ### Basic Usage ###
%
% 1. User creates a Documentor instance with the list of .m file paths to
% document:
%
% Dr = Documentor( mFiles ) ;
%
% 2. To create the .md documentation, the user calls: 
% 
% Dr.write ; 
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
    Info Informer ; % = Informer( which( "Documentor.m" ) ) ;

    % List of .m files to document (string scalar or vector of full file paths)
    mFiles {mustBeFile} = string( [ mfilename('fullpath') '.m' ] ) ;

    % Index of next .m file in mFiles list to document
    iM(1,1) uint64 {mustBePositive, mustBeInteger} = uint64(1) ;

    % Toggle whether to overwrite existing documentation files
    isOverwriting(1,1) {mustBeNumericOrLogical} = false ;
    
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
    syntax(1,1) string { mustBeMember( syntax, ["mat" "mkd"] ) } = "mkd" ;


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
    Dr.iM = iM ;

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
    
    % [ Dr.dirIn, Dr.nameIn, Dr.extIn ] = arrayfun( @fileparts, mFiles ) ;
    % if unassigned, set .dirOutTop
    % if strcmp( Dr.dirOutTop, "" )
    %    Dr.dirOutTop = Dr.dirInTop ;
    % end 

end
% =========================================================================    
function [mdDoc] = get.mdDoc( Dr )  
    
    % TODO: Choose desired formating!
   
    Info  = Dr.Info.Attributes ;

    switch Info.mType{:}
        case 'script'
            mdDoc = documentbasic( Info ) ; 
        case 'function'
            mdDoc = documentfunction( Info ) ;
        case 'classdef'
            mdDoc = documentclassdef( Info ) ;
    end

    if Dr.syntax == "mkd"
        mdDoc = reformatasmarkdown( mdDoc ) ;
    end
    
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

    % if strcmp( mfiletype( Dr.pathIn ), 'classdef' )
    %
    % mdDoc = Dr.mHelp ;
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

    
    function [mdDoc] = reformatasmarkdown( mdDoc ) 
    %% Mkdocs Markdown
    %
    % Reformat links and link-text to Markdown syntax:
    % i.e. MATLAB markup uses: <https://www.thisSite.com text to display>
    % whereas Markdown uses: [text to display](https://www.thisSite.com)
    %
    % NOTE the following works for weblinks, but will need to be elaborated
    % for local links to custom functions & classes: either tags or relative
    % paths could be used for Mkdocs build.
    try 
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
    catch Me
        warning( 'Link replacement fail: TODO fix' ) ;
    end    

    mdDoc = strip( splitlines( mdDoc ) ) ;
    
    % trim any empty terminating lines 
     while( strcmp( mdDoc(end), "" ) )
         mdDoc(end) = [] ;
     end

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
    
    %% Local functions (could make these Documentor methods instead)
    % 
    function [mdDoc] = documentbasic( Info )
    %% Basic documentation: will apply to all .m file types
        mdDoc    = strings( 6,1 ) ;
        mdDoc(1) = strcat( "#", Info.Name, "#", " (a Matlab ", Info.mType, ")" ) ;
        mdDoc(2) = "" ;
        mdDoc(3) = strcat( "_", Info.Description, "_" ) ;
        mdDoc(4) = "" ;
        mdDoc(5) = "### Description ###"
        mdDoc    = [mdDoc ; Info.DetailedDescription ] ;
    end 
    
    function [mdDoc] = documentfunction( Info )
    %DOCUMENTFUNCTION adds function-specific info to documentation  
    % (NOTE: for now, this is just nArgin/nArgout but this should be elaborated
    % in Informer.m -- e.g. by parsing the function arguments block when it exists)
        mdDoc = documentbasic( Info ) ;
        
        fields = string( fieldnames( Info ) ) ;
    
        % remove fields included in documentbasic 
        fields( fields=="mType" )               = [] ;
        fields( fields=="Name" )                = [] ;
        fields( fields=="Description" )         = [] ;
        fields( fields=="DetailedDescription" ) = [] ;

        mdDoc = [mdDoc ; "" ; "### Attributes ###"] ;
   
        for iField = 1 : numel(fields)
            field = char( fields( iField ) ) ;
            mdDoc(end+1) = "" ;
            mdDoc(end+1) = strcat( "- ", fields(iField), " : ", string( Info.( field ) ) ) ;
        end
    end
    
    function [mdDoc] = documentclassdef( Info )
    %DOCUMENTCLASSDEF adds class-specific info to documentation  
        mdDoc = documentbasic( Info ) ;
        
        fields = string( fieldnames( Info ) ) ;
    
        % remove fields included in documentbasic 
        fields( fields=="mType" )               = [] ;
        fields( fields=="Name" )                = [] ;
        fields( fields=="Description" )         = [] ;
        fields( fields=="DetailedDescription" ) = [] ;

        mdDoc = [mdDoc ; "" ; "### Attributes ###"] ;
        
        for iField = 1 : numel(fields)
            
            field = char( fields( iField ) ) ;
            
            if isempty( Info.(field) )
                mdDoc(end+1) = strcat( "- ", field, " : [N/A] " ) ;
            elseif strcmp( field, 'SuperclassList' )
                mdDoc(end+1) = strcat( "- Superclasses: ", strjoin( Info.SuperclassList, ", " ) ) ;
            elseif ~isstruct( Info.(field) ) % Property + MethodList etc. will be structs (handle them separately)
                mdDoc(end+1) = strcat( "- ", fields(iField), " : ", string( Info.( field ) ) ) ;
            end
        end
        
        mdDoc = [ mdDoc ; documentclassproperties( Info ) ] ;    
        mdDoc = [ mdDoc ; documentclassmethods( Info ) ] ;    
    end
    
    function methDoc = documentclassmethods( Info )
    %DOCUMENTCLASSMETHODS
        methDoc = [ "" ; "### Methods ###" ; "" ] ;

        if isempty( Info.MethodList )
            methDoc(end+1) = "[No Methods]" ;
        else
            methFields = fieldnames( Info.MethodList ) ;
             
            for iMeth = 1 : numel( Info.MethodList ) 
                methDoc = [ methDoc ; "" ] ;

                Meth = Info.MethodList( iMeth ) ; 
                
                if ~isempty( Meth.Description )
                    nameAndDescription = strcat( "####", Meth.Name, "####", " : ", "_", Meth.Description, "_" )
                    methDoc(end+1) = nameAndDescription;
                else % name only
                    methDoc(end+1) = strcat( "####", Meth.Name, "####" );
                end

                if ~isempty( Meth.DetailedDescription )
                    % NOTE: could change Informer so it doesn't fill DetailedDescription with the copy of Description
                    if ~strcmp( Meth.DetailedDescription, Meth.Description )                         
                        methDoc = [ methDoc ; "Description: " ; Meth.DetailedDescription ] ;
                    end
                end

                for iField = 4 : length( methFields ) % start at 4 to skip name, description, detailed description:
                
                    field = string( methFields( iField ) ) ;

                    if isempty( Meth.(field) )
                        methDoc(end+1) = strcat( "- ", field, " : [N/A] " ) ;
                    else 
                        methDoc(end+1) = strcat( "- ", field, " : ", strjoin( string( Meth.( field ) ), ", " ) ) ;
                    end
                end
            end
        end
    end

    function propsDoc = documentclassproperties( Info )
    %DOCUMENTCLASSPROPERTIES
        propsDoc = [ "" ; "### Properties ###" ; "" ] ;

        if isempty( Info.PropertyList )
            propsDoc(end+1) = "[No Properties]" ;
        else
            propFields = fieldnames( Info.PropertyList ) ;
             
            for iProp = 1 : numel( Info.PropertyList ) 
                propsDoc = [ propsDoc ; "" ] ;

                Prop = Info.PropertyList( iProp ) ; 
                
                if ~isempty( Prop.Description )
                    nameAndDescription = strcat( "*", Prop.Name, "*", " : ", "_", Prop.Description, "_" )
                    propsDoc(end+1) = nameAndDescription;
                else % name only
                    propsDoc(end+1) = strcat( "*", Prop.Name, "*" );
                end

                if ~isempty( Prop.DetailedDescription )
                    % NOTE: could change Informer so it doesn't fill DetailedDescription with the copy of Description
                    if ~strcmp( Prop.DetailedDescription, Prop.Description )                         
                        propsDoc = [ propsDoc ; "Description: " ; Prop.DetailedDescription ] ;
                    end
                end
                
                for iField = 4 : length( propFields ) % start at 4 to skip name, description, detailed description:
                    
                    field = string( propFields( iField ) ) ;

                    if isempty( Prop.(field) )
                        propsDoc(end+1) = strcat( "- ", field, " : [N/A] " ) ;
                    
                    elseif strcmp( field, "Validation" )

                        propsDoc(end+1) = "- Validation: " ;

                        if Prop.Validation.Class ~= "" 
                            propsDoc(end+1) = strcat( "Class: ", Prop.Validation.Class ) ;
                        end

                        propsDoc(end+1) = strcat( "Validator functions: ", strjoin( Prop.Validation.ValidatorFunctions, "," ) ) ;
                        
                        %TODO add size constraints
                    else 
                        propsDoc(end+1) = strcat( "- ", field, " : ", string( Prop.( field ) ) ) ;
                    end
                end 
            end
        end
    end



end
% =========================================================================    
function [membersDoc] = docclassmethods( Dr )
%DOCCLASSMETHODS

end
% =========================================================================    
function [pathOut] = write( Dr )  
%WRITE Write doc contents to file
% 
% [pathOut] = write( Self ) 

    pathOut = [Dr.dirOutTop+filesep+Dr.nameIn+Dr.extOut] ;
    
    assert( Dr.isOverwriting || ~exist( pathOut ), ...
        ['Doc file already exists. Assign a different file path for the output,' ...
         'or set ' inputname(1) '.isOverwriting =true to force overwrite'], '%s' );
    
    [fid, errMsg] = fopen( pathOut, 'w+') ;
    assert( fid~=-1, ['Write failed:' errMsg] ) ;

    fprintf( strcat("Writing doc: ", pathOut, "\n") ) ;
    fprintf( fid, '%s\n', Dr.mdDoc);

    fclose(fid);

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

end
% =========================================================================    
% =========================================================================    
methods
    %.....
    [mFiles] = findfilestodocument( Dr, pathIn ) ;
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
