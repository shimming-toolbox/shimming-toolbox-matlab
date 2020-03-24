function [txt] = readme( isCopying, S )
%README Create [return|copy|print] README text with title, sections, & index/ToC 
%   
%   []    = readme( [false] ) 
%   [txt] = readme(  ) 
%
% When called without any input or output arguments, README prints the
% character vector of default 'readme' content `txt` to the standard output. 
% To suppress the print display, call the function with a return argument, 
% e.g. using the dummy argument `[~] = readme();`.
%
% When the first input is `[ 1 | true ]`, the `txt` is copied to the clipboard.
% This should be adequate for most cases: Paste to your file and begin editing.


%TODO test with options...
% __OPTIONS__
%
%   [txt] = readme( Sections, ToC )
%   [txt] = readme( Sections, ToC, isCopying )
%   [txt] = readme( Sections, ToC, filename )
%   [txt] = readme( Sections, ToC, filename, isOverwriting )
% 
% __DESCRIPTION__
%
% S is string cell array with columns 1-3 corresponding to:
%
% 1. Section names. The first element is the title.
%
% 2. Subsection depth-levels. All *must* contain non-negative integers.
% The first/title element should be 0 and other sections should be <=6 for
% proper Markdown formatting. Reminder warnings are issued if these conditions
% are not met, however they can be suppressed by issuing these commands before
% calling this function:
%
% `warning('off':'helpDocMd:titleSectionDepth')`
% `warning('off':helpDocMd:subsectionDepth')`
%
% 3. Section bullets (character vectors preceding each section heading 
% [default: '#' repeated n=depth-level times]
    arguments
        isCopying(1,1) {mustBeBoolean}     = false ;
        S(:,:) { mustBeA( S, "string" ) } = DEFAULT_S ;
    end

% Total number of sections incl. title
nS = size(S, 1) ; 

%------------------------------ 
% 2. Subsection depth-levels

% convert to numeric
try
    Sn.depth = arrayfun( @str2num, S(:, 1) ) ;
catch Me
    Me.rethrow ;
end

mustBeInteger(Sn.depth) ;
mustBeNonnegative(Sn.depth) ;

if Sn.depth(1) ~= 0
    warning( ['The first column of input S assigns the subsection depth.' ...
    'S(1,1) value applies to the main title and should contain "0" for proper formating.', 'helpDocMd:titleSectionDepth' ] ) ;

elseif any( ( Sn.depth ) ) > 6
    warning('Subsections depths should be less than 7 for proper Markdown formatting.', 'helpDocMd:subsectionDepth') ;
end

%------------------------------ 
% 3. Section bulets 
if ( size( S, 2 ) < 3 ) || all( strcmp( S(:,3), " " ) )

    S(:, 3) = strings( 1,  size( S, 1 ) ) ;
    for iS = 1 : size( S, 1 )
        S(iS, 3) =  "#" + repmat( '#', [1 str2num( S(iS,1) ) ]  ) ;
    end
end

%------------------------------ 
%% Populate output text

% Title line + empty line
txt = [ S(1,3) + " " + S(1,2) ; " " ] ;

% extract depth=1 subsection names for t.o.c 
txt = [ txt ; maketableofcontents( S( Sn.depth==1, 2 ) ) ; " " ] ;

% Add subsection headers, each followed by an empty line
for iS = 2 : nS
    txt = [ txt ; S(iS,3) + " " + S(iS,2) ; " " ] ;
end

%------------------------------ 
%% Output
txt = sprintf( '%s', compose( txt + '\n' ) ) ;

if isCopying
    clipboard( 'copy', txt ) ;
end

if nargout == 0
    fprintf(txt) ;
    clear txt ; 
    return
end

end %readme()

%---------------------------------------------------------------------------  
function [D_S] =DEFAULT_S( )
% Default content:

% 1. Section depth-levels
D_S(:,1) = string( [ 0 ; 1 ; 1 ; 2 ; 1 ; 2 ;1 ;1 ] ) ;

% 2. Names
D_S(:,2) = [ "My Project" ; "Purpose" ; "Getting started" ; "Dependencies" ; ...
             "Basics" ; "Tips" ; "References" ; "License" ]  ;

end

%---------------------------------------------------------------------------  
function [cTable] = maketableofcontents( S1 )
% Return table of contents string from (depth=1) subsection names
% (input/output are string column-vectors)

assert( size(S1,2) == 1 )

for iS = 1 : length( S1 )
    cTable(iS) = [ "- " + ... % bullet 
        "[" + S1(iS) + "]" + ... % title
        "(#" + lower( replace( S1(iS)," ","" ) ) + ")" ] ; % link
end

cTable = cTable(:) ;

end
     

% % 2. Section depth-levels
% Sx(:,2) = { 1 ; 2 ; 2 ; 3 ; 2 ;3 ;2 ;2 } ;
% % Sx(:,2) = cellstr( string( [S{:,2}] ) ) ; % initialize as characters
%
% if iscellstr( Sx(:,2) ) % If input as characters, convert to numeric
%     Sx(:,2) = mat2cellRows( cellfun( @str2num, Sx(:,2) ) ) ;
% elseif ~all( cellfun( @isnumeric, ( Sx(:,2) ) ) ) 
%     error('Column 2 of *Sections* cell must contain positive integers.' )
% end
%
% %% nSubsections per depth-level
% errMsg = [ 'The first section element refers to the title and its depth must be set to 1.\n'... 
%            'Subsections depth-levels *must* be integers greater than 1 and should be less than 7 for proper Markdown formatting.\n' ] ;
%
% end
%
% if ~iscellstr( Sx(:,3) ) % if input as numeric [1, 2, ...] convert to char
%     Sx(:,3) = cellstr( cellfun( @num2str, Sx(:,3) ) ) ;
% end
%
% % 3. Section bullets
% Sx(:,3) = cell( nSections, 1 ) ;
% for iS = 1 : nSections
%     Sx(iS,3) = cellstr( repmat( '#', [1 Sx{iS,2} ] ) ) ;
% end
%
%
% % % Count 
% % for iS = 1 : size( S, 1 ) 
% %     nSubs(iS) = nnz( S.depth == iS ) ;
% % end
%
% if size(S,2) < 2
%     S(:,2) = [ "1" repmat("2",size(
% if size(S,2) < 3
%     
%     Sections.i = strings( size( Sections.name ) ) ;
%
%     for iS = 1 : numel(Sections.name) 
%         Sections.i(iS) = repmat( "#", [1 Sections.depth(i)] ) ;
%     end
% else
%     assert( isequal( size( Sections.i ), size( Sections.name ) ) )
% end
%
% % Title line
% txt = [ Sections.depth(1) + " " + Sections.name(1) ] ;
%
%
%
% nSections = [1 length( Sections.name ) ;
%
% for iS = 1 : numel( Sections.name ) 
%     Sections.i = repmat( "#", [1 Sections.depth(iS) ] ) ;
% end
%
% txt = [ txt ; Sections.i(1) + " " + Sections.name(1) ] ;
% nSections
% Toc.i = repmat( "-", [1 nSections-1] ) ;
%
% for iS = nSections-1
%
%
% if exist('Toc', 'var') 
%    
%     for iS = 2 : nSections 
%    
%         txt(end+1) = "
% for iS = 1 : nSections 
%
%     txt = [ txt ; Sections.i(iS) + " " + Sections.name(iS) ] ;
%
% S.i     = strings( size( S.hdrs ) ) ;
% S.i(1) = "#"
% s.i(2:end) = "##" ;
% S.hdrs = S.hdrs(:) ;
%
% nS     = numel( S.hdrs ) ;
% S.i    = [ "#" repmat( "##",   ) )
%
% for iHdr = 1 : numel( Hdrs )
%     txt  = [ bullers(1) + " " + Hdrs(1) ; " " ]
%
%
% % arguments
% %         title(1,1)
% %         sections(:,:) 
% %         Specs(1,1) struct                         = getdefaultcontents ; 
% %         isLinked ...
% %         isCopying(1,1) {mustBeBoolean}            = true ;
% %         filename {mustBeStringScalarOrCharVector} = "README.md" ;
% %         isOverwriting(1,1)                        = false ;
% %     end
