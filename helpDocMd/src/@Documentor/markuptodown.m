function [ mdDocStr ] = markuptodown( muDocStr ) 
%MARKUPTODOWN Replace MATLAB Markup elements with corresponding Mkdocs-style Markdown
%
% Reformat links and link-text to Markdown syntax:
% i.e. MATLAB markup uses: <https://www.thisSite.com text to display>
% whereas Markdown uses: [text to display](https://www.thisSite.com)
%
% [ mdDocStr ] = MARKUPTODOWN( muDocStr )
% 
% TODO
%
% + currrent implementation is simplistic: Basically works for weblinks, but needs
% to be elaborated for local links to custom functions & classes: either tags
% or relative paths could be used for Mkdocs build.
%
% + it doesn't distinguish between links and embedded HTML, so it
% messes up the latter (see: Documentor.tableattributes)
% 
% + replace instances of 'MATLAB' and/or 'MATLAB(R)' with 'MATLAB&reg;'

assert( isstring( muDocStr ), 'Input argument must be a string' ) ;

mdDocStr = muDocStr ;

for iLine = 1 : length( mdDocStr )
   docLine = mdDocStr(iLine) ; 
    try 
        links = extractBetween( docLine, "<", ">" ) ;

        for iLink = 1 : numel(links)
        
            substrings = split( links(iLink) ) ;

            linkUrl = substrings{1} ;
            
            if numel(substrings) == 1 %URL-only
                linkText=substrings{1} ;
            else
                linkText = strcat( substrings{2:end} ) ;
            end
        
            mdDocStr(iLine) = replace( docLine, strcat("<", links{iLink}, ">"), ...
                        strcat( "[", linkText, "](", linkUrl,  ")" ) ) ;
            
        end

    catch Me
        warning( 'Link replacement fail: TODO fix' ) ;

    end    
end
% mdDocStr = strip( splitlines( mdDocStr ) ) ;

% trim any empty terminating lines 
while( strcmp( mdDocStr(end), "" ) )
    mdDocStr(end) = [] ;
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
