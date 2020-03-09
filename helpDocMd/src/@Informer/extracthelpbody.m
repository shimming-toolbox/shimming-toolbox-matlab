function [mHelpBody] = extracthelpbody( mHelp )  
%EXTRACTHELPBODY Return the body of help-text, trimmed of the leading line
%
% EXTRACTHELPBODY returns a body of help-text (i.e. from Informer.gethelptext) 
% and trims it of its leading line of text. (If the text consists solely of the
% leading line, then this line is returned.)
%
% ### Syntax
%      
%     mHelpBody = EXTRACTHELPBODY( mHelp ) ;
% 
% ### Example
% 
% % To display the current section of text, without the title line:
% mHelpBody = Informer.extracthelpbody( Informer.gethelptext( 'Informer.extracthelpbody' ) )
%
% ### References 
%
% See also
%
% -Informer.gethelptext
% -Informer.extracthelpheader 
    arguments
        mHelp {mustBeStringOrCharOrCellstr} ;
    end
    
if iscellstr( mHelp ) || ischar( mHelp )
   mHelp = string( mHelp ) ;
end 

mHelp = splitlines( mHelp ) ;

%% Remove any leading textless lines:
while( mHelp(1) == "" )
    mHelp(1) = [] ;
end

% Remove first line of text (unless that's all there is) 
if length(mHelp)<2
    mHelpBody = mHelp ;
    return ;
else
    % if a line begins 4 or more empty spaces, assume it is a code block.
    % Otherwise, strip the empty spaces.
    for iLine = 1 : length( mHelp )
        if ~startsWith( mHelp(iLine), "    " )
            mHelp(iLine) = strip( mHelp(iLine) ) ;
        end
    end

    mHelpBody = mHelp(2:end) ;
end

if length(mHelpBody)>2 && (mHelpBody(1) == "")
    % Remove any leading textless lines:
    while( mHelpBody(1) == "" )
        mHelpBody(1) = [] ;
    end
end

end
