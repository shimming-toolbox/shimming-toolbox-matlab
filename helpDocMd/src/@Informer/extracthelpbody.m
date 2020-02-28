function [mHelpBody] = extracthelpbody( mHelp )  
%EXTRACTHELPBODY Return the body of help-text, trimmed of the leading line
%
% EXTRACTHELPHEADER returns a body of help-text (i.e. from Informer.gethelptext) 
% and trims it of its leading line of text. (If the text consists solely of the
% leading line, then this line is returned.)
%
% mHelpBody = EXTRACTHELPHEADER( mHelp )
% 
% ### Example ###
% 
% % To display the current section of text, without the title line:
% mHelpBody = Informer.extracthelpbody( Informer.gethelptext( 'Informer.extracthelpbody' ) )
%
% See also
%
% Informer.gethelptext
% Informer.extracthelpheader 
    arguments
        mHelp {mustBeStringOrCharOrCellstr} ;
    end
    
    if iscellstr( mHelp ) || ischar( mHelp )
       mHelp = string( mHelp ) ;
    end 

    % Avoid any leading textless lines 
    mHelp     = strip( splitlines( mHelp ) ) ;
    mHelpBody = mHelp(mHelp ~= "") ;
    % Remove first line of text (unless that's all there is) 
    if length(mHelp)>1
        mHelpBody = mHelp(2:end) ;
    end

end

