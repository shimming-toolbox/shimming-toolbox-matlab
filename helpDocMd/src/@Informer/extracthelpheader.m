function [mHelpHeader] = extracthelpheader( mHelp, name )  
%EXTRACTHELPHEADER Return the leading line of help-text
%    
%     mHelpHeader = EXTRACTHELPHEADER( mHelp )
%     mHelpHeader = EXTRACTHELPHEADER( mHelp, name )
% 
% If 'name' is provided as a second input argument and it appears as the first
% word of the header line (irrespective of case) it will be removed from the
% returned string.
%
% __EXAMPLE__
%
%    mHelp = Informer.gethelptext( 'Informer.gethelptext' ) 
% 
% % diplays:
% %
% % "GETHELPTEXT Return help-text of script|function|method|property as string-vector"
% %  ""
% % "mHelp = GETHELPTEXT( name )"
% 
%    Informer.extracthelpheader( mHelp, "GETHELPTEXT" )
%
% % diplays:
% %
% % "Return help-text of script|function|method|property as string-vector"
%
% __ETC__
%
% See also
% Informer.gethelptext
    arguments
        mHelp {mustBeStringOrCharOrCellstr} ;
        name {mustBeStringOrChar} = "" ;
    end
    
if iscellstr( mHelp ) || ischar( mHelp )
   mHelp = string( mHelp ) ;
end 

% Avoid any leading textless lines 
mHelp = strip( splitlines( mHelp ) ) ;
mHelp = mHelp(mHelp ~= "") ;
mHelpHeader = mHelp(1) ;

% If present remove the name/title (unless it's the only word there)
words = split( mHelpHeader ) ;

if ( numel( words ) > 1 ) && strcmpi( words(1), name ) 
    mHelpHeader = join( words(2:end) ) ;
end

end
