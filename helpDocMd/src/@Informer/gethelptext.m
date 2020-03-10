function [mHelp] = gethelptext( name )  
%GETHELPTEXT Return help-text of script|function|method|property|class as string-vector
%
% mHelp = GETHELPTEXT( name )
    arguments
        name {mustBeStringOrChar} = "Informer.gethelptext" ;
    end

mHelp = help( name ) ;
assert( ~isempty(mHelp), 'Nothing found. Verify item of interest is on the MATLAB path.' ) ;

mHelp = splitlines( string( mHelp ) ) ;
     
while( strcmp( mHelp(end), "" ) ) % trim blank concluding lines
     mHelp(end) = [] ;
end

end
