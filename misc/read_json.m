function [info] = read_json( filename )
%READ_JSON Return struct from json file (system path or URL)
%    
%    info = read_json( filename )
%
% See also  
% [ jsondecode ](https://www.mathworks.com/help/matlab/ref/jsondecode.html)
    arguments
        filename(1,1) string; 
    end

    try
        if isfile( filename )
            jsonTxt = fileread( filename );
        else
            jsonTxt = webread( filename, weboptions( 'ContentType', 'text') );
        end

        info = jsondecode( jsonTxt );
    
    catch Me
        disp('Failed to read file.')
        Me.rethrow;
    end

end
