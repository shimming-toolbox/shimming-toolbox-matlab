function [txtOut] = typecast( txtIn, castAs )
% TYPECASTTEXT Convert between string, char, cellstr
% 
%     [txtOut] = typecasttext( txtIn ) 
%     [txtOut] = typecasttext( txtIn, castAs ) 
%
% Converts an array of text `txtIn` to the class specified in `castAs`, which
% may be `"string"`, `"char"` or `"cell"` (or, equivalently, `"cellstr"`). 
% When called with a single input, the function returns the text as a string.
% As with [cellstr], trailing white space is deleted when converting to
% string from char.
%
% __ETC__
% 
% Matlab has three options for text-containers:
%
% 1. `char`, which is the most primitive and least flexible.
%
% 2. `cell`, which is a general container for any data type.
%
% 3. `string` which has the virtue of being much more flexible than `char`
% while still being specialized for text, unlike `cell`. However, as the most
% recently introduced (2016b) of the three, it is often necessary to convert
% to one of the other two.
%
% [string](https://www.mathworks.com/help/matlab/ref/string.html)
% [cellstr](https://www.mathworks.com/help/matlab/ref/cellstr.html)
%
% See also
% CHAR,CELL,CELLSTR,STRING
    arguments
        txtIn {mustBeStringOrCharOrCellstr} ;
        castAs(1,:) {mustBeStringScalarOrCharVector, ...
            mustBeMember( castAs, ["string" "char" "cell" "cellstr"] ) } = "string" ;
    end

    if strcmp( castAs, 'cellstr' )
        castAs = 'cell' ;
    end

    castFrom = class( txtIn ) ;

    if strcmp( castFrom, castAs )
        txtOut = txtIn ;
        return ;
    end
   
    switch castAs
        case 'cell'
            funcast = @cellstr ;
        otherwise
            funcast = str2func( castAs ) ;
    end

    txtOut = funcast( txtIn ) ;
    
    if strcmp( castFrom, 'char' ) && strcmp( castAs, 'string' )
        txtOut = deblank( txtOut ) ;
    end

end
