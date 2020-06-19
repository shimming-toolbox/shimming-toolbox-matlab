function [] = write_json( obj, filename, Options )
%WRITE_JSON Encodes an object instance as json and writes to file
%    
%    save_json( obj, filename )
%    save_json( obj, filename, "isOverwriting", true )
% 
% Encodes the scalar struct or object `obj` as json and writes to the path
% string `filename`.
%
% To overwrite an existing file, use the `(..., "isOverwriting", true)`
% name-value pair. 
%
% __ETC__
%
% `save_json` wraps to Matlab's `jsonencode` which generates content in "compact"
% form. The file can be "prettified" for readability using tools such as:
%
% 1. [webbrowser](https://jsonformatter.org/json-pretty-print)  
% 2. python CLI: `python -m json.tool filename`
% 3. [vscode](https://marketplace.visualstudio.com/items?itemName=vthiery.prettify-selected-json)
%
% See also  
% [ jsonencode ](https://www.mathworks.com/help/matlab/ref/jsonencode.html)
    arguments
        obj(1,1) {mustBeStructOrObject};
        filename(1,1) string;
        Options.isOverwriting {mustBeScalarOrEmpty, valid.mustBeBoolean} = false;
    end
    
    if ~contains(filename, '.') % Add extension
        filename = strcat( filename, '.json') ;
    end

    assert( [~isfile(filename) | Options.isOverwriting], ['File already exists.' ...
         'Use the name-value pair (...,"isOverwriting", true) to force overwrite.'])

    json = jsonencode( obj ) ;
    
    fid = fopen( filename, 'w' ) ;
    fwrite( fid, json ) ; 
    fclose(fid) ;
    
end

function [] = mustBeStructOrObject( x )
    assert( [isstruct(x) | isobject(x)], 'Input must be a scalar struct or object' )
end
