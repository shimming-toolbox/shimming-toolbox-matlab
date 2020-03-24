function [parentDir, baseDir] = parent( pathIn )
% PARENT Return immediate and common (ap-)parent directories 
%    
%    [ parentDir, baseDir ] = parent( pathIn ) ;
%
% Wraps to Matlab function [fileparts] to return the immediate (single level)
% parent folder(s) of `pathIn` as `parentDir`. If `pathIn` contains multiple
% paths, then the path to the first common base folder (the "fork") is returned
% as `baseDir`.
% 
% **NOTE** Like `fileparts`, the input does not need to refer to existing
% file system paths (consider the mnemonic "ap-*parent*" for this function).
% Unlike `fileparts`, the input to `parent` can contain multiple path elements.
% 
% __INPUTS__
%     
%   pathIn
%     A string-, char-, or cellstr-array of paths.
%
% __OUTPUTS__
%
%   parentDir
%     Paths to the immediate parent folders of the elements of `pathIn`, formed
%     by calling `fileparts( )` for each element.
%  
%   baseDir
%     The first common parent directory the elements of `pathIn`. (The
%     "fork").  If pathIn consists of a single path, `baseDir` is the same as
%     `parentDir`. If the input paths have nothing in common (even the leading
%     character), then baseDir is assigned '?'.
%
% __Example 1__
% ```
% % pathIn is a triple-row string vector
% >> pathIn = [ "/A/B/C" ; "/A/B/C2/D2/E2" ; "/A/B/C3/D3/E3" ; ] ;
% 
% % omiting the terminating semi-colon for display produces:
% >> [parentDir, baseDir ] = parent( pathIn )
%
% parentDir =
%
%   3x1 string array
%
%     "/A/B"
%     "/A/B/C2/D2"
%     "/A/B/C3/D3"
%
% baseDir =
%
%     "/A/B"
%
% ```
% __ETC__
% - [fileparts](https://www.mathworks.com/help/matlab/ref/fileparts.html)
%
% See also
% FILEPARTS
    arguments
        pathIn {mustBeStringOrCharOrCellstr}
    end

P         = Pathologist( pathIn ) ;
parentDir = P.parentDir ;
baseDir   = P.baseDir ;

end
