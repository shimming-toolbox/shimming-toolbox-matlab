function [parentDir] = parent( pathIn, level )
% PARENT Return immediate and common (ap-)parent directories 
%    
%    [ parentDir ] = parent( pathIn, 'single' ) ;
%    [ parentDir ] = parent( pathIn, 'common' ) ;
%
% When called with a single argument, or with the second argument 'single',
% `parent` wraps to Matlab function [fileparts] to return the immediate (single level)
% parent folder(s) of `pathIn` as `parentDir`. 
%
% With a second argument of 'common', the unique path to the common base
% directory among `pathIn` (i.e. the "fork") is returned. If the input paths
% have nothing in common (even the leading character), then baseDir is assigned
% '?'. 
%
% If `pathIn` consists of a single entry the result will be the same in both cases.
%
% **NOTE** Like `fileparts`, the input does not need to refer to existing
% file system paths (consider the mnemonic "ap-*parent*" for this function).
% Unlike `fileparts`, the input to `parent` can contain multiple path elements.
%
% __EXAMPLE__
% ```
% % pathIn is a triple-row string vector
% >> pathIn = [ "/A/B/C" ; "/A/B/C2/D2/E2" ; "/A/B/C3/D3/E3" ; ] ;
% 
% % Omiting the terminating semi-colon for display:
% >> parent( pathIn )
%
% ans =
%
%   3x1 string array
%
%     "/A/B"
%     "/A/B/C2/D2"
%     "/A/B/C3/D3"
%
% >> parent( pathIn, 'common' )
%
% ans =
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
        level {mustBeMember(level,["single" "common"])} = "single" ;
    end

P = Pathologist( pathIn ) ;

if strcmp( level, 'single' ) || ( P.nPaths == 1 )
    parentDir = P.parentDir ;
else
    parentDir = P.baseDir ;
end

end
