function [ Info ] = whatr( baseDir, isExcludingHidden )
%WHATR Returns struct array from recursive calls to WHAT across subdirectories
%
% WHATR denotes _recursive_ WHAT ("what -r"?)
% 
% ### Usage ###
% 
% Info = WHATR( baseDir, isExcludingHidden )
%
% #### Inputs ####
% 
% baseDirectory 
%
% - parent/top directory path as a string scalar or character vector.
%   [default: "." ]
%
% - isExcludingHidden     
%   [default = 1]
%   Boolean toggle includes hidden folders in the returned list when set to 0 (false).
%
% ### References ###
%
% See also:
%
% - [mapdirectorytree](local link?) 
%
% - <https://www.mathworks.com/help/matlab/ref/what.html>
    arguments
        baseDir(1,:) { mustBeStringScalarOrCharVector, mustBeFolder } = "." ;
        isExcludingHidden(1,1) { mustBeBoolean } = true ;
    end

baseDir = string( baseDir ) ;

%% recursive what() 
Info = what( baseDir ) ;
dirs = mapdirectorytree( baseDir, false, isExcludingHidden ) ;   

for iSubDir = 1 : numel( dirs )
    Info = [ Info ; what( dirs( iSubDir ) ) ] ;
end

end
