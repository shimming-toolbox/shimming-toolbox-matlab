function [ Info ] = whatr( baseDir, isExcludingHidden )
% WHATR Recursive` what()` 
%     
%     Info = whatr( baseDir, isExcludingHidden )
%
%  Calls MATLAB function [what] for `baseDir` and each of its subdirectories,
%  and arrays the results to return `Info`.
%
% __INPUTS__ 
%     
%   `baseDir=["."]`
%     The parent/top directory path as a string scalar or character vector.
%
%   `isExcludingHidden=[true|1]`
%     A Boolean toggle: when set to false, hidden folders are included in the returned list.
%
% __ETC__ 
%
% To help remember the function name, think "what -r" on the commandline, and
% "what files are there?" in English.
%
% - [what](https://www.mathworks.com/help/matlab/ref/what.html)
%
% See also
% WHAT
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
