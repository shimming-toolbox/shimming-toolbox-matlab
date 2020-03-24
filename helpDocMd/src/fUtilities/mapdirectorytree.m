function [dirs] = mapdirectorytree( baseDir, isReturnRelative, isExcludingHidden )
% MAPDIRECTORYTREE Returns list of subdirectories
%
% ### Usage ### 
%
% dirs = MAPDIRECTORYTREE( baseDir, isReturnRelative, isExcludingHidden )
% 
% Returns a string column vector wherein each element is a path to a subdirectory.
%
% ### Inputs ###
%
% - baseDir
%   [default = "."]
%   Parent/top directory as string scalar
%
% - isReturnRelative 
%   [default = 1]
%   Boolean toggle to return the subdirectory paths as relative (true) or absolute (false). 
%
% - isExcludingHidden     
%   [default = 1]
%   Boolean toggle includes hidden folders in the returned list when set to 0 (false).
%
% ### References ###
%
% See also fileattrib, dir
    arguments
        baseDir(1,:) { mustBeStringScalarOrCharVector, mustBeFolder } = "." ;
        isReturnRelative(1,1) { mustBeBoolean }  = true ;
        isExcludingHidden(1,1) { mustBeBoolean } = true ;
    end

%% Map subdirectories 
[~,Info,~] = fileattrib( strcat( baseDir, filesep, "*") ) ;
Info       = Info( [ Info.directory ] ) ;

if isExcludingHidden
    if ispc 
        Info    = Info( ~[ Info.hidden ] ) ;
    else
        iHidden = contains( string( { Info.Name } ), strcat(filesep,".") ) ;
        Info    = Info( ~iHidden ) ;
    end
end

dirs = [ string( { Info.Name } )' ] ;

%% Trim out the abs. parent dir to return the relative paths
if isReturnRelative
    dirs = replace( dirs, baseDir, "." ) ;
end

end
