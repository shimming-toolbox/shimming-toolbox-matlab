function [dirs] = mapdirectorytree( baseDir, isReturnRelative, isExcludingHidden )
% MAPDIRECTORYTREE Returns list of subdirectories
%    
%     dirs = MAPDIRECTORYTREE( baseDir, isReturnRelative, isExcludingHidden )
% 
% Returns paths to the subdirectories of `baseDir` as elements of the string
% column vector `dirs`.
%
% __INPUTS__
%
%   baseDir=["."]
%     Parent/top directory as string scalar
%
%   isReturnRelative=[true|1]
%     Boolean toggle: Returns the subdirectory paths as relative to `baseDir` 
%     when true, or in absolute terms when false.
%
%   isExcludingHidden=[true|1]
%     Boolean toggle includes hidden folders in the returned list when set to 0
%     (false).
%
% __ETC__
%
% See also 
% FILEATTRIB, DIR
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
