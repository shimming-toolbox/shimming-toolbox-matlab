function list = findfiles( searchDir, poi, isRecursive )
%FINDFILES List files matching a given naming pattern
%
% FINDFILES searches a directory (recursively by default) for files matching a
% given naming pattern of interest (all-files wildcard '*.*' by default) and
% returns the list of matches as a 1-D cell array with the full file paths
% along the rows.
%
% Usage [default arguments in square brackets]
% 
%   list = FINDFILES( searchDir ['./'], patternOfInterest ['*.*'], isRecursive [true] ) 
 
% NOTE A supposedly equivalent function: 
% <a href="matlab: web('https://www.mathworks.com/matlabcentral/fileexchange/57298-recursively-search-for-files')">Matlabcentral</a>
% However it does not work for me (RT), running:
% MATLAB Version: 9.7.0.1216025 (R2019b) Update 1
% MATLAB License Number: STUDENT
% Operating System: Mac OS X  Version: 10.15.1 Build: 19B88

DEFAULT_SEARCHDIR = ' ./' ;
DEFAULT_POI       = '*.*' ;

if ( nargin < 1 ) || isempty( searchDir ) 
    searchDir  = DEFAULT_SEARCHDIR ;
elseif ~( ischar( searchDir ) || isstring( searchDir ) )
    error('Input searchDir must be a character array or string.');
end

list = {} ;

if ~exist( searchDir ) 
    
    return ;

elseif isfile( searchDir ) 
    
    list{1} = searchDir ; 
    
    return ;

elseif isfolder( searchDir )

    if ~strcmp( searchDir(end), filesep )
        searchDir(end+1) = filesep ;    
    end

    if isRecursive
        Dirs = dir( [ searchDir '**' filesep poi ] ) ;
    else
        Dirs = dir( [ searchDir poi ] ) ;
    end

    nDirs = length(Dirs) ;
    list  = cell( nDirs, 1 ) ;

    for iDir = 1 : nDirs 
         list{ iDir } = [ Dirs(iDir).folder filesep Dirs(iDir).name ] ;
    end
    
end

end

function [] = mustBeChar( c )
%MUSTBECHAR Function validator: Error if input is not character array/string
    if ~( ischar( c ) || isstring(c) )
        error('Must be a character array.');
    end
end
