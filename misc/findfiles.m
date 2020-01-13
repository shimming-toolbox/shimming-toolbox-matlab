function list = findfiles( searchPath, poi, isRecursive )
%FINDFILES Return list of files matching given naming pattern
%
% Searches directory (recursively by default) for files matching a
% given naming pattern of interest (all-files wildcard '*.*' by default)
% to return the list of matches as a cell-of-cells, whereby outer-cells
% correspond to subdirectories containing file matches, and their corresponding
% inner-cell entries contain the full file paths.
%
% Usage [default arguments in square brackets]
% 
%   list = FINDFILES( searchPath, patternOfInterest ['.*'], isRecursive [true] ) 
% 
% NOTE A supposedly equivalent function exists on <a href="matlab: web('https://www.mathworks.com/matlabcentral/fileexchange/57298-recursively-search-for-files')">Matlabcentral</a>
% However it does not work for me (RT), running:
% MATLAB Version: 9.7.0.1216025 (R2019b) Update 1
% MATLAB License Number: STUDENT
% Operating System: Mac OS X  Version: 10.15.1 Build: 19B88

DEFAULT_POI         = '*.*' ;
DEFAULT_ISRECURSIVE = true ;

if ( nargin < 1 ) || isempty( searchPath ) 
    error('Search path not provided: Not enough input arguments.') ;
end

if ( nargin < 2 ) || isempty( poi ) 
    poi = DEFAULT_POI ;
end

if ( nargin < 3 ) || isempty( isRecursive ) 
    isRecursive = DEFAULT_ISRECURSIVE ;
end

list = {} ;

if ~exist( searchPath ) 
    
    return ;

elseif isfile( searchPath ) 
    
    [~, ~, ext] = fileparts( searchPath ) ;
   
    if strcmp( ext, poi )   
        list{1} = { searchPath } ; % return cell-containing-cell to be consistent with folder/directory case
    end
    
    return ;

elseif isfolder( searchPath )

    if ~strcmp( searchPath(end), filesep )
        searchPath(end+1) = filesep ;    
    end

    if isRecursive
        Dirs = dir( [ searchPath '**' filesep poi ] ) ;
    else
        Dirs = dir( [ searchPath poi ] ) ;
    end

    % find unique subdirectory names and indices of the corresponding Dirs entries 
    DirsC = struct2cell( Dirs ) ;
    [ ~, ~, ithSubDir ] = unique( DirsC(2,:) ) ;
    nSubDirs = numel( unique( ithSubDir ) ) ;

    list = cell( 1, nSubDirs ) ;

    for i = 1 : nSubDirs 
    
        % indices of Dirs entries corresponding to ithSubDir 
        ind_ithSubDir = find( ithSubDir == i ) ;
        
        % number of [.poi] files within ithSubDir 
        nFiles_ithSubDir = length( ind_ithSubDir ) ;
        
        % file list cell-array for the ithSubDir 
        list_ithSubDir = cell( nFiles_ithSubDir, 1 ) ;
        
        for iFile = 1 : nFiles_ithSubDir
            list_ithSubDir{iFile, 1} = [ Dirs( ind_ithSubDir(iFile) ).folder filesep Dirs( ind_ithSubDir(iFile) ).name ] ;
        end
        
        list{i} = list_ithSubDir ; 

    end 
    
end

end
