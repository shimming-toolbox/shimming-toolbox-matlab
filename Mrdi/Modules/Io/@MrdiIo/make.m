function [Imgs] = make( varargin ) 
%MAKE  Return cell array of Mrdi objects 
%     
%     [Imgs] = make( searchDir ) 
%     [Imgs] = make( Params ) 
%
% Searches for images files and initializes Mrdi objects with the image data. 
%
% __INPUTS__
%   
%   `searchDir` 
%     directory to search for images (string or char array)

%   Options: [default: = MrdiIo.defaults ]
%   
%   Struct of parameters for which the following fields are supported:
%
%   .fileExt [default: { '.dcm' ; '.IMA' } ( i.e. = MrdiIo.SUPPORTED_INPUTS ) ]
%       file extensions of images (char or char-containing cell array) 
%       (fileExt must be a member of the given default values)
%
%   .isSearchRecursive [default: TRUE]
%       if true, subdirectories are included in the image search.
%
%   .nBytesMax [default: 2 * 1E9]
%       byte limit on load. Function aborts if limit is exceeded.
%       (NOTE: matlab 'memory' function exists to query available memory but
%       currently only works on Windows)
%
% When the image type is known to correspond to a specific Mrdi
% child class, the object will be typecast accordingly.
 
    % arguments
    %     searchDir(1,:) ;
    %     Options(:,1) struct = MrdiIo.defaults() ;
    % end

DEFAULTS = MrdiIo.makedefaults ;

%% Check inputs
narginchk(1,1) ;

if isstring( varargin{1} ) || ischar( varargin{1} )

    searchDir = varargin{1} ;
    validateattributes( searchDir, {'string', 'char'}, {'scalartext'}, mfilename, 'searchDir', 1 ) ;
    Params = DEFAULTS ;

elseif isstruct( varargin{1} )

    Params = varargin{1} ;
    Params = assignifempty( Params, DEFAULTS ) ;

    Mc            = metaclass( MrdiIo ) ;
    propertyNames = string( { Mc.PropertyList.Name } ) ;

    for iParam = 1 : length( DEFAULTS )

        if
        Mc.PropertyList( iParam ).

    
%% Find image files 
List = MrdiIo.findimagefiles( searchDir, fileExt, isSearchRecursive ) ;

%% Read in files 
[imgs, Hdrs] = loadandsortimages( List, nBytesMax )

%% ----- 
% Initialize Mrdi objects from image data
for iSeries = 1 : numel( imgs )
    Imgs{iImg} = Mrdi( imgs{iSeries}, Hdrs{iSeries} ) ;
end

end

