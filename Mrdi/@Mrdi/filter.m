function [] = filter( Img, weights, Params )
%FILTER  3D low-pass filtering (weighted, unweighted, or median)
%
% Wraps to smooth3() or medfilt3().
%
% Usage:
% -----
%
% (where Img is a Mrdi-type image object...)
%
% [] = Img.filter( )
% [] = Img.filter( weights )
% [] = Img.filter( weights, Params )
% 
% Inputs:
% ------
%
% weights
%   an array of data weights (>=0) to penalize (for smooth3) or exclude (for
%   medfilt3()) identifiable outliers.  Dimensions of weights must be the
%   same as the image volume (i.e. Img.size )
%
% Params
%   an optional struct for which the following Params.fields are supported
%
%   .kernelSize
%       in number of voxels
%       default = [3 3 3] 
%
%   .method
%       'gaussian' OR 'box' OR 'median'
%       default = 'gaussian'
% 
% TODO
%   Add support for 2d (single slice) images 

assert( Img.size(3) > 1, 'Method does not currently support 2d images. TODO' )

%% -----
% Validate inputs, assign defaults
DEFAULTS.kernelSize = [3 3 3] ;
DEFAULTS.method     = 'gaussian' ;

if nargin < 2 || isempty( weights )
    weights = ones( Img.size ) ;
else
    assert( isequal( size( weights ), Img.size ), ...
        'Filter weights and image volume must possess the same dimensions' ) ;
end

if nargin < 3 || isempty(Params)
    Params.dummy = [] ;
end

Params = assignifempty( Params, DEFAULTS ) ;

%% -----
% Filtering

for iImg = 1 : Img.nVolumes

    img = Img.img(:,:,:, iImg ) ;
    w   = weights(:,:,:, iImg ) ;

    switch Params.method
        case 'median'
            img( ~w ) = NaN ; % medfilt3() will ignore these values
            img       = medfilt3( img, Params.kernelSize ) ; 
            img( ~w ) = 0 ;

        otherwise 
            weightsSmoothed = smooth3( w, Params.method, Params.kernelSize ) ;
            img             = smooth3( w .* img, Params.method, Params.kernelSize ) ./ weightsSmoothed ; 
    end
    
    Img.img(:,:,:, iImg ) = img ;

end

end
