function [] = filter( Img, weights, Params )
%FILTER  3D low-pass filtering (weighted, unweighted, or median)
%
% Wraps to smooth3() or medfilt3().
%
% [] = FILTER( Img )
% [] = FILTER( Img, weights )
% [] = FILTER( Img, weights, Params )
% 
% Inputs:
%
% Img 
%   the Mrdi-type image volume.
%
% weights
%   an array of data weights (>=0) to penalize (for smooth3) or exclude (for
%   medfilt3()) identifiable outliers.  Dimensions of weights must be the
%   same as the image volume (i.e. Img.getgridsize() )
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

DEFAULTS.kernelSize = [3 3 3] ;
DEFAULTS.method     = 'gaussian' ;

if nargin < 3 || isempty(Params)
    Params.dummy = [] ;
end

Params = assignifempty( Params, DEFAULTS ) ;

nEchoes  = size( Img.img, 4 ) ;
nVolumes = size( Img.img, 5 ) ;

if nargin < 2 || isempty( weights )
    weights = [ ones( size(Img.img) ) ] ;
else
    assert( all( size(weights) == size(Img.img) ), ...
        'Filter weights and image volume must possess the same dimensions' ) ;
end

for iVolume = 1 : nVolumes
    for iEcho = 1 : nEchoes

        img = Img.img(:,:,:, iEcho, iVolume ) ;
        w   = weights(:,:,:, iEcho, iVolume ) ;

        switch Params.method
            case 'median'
                img( ~w ) = NaN ; % medfilt3() will ignore these values
                img = medfilt3( img, Params.kernelSize ) ; 
                img( ~w ) = 0 ;
            otherwise 
                weightsSmoothed = smooth3( w, Params.method, Params.kernelSize ) ;
                Img.img = smooth3( w .* img, Params.method, Params.kernelSize ) ./ weightsSmoothed ; 
        end
        
        Img.img(:,:,:, iEcho, iVolume) = img ;
    end
end

end
