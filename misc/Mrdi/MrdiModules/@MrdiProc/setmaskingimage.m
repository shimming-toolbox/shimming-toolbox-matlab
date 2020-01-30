function [] = setmaskingimage( Img, mask )
%SETMASKINGIMAGE Copies valid mask (a logical array of 1's and 0's) to Img.Hdr.MaskingImage
% 
% [] = SETMASKINGIMAGE( Img, mask ) 
%
% The purpose of this function is to specify the signal spatial support (e.g. 
% of mag, phase, field data) within the image grid. e.g. it is called during 
% phase unwrapping/field mapping, but it might also be called prior to regridding
% if the interpolation should exclude certain voxels.
% 
% To be valid, mask must either be the same size as Img.img OR the same size as
% Img.getgridsize() (i.e. the size of a single image volume of a
% multi-echo/multi-measurement stack), in which case, the single mask is simply
% copied such that the assigned Img.Hdr.MaskingImage always possesses the same
% dimensions as Img.img

assert( nargin == 2, 'Function requires at least 2 input arguments. See HELP MrdiProc.setmaskingimage' ) ;

if ~islogical( mask )
    assert( numel(mask) == ( nnz( mask(:) == 0 ) + nnz( mask(:) == 1 ) ), ...
        'Input mask must consist exclusively of zeros and ones. See HELP MrdiProc.setmaskingimage' ) ;
    mask = logical( mask ) ;
end

maskSize = size( mask ) ;

%% -----
% assign masking image 
if ndims( Img.img ) == ndims( mask ) && all( size(Img.img) == maskSize )
    Img.Hdr.MaskingImage = mask ;

elseif ndims( Img.img ) > ndims( mask ) 
    
    if numel( maskSize ) == 2
        maskSize = [ maskSize 1 ] ;
    end

    if ( numel( Img.getgridsize() ) == numel( maskSize ) ) && ... 
       ( all( Img.getgridsize() == maskSize ) ) % single mask provided
        Img.Hdr.MaskingImage = repmat( mask, [1 1 1 size(Img.img, 4) size(Img.img, 5)] ) ;
    else
        error('Input mask and Img.img should possess the same dimensions' ) ;
    end
else 
    error('Input mask and Img.img should possess the same dimensions' ) ;
end

end

