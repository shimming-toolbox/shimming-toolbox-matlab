function [x, y, z] = get_nii_coordinates( nii )
%GET_NII_COORDINATES Return voxel position arrays in mm
%     
%     [x, y, z] = get_nii_coordinates( niiFile )
%     [x, y, z] = get_nii_coordinates( niiInfo )
% 
% Returns `x, y, z`: three 3-D arrays of voxel positions in mm.
%
% The single input can be given either as:
%
% 1. `niiFile`—A path string to the NIfTI file of interest; or, 
%
% 2. `niiInfo`—A struct of the form `niiInfo = niftiinfo( niiFile )`
%
% __NOTE__
% The function is implemented as
% ```
% [i, j, k] = ndgrid( [0:niiInfo.ImageSize(1)-1], [0:niiInfo.ImageSize(2)-1], [0:niiInfo.ImageSize(3)-1] ) ;
% [x, y, z] = niiInfo.Transform.transformPointsForward( i, j, k ) ;  
% ```
% * Output looks reasonable but the function should be further tested.
% * TODO: support for 2d case. Current implementation support restricted to 3d image volumes.

% -----------
%% Check input
narginchk(1,1);

if isfile( nii )
    
    [~, info] = img.read_nii( nii );

elseif isstruct( nii ) && isfield( nii, 'ImageSize' ) && ...
        isfield( nii, 'Transform' ) && isa( nii.Transform, 'affine3d' ) 

    info = nii ;

else
    error('Input must be a path to a .nii file, or a struct returned from `niftiinfo`')
end

assert( all(info.ImageSize(1:3)>1), 'Current implementation restricted to 3d image volumes. TODO' )
% NOTE: if a single element of info.ImageSize(1:3) == 1, the following call to
% `ndgrid` will return 2 arrays; however, even for a 2d image we still need the
% coordinates of the remaining 3rd dimension (e.g. slice), which can also vary
% between pixels (e.g. for an oblique slice)

% ---------------
%% Get coordinates

% voxel indices
[i, j, k] = ndgrid( [0:info.ImageSize(1)-1], [0:info.ImageSize(2)-1], [0:info.ImageSize(3)-1] ) ;
[x, y, z] = info.Transform.transformPointsForward( i, j, k ) ;

end
