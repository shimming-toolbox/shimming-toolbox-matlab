function [x, y, z] = get_nii_coordinates( nii )
%GET_NII_COORDINATES Return voxel position "world-coordinates" in mm
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
% __TODO__
% * function should be further tested.

% -----------
%% Check input
narginchk(1,1);

if [ isstring( nii ) | ischar( nii ) ] & isfile( nii )
    
    [~, info] = img.read_nii( nii );

elseif isstruct( nii ) && isfield( nii, 'ImageSize' ) && ...
        isfield( nii, 'Transform' ) && isa( nii.Transform, 'affine3d' ) 

    info = nii ;

else
    error('Input must be a path to a .nii file, or a struct returned from `niftiinfo`')
end

% ---------------
%% Get coordinates

% voxel indices
[i, j, k] = ndgrid( [0:info.ImageSize(1)-1], [0:info.ImageSize(2)-1], [0:info.ImageSize(3)-1] ) ;
[x, y, z] = info.Transform.transformPointsForward( i, j, k ) ;

end
