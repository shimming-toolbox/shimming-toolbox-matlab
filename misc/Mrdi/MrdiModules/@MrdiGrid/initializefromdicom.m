function [Grid] = initializefromdicom( Grid, Hdrs ) 

% "If Anatomical Orientation Type (0010,2210) is absent or has a value of
% BIPED, the x-axis is increasing to the left hand side of the patient. The
% y-axis is increasing to the posterior side of the patient. The z-axis is
% increasing toward the head of the patient."
%
% from DICOM standard: https://www.dabsoft.ch/dicom/3/C.7.6.2.1.1/
% assert( ~myisfield( Grid.Img.Hdr, 'AnatomicalOrientationType' ) || ...
%             strcmp( Grid.Img.Hdr.AnatomicalOrientationType, 'BIPED' ), ...
%             'Error: AnatomicalOrientationType not supported.' ) ;

Grid.size = uint64( [Hdrs(1).Rows Hdrs(2).Columns size(Hdrs,3) ] ) ;

Grid.imageOrientationPatient = Hdrs(1).ImageOrientationPatient ;

if isfield( Hdrs{1}, 'SpacingBetweenSlices' )
    Grid.spacing = [ Hdrs{1}.PixelSpacing(1) Hdrs{1}.PixelSpacing(2) Hdrs{1}.SpacingBetweenSlices ] ;
else
    Grid.spacing = [ Hdrs{1}.PixelSpacing(1) Hdrs{1}.PixelSpacing(2) Hdrs{1}.SliceThickness ] ;
end

for iSlice = 1 : Grid.size(3)
    Grid.imagePositionPatient(:,iSlice) = Hdrs{iSlice,1}.ImagePositionPatient ;
end


% if size( Grid.Img.img, 3 ) > 1
% % Determine: ascending or descending slices?
%
%     % Estimate positions of last slice using the position of the 1st image 
%     % and the 2 possibilities for the SliceNormalVector:
%
%     % 1. using cross( c, r ) ;
%     % [X1,Y1,Z1] = [Grid.voxelPositions.X, Grid.voxelPositions.Y, Grid.voxelPositions.Z] ;     
%     [X1,Y1,Z1] = deal( Grid.voxelPositions.X,  Grid.voxelPositions.Y, Grid.voxelPositions.Z ) ;
%     
%     % 2. using the reverse
%     Grid.SliceNormalVector = cross( r, c ) ;
%     [X2,Y2,Z2] = deal( Grid.voxelPositions.X,  Grid.voxelPositions.Y, Grid.voxelPositions.Z ) ;
%     
%     % Actual position corresponding to the slice direction can be increasing or
%     % decreasing with slice/image number. So, which estimate is closer: 1 or 2? 
%     if norm( Grid.Img.Hdrs{end,1,1}.ImagePositionPatient' - [ X1(1,1,end) Y1(1,1,end) Z1(1,1,end) ] ) < ...
%             norm( Grid.Img.Hdrs{end,1,1}.ImagePositionPatient' - [ X2(1,1,end) Y2(1,1,end) Z2(1,1,end) ] ) 
%         % if true, then 1. corresponds to the correct orientation
%         Grid.SliceNormalVector = cross( c, r ) ;
%     end
% end
%     
% function [X,Y,Z] = calculategridpositions( size, spacing, rotationMatrix, imagePositionPatient ) 
%
%     % Arrays of voxel row, column, and slice indices
%     [iRows,iColumns,iSlices] = ndgrid( [0:1:Grid.size(1)-1], ...
%                       [0:1:Grid.size(2)-1], ...
%                       [0:1:Grid.size(3)-1] ) ; 
%
%     % Rotation and Scaling matrix: RS  
%     RS = Grid.rotationMatrix * diag(Grid.spacing) ;
%
%     % Scale and rotate to align row direction with x-axis, column direction
%     % with y-axis, slice with z-axis; then translate w.r.t origin (via ImagePositionPatient)
%     X = ( RS(1,1)*iRows + RS(1,2)*iColumns + RS(1,3)*iSlices ) + Grid.Img.Hdr.ImagePositionPatient(1) ;
%     Y = ( RS(2,1)*iRows + RS(2,2)*iColumns + RS(2,3)*iSlices ) + Grid.Img.Hdr.ImagePositionPatient(2) ;
%     Z = ( RS(3,1)*iRows + RS(3,2)*iColumns + RS(3,3)*iSlices ) + Grid.Img.Hdr.ImagePositionPatient(3) ;
%
% end

end

