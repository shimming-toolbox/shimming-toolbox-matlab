function [self] = initializefromdicom( self, Hdrs ) 
%INITIALIZEFROMDICOM

% "If Anatomical Orientation Type (0010,2210) is absent or has a value of
% BIPED, the x-axis is increasing to the left hand side of the patient. The
% y-axis is increasing to the posterior side of the patient. The z-axis is
% increasing toward the head of the patient."
%
% from DICOM standard: https://www.dabsoft.ch/dicom/3/C.7.6.2.1.1/
% assert( ~myisfield( self.Img.Hdr, 'AnatomicalOrientationType' ) || ...
%             strcmp( self.Img.Hdr.AnatomicalOrientationType, 'BIPED' ), ...
%             'Error: AnatomicalOrientationType not supported.' ) ;

self.size = double( [Hdrs(1).Rows Hdrs(2).Columns size(Hdrs, 1) ] ) ;

self.imageOrientationPatient = Hdrs(1).ImageOrientationPatient ;

if isfield( Hdrs, 'SpacingBetweenSlices' )
    self.spacing = ...
        [ Hdrs(1).PixelSpacing(1) Hdrs(1).PixelSpacing(2) Hdrs(1).SpacingBetweenSlices ] ;
else
    self.spacing = ...
        [ Hdrs(1).PixelSpacing(1) Hdrs(1).PixelSpacing(2) Hdrs(1).SliceThickness ] ;
end

for iSlice = 1 : self.size(3)
    self.imagePositionPatient(:, iSlice) = Hdrs(iSlice, 1).ImagePositionPatient ;
end


% if size( self.Img.img, 3 ) > 1
% % Determine: ascending or descending slices?
%
%     % Estimate positions of last slice using the position of the 1st image 
%     % and the 2 possibilities for the SliceNormalVector:
%
%     % 1. using cross( c, r ) ;
%     % [X1,Y1,Z1] = [self.voxelPositions.X, self.voxelPositions.Y, self.voxelPositions.Z] ;     
%     [X1,Y1,Z1] = deal( self.voxelPositions.X,  self.voxelPositions.Y, self.voxelPositions.Z ) ;
%     
%     % 2. using the reverse
%     self.SliceNormalVector = cross( r, c ) ;
%     [X2,Y2,Z2] = deal( self.voxelPositions.X,  self.voxelPositions.Y, self.voxelPositions.Z ) ;
%     
%     % Actual position corresponding to the slice direction can be increasing or
%     % decreasing with slice/image number. So, which estimate is closer: 1 or 2? 
%     if norm( self.Img.Hdrs{end,1,1}.ImagePositionPatient' - [ X1(1,1,end) Y1(1,1,end) Z1(1,1,end) ] ) < ...
%             norm( self.Img.Hdrs{end,1,1}.ImagePositionPatient' - [ X2(1,1,end) Y2(1,1,end) Z2(1,1,end) ] ) 
%         % if true, then 1. corresponds to the correct orientation
%         self.SliceNormalVector = cross( c, r ) ;
%     end
% end
%     

end
