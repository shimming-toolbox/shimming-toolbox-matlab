function [Self] = initializefromdicom( Self, Hdrs ) 
%INITIALIZEFROMDICOM

% "If Anatomical Orientation Type (0010,2210) is absent or has a value of
% BIPED, the x-axis is increasing to the left hand side of the patient. The
% y-axis is increasing to the posterior side of the patient. The z-axis is
% increasing toward the head of the patient."
%
% from DICOM standard: https://www.dabsoft.ch/dicom/3/C.7.6.2.1.1/
% assert( ~myisfield( Self.Img.Hdr, 'AnatomicalOrientationType' ) || ...
%             strcmp( Self.Img.Hdr.AnatomicalOrientationType, 'BIPED' ), ...
%             'Error: AnatomicalOrientationType not supported.' ) ;

Self.size = double( [Hdrs(1).Rows Hdrs(2).Columns size(Hdrs, 1) ] ) ;

Self.imageOrientationPatient = Hdrs(1).ImageOrientationPatient ;

if isfield( Hdrs, 'SpacingBetweenSlices' )
    Self.spacing = [ Hdrs(1).PixelSpacing(1) Hdrs(1).PixelSpacing(2) Hdrs(1).SpacingBetweenSlices ] ;
else
    Self.spacing = [ Hdrs(1).PixelSpacing(1) Hdrs(1).PixelSpacing(2) Hdrs(1).SliceThickness ] ;
end

for iSlice = 1 : Self.size(3)
    Self.imagePositionPatient(:, iSlice) = Hdrs(iSlice, 1).ImagePositionPatient ;
end


% if size( Self.Img.img, 3 ) > 1
% % Determine: ascending or descending slices?
%
%     % Estimate positions of last slice using the position of the 1st image 
%     % and the 2 possibilities for the SliceNormalVector:
%
%     % 1. using cross( c, r ) ;
%     % [X1,Y1,Z1] = [Self.voxelPositions.X, Self.voxelPositions.Y, Self.voxelPositions.Z] ;     
%     [X1,Y1,Z1] = deal( Self.voxelPositions.X,  Self.voxelPositions.Y, Self.voxelPositions.Z ) ;
%     
%     % 2. using the reverse
%     Self.SliceNormalVector = cross( r, c ) ;
%     [X2,Y2,Z2] = deal( Self.voxelPositions.X,  Self.voxelPositions.Y, Self.voxelPositions.Z ) ;
%     
%     % Actual position corresponding to the slice direction can be increasing or
%     % decreasing with slice/image number. So, which estimate is closer: 1 or 2? 
%     if norm( Self.Img.Hdrs{end,1,1}.ImagePositionPatient' - [ X1(1,1,end) Y1(1,1,end) Z1(1,1,end) ] ) < ...
%             norm( Self.Img.Hdrs{end,1,1}.ImagePositionPatient' - [ X2(1,1,end) Y2(1,1,end) Z2(1,1,end) ] ) 
%         % if true, then 1. corresponds to the correct orientation
%         Self.SliceNormalVector = cross( c, r ) ;
%     end
% end
%     

end
