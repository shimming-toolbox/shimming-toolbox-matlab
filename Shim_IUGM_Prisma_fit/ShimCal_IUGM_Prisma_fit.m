classdef ShimCal_IUGM_Prisma_fit < ShimCal 
%SHIMCAL_IUGM_PRISMA_FIT - Shim Calibration 
%
% ShimCal_IUGM_Prisma_fit is a ShimCal subclass 
%
% =========================================================================
% Updated::20180325::ryan.topfer@polymtl.ca
% =========================================================================

properties
    % img ;
    % Hdr ;
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimCal_IUGM_Prisma_fit( Params )
%SHIMCAL - Shim calibration (i.e. reference maps) 

Shim.img = [] ;
Shim.Hdr = [] ;

if nargin < 1
    Params = ShimCal_IUGM_Prisma_fit.declarecalibrationparameters20180319( ) ;
end

% -------
% map dB/dI
[Shim.img, Shim.Hdr] = ShimCal.mapdbdi_allchannels( Params ) ;

% insert channel (not measured) for "0th order" (Larmor transmit frequency):
tmp = zeros( size(Shim.img) + [0 0 0 1] ) ;
tmp(:,:,:,1)     = double( Shim.Hdr.MaskingImage );
tmp(:,:,:,2:end) = Shim.img ;

Shim.img = tmp; 

save( Params.filenameSave ,'Shim' );

end
% =========================================================================
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================
function Params = declarecalibrationparameters20180319( )
%DECLARECALIBRATIONPARAMETERS20180319
% 
% Initializes parameters for shim reference map construction (i.e. shim calibration)

fprintf('##### \n WARNING: Ignoring EchoTime DICOM header when normalizing phase to field. \n\n')
Params.isHardCodingEchoTime = true ;
Params.EchoTimes  = [4.92 7.38] ;
Params.echoTimeDifference = ( Params.EchoTimes(2) - Params.EchoTimes(1) ) ; 
Params.nChannels  = 8 ;
Params.nCurrents  = 2 ;

Params.currents = [ -30 30; % A11
                    -30 30; % B11
                    -30 30; % A10
                    -600 600; % A20
                    -600 600; % A21
                    -600 600; % B21
                    -800 800; % A22
                    -800 800;] ; % B22

fileIn = fopen('/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_047p/20170623UnfPrismaShimCalibrationDataDirList.txt') ;
Params.dataLoadDirectories = textscan( fileIn, '%s' ) ;
Params.dataLoadDirectories = Params.dataLoadDirectories{1}; % not sure why this is necessary?
fclose(fileIn);
% NB: zero-current acquisition assumed to be the first 2 directories (mag, phase)

Params.filenameSave = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/ShimReferenceMaps_IUGM_Prisma_fit_20180325' ;

Params.Filtering.isFiltering  = true ;
Mag                           = MaRdI( Params.dataLoadDirectories{1} ) ;
voxelSize                     = Mag.getvoxelsize() ;
Params.Filtering.filterRadius = 2*voxelSize(1) ;

Params.reliabilityMask = Mag.img > 0.01 ; % region of reliable SNR for unwrapping
Params.reliabilityMask(1:16,:,:) = 0; % artifacts observed around the 'posterior' edge of the phantom (partial volume?)
Params.reliabilityMask(end-15:end,:,:) = 0; % artifacts observed around the 'posterior' edge of the phantom (partial volume?)
Params.reliabilityMask(:,1:4,:) = 0; 
Params.reliabilityMask(:,30:end,:) = 0; 

Params.Extension.isExtending = true ; % harmonic field extrapolation 
Params.Extension.voxelSize = voxelSize ;
Params.Extension.radius     = 6 ;

Params.unwrapper = 'AbdulRahman_2007' ;        

end
% =========================================================================
end

% =========================================================================
% =========================================================================
% =========================================================================
end




% % =========================================================================
% % function [] = acquirecalibrationdata20170323( )
%
% baselineDacValues.x    = 98.13 ;
% baselineDacValues.y    = -1172.15 ;
% baselineDacValues.z    = -1495.16 ; % note: dicom filenames were mislabeled as  _A01 instead of _A10
% baselineDacValues.z2   = 4.61 ;
% baselineDacValues.zx   = 0.47 ;
% baselineDacValues.zy   = 5.15 ;
% baselineDacValues.x2y2 = -63.91 ;
% baselineDacValues.xy   = -21.38 ;
%
%
%
% dacValues = [ baselineDacValues.x, (baselineDacValues.x - 30.0)   , (baselineDacValues.x + 30.0) ;
%               baselineDacValues.y, (baselineDacValues.y - 30.0)   , (baselineDacValues.y + 30.0) ;
%               baselineDacValues.z, (baselineDacValues.z - 30.0)   , (baselineDacValues.z + 30.0) ;
%               baselineDacValues.z2, (baselineDacValues.z2 - 600.0) , (baselineDacValues.z2 + 600.0) ;
%               baselineDacValues.zx, (baselineDacValues.zx - 600.0)  , (baselineDacValues.zx + 600.0) ;
%               baselineDacValues.zy, (baselineDacValues.zy - 600.0)  , (baselineDacValues.zy + 600.0) ;
%               baselineDacValues.x2y2, (baselineDacValues.x2y2 - 800.0)  , (baselineDacValues.x2y2 + 800.0) ;
%               baselineDacValues.xy, (baselineDacValues.xy - 800.0)  , (baselineDacValues.xy + 800.0) ; ] ;
