classdef ShimCal_Acdc < ShimCal 
%SHIMCAL_ACDC - Shim Calibration for the 8-channel "Ac/Dc" coil.
% 
% ShimCal_Acdc is a ShimCal subclass [ShimCalAcdc < ShimCal]
%
% =========================================================================
% Updated::20180328::ryan.topfer@polymtl.ca
% =========================================================================

properties
    % img ;
    % Hdr ;
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimCal_Acdc( Params )
%SHIMCAL - Shim calibration (i.e. reference maps) 

Shim.img = [] ;
Shim.Hdr = [] ;

if nargin < 1
    % Params = ShimCal_Acdc.declarecalibrationparameters20171018( ) ;
    % Params = ShimCal_Acdc.declarecalibrationparameters20171024( ) ;
    % Params = ShimCal_Acdc.declarecalibrationparameters20171107( ) ;
    % Params = ShimCal_Acdc.declarecalibrationparameters20171209( ) ;
    % Params = ShimCal_Acdc.declarecalibrationparameters20180322( ) ;
    Params = ShimCal_Acdc.declarecalibrationparameters20180326( ) ;
end

% -------
% map dB/dV
% [Shim.img, Shim.Hdr] = Shim.mapdbdv( Params ) ;

% -------
% map dB/dI
[Shim.img, Shim.Hdr] = ShimCal.mapdbdi_allchannels( Params ) ;

save( Params.filenameSave ,'Shim' );

end
% =========================================================================
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================
function Params = declarecalibrationparameters20171018( )
%DECLARECALIBRATIONPARAMETERS20171018 
% 
% Initializes parameters for shim reference map construction (i.e. shim calibration)

fprintf('##### \n WARNING: Ignoring EchoTime DICOM header when normalizing phase to field. \n\n')
Params.isHardCodingEchoTime = true ;
Params.EchoTimes  = [4.92 2.46] ;
Params.echoTimeDifference = ( Params.EchoTimes(2) - Params.EchoTimes(1) ) ; 
Params.nChannels  = 8 ;
Params.nVoltages  = 2 ;

Params.voltages = [100 300; %ch1
                    50 100; %
                    50 100; %
                    50 100; %...
                    50 100; %
                    50 100; %
                    50 100; % 
                    50 100;] ; %ch8

fileIn = fopen('/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_07p/20171018AcdcShimCalibrationDataDirList.txt') ;
Params.dataLoadDirectories = textscan( fileIn, '%s' ) ;
Params.dataLoadDirectories = Params.dataLoadDirectories{1}; % not sure why this is necessary?
fclose(fileIn);
% NB: zero-current acquisition assumed to be the first 2 directories (mag, phase)

Params.filenameSave = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/AcdcReferenceMaps20171018' ;

Params.Filtering.isFiltering  = true ;
Mag                           = MaRdI( Params.dataLoadDirectories{1} ) ;
voxelSize                     = Mag.getvoxelsize() ;
Params.Filtering.filterRadius = 2*voxelSize(1) ;

Params.reliabilityMask = Mag.img > 0.01 ; % region of reliable SNR for unwrapping
Params.reliabilityMask(:,48:end,:) = 0; % artifacts observed around the 'posterior' edge of the phantom (partial volume?)

Params.Extension.isExtending = false ; % harmonic field extrapolation 

end
% =========================================================================
function Params = declarecalibrationparameters20171024( )
%DECLARECALIBRATIONPARAMETERS20171024
% 
% Initializes parameters for shim reference map construction (i.e. shim calibration)
%  
%  NOTE
% Difference from declarecalibrationparameters20171018 :
%   replacing the 100 mV map for channel 8 with the 0 mV-all channel map
%   since the 100 mV ch8 map was observed to have some strange artifacts.

fprintf('##### \n WARNING: Ignoring EchoTime DICOM header when normalizing phase to field. \n\n')
Params.isHardCodingEchoTime = true ;
Params.EchoTimes  = [4.92 2.46] ;
Params.echoTimeDifference = ( Params.EchoTimes(2) - Params.EchoTimes(1) ) ; 
Params.nChannels  = 8 ;
Params.nVoltages  = 2 ;

Params.voltages = [100 300; %ch1
                    50 100; %
                    50 100; %
                    50 100; %...
                    50 100; %
                    50 100; %
                    50 100; % 
                    50 0;] ; %ch8 -- RT::20171024 see above note 

fileIn = fopen('/Users/ryan/Projects/Shimming/Static/Calibration/Data/acdc_07p/20171018AcdcShimCalibrationDataDirList.txt') ;
Params.dataLoadDirectories = textscan( fileIn, '%s' ) ;
Params.dataLoadDirectories = Params.dataLoadDirectories{1}; % not sure why this is necessary?
fclose(fileIn);
% NB: zero-current acquisition = the first 2 directories (mag, phase)

% RT::20171024 see above note 
Params.dataLoadDirectories(end-1) = Params.dataLoadDirectories(1) ;
Params.dataLoadDirectories(end)   = Params.dataLoadDirectories(2) ;


Params.filenameSave = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/AcdcReferenceMaps20171024' ;

Params.Filtering.isFiltering  = true ;
Mag                           = MaRdI( Params.dataLoadDirectories{1} ) ;
voxelSize                     = Mag.getvoxelsize() ;
Params.Filtering.filterRadius = 2*voxelSize(1) ;

Params.reliabilityMask = Mag.img > 0.01 ; % region of reliable SNR for unwrapping
Params.reliabilityMask(:,48:end,:) = 0; % artifacts observed around the 'posterior' edge of the phantom (partial volume?)

Params.Extension.isExtending = true ; % harmonic field extrapolation 
Params.Extension.voxelSize = voxelSize ;
Params.Extension.radius     = 6 ;

end
% =========================================================================
function Params = declarecalibrationparameters20171107( )
%DECLARECALIBRATIONPARAMETERS20171107
% 
% Initializes parameters for shim reference map construction (i.e. shim calibration)
%  
%  NOTE
% Difference from declarecalibrationparameters20171025 :
%
%   -DAC reference voltage has been changed from 0 to 1.25 V
%    
%   -Chokes added to cabling (board end only) to prevent gradient interaction
%   
%   -Nibardo changed chokes (lower inductance, lower capacitance) after 
%   Skype conversation with Jason (20171102)

fprintf('##### \n WARNING: Ignoring EchoTime DICOM header when normalizing phase to field. \n\n')
Params.isHardCodingEchoTime = true ;
Params.EchoTimes  = [4.92 2.46] ;
Params.echoTimeDifference = ( Params.EchoTimes(2) - Params.EchoTimes(1) ) ; 
Params.nChannels  = 8 ;
Params.nCurrents  = 2 ;

Params.currents = [-0.4 0.4; %ch1, [units: A]
                   -0.4 0.4; %
                   -0.4 0.4; %
                   -0.4 0.4; %...
                   -0.4 0.4; %
                   -0.4 0.4; %
                   -0.4 0.4; % 
                   -0.4 0.4;] ; %ch8

% data from
fileIn = fopen('/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_11p/20171106AcdcShimCalibrationDataDirList.txt') ;
Params.dataLoadDirectories = textscan( fileIn, '%s' ) ;
Params.dataLoadDirectories = Params.dataLoadDirectories{1}; % not sure why this is necessary?
fclose(fileIn);
% NB: zero-current acquisition = the first 2 directories (mag, phase)

Params.filenameSave = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/AcdcReferenceMaps20171107' ;

Params.Filtering.isFiltering  = true ;
Mag                           = MaRdI( Params.dataLoadDirectories{1} ) ;
voxelSize                     = Mag.getvoxelsize() ;
Params.Filtering.filterRadius = 2*voxelSize(1) ;

Params.reliabilityMask = Mag.img > 0.01 ; % region of reliable SNR for unwrapping
Params.reliabilityMask(:,47:end,:) = 0; % artifacts observed around the 'posterior' edge of the phantom (partial volume?)

Params.Extension.isExtending = true ; % harmonic field extrapolation 
Params.Extension.voxelSize   = voxelSize ;
Params.Extension.expansionOrder = 2 ;
Params.Extension.radius     = 6 ;

end
% =========================================================================
function Params = declarecalibrationparameters20171209( )
%DECLARECALIBRATIONPARAMETERS20171209
% 
% Initializes parameters for shim reference map construction (i.e. shim calibration)
%  
%  NOTE
% Difference from declarecalibrationparameters20171107 :
%
% -resolving issue of missing negative sign in ref. maps?
% -using complex difference between +/- current field maps


fprintf('##### \n WARNING: Ignoring EchoTime DICOM header when normalizing phase to field. \n\n')
Params.isHardCodingEchoTime = true ;
Params.EchoTimes  = [4.92 2.46] ;
Params.echoTimeDifference = ( Params.EchoTimes(2) - Params.EchoTimes(1) ) ; 
Params.nChannels  = 8 ;
Params.nCurrents  = 2 ;

Params.currents = [-0.4 0.4; %ch1, [units: A]
                   -0.4 0.4; %
                   -0.4 0.4; %
                   -0.4 0.4; %...
                   -0.4 0.4; %
                   -0.4 0.4; %
                   -0.4 0.4; % 
                   -0.4 0.4;] ; %ch8

% data from
fileIn = fopen('/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_11p/20171106AcdcShimCalibrationDataDirList.txt') ;
Params.dataLoadDirectories = textscan( fileIn, '%s' ) ;
Params.dataLoadDirectories = Params.dataLoadDirectories{1}; % not sure why this is necessary?
fclose(fileIn);
% NB: zero-current acquisition = the first 2 directories (mag, phase)

Params.filenameSave = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/AcdcReferenceMaps20171209' ;

Params.Filtering.isFiltering  = true ;
Mag                           = MaRdI( Params.dataLoadDirectories{1} ) ;
voxelSize                     = Mag.getvoxelsize() ;
Params.Filtering.filterRadius = 2*voxelSize(1) ;

Params.reliabilityMask = Mag.img/max(Mag.img(:)) > 0.01 ; % region of reliable SNR for unwrapping
% Params.reliabilityMask(:,47:end,:) = 0; % artifacts observed around the 'posterior' edge of the phantom (partial volume?)

Params.Extension.isExtending = true ; % harmonic field extrapolation 
Params.Extension.voxelSize   = voxelSize ;
Params.Extension.expansionOrder = 2 ;
Params.Extension.radius     = 6 ;

end
% =========================================================================
function Params = declarecalibrationparameters20180322( )
%DECLARECALIBRATIONPARAMETERS20180322
% 
% Initializes parameters for shim reference map construction (i.e. shim calibration)
%  
%  NOTE
% Difference from declarecalibrationparameters20171209 :
% 
% -Numerous changes to the coil electronics (e.g. new power source)
% -Different GRE acquisition parameters 
% -Different phantom used (long ~40+ cm Siemens torpedo)


fprintf('##### \n WARNING: Ignoring EchoTime DICOM header when normalizing phase to field. \n\n')
Params.isHardCodingEchoTime = true ;
Params.EchoTimes  = [ 4.92 7.38 ] ;
Params.echoTimeDifference = ( Params.EchoTimes(2) - Params.EchoTimes(1) ) ; 
Params.nChannels  = 8 ;
Params.nCurrents  = 2 ;

Params.currents = [-0.4 0.4; %ch1, [units: A]
                   -0.4 0.4; %
                   -0.4 0.4; %
                   -0.4 0.4; %...
                   -0.4 0.4; %
                   -0.4 0.4; %
                   -0.4 0.4; % 
                   -0.4 0.4;] ; %ch8

% data from
fileIn = fopen('/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/20180320AcdcShimCalibrationDataDirList.txt') ;
Params.dataLoadDirectories = textscan( fileIn, '%s' ) ;
Params.dataLoadDirectories = Params.dataLoadDirectories{1}; % not sure why this is necessary?
fclose(fileIn);
% NB: zero-current acquisition = the first 2 directories (mag, phase)

Params.filenameSave = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/ShimReferenceMaps_Acdc_20180322' ;

Params.Filtering.isFiltering  = true ;
Mag                           = MaRdI( Params.dataLoadDirectories{1} ) ;
voxelSize                     = Mag.getvoxelsize() ;
Params.Filtering.filterRadius = 2*voxelSize(1) ;

Params.reliabilityMask = Mag.img/max(Mag.img(:)) > 0.01 ; % region of reliable SNR for unwrapping

Params.reliabilityMask(:,1:5,:)    = 0 ; % based on visual inspection of magnitude (there is aliasing in phase)
Params.reliabilityMask(:,71:end,:) = 0 ; 

Params.Extension.isExtending = true ; % harmonic field extrapolation 
Params.Extension.voxelSize   = voxelSize ;
Params.Extension.expansionOrder = 2 ;
Params.Extension.radius     = 6 ;

end
% =========================================================================
function Params = declarecalibrationparameters20180326( )
%DECLARECALIBRATIONPARAMETERS20180326
% 
% Initializes parameters for shim reference map construction (i.e. shim calibration)
%  
%  NOTE
% Difference from declarecalibrationparameters20180322 :
% 
% Same acquistion data but DIS3D (distortion correction) has been applied (+before unwrapping)


fprintf('##### \n WARNING: Ignoring EchoTime DICOM header when normalizing phase to field. \n\n')
Params.isHardCodingEchoTime = true ;
Params.EchoTimes  = [ 4.92 7.38 ] ;
Params.echoTimeDifference = ( Params.EchoTimes(2) - Params.EchoTimes(1) ) ; 
Params.nChannels  = 8 ;
Params.nCurrents  = 2 ;

Params.currents = [-0.4 0.4; %ch1, [units: A]
                   -0.4 0.4; %
                   -0.4 0.4; %
                   -0.4 0.4; %...
                   -0.4 0.4; %
                   -0.4 0.4; %
                   -0.4 0.4; % 
                   -0.4 0.4;] ; %ch8

% data from
fileIn = fopen('/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/20180326AcdcShimCalibrationDataDirList.txt') ;
Params.dataLoadDirectories = textscan( fileIn, '%s' ) ;
Params.dataLoadDirectories = Params.dataLoadDirectories{1}; % not sure why this is necessary?
fclose(fileIn);
% NB: zero-current acquisition = the first 2 directories (mag, phase)

% load the uncorrected images as well
fileIn = fopen('/Users/ryan/Projects/Shimming/Acdc/Calibration/data/acdc_19pb/20180320AcdcShimCalibrationDataDirList.txt') ;
Params.dataLoadDirectoriesNd = textscan( fileIn, '%s' ) ;
Params.dataLoadDirectoriesNd = Params.dataLoadDirectoriesNd{1};
fclose(fileIn);

Params.filenameSave = '/Users/ryan/Projects/Shimming/Acdc/Calibration/data/ShimReferenceMaps_Acdc_20180326' ;

Params.Filtering.isFiltering  = true ;
Mag                           = MaRdI( Params.dataLoadDirectories{1} ) ;
voxelSize                     = Mag.getvoxelsize() ;
Params.Filtering.filterRadius = 2*voxelSize(1) ;

Params.reliabilityMask = Mag.img/max(Mag.img(:)) > 0.01 ; % region of reliable SNR for unwrapping

Params.reliabilityMask(:,1:5,:)    = 0 ; % based on visual inspection of magnitude (there is aliasing in phase-encore dir)
Params.reliabilityMask(:,71:end,:) = 0 ; 

Params.Extension.isExtending = true ; % harmonic field extrapolation 
Params.Extension.voxelSize   = voxelSize ;
Params.Extension.expansionOrder = 2 ;
Params.Extension.radius     = 6 ;

end
% =========================================================================
end

% =========================================================================
% =========================================================================
% =========================================================================
end
