classdef ShimCalAcdc < ShimCal 
%SHIMCALC - Shim Calibration 
%
% .......
% 
% Usage
%
% Shim = ShimCal( )
%
%   Shim contains fields
%
%       .img
%           Shim reference maps (4d array) [units: Hz/mV]
%
%       .Hdr
%           Info re: calibration data
%           (e.g. Hdr.MaskingImage defines the spatial support of the ref maps)
%
%
% =========================================================================
% Notes
%
% Part of series of classes pertaining to shimming:
%
%    Tracking
%    ShimCal
%    ShimCom
%    ShimOpt
%    ShimSpecs
%    ShimUse
%    ShimTest 
%     
% 
% ShimCal is a MaRdI subclass [ShimCal < MaRdI]
%
% =========================================================================
% Updated::20171101::ryan.topfer@polymtl.ca
% =========================================================================

properties
    % img ;
    % Hdr ;
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimCalAcdc( Params )
%SHIMCAL - Shim calibration (i.e. reference maps) 

Shim.img = [] ;
Shim.Hdr = [] ;

if nargin < 1
    % Params = ShimCalAcdc.declarecalibrationparameters20171018( ) ;
    % Params = ShimCalAcdc.declarecalibrationparameters20171024( ) ;
    Params = ShimCalAcdc.declarecalibrationparameters20171107( ) ;
end

% -------
% map dB/dV
% [Shim.img, Shim.Hdr] = Shim.mapdbdv( Params ) ;

% -------
% map dB/dI
[Shim.img, Shim.Hdr] = Shim.mapdbdi( Params ) ;

save( Params.filenameSave ,'Shim' );

end
% =========================================================================
function [ img, Hdr ] = mapdbdv( Shim, Params )
%MAPDBDV
% 
% [ img, Hdr ] = mapdbdv( Shim, Params  ) 
% map dB/dV --- change in field [Hz] per unit voltage [mV]
% 
% Params
%   .reliabilityMask 
%       binary array indicating where phase unwrapping should take place 
img = [] ;
Hdr = [] ;

Params.DataDir = [] ;

Mag     = MaRdI( Params.dataLoadDirectories{1} ) ;
dBdVRaw = zeros( [Mag.getgridsize Params.nChannels] ) ;

Params.isFilteringField = false ;
Params.mask = Params.reliabilityMask ;

for iChannel = 1 : Params.nChannels 

    for iVoltage = 1 : (Params.nVoltages)

        iImg  = 4*iChannel + 2*(iVoltage-1) ;

        Params.DataDir.mag( iVoltage )   = Params.dataLoadDirectories(iImg-1) ;
        Params.DataDir.phase( iVoltage ) = Params.dataLoadDirectories(iImg) ;

        % -------
        % PROCESS GRE FIELD MAPS
        % -------
        ImgArray = cell( 1, 1 ) ;

        ImgArray{1,1}  = MaRdI( Params.DataDir.mag{ iVoltage }  ) ;
        ImgArray{1,2}  = MaRdI( Params.DataDir.phase{ iVoltage }  ) ;

        [Field,Extras] = FieldEval.mapfield( ImgArray, Params ) ;

        fieldMaps( :,:,:, iVoltage ) = Field.img ;

    end

    disp(['Channel ' num2str(iChannel) ' of ' num2str(Params.nChannels) ] )        
    
    % sends voltages tested to ShimCal.mapdbdi instead of currents but this doesn't matter,
    % we just end up with reference maps in terms of field shift per millivolt instead of per ampere:
    dBdVRaw(:,:,:, iChannel) = ShimCal.mapdbdi( fieldMaps, Params.voltages(iChannel,:), Params.reliabilityMask, Params ) ;  

end 

% -------
% filter the dB/dI maps
if Params.Filtering.isFiltering

    disp(['Filtering dB/dV maps...'] )
    Params.filterRadius = Params.Filtering.filterRadius ;

    dBdVFiltered = zeros( size( dBdIRaw ) ) ;

    Tmp = Field.copy();

    for iChannel = 1 : Params.nChannels 

        disp(['iChannel ' num2str(iChannel)] )

        Tmp.img = dBdVRaw(:,:,:, iChannel) ;
        [~,TmpFiltered] = Tmp.extractharmonicfield( Params ) ;
        dBdVFiltered(:,:,:,iChannel) = TmpFiltered.img ;

    end

    Hdr = TmpFiltered.Hdr ;
    img = dBdVFiltered ;

    if Params.Extension.isExtending

        tic
        disp(['Performing harmonic extension...']) ;
        disp(['Channel 1 of ' num2str(Params.nChannels) ] ) ;

        gridSizeImg = size( Params.reliabilityMask ) ;
        maskFov     = ones( gridSizeImg )  ; % extended spatial support

        [ tmp, A, M ] = extendharmonicfield( img(:,:,:, 1) , maskFov, Hdr.MaskingImage, Params.Extension ) ;

        for iChannel = 1 : Params.nChannels 
            disp(['Channel ' num2str(iChannel) ' of ' num2str(Params.nChannels) ] ) ;

            reducedField         = Hdr.MaskingImage .* img(:,:,:, iChannel) ;
            extendedField        = reshape( M'*A*reducedField(:), gridSizeImg ) ;
            img(:,:,:, iChannel) = extendedField + reducedField ;
        end

        Hdr.MaskingImage = (sum(abs(img),4))~=0 ; % extended spatial support

        toc

    end

else
    Hdr = ImgArray{1,2}.Hdr ;
    img = dBdVRaw ;
end

end
% =========================================================================
function [ img, Hdr ] = mapdbdi( Shim, Params )
%MAPDBDV
% 
% [ img, Hdr ] = mapdbdv( Shim, Params  ) 
% map dB/dI --- change in field [Hz] per unit current (A)
% 
% Params
%   .reliabilityMask 
%       binary array indicating where phase unwrapping should take place 
img = [] ;
Hdr = [] ;

Params.DataDir = [] ;

Mag     = MaRdI( Params.dataLoadDirectories{1} ) ;
dBdIRaw = zeros( [Mag.getgridsize Params.nChannels] ) ;

Params.isFilteringField = false ;
Params.mask = Params.reliabilityMask ;
dbstop in ShimCalAcdc at 200
for iChannel = 1 : Params.nChannels 

    for iCurrent = 1 : (Params.nCurrents)

        iImg  = 4*iChannel + 2*(iCurrent-1) ;

        Params.DataDir.mag( iCurrent )   = Params.dataLoadDirectories(iImg-1) ;
        Params.DataDir.phase( iCurrent ) = Params.dataLoadDirectories(iImg) ;

        % -------
        % PROCESS GRE FIELD MAPS
        % -------
        ImgArray = cell( 1, 1 ) ;

        ImgArray{1,1}  = MaRdI( Params.DataDir.mag{ iCurrent }  ) ;
        ImgArray{1,2}  = MaRdI( Params.DataDir.phase{ iCurrent }  ) ;
        
        [Field,Extras] = FieldEval.mapfield( ImgArray, Params ) ;

        fieldMaps( :,:,:, iCurrent ) = Field.img ;

    end

    disp(['Channel ' num2str(iChannel) ' of ' num2str(Params.nChannels) ] )        
    
    dBdIRaw(:,:,:, iChannel) = ShimCal.mapdbdi( fieldMaps, Params.currents(iChannel,:), Params.reliabilityMask, Params ) ;  

end 

% -------
% filter the dB/dI maps
if Params.Filtering.isFiltering

    disp(['Filtering dB/dI maps...'] )
    Params.filterRadius = Params.Filtering.filterRadius ;

    dBdIFiltered = zeros( size( dBdIRaw ) ) ;

    Tmp = Field.copy();

    for iChannel = 1 : Params.nChannels 

        disp(['iChannel ' num2str(iChannel)] )

        Tmp.img = dBdIRaw(:,:,:, iChannel) ;
        [~,TmpFiltered] = Tmp.extractharmonicfield( Params ) ;
        dBdIFiltered(:,:,:,iChannel) = TmpFiltered.img ;
        

    end

    Hdr = TmpFiltered.Hdr ;
    img = dBdIFiltered ;

    if Params.Extension.isExtending

        tic
        disp(['Performing harmonic extension...']) ;
        disp(['Channel 1 of ' num2str(Params.nChannels) ] ) ;

        gridSizeImg = size( Params.reliabilityMask ) ;
        maskFov     = ones( gridSizeImg )  ; % extended spatial support

        [ ~, A, M ] = extendharmonicfield( img(:,:,:, 1) , maskFov, Hdr.MaskingImage, Params.Extension ) ;

        for iChannel = 1 : Params.nChannels 
            disp(['Channel ' num2str(iChannel) ' of ' num2str(Params.nChannels) ] ) ;

            reducedField         = Hdr.MaskingImage .* img(:,:,:, iChannel) ;
            extendedField        = reshape( M'*A*reducedField(:), gridSizeImg ) ;
            img(:,:,:, iChannel) = extendedField + reducedField ;
        end

        Hdr.MaskingImage = (sum(abs(img),4))~=0 ; % extended spatial support

        toc

    end

else
    Hdr = ImgArray{1,2}.Hdr ;
    img = dBdIRaw ;
end

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
end

% =========================================================================
% =========================================================================
% =========================================================================
end
