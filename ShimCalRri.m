classdef ShimCalRri < ShimCal 
%SHIMCAL - Shim Calibration 
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
%           Shim reference maps (4d array) [units: Hz/ampere]
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
% Updated::20170412::ryan.topfer@polymtl.ca
% =========================================================================

properties
    % img ;
    % Hdr ;
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimCalRri( Params )
%SHIMCAL - Shim calibration (i.e. reference maps) 

Shim.img = [] ;
Shim.Hdr = [] ;

if nargin < 1
    % Params = ShimCal.declarecalibrationparameters20160809() ;
    % Params = ShimCal.declarecalibrationparameters20161003() ;
    % Params = ShimCal.declarecalibrationparameters20161004() ;
    % Params = ShimCal.declarecalibrationparameters20161007() ;
    % Params = ShimCalRri.declarecalibrationparameters20170324() ;
    % Params = ShimCalRri.declarecalibrationparameters20170410() ;
    % Params = ShimCalRri.declarecalibrationparameters20170418() ;
    Params = ShimCalRri.declarecalibrationparameters20170706() ;
end

% -------
% map dB/dI
% [Shim.img, Shim.Hdr] = ShimCalRri.mapdbdi20170410( Params ) ;
[Shim.img, Shim.Hdr] = ShimCalRri.mapdbdi20170706(  ) ;

save( Params.filenameSave ,'Shim');


end
% =========================================================================
% =========================================================================
% =========================================================================
end
% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================

function [] = acquirecalibrationdata( )
%ACQUIRECALIBRATIONDATA 
%
%  Re-acquiring calibration data on "waterproof phantom" (28 L water) 
%  Experiment performed: 20160809
%
%  Notes
%
%   **This is a demo script to accompany reference field map acquisition 
%   ( i.e. it's not meant to be run as a standalone function).**
%
%   Halfway through the acquisitions I realised I should record the channel
%   output currents rather than simply the assigned test currents (ideally,
%   these would be identical, but owing to the finite (16 bit) precision of 
%   the DACs, and perhaps other factors, they are close, but not identical).
%
%   Therefore I created a variable called measuredCurrents to keep track of 
%   subsequent recordings, and printed the terminal display output to a file
%   ( ~/Projects/Shimming/Static/Calibration/Data/20160809_currentsLog.pdf)
%   Unfortunately, "infinte scrollback" was disabled on iterm2, so, the 
%   measured currents for the first 7.5 channels were lost. (7.5 because 
%   the *positive* measured current for the 8th channel was recovered.)
%
%   To correct for this, I'm substituting the average measured currents for
%   the remaining channels into the absent values for these first 7.5 
%   channels. (See the last section below.)
%

% -------
% Connect to shim amplifier 
ShimsC = ShimComRri;
ShimsC.getsystemheartbeat 

channelsToBankKey = ShimsC.getchanneltobankkey ;

testCurrents = [-0.5 0.5]; % [units : A]

iChannel = 0;

measuredCurrents = zeros(24, length(testCurrents) ) ;
% -------
% Perform routine for iChannel = 1 : 24
iChannel = iChannel + 1 

ShimsC.resetallshims()
iTestCurrent = 1 ;

ShimsC.setandloadshim( channelsToBankKey(iChannel, 2), channelsToBankKey(iChannel,3), testCurrents(iTestCurrent) ) ;

pause(1)
ChannelOutput = ShimsC.getchanneloutput(channelsToBankKey(iChannel,2), channelsToBankKey(iChannel,3) )

measuredCurrents( iChannel, iTestCurrent ) = ChannelOutput.current ;


ShimsC.resetallshims()
iTestCurrent = 2 ;

ShimsC.setandloadshim( channelsToBankKey(iChannel, 2), channelsToBankKey(iChannel,3), testCurrents(iTestCurrent) ) ;

pause(1)
ChannelOutput = ShimsC.getchanneloutput(channelsToBankKey(iChannel,2), channelsToBankKey(iChannel,3) )

measuredCurrents( iChannel, iTestCurrent ) = ChannelOutput.current ;
% ------- 
% Account for missing measured currents & save as .txt (see 'Notes' above)

measuredCurrents(1:8,1) = mean( measuredCurrents(9:end,1) ) ;
measuredCurrents(1:7,2) = mean( measuredCurrents(8:end,2) ) ;

filenameSave = ...
    '~/Projects/Shimming/Static/Calibration/Data/20160809_measuredCurrents.txt' ;

fid = fopen( filenameSave, 'w' ) ;
fprintf( fid, '%f %f \n', measuredCurrents' ) ;
fclose( fid ) ;

end
% =========================================================================
function [] = acquirecalibrationdata20170410( )
%ACQUIRECALIBRATIONDATA20170410 
%
%  Re-acquiring calibration data on "waterproof phantom" (28 L water) 
%  (i.e. shim reference maps post-Prisma upgrade).
%
%  Experiment performed: 20170323, and again: 20170410
%
%  Notes
%
%   **This is a demo script to accompany reference field map acquisition 
%   ( i.e. it's not meant to be run as a standalone function).**
%

% filenameSave = ...
%     '~/Projects/Shimming/Static/Calibration/Data/20170323_measuredCurrents.txt' ;

filenameSave = ...
    '~/Projects/Shimming/Static/Calibration/Data/20170410_measuredCurrents.txt';



% -------
% Connect to shim amplifier 
ShimsC = ShimComRri();
ShimsC.getsystemheartbeat() ;

channelsToBankKey = ShimsC.getchanneltobankkey ;

testCurrents  = [-0.5 0.5]; % [units : A]
nTestCurrents = length(testCurrents) ;

iChannel  = 0;
nChannels = size(channelsToBankKey, 1);
measuredCurrents = zeros( nChannels, length(testCurrents) ) ;

% -------
% Perform routine for iChannel = 1 : 24
iChannel = iChannel + 1 

ShimsC.resetallshims()
iTestCurrent = 1 ;

ShimsC.setandloadshim( channelsToBankKey(iChannel, 2), channelsToBankKey(iChannel,3), testCurrents(iTestCurrent) ) ;

pause(0.5)
ChannelOutput = ShimsC.getchanneloutput(channelsToBankKey(iChannel,2), channelsToBankKey(iChannel,3) )

measuredCurrents( iChannel, iTestCurrent ) = ChannelOutput.current ;


ShimsC.resetallshims()
iTestCurrent = 2 ;

ShimsC.setandloadshim( channelsToBankKey(iChannel, 2), channelsToBankKey(iChannel,3), testCurrents(iTestCurrent) ) ;

pause(0.5)
ChannelOutput = ShimsC.getchanneloutput(channelsToBankKey(iChannel,2), channelsToBankKey(iChannel,3) )

measuredCurrents( iChannel, iTestCurrent ) = ChannelOutput.current ;


% -------
% save measured currents
ShimsC.resetallshims()

fid = fopen( filenameSave, 'w' ) ;
fprintf( fid, '%f %f \n', measuredCurrents' ) ;
fclose( fid ) ;

end
% =========================================================================
function [] = comparereferencemaps( )
%COMPAREREFERENCEMAPS 
%

Params.pathToShimReferenceMaps = '~/Projects/Shimming/RRI/data/SpineShimReferenceMaps.mat' ; 
ShimsOld = ShimOpt(Params) ;

Params.pathToShimReferenceMaps = '~/Projects/Shimming/RRI/data/SpineShimReferenceMaps20160809.mat' ;
ShimsNew = ShimOpt(Params) ;

[X,Y,Z] = ShimsOld.getvoxelpositions() ;
ShimsNew.resliceimg( X, Y, Z ) ;

NiiOptions.filename = './old_ref_maps' ;
nii( ShimsOld.img, NiiOptions ) ;

NiiOptions.filename = './new_ref_maps' ;
nii( ShimsNew.img, NiiOptions ) ;
end
% =========================================================================
function Params = declarecalibrationparameters20160809( )
%DECLARECALIBRATIONPARAMETERS20160809 
% 
% Initializes parameters for shim reference map construction (i.e. shim calibration)


fprintf('##### \n WARNING: Ignoring EchoTime DICOM header when normalizing phase to field. \n\n')
Params.isHardCodingEchoTime = true ;
Params.EchoTimes  = [4.92 7.38] ;
Params.nChannels  = 24 ;
Params.nCurrents  = 3 ;

fileIn = fopen('/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_012p/20160809_measuredCurrents.txt') ;
Params.currents = fscanf( fileIn, '%f', [2 Params.nChannels ] )' ;
fclose(fileIn);

fileIn = fopen('/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_012p/shimcalibrationdata_dirlist_20160809.txt') ;
Params.dataLoadDirectories = textscan( fileIn, '%s' ) ;
Params.dataLoadDirectories = Params.dataLoadDirectories{1}; % not sure why this is necessary?
fclose(fileIn);
% NB: zero-current acquisition assumed to be the first 2 directories (mag, phase)

Params.filenameSave = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/SpineShimReferenceMaps20160809' ;

Params.Filtering.isFiltering  = true ;
Mag                           = MaRdI( Params.dataLoadDirectories{1} ) ;
voxelSize                     = Mag.getvoxelsize() ;
Params.Filtering.filterRadius = 2*voxelSize(1) ;

Params.Extension.isExtending = false ; % harmonic field extrapolation 
% Define data support (used for all field maps)
Params.limitsForUnwrapping = [2 90; 2 28; 1 80] ; 

end
% =========================================================================
function Params = declarecalibrationparameters20161003( )
    
Params = ShimCal.declarecalibrationparameters20160809( );
Params.Filtering.isFiltering  = false ;

Params.filenameSave = ...
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/SpineShimReferenceMaps20161003_nofiltering' ;

end
% =========================================================================
function Params = declarecalibrationparameters20161004( )
    
Params = ShimCal.declarecalibrationparameters20160809( );

% Define data support (used for all field maps)
Params.limitsForUnwrapping = [11 80; 3 26; 12 69] ; 

Params.filenameSave = ...
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/SpineShimReferenceMaps20161004' ;

end
% =========================================================================
function Params = declarecalibrationparameters20161007( )
    
Params = ShimCal.declarecalibrationparameters20160809( );

% Define data support (used for all field maps)
Params.limitsForUnwrapping = [11 80; 3 26; 12 69] ; 

Params.Extension.isExtending = true ; % harmonic field extrapolation 
Params.Extension.radius      = 10 ;
Params.Extension.isDisplayingProgress = true ;
Params.Extension.limitsForExtension = [10 80; 1 28; 13 70] ;

Params.filenameSave = ...
    '/Users/ryan/Projects/Shimming/Static/Calibration/Data/SpineShimReferenceMaps20161007' ;

end
% =========================================================================
function Params = declarecalibrationparameters20170324( )

fprintf('##### \n WARNING: Ignoring EchoTime DICOM header when normalizing phase to field. \n\n')
Params.isHardCodingEchoTime = true ;
Params.EchoTimes  = [4.92 7.38] ;
Params.nChannels  = 24 ;
Params.nCurrents  = 2 ;

fileIn = fopen('/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_031p/20170323_measuredCurrents.txt') ;
Params.currents = fscanf( fileIn, '%f', [2 Params.nChannels ] )' ;
fclose(fileIn);

fileIn = fopen('/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_031p/20170323_shimCalibrationDataDirList.txt') ;
Params.dataLoadDirectories = textscan( fileIn, '%s' ) ;
Params.dataLoadDirectories = Params.dataLoadDirectories{1}; % not sure why this is necessary?
fclose(fileIn);

% NB: zero-current acquisition assumed to be the first 2 directories (mag, phase)

Params.filenameSave = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/ShimReferenceMapsRri20170324' ;

Params.Filtering.isFiltering  = false ;
Mag                           = MaRdI( Params.dataLoadDirectories{1} ) ;
voxelSize                     = Mag.getvoxelsize() ;
Params.Filtering.filterRadius = 2*voxelSize(1) ;

Params.Extension.isExtending = false ; % harmonic field extrapolation 
% Define data support (used for all field maps)
Params.limitsForUnwrapping = [2 90; 2 28; 1 80] ; 

end
% =========================================================================
function Params = declarecalibrationparameters20170410( )

Params.filenameSave = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/ShimReferenceMapsRri20170410' ;

fprintf('##### \n WARNING: Ignoring EchoTime DICOM header when normalizing phase to field. \n\n')
Params.isHardCodingEchoTime = true ;
Params.EchoTimes  = [4.92 7.38] ;
Params.nChannels  = 24 ;
Params.nCurrents  = 2 ;

fileIn = fopen('/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/20170410_measuredCurrents.txt') ;
Params.currents = fscanf( fileIn, '%f', [2 Params.nChannels ] )' ;
fclose(fileIn);

fileIn = fopen('/Users/ryan/Projects/Shimming/Static/Calibration/Data/shim_034p/20170410_shimCalibrationDataDirList.txt') ;
Params.dataLoadDirectories = textscan( fileIn, '%s' ) ;
Params.dataLoadDirectories = Params.dataLoadDirectories{1}; % not sure why this is necessary?
Params.dataLoadDirectories = flipud( Params.dataLoadDirectories ) ; % so ch1 comes first and ch24 last...
fclose(fileIn);

Params.Filtering.isFiltering  = true ;
Mag                           = MaRdI( Params.dataLoadDirectories{1} ) ;
voxelSize                     = Mag.getvoxelsize() ;
Params.Filtering.filterRadius = voxelSize(1) ;

Params.Extension.isExtending = false ; % harmonic field extrapolation 

% Define data support (used for unwrapping all phase images)
Params.reliabilityMask = ones( Mag.getgridsize ) ;

Params.reliabilityMask(1:9,:,:)       = 0 ;
Params.reliabilityMask(end-8:end,:,:) = 0;

Params.reliabilityMask(:,1:4,:)       = 0 ;
Params.reliabilityMask(:,end,:)       = 0 ;

Params.reliabilityMask(:,:,1:10)       = 0 ;
Params.reliabilityMask(:,:,end-9:end)  = 0 ;


end
% =========================================================================
function Params = declarecalibrationparameters20170418( )

Params = ShimCalRri.declarecalibrationparameters20170410() ;

Params.filenameSave = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/ShimReferenceMapsRri20170418' ;

Mag = MaRdI( Params.dataLoadDirectories{1} ) ;

Params.Filtering.isFiltering  = true ;
Mag                           = MaRdI( Params.dataLoadDirectories{1} ) ;
voxelSize                     = Mag.getvoxelsize() ;
Params.Filtering.filterRadius = 2*voxelSize(1) ;

% same as 2017/04/10 but with extrapolation
Params.Extension.isExtending = true ; % harmonic field extrapolation 
Params.Extension.radius      = 12 ;
Params.Extension.isDisplayingProgress = true ;
Params.Extension.voxelSize = Mag.getvoxelsize() ;

% and more limited initial data support (used for unwrapping all phase images)
% since the remaining artifacts in the region defined for 20170410 corrupt 
% some of the extrapolation
Params.reliabilityMask = ones( Mag.getgridsize ) ;

Params.reliabilityMask(1:15,:,:)      = 0 ;
Params.reliabilityMask(end-8:end,:,:) = 0;

Params.reliabilityMask(:,1:4,:)       = 0 ;
Params.reliabilityMask(:,end,:)       = 0 ;

Params.reliabilityMask(:,:,1:10)      = 0 ;
Params.reliabilityMask(:,:,end-9:end) = 0 ;


end
% =========================================================================
function Params = declarecalibrationparameters20170706( )

Params.filenameSave = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/ShimReferenceMapsRri20170706' ;

% NOTE 
% Params are otherwise the same as 2017/04/10 but with the extrapolated support
% of from 2017/04/18 
%
% It seems the two in vivo experiments using the 2017/04/10 data (shim037 &
% shim038) were actually quite accurate. The problem was that the scanner
% updated the Tx-freq. for the EPI once the shims were ON. So, the shimmed
% field maps looked pretty good, whereas the EPI were garbage.
% 
% This hybrid set of reference maps is a temporary solution for today
% 2017/07/06 : shim052, 053, 054


end
% =========================================================================
% =========================================================================
function [img, Hdr] = mapdbdi20170410( Params ) 
%MAPDBDI
% 
% map dB/dI --- change in field [Hz] per unit current [A]
%
% dBdI = mapdbdi20170410( Params )
%
% Params
%   .reliabilityMask 
%       binary array indicating where phase unwrapping should take place 


Params.DataDir = [] ;

Params.echoTimeDifference = ( Params.EchoTimes(2) - Params.EchoTimes(1) ) ; 
Params.isFilteringField = false ;
Params.mask = Params.reliabilityMask ;

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

        [Field,Extras] = ShimOpt.mapfield( ImgArray, Params ) ;

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

    Tmp = ImgArray{1,2}.copy();

    tic
    for iChannel = 1 : 24

        disp(['iChannel ' num2str(iChannel)] )

        Tmp.img = dBdIRaw(:,:,:, iChannel) ;
        [~,TmpFiltered] = Tmp.extractharmonicfield( Params ) ;
        dBdIFiltered(:,:,:,iChannel) = TmpFiltered.img ;

    end
    toc

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
function [img, Hdr] = mapdbdi20170706(  ) 
%MAPDBDI20170706
%
% Combining dBdI from 20170410( ) and 20170706

Params.pathToShimReferenceMaps = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/ShimReferenceMapsRri20170410.mat' ;
Shims20170410 = ShimOptRri( Params ) ;

Params.pathToShimReferenceMaps = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/ShimReferenceMapsRri20170418.mat' ;
Shims20170418 = ShimOptRri( Params ) ;

shimSupport20170410 = Shims20170410.getshimsupport() ;
shimSupport4d       = repmat( shimSupport20170410, 1, 1, 1, 24 ) ;

% extended region courtesy of the 20170418 ref maps
img = Shims20170418.img ;
% interior region via the earlier 20170410 ref maps
img( shimSupport4d ) = Shims20170410.img( shimSupport4d ) ;

Hdr = Shims20170418.Hdr ;

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
