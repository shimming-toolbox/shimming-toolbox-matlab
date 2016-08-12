classdef ShimCal < MaRdI
%SHIMCAL
%
% Shim Calibration, a MaRdI subclass (i.e. ShimCal < MaRdI)
%
% .......
%
% Shim = ShimCal( )
%
%   Shim contains fields
%
%       .img
%           Shim reference maps
%
%       .Hdr
%           Info re: calibration data
%           (e.g. Hdr.MaskingImage defines the spatial support of the ref maps)
%
%
% =========================================================================
% Part of series of classes pertaining to shimming:
%
%     ProbeTracking
%     ShimCal
%     ShimCom
%     ShimOpt
%     ShimSpecs
%     ShimUse
%     
% =========================================================================
% Updated::20160810::ryan.topfer@polymtl.ca
% =========================================================================

properties

end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimCal(  )
%SHIMCAL - Shim calibration (i.e. reference maps) 
%

end
% =========================================================================
function Shim = calibratespineshim(  )
% -------
fprintf('##### \n WARNING: This script ignores .EchoTime DICOM header when normalizing phase to field. \n\n')

isHardCodingEchoTime = true ;

% -------
Params.EchoTimes = [4.92 7.38] ;

Params.projectDirectory  = '~/Projects/Shimming/Static/Calibration/' ;
Params.dataLoadDirectory = [Params.projectDirectory '/Data/shim_012p/'] ;
Params.tmpSaveDirectory  = [Params.projectDirectory 'Tmp/'] ;

Params.nChannels  = 24 ;
Params.nCurrents  = 3 ;

Params.limitsForUnwrapping = [2 90; 2 28; 1 80] ; 

% Params.Filtering.isFiltering = false ; % harmonVic field filtering
% Params.Filtering.filterType  = 'pdf' ;
% Params.Extension.isExtending = false ; % harmonic field extrapolation 

field = [] ; % the field maps to be fitted to current
dBdI  = [] ; % change in field per change in shim current

% -------
% Process Zero-current data
Mag   = MaRdI( [Params.dataLoadDirectory '118-gre_field_mapping_0A/echo_7.38/'] ) ;
Phase = MaRdI( [Params.dataLoadDirectory '119-gre_field_mapping_0A/echo_7.38/'] ) ;

%-----
% Define data support
maskForUnwrapping  = zeros( Phase.getgridsize( ) ) ;
maskForUnwrapping( Params.limitsForUnwrapping(1,1) : (Params.limitsForUnwrapping(1,2)), ...
                   Params.limitsForUnwrapping(2,1) : (Params.limitsForUnwrapping(2,2)), ...
                   Params.limitsForUnwrapping(3,1) : (Params.limitsForUnwrapping(3,2)) ) = 1 ;

%-----
% Unwrap phase
Phase.Hdr.MaskingImage = logical( maskForUnwrapping ) ;
Phase = Phase.unwrapphase(  ) ;

%-----
if isHardCodingEchoTime
    Phase.Hdr.EchoTime = Params.EchoTimes(2) - Params.EchoTimes(1) ;
end

Field1 = Phase.scalephasetofrequency( ) ;

Params.voxelSize = Field.getvoxelsize();    
Params.maxIterations = 4000 ;
[ localField, reducedMask, bkgrField ] = projectontodipolefields( ...
    permute( Field.img, [1 3 2] ) , ...
    permute( maskForUnwrapping, [1 3 2] ) , ...
    permute( maskForUnwrapping .* Mag.img, [1 3 2] ) , ...
    Params ) ; 

%     NiiOptions.filename = ['../Tmp/' 'total0' ];
%     nii( Field.img, NiiOptions ) ;
%     
%     NiiOptions.filename = ['../Tmp/' 'local0' ];
%     nii( reducedMask.*localField, NiiOptions ) ;
%     
%     NiiOptions.filename = ['../Tmp/' 'bkgr0' ];
%     nii( bkgrField, NiiOptions ) ;
%     
%     NiiOptions.filename = ['../Tmp/' 'percentLocal0' ];
%     nii( 100*reducedMask.*abs(localField./Field.img), NiiOptions ) ;
%
%     dbstop in calibratespinecoil2016 at 92;
%     Field.img = reducedMask .* ( Field.img - localField ) ;


end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================

function [] = acquirecalibrationdata( )
% Acquire calibration data 
%  Re-acquiring calibration data on "waterproof phantom" (28 L water) 
%  Experiment performed: 20160809
%
%     Notes
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
ShimsC = ShimCom;
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

end
% =========================================================================
% =========================================================================

end

