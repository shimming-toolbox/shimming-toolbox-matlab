classdef ShimCalUnfPrisma < ShimCal 
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
% Updated::20171116::ryan.topfer@polymtl.ca
% =========================================================================

properties
    % img ;
    % Hdr ;
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimCalUnfPrisma( Params )
%SHIMCAL - Shim calibration (i.e. reference maps) 

Shim.img = [] ;
Shim.Hdr = [] ;

if nargin < 1
    Params = ShimCalUnfPrisma.declarecalibrationparameters20171116( ) ;
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
function Params = declarecalibrationparameters20171116( )
%DECLARECALIBRATIONPARAMETERS20171018 
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

Params.filenameSave = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/UnfPrismaShimReferenceMaps20171116' ;

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
