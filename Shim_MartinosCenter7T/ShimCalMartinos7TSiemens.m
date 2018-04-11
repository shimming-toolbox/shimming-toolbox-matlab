classdef ShimCalMartinos7TSiemens < ShimCal
%SHIMCAL - Martinos 7T, Siemens Shim Calibration 
%
% .......
% 
% Usage
%
% Shim = ShimCalMartinos7TSiemens( )
%
%   Shim contains fields
%
%       .img
%           Shim reference maps (4d array) [units: Hz/DacValue]
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
%    ProbeTracking
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
% Updated::20170206::ryan.topfer@polymtl.ca
% =========================================================================


properties

end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Shim = ShimCalMartinos7TSiemens( Params )
%SHIMCALMARTINOS7TSIEMENS - Shim calibration (i.e. reference maps) 

Shim.img = [] ;
Shim.Hdr = [] ;

if nargin < 1
    Params = ShimCalMartinos7TSiemens.declarecalibrationparameters20170216() ;
end

% -------
% Process Zero-current data
%   zero-current acquisition assumed to be the first 2 directories (mag, phase)
Mag   = MaRdI( Params.dataLoadDir{1} ) ;
Mag.img = Mag.img(:,:,1:2:end) ; % ignore 2nd echo.
Mag.Hdr.NumberOfSlices = size( Mag.img, 3 ) ;

Phase = MaRdI( Params.dataLoadDir{2} ) ;

%-----
% Define data support (used for all field maps)
mask = (Mag.img > Params.magnitudeThreshold) ;
mask = dilater( mask, 1 ) ; % fills in internal 'holes'
mask = shaver( mask, 3 ) ; % erode voxels off edge since it exhibits artifacts

Phase.Hdr.MaskingImage = mask ;


BaselineField = ShimCal.preprocessfielddata( Phase, Params ) ;

% %----
% % crop to speed up processing?
% BaselineField = BaselineField.cropimg([94 94 80], [64 64 40] ) ;

dBdIRaw  = BaselineField.img;

Shim.Hdr = BaselineField.Hdr ;
Shim.img = zeros( BaselineField.getgridsize ) ;

Params.DataDir.Phase = Params.dataLoadDir(2) ; % just to fill field/declare variable

fieldMaps(:,:,:,1) = BaselineField.img ;

for iChannel = 1 : Params.nChannels 

    for iCurrent = 1 : (Params.nCurrents-1)

        iImg  = 4*iChannel + 2*(iCurrent-1) ;

        Params.DataDir.Phase( iCurrent ) = Params.dataLoadDir(iImg) ;

        Phase = MaRdI( Params.DataDir.Phase{ iCurrent }  ) ;
        Phase.Hdr.MaskingImage = mask ;    

        Field = ShimCal.preprocessfielddata( Phase, Params ) ;

        fieldMaps( :,:,:, iCurrent ) = Field.img ;

    end

    disp(['Channel ' num2str(iChannel) ' of ' num2str(Params.nChannels) ] )        

    dBdIRaw(:,:,:, iChannel) = ShimCal.mapdbdi( fieldMaps, Params.currents(iChannel,2:3), BaselineField.Hdr.MaskingImage, Params ) ;  

end 



% ShModel = calc_spherical_harmonics_arb_points_cz( (0:2), voxelPositions );
%
% shFields = X ;
% for iChannel = 1 : 9
%     shFields(:,:,:, iChannel)

% fm_dim = [128 128 80];  % cm
% FOV = [25.6 25.6 16.0];  % use 2mm iso. grid
% offset = FOV./fm_dim/2;
% [XX,YY,ZZ] = ndgrid(linspace(-FOV(1)/2+offset(1),FOV(1)/2-offset(1),fm_dim(1)),...
%     linspace(-FOV(2)/2+offset(2),FOV(2)/2-offset(2),fm_dim(2)),...
%     linspace(-FOV(3)/2+offset(3), FOV(3)/2-offset(3), fm_dim(3)));
%
% Params.Fitting.order = 9 ;  % polynomial fit order
%
% dBdIFit = [] ;
%
% if Params.Fitting.isFitting
%     disp(['Polynomial fitting... '] )        
%
%     % Polynomial fitting
%     for iChannel = 1 : Params.nChannels 
%
%         disp(['Channel ' num2str(iChannel) ' of ' num2str(Params.nChannels) ] )        
%
%         tmp       = dBdIRaw(:,:,:, iChannel);
%         tmp       = tmp(mask) ;
%
%         PolyModel = polyfitn( [X(mask), Y(mask), Z(mask)], tmp, Params.Fitting.isFitting );
%
%         tmp     = polyvaln( PolyModel, [X(:), Y(:), Z(:)]);
%
%         dBdIFit(:,:,:,iChannel) = reshape(tmp,[BaselineField.getgridsize()]);
%
%     end 
%
% end


if Params.Fitting.isFitting

    Params.ordersToGenerate = [0:2] ;
    ModelShim = ShimOptSHarmonics( Params, BaselineField ) ;
    ModelShim.setshimvolumeofinterest( BaselineField.Hdr.MaskingImage ) ;
    %       ModelShim.img(:,:,:,1) corresponds to 0th order, 
    %
    %       ModelShim.img(:,:,:,2:4) to 1st orders
    %         2 -> (y)  
    %         3 -> (z)   
    %         4 -> (x)   
    %
    %       ModelShim.img(:,:,:,5:9) to 2nd orders
    %         5 -> (xy)  
    %         6 -> (zy)   
    %         7 -> (z2)   
    %         8 -> (zx) 
    %         9 -> (x2y2)  

    dBdIFit = [] ;

    M         = ModelShim.gettruncationoperator() ;

    % -------
    % 4 -> (x)  
    iChannel  = 1; 
    dBdIModel = ModelShim.img(:,:,:,4) ;
    dBdIData  = dBdIRaw(:,:,:,iChannel) ;

    b = dBdIModel(:)'*M'*M*dBdIData(:) ;
    a(iChannel) = inv( dBdIModel(:)'*M'*M*dBdIModel(:) )*b ; % ideal scaling param

    dBdIFit(:,:,:, iChannel) = a(iChannel)*dBdIModel ;

    % -------
    % 2 -> (y)   
    iChannel  = iChannel + 1;
    dBdIModel = ModelShim.img(:,:,:,2) ;
    dBdIData  = dBdIRaw(:,:,:, iChannel) ;

    b    = dBdIModel(:)'*M'*M*dBdIData(:) ;
    a(iChannel) = inv( dBdIModel(:)'*M'*M*dBdIModel(:) )*b ; % ideal scaling param

    dBdIFit(:,:,:, iChannel) = a(iChannel)*dBdIModel ;

    % -------
    % 3 -> (z)   
    iChannel  = iChannel + 1;
    dBdIModel = ModelShim.img(:,:,:,3) ;
    dBdIData  = dBdIRaw(:,:,:, iChannel) ;

    b    = dBdIModel(:)'*M'*M*dBdIData(:) ;
    a(iChannel) = inv( dBdIModel(:)'*M'*M*dBdIModel(:) )*b ; % ideal scaling param

    dBdIFit(:,:,:, iChannel) = a(iChannel)*dBdIModel ;

    % -------
    % 7 -> (z2)   
    iChannel  = iChannel + 1;
    dBdIModel = ModelShim.img(:,:,:,7) ;
    dBdIData  = dBdIRaw(:,:,:, iChannel) ;

    b    = dBdIModel(:)'*M'*M*dBdIData(:) ;
    a(iChannel) = inv( dBdIModel(:)'*M'*M*dBdIModel(:) )*b ; % ideal scaling param

    dBdIFit(:,:,:, iChannel) = a(iChannel)*dBdIModel ;


    % -------
    % 8 -> (zy)   
    iChannel  = iChannel + 1;
    dBdIModel = ModelShim.img(:,:,:,8) ;
    dBdIData  = dBdIRaw(:,:,:, iChannel) ;

    b    = dBdIModel(:)'*M'*M*dBdIData(:) ;
    a(iChannel) = inv( dBdIModel(:)'*M'*M*dBdIModel(:) )*b ; % ideal scaling param

    dBdIFit(:,:,:, iChannel) = a(iChannel)*dBdIModel ;

    % -------
    % 6 -> (zx)   
    iChannel  = iChannel + 1;
    dBdIModel = ModelShim.img(:,:,:,6) ;
    dBdIData  = dBdIRaw(:,:,:, iChannel) ;

    b    = dBdIModel(:)'*M'*M*dBdIData(:) ;
    a(iChannel) = inv( dBdIModel(:)'*M'*M*dBdIModel(:) )*b ; % ideal scaling param

    dBdIFit(:,:,:, iChannel) = a(iChannel)*dBdIModel ;

    % -------
    % 9 -> (x2y2)   
    iChannel  = iChannel + 1;
    dBdIModel = ModelShim.img(:,:,:,9) ;
    dBdIData  = dBdIRaw(:,:,:, iChannel) ;

    b    = dBdIModel(:)'*M'*M*dBdIData(:) ;
    a(iChannel) = inv( dBdIModel(:)'*M'*M*dBdIModel(:) )*b ; % ideal scaling param

    dBdIFit(:,:,:, iChannel) = a(iChannel)*dBdIModel ;

    % -------
    % 5 -> (xy)   
    iChannel  = iChannel + 1;
    dBdIModel = ModelShim.img(:,:,:,5) ;
    dBdIData  = dBdIRaw(:,:,:, iChannel) ;

    b    = dBdIModel(:)'*M'*M*dBdIData(:) ;
    a(iChannel) = inv( dBdIModel(:)'*M'*M*dBdIModel(:) )*b ; % ideal scaling param

    dBdIFit(:,:,:, iChannel) = a(iChannel)*dBdIModel ;

    for iChannel = 1 : 8
        err(:,:,:,iChannel) = 100*mask.* abs( dBdIRaw(:,:,:,iChannel) - dBdIFit(:,:,:,iChannel) )./ abs(dBdIFit(:,:,:,iChannel))  ;
    end

    Shim.img = dBdIFit ;
else
    Shim.img = dBdIRaw ;
end

Shim.Hdr    = BaselineField.Hdr ;

save( Params.filenameSave ,'Shim');

%
% for iChannel = 1 : 8
%     figure
%     subplot(2,3,1),imagesc(mask(:,:,end/2).*dBdIRaw(:,:,end/2,iChannel)),axis image, axis off,colormap(jet),title([num2str(iChannel),': raw axial']),colorbar
%     subplot(2,3,2),imagesc(squeeze(mask(:,end/2,:).*dBdIRaw(:,end/2,:,iChannel))),axis image, axis off,colormap(jet),title('raw sagittal'),colorbar
%     subplot(2,3,3),imagesc(squeeze(mask(end/2,:,:).*dBdIRaw(end/2,:,:,iChannel))),axis image, axis off,colormap(jet),title('raw coronal'),colorbar
%     
%     subplot(2,3,4),imagesc(mask(:,:,end/2).*dBdIFit(:,:,end/2,iChannel)),axis image, axis off,colormap(jet),colorbar,title('polynomial interp. axial')
%     subplot(2,3,5),imagesc(squeeze(mask(:,end/2,:).*dBdIFit(:,end/2,:,iChannel))),axis image, axis off,colormap(jet),colorbar,title('polynomial interp. sagittal')
%     subplot(2,3,6),imagesc(squeeze(mask(end/2,:,:).*dBdIFit(end/2,:,:,iChannel))),axis image, axis off,colormap(jet),colorbar,title('polynomial interp. coronal')
%
%     title(['Channel ' num2str(iChannel)])
% end
%
% nii( 100 * (abs(dBdIFit - dBdIRaw)./abs(dBdIFit)) ) ;
%     
    
% if Params.Extension.isExtending
%     
%     disp(['Performing harmonic extension...']) ;
%     disp(['Channel 1 of ' num2str(Params.nChannels) ] ) ;
%
%     gridSizeImg = size( Shim.img(:,:,:,1) ) ;
%     Params.Extension.voxelSize = Shim.getvoxelsize() ;
%     
%     maskReduced = (sum(abs(Shim.img),4))~=0 ; % initial spatial support
%     mask        = maskReduced ; % extended spatial support
%     mask( Params.Extension.limitsForExtension(1,1) : (Params.Extension.limitsForExtension(1,2)), ...
%         Params.Extension.limitsForExtension(2,1) : (Params.Extension.limitsForExtension(2,2)), ...
%         Params.Extension.limitsForExtension(3,1) : (Params.Extension.limitsForExtension(3,2)) ) = 1 ;
%     
%     BackgroundField.Hdr.MaskingImage = mask ;
%
%     [ ~, A, M ] = extendharmonicfield( Shim.img(:,:,:, 1) , mask, maskReduced, Params.Extension ) ;
%
%     for iChannel = 1 : Params.nChannels 
%         disp(['Channel ' num2str(iChannel) ' of ' num2str(Params.nChannels) ] ) ;
%         
%         reducedField              = maskReduced .* Shim.img(:,:,:, iChannel) ;
%         extendedField             = reshape( M'*A*reducedField(:), gridSizeImg ) ;
%         Shim.img(:,:,:, iChannel) = extendedField + reducedField ;
%     end
%     
% end
%


end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================
function [dacValues] = getshimdacvalues20170130( )
%GETSHIMDACVALUES20170130
%
% dacValues = GETSHIMDACVALUES20170130()
%
% Returns full set of shim DAC values entered into Siemens console for experiment 
% on Martinos 7T scanner evening of 2017.01.30 (i.e. "boris_testing_32ch")
%
% dacValues is a matrix containing 8 shim terms (1st through 2nd order).
% The 3 columns correspond, from left to right, to the 
%   [ baseline_shim_offset, negative_shim_offset, positive_shim_offset]
% 
% And the 8 rows are as follows
% .x    
% .y    
% .z    
% .z2  
% .zx   
% .zy   
% .x2y2 
% .xy  

baselineDacValues.x    = 17.10 ;
baselineDacValues.y    = 9.00 ;
baselineDacValues.z    = -4.31 ;
baselineDacValues.z2   = -278.06 ;
baselineDacValues.zx   = 109.01 ;
baselineDacValues.zy   = -44.63 ;
baselineDacValues.x2y2 = -14.26 ;
baselineDacValues.xy   = 72.21 ;

dacValues = [ baselineDacValues.x, (baselineDacValues.x - 30.0)   , (baselineDacValues.x + 30.0) ;
              baselineDacValues.y, (baselineDacValues.y - 30.0)   , (baselineDacValues.y + 30.0) ;
              baselineDacValues.z, (baselineDacValues.z - 30.0)   , (baselineDacValues.z + 30.0) ;
              baselineDacValues.z2, (baselineDacValues.z2 - 600.0) , (baselineDacValues.z2 + 600.0) ;
              baselineDacValues.zx, (baselineDacValues.zx - 600.0)  , (baselineDacValues.zx + 600.0) ;
              baselineDacValues.zy, (baselineDacValues.zy - 600.0)  , (baselineDacValues.zy + 600.0) ;
              baselineDacValues.x2y2, (baselineDacValues.x2y2 - 800.0)  , (baselineDacValues.x2y2 + 800.0) ;
              baselineDacValues.xy, (baselineDacValues.xy - 800.0)  , (baselineDacValues.xy + 800.0) ; ] ;

end
% =========================================================================
function Params = declarecalibrationparameters20170216( )
%DECLARECALIBRATIONPARAMETERS20160130 
% 
% Initializes parameters for shim reference map construction (i.e. shim calibration)

fprintf('##### \n WARNING: Ignoring EchoTime DICOM header when normalizing phase to field. \n\n')
Params.isHardCodingEchoTime = true ;
Params.EchoTimes  = [3.16 4.18] ;

Params.nChannels  = 8 ;
Params.nCurrents  = 3 ;

Params.currents = ShimCalMartinos7TSiemens.getshimdacvalues20170130( ) ; % To do? resolve misnomer - these are not currents. 

fileIn = fopen('/Users/ryan/Projects/Shimming/MGH/20170130_referenceMaps7T/data/20170130_referenceMaps7Tdirlist.txt') ;
Params.dataLoadDir = textscan( fileIn, '%s' ) ;
Params.dataLoadDir = Params.dataLoadDir{1}; % not sure why this is necessary?
fclose(fileIn);
% NB: baseline shim acquisition assumed to be the first 2 directories (mag, phase)

Params.Fitting.isFitting      = false ;
if Params.Fitting.isFitting
    Params.filenameSave = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/Martinos7TReferenceMaps20160216_fit' ;
else
    Params.filenameSave = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/Martinos7TReferenceMaps20160216_raw' ;
end

Params.Filtering.isFiltering  = false ;
Mag                           = MaRdI( Params.dataLoadDir{1} ) ;
voxelSize                     = Mag.getvoxelsize() ;
Params.Filtering.filterRadius = 2*voxelSize(1) ;


Params.Extension.isExtending = false ; % harmonic field extrapolation 

% % Define data support (used for all field maps)
% Params.limitsForUnwrapping = [2 90; 2 28; 1 80] ; 

Params.magnitudeThreshold = 0.05 ; % as fraction of max intensity

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
