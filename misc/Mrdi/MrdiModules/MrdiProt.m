classdef (Abstract) MrdiProt < handle
%MrdiProt   MR DICOM Image Protocol 
%
% Member methods for retrieving information about the imaging/acquisition protocol.
% 
% e.g.
%   adjvalidateshim()
%   getimagingfrequency()
%   getacquisitiontime()
%   ...etc.
% 
% For documentation, type doc MrdiProt
%
% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Prot = MrdiProt( Img )

end
% =========================================================================
end

% =========================================================================
% =========================================================================    
methods(Sealed=true)
% =========================================================================
function [f0, g0, s0] = adjvalidateshim( Img )
%ADJVALIDATESHIM    Return protocol shim settings from Siemens private header
% 
% [f0, g0, s0] = ADJVALIDATESHIM( Img )
% 
% ADJVALIDATESHIM is not a particularly revealing name for a function; however,
% it is based on the Siemens AdjValidate commandline tool, with a "-shim" input
% argument.
%
% f0 = scalar Larmor (Tx) freq. [ units : Hz ]
% g0 = 3-element vector of the linear gradient offsets (gX, gY, gZ) [units : DAC bits]
% s0 = 5-element vector of 2nd order shim currents [units : mA] 
%
% To convert to the units of the 3D Shim card on the Siemens (Prisma) console,
% see 
% ShimSpecs_IUGM_Prisma_fit.converttomultipole( )
% ShimSpecs_HGM_Prisma_fit.converttomultipole( )

% Larmor (Tx) freq.
f0 = Img.getimagingfrequency() ; 

% linear gradients
g0 = [ Img.Hdr.MrProt.sGRADSPEC.asGPAData.lOffsetX ; ...
       Img.Hdr.MrProt.sGRADSPEC.asGPAData.lOffsetY ; ...
       Img.Hdr.MrProt.sGRADSPEC.asGPAData.lOffsetZ ] ;

% -------
% 2nd order shims (5-element vector)
s0 = Img.Hdr.MrProt.sGRADSPEC.alShimCurrent ;

end
% =========================================================================
function [f0] = getimagingfrequency( Img ) 
%GETIMAGINGFREQUENCY    Return Larmor frequency [units: Hz]
% 
% f0 = GETIMAGINGFREQUENCY( Img ) ;  

f0 = Img.Hdr.MrProt.sTXSPEC.asNucleusInfo.lFrequency ;

end
% =========================================================================
function [nVolumes] = getnumberofmeasurements( Img )
%GETNUMBEROFMEASUREMENTS    Return the number of measurements
%
% n = GETNUMBEROFMEASUREMENTS( Img )
%
% NOTE
%
% GETNUMBEROFMEASUREMENTS( Img ) is equivalent to n = size( Img.img, 5 ),
% however GETNUMBEROFMEASUREMENTS also checks the DICOM header in Img.Hdr and
% issues a warning if n differs from the expected value
% (Img.Hdr.MrProt.lRepetitions +1).

nVolumes = size( Img.img, 5 ) ;

if myisfield( Img.Hdr.MrProt, 'lRepetitions' ) 
    % number of measurements according to hdr:
    nVolumesHdr = Img.Hdr.MrProt.lRepetitions + 1 ;
    if nVolumes ~= nVolumesHdr
        warning([ num2str(nVolumes) ' image volumes have been loaded, however the DICOM Hdr indicates ' num2str(nVolumesHdr) ' were acquired.'] ) ;
    end 
end

end
% =========================================================================
function [nSlices] = getnumberofslices( Img ) 
%GETNUMBEROFSLICES  Return number of acquired slices
% 
% nSlices = GETNUMBEROFSLICES( Img ) 
%
% NOTE: nSlices is not necessarily equal to size( Img.img, 3).
% e.g. For a 3d (slab) encoding, GETNUMBEROFSLICES returns 1. 

nSlices = Img.Hdr.MrProt.sSliceArray.lSize ;

end
% =========================================================================
function [pff] = getpartialfourierfactors( Img ) 
%GETPARTIALFOURIERFACTORS  Return fraction of k-space coverage in each dim
% 
% pff = GETPARTIALFOURIERFACTORS( Img ) 
% 
% Returns 3-element vector of partial Fourier factors in read, phase (in-plane),
% and partition (slice) encoding directions.
%
% NOTE:
% Siemens uses an enumeration scheme to store Partial Fourier info in the DICOM header:
%
% pff --> Siemens DICOM Hdr entry 
% 4/8 --> 0x1  = 1
% 5/8 --> 0x2  = 2
% 6/8 --> 0x4  = 4
% 7/8 --> 0x8  = 8
% 8/8 --> 0x10 = 16 (i.e. no partial Fourier)
% 
% See: https://github.com/malaterre/GDCM/blob/master/Source/DataDictionary/CSAHeader.xml

dcmFields = { 'Img.Hdr.MrProt.sKSpace.ucReadoutPartialFourier' ;
              'Img.Hdr.MrProt.sKSpace.ucPhasePartialFourier' ;
              'Img.Hdr.MrProt.sKSpace.ucSlicePartialFourier' ; } ;

pfAsInt = [ Img.Hdr.MrProt.sKSpace.ucReadoutPartialFourier
            Img.Hdr.MrProt.sKSpace.ucPhasePartialFourier
            Img.Hdr.MrProt.sKSpace.ucSlicePartialFourier ]' ;

pff     = [ 1 1 1 ] ;

for iDim = 1 : 3
    switch pfAsInt(iDim)
        case 1
            pff(iDim) = 4/8 ;
        case 2
            pff(iDim) = 5/8 ;
        case 4
            pff(iDim) = 6/8 ;
        case 8
            pff(iDim) = 7/8 ;
        case 16
            pff(iDim) = 1 ;
        otherwise
            error( [ 'Unexpected value of ' num2str( pfAsInt(iDim) ) ' in ' dicomFields{iDim} ] ) ;
    end
end

end
% =========================================================================
function [echoTime] = getechotime( Img, iEcho ) 
%GETECHOTIME    Return vector of echo times [units: ms]
% 
% TE = GETECHOTIME( Img )
% TE = GETECHOTIME( Img, iEcho )
%
% If 2nd argument (echo index iEcho) is provided, GETECHOTIME returns the TE of
% the corresponding echo.

if nargin == 1
    if strcmp( Img.Hdr.SequenceName, '*fm2d2' ) && Img.isphase() 
        echoTime = ( Img.Hdr.MrProt.alTE(2) - Img.Hdr.MrProt.alTE(1) )/1000  ;
    else
        nEchoes  = size( Img.img, 4 ) ;
        echoTime = Img.Hdr.MrProt.alTE(1:nEchoes)/1000 ;
    end
else
    echoTime = Img.Hdr.MrProt.alTE(iEcho)/1000 ;
end 

end
% =========================================================================
function [t0] = estimatekorigintime( Img ) 
%ESTIMATEKORIGINTIME    Return estimate of time when k-space origin was sampled
% 
% t0 = ESTIMATEKORIGINTIME( Img )
%
% Returns an estimate of when the k-space origin of an image was sampled
% relative to the AcquisitionTime (field in Siemens DICOM header) as a double
% in units of milliseconds. 
% 
% See also: MrdiProt.getacquisitiontime()
% 
% NOTE: This is a crude estimate and only the case of Cartesian k-space
% sampling, beginning at the k_min periphery, has been considered in the
% current implementation!
 
nSlices = Img.getnumberofslices() ;
nEchoes = length( Img.getechotime() ) ;

if myisfield( Img.Hdr.MrProt, 'lRepetitions' ) 
    nMeasurements = ( Img.Hdr.MrProt.lRepetitions + 1) ;
else
    nMeasurements = 1;
end

assert( nMeasurements == size( Img.img, 5 ), 'Invalid .Hdr' ) ;

tAcq = Img.getacquisitiontime() ;

% Estimate time from excitation to k-space origin:
dt = 0; % [units: ms]

if Img.Hdr.MrProt.sKSpace.ucTrajectory == 1
    
    if ~strcmp( Img.Hdr.ScanningSequence, 'EPI' )
    % NOTE: The Siemens DICOM field SliceMeasurementDuration appears to be a misnomer:
    % Rather, it (usually?) refers to the duration of an entire volume. 
        dt = Img.Hdr.Img.SliceMeasurementDuration ; % [units: ms]
    else
    % *Exception*: Apparently EPI are handled differently as said DICOM field
    % is oddly large (probably includes dummy volumes?)
        dt = Img.Hdr.RepetitionTime/nSlices ; % [units: ms]
    end
    
    pf = Img.getpartialfourierfactors() ;
    
    if all( pf == 1 )
        dt = dt/2 ;
    else
        switch Img.Hdr.MRAcquisitionType 
            case '2D'
                dt = dt*( 3/2 - pf(2) ) ; 
            case '3D' % not sure what the function is here...
                error('Not implemented!') ;
        end
    end

else    
    display('Estimating time at which the k-space origin was sampled')
    warning(' Current strategy is simplistic and likely wrong for certain encoding trajectories.')
    dt = 0 ;
end

assert( dt >= 0, 'Unexpected value for time-to-origin since excitation (AcquisitionTime). See code + check DICOM entries are valid.' )

t0 = tAcq + dt ;

end
% =========================================================================
function [t] = getacquisitiontime( Img ) 
% GETACQUISITIONTIME  Return array of image acquisition times [ms elapsed since midnight]
%
% Derives from the AcquisitionTime field in Siemens DICOM header.
%
% t = GETACQUISITIONTIME( Img )
%
% dimensions of t: [ nSlices x nEchoes x nMeasurements ]
%
% t - t(1) yields the elapsed time since acquisition of the 1st k-space point
% in the series.
%
% For EPI MOSAIC (which has a single AcquisitionTime value for each volume),
% t( iSlice ) = AcquisitionTime (first slice) + iSlice*(volumeTR/nSlices)

nSlices = Img.getnumberofslices() ;
nEchoes = length( Img.getechotime() ) ;

if myisfield( Img.Hdr.MrProt, 'lRepetitions' ) 
    nMeasurements = ( Img.Hdr.MrProt.lRepetitions + 1) ;
else
    nMeasurements = 1;
end

assert( nMeasurements == size( Img.img, 5 ), 'Invalid .Hdr' ) ;

t = zeros( nSlices, nEchoes, nMeasurements ) ;

if isempty(strfind( Img.Hdrs{1}.ImageType, 'MOSAIC' ) )
    
    for iMeasurement = 1 : nMeasurements
        for iEcho = 1 : nEchoes 
            for iSlice = 1 : nSlices
                % RT: commented 20191013, added convertacquisitiontime() -- which is probably more correct/accurate?
                %
                % t_iAcq = Img.Hdrs{iSlice,iEcho,iMeasurement}.AcquisitionTime ;
                % t_iAcq = str2double( t_iAcq(1:2) )*60*60 + str2double( t_iAcq(3:4) )*60 + str2double( t_iAcq(5:end) ) ;
                % t(iSlice, iEcho, iMeasurement) = t_iAcq ;
                
                t(iSlice, iEcho, iMeasurement) = convertacquisitiontime( Img.Hdrs{iSlice,iEcho,iMeasurement}.AcquisitionTime ) ;
            end
        end
    end

else % special case for EPI MOSAIC (single DICOM header for multiple slices)
   
    sliceOrder = Img.Hdr.MrProt.sSliceArray.alSliceAcqOrder ;
    sliceMeasurementDuration = Img.Hdr.RepetitionTime/nSlices ; % [units: ms]

    for iMeasurement = 1 : nMeasurements
        for iEcho = 1 : nEchoes 
            for iSlice = 1 : nSlices
                % See above comment 20191013
                % t_iAcq = Img.Hdrs{1,iEcho,iMeasurement}.AcquisitionTime  ;
                % t_iAcq = str2double( t_iAcq(1:2) )*60*60 + str2double( t_iAcq(3:4) )*60 + str2double( t_iAcq(5:end) ) ;
                t_iAcq = convertacquisitiontime( Img.Hdrs{1,iEcho,iMeasurement}.AcquisitionTime ) ;
                t_iAcq = t_iAcq + sliceOrder(iSlice)*sliceMeasurementDuration ;
                t(iSlice, iEcho, iMeasurement) = t_iAcq ;
            end
        end
    end

end

function t_s = convertacquisitiontime( AcquisitionTime )
%CONVERTACQUISITIONTIME
% 
% Converts the time-stamp string in the .AcquisitionTime DICOM header
%
% For details, see https://cfn.upenn.edu/aguirre/wiki/public:pulse-oximetry_during_fmri_scanning

% hours + minutes + seconds + fractional seconds :
t_s = str2double( AcquisitionTime(1:2) )*60*60 + str2double( AcquisitionTime(3:4) )*60 + ...
      str2double( AcquisitionTime(5:6) ) + (str2double( AcquisitionTime(7:10) )/10)/1000 ;
t_s = t_s * 1000 ; % [units: ms]

end

end
% =========================================================================
function GYRO = getgyromagneticratio( Img )
%GETGYROMAGNETICRATIO  Return gyromagnetic ratio of imaged nucleus in units of rad/s/T.
%
% Gyro = getgyromagneticratio( Img )
%
% NOTE: Only supports 1H protons!

switch Img.Hdr.ImagedNucleus 
    case '1H' 
        GYRO = 267.513E6 ; 
    otherwise
        error('Not implemented.') ;
end

end
% =========================================================================
function xyzIso = isocenter( Img )
%ISOCENTER
% 
% xyzIso = ISOCENTER( Img ) 
%
% Returns the 3-element vector of the x, y and z coordinates of the magnet
% isocenter in the patient coordinate system

xyzIso = Img.Hdr.Img.ImaRelTablePosition()' ;

assert( xyzIso(1) == 0, 'Table shifted in L/R direction?' ) ;
assert( xyzIso(2) == 0, 'Table shifted in A/P direction?' ) ;

end
% =========================================================================

end
% =========================================================================
% =========================================================================



end
