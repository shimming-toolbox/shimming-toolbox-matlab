classdef (Sealed=true) MrdiProt 
%MrdiProt  MR DICOM Image Protocol: Image acquisition protocol properties

% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% =========================================================================
properties( Transient=true )
    Img Mrdi ;
end

properties( Dependent, SetAccess=private )

    % Image acquisition times [ms elapsed since midnight]. Array dims [ nSlices x nEchoes x nMeasurements ]
    % 
    % Derives from the AcquisitionTime field in Siemens DICOM header.
    %
    % acquisitionTime - acquisitionTime(1) yields the elapsed time since
    % acquisition of the 1st k-space point in the series.
    %
    % For EPI MOSAIC (which has a single AcquisitionTime value in the DICOM
    % header for each volume), 
    % acquisitionTime( iSlice ) = AcquisitionTime (first slice) + iSlice*(volumeTR/nSlices)
    acquisitionTime {mustBeNonnegative} = 0 ; 

    % Gyromagnetic ratio of imaged nucleus [units: rad/s/T] (Note: Only implemented for 1H protons!)
    gyro(1,1) {mustBeReal} = 0 ; 
    
    % Larmor Tx. frequency [units: Hz]
    imagingFrequency(1,1) {mustBeNonnegative} = 0 ; 
    
    % [x,y,z] coordinates of magnet isocenter in the patient coordinate system
    isocenter(1,3) {mustBeReal} = [ 0 0 0 ] ; 
    
    % Number of image measurements
    nMeasurements(1,1) {mustBeInteger} = 1 ; 
    
    % Number of acquired slices (=1 for 3d slab encoding)
    nSlices(1,1) {mustBeInteger} = 1 ;  

    % Fraction of k-space coverage in each dim [read phase partition]
    partialFourierFactors {mustBePositive} = [ 1 1 1 ] ; 
    
    % Slice [units: mm] : DICOM Hdr tag ( )
    sliceThickness(1,1) {mustBePositive} = 1 ; 
    
    % spacingBetweenSlices [units: mm] : DICOM Hdr tag ( )
    spacingBetweenSlices(1,1) {mustBePositive} = 1 ; 
    
    % Vector of echo times [units: ms]
    echoTime {mustBeNonnegative} = 0 ; 
    
    % Repetition time [units: ms]
    repetitionTime(1,1) {mustBeNonnegative} = 0 ;

end
% =========================================================================
% =========================================================================

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Prot = MrdiProt( varargin )
    
    if nargin == 0
        return ;
    elseif nargin == 1 && isa( varargin{1}, 'Mrdi' ) ;
        Prot.Img  = vararargin{1} ; 
    else
        error('Invalid input') ;
    end

end
% =========================================================================
function [t] = get.acquisitionTime( Prot ) 
% GETACQUISITIONTIME  Return array of image acquisition times [ms elapsed since midnight]
%
% Derives from the AcquisitionTime field in Siemens DICOM header.
%
% t = GETACQUISITIONTIME( Prot )
%
% dimensions of t: [ nSlices x nEchoes x nMeasurements ]
%
% t - t(1) yields the elapsed time since acquisition of the 1st k-space point
% in the series.
%
% For EPI MOSAIC (which has a single AcquisitionTime value for each volume),
% t( iSlice ) = AcquisitionTime (first slice) + iSlice*(volumeTR/nSlices)

nSlices       = Prot.nSlices ;
nEchoes       = length( Prot.echoTime ) ;
nMeasurements = Prot.nMeasurements ;

t = zeros( nSlices, nEchoes, nMeasurements ) ;

if isempty( strfind( Prot.Hdrs(1).ImageType, 'MOSAIC' ) )
    
    for iMeasurement = 1 : nMeasurements
        for iEcho = 1 : nEchoes 
            for iSlice = 1 : nSlices
                t(iSlice, iEcho, iMeasurement) = convertacquisitiontime( Prot.Hdrs(iSlice,iEcho,iMeasurement,1).AcquisitionTime ) ;
            end
        end
    end

else % special case for EPI MOSAIC (single DICOM header for multiple slices)
   
    sliceOrder = Prot.Hdrs(1).MrProt.sSliceArray.alSliceAcqOrder ;
    sliceMeasurementDuration = Prot.repetitionTime/nSlices ; % [units: ms]

    for iMeasurement = 1 : nMeasurements
        for iEcho = 1 : nEchoes 
            for iSlice = 1 : nSlices
                t_iAcq = convertacquisitiontime( Prot.Hdrs(1,iEcho,iMeasurement,1).AcquisitionTime ) ;
                t_iAcq = t_iAcq + sliceOrder(iSlice)*sliceMeasurementDuration ;
                t(iSlice, iEcho, iMeasurement) = t_iAcq ;
            end
        end
    end

end

function t_s = convertacquisitiontime( AcquisitionTime )
%CONVERTACQUISITIONTIME  Converts the time-stamp string in the .AcquisitionTime DICOM header
% For details, see https://cfn.upenn.edu/aguirre/wiki/public:pulse-oximetry_during_fmri_scanning

    % hours + minutes + seconds + fractional seconds :
    t_s = str2double( AcquisitionTime(1:2) )*60*60 + str2double( AcquisitionTime(3:4) )*60 + ...
          str2double( AcquisitionTime(5:6) ) + (str2double( AcquisitionTime(7:10) )/10)/1000 ;
    t_s = t_s * 1000 ; % [units: ms]

end

end
% =========================================================================
function GYRO = get.gyro( Prot )
%GETGYRO  Return gyromagnetic ratio of imaged nucleus [units: rad/s/T]
%
% NOTE: Only supports 1H protons!

    switch Prot.Hdrs(1).ImagedNucleus 
        case '1H' 
            GYRO = 267.513E6 ; 
        otherwise
            error('Not implemented.') ;
    end

end
% =========================================================================
function xyzIso = get.isocenter( Prot )
%ISOCENTER  Return magnet isocenter [x,y,z] coordinates w.r.t. patient coordinate system
% 
% xyzIso = ISOCENTER( Prot ) 

    xyzIso = Prot.Hdrs(1).Prot.ImaRelTablePosition()' ;

    assert( xyzIso(1) == 0, 'Table shifted in L/R direction?' ) ;
    assert( xyzIso(2) == 0, 'Table shifted in A/P direction?' ) ;

end
% =========================================================================
function [f0] = get.imagingFrequency( Prot ) 
%GETIMAGINGFREQUENCY  Return Larmor frequency [units: Hz]

    f0 = Prot.Hdr.MrProt.sTXSPEC.asNucleusInfo.lFrequency ;

end
% =========================================================================
function [nVolumes] = get.nMeasurements( Prot )
%GETNMEASUREMENTS    Return the number of image measurements
%
% n = GETNMEASUREMENTS( Prot )
%
% NOTE: GETNMEASUREMENTS( Prot ) is equivalent to n = size( Prot.img, 5 ),
% however GETNUMBEROFMEASUREMENTS also checks the DICOM header in Prot.Hdr and
% issues a warning if n differs from the expected value
% (Prot.Hdr.MrProt.lRepetitions +1).

    nVolumes = size( Prot.img, 5 ) ;

    if myisfield( Prot.Hdr.MrProt, 'lRepetitions' ) 
        % number of measurements according to Hdr:
        nVolumesHdr = Prot.Hdr.MrProt.lRepetitions + 1 ;
        if nVolumes ~= nVolumesHdr
            warning([ num2str(nVolumes) ...
                ' image volumes have been loaded, however the DICOM Hdr indicates ' num2str(nVolumesHdr) ' were acquired.'] ) ;
        end 
    end

end
% =========================================================================
function [nSlices] = get.nSlices( Prot ) 
%GETNSLICES  Return number of acquired slices
%
% NOTE: nSlices is not necessarily equal to size( Prot.img, 3).
% e.g. For a 3d (slab) encoding, GETNUMBEROFSLICES returns 1. 

    nSlices = Prot.Hdr.MrProt.sSliceArray.lSize ;

end
% =========================================================================
function [pff] = get.partialFourierFactors( Prot ) 
%GETPARTIALFOURIERFACTORS  Return fraction of k-space coverage in each dim
% 
% pff = GETPARTIALFOURIERFACTORS( Prot ) 
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

    dcmFields = { 'Prot.Hdr.MrProt.sKSpace.ucReadoutPartialFourier' ;
                  'Prot.Hdr.MrProt.sKSpace.ucPhasePartialFourier' ;
                  'Prot.Hdr.MrProt.sKSpace.ucSlicePartialFourier' ; } ;

    pfAsInt = [ Prot.Hdr.MrProt.sKSpace.ucReadoutPartialFourier
                Prot.Hdr.MrProt.sKSpace.ucPhasePartialFourier
                Prot.Hdr.MrProt.sKSpace.ucSlicePartialFourier ]' ;

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
                error( [ 'Unexpected value of ' num2str( pfAsInt(iDim) ) ' in ' dcmFields{iDim} ] ) ;
        end
    end

end
% =========================================================================
function [te] = get.te( Prot ) 
%GETTE    Return vector of echo times [units: ms]
% 
% TE = GETECHOTIME( Prot )

    if strcmp( Prot.Hdr.SequenceName, '*fm2d2' ) && Prot.isphase() 
        te = ( Prot.Hdr.MrProt.alTE(2) - Prot.Hdr.MrProt.alTE(1) )/1000  ;
    else
        nEchoes  = size( Prot.img, 4 ) ;
        te = Prot.Hdr.MrProt.alTE(1:nEchoes)/1000 ;
    end

end
% =========================================================================
function tr = get.tr( Prot ) 
%GETTR repetition time [units: ms]

    tr = Prot.Hdrs(1).RepetitionTime ;

end
% % =========================================================================
% function [Hdrs] = set.Hdrs( Prot )
% %SETHDRS Update protocol with new cell array of DICOM Hdrs. 
%     % Prot.Hdrs = Hdrs ; 
% end
% =========================================================================

end
% =========================================================================
% =========================================================================    
methods(Sealed=true)
% =========================================================================
function [f0, g0, s0] = adjvalidateshim( Prot )
%ADJVALIDATESHIM    Return protocol shim settings from Siemens private header
% 
% [f0, g0, s0] = ADJVALIDATESHIM( Prot )
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
f0 = Prot.imagingFrequency ; 

% linear gradients
g0 = [ Prot.Hdr.MrProt.sGRADSPEC.asGPAData.lOffsetX ; ...
       Prot.Hdr.MrProt.sGRADSPEC.asGPAData.lOffsetY ; ...
       Prot.Hdr.MrProt.sGRADSPEC.asGPAData.lOffsetZ ] ;

% -------
% 2nd order shims (5-element vector)
s0 = Prot.Hdr.MrProt.sGRADSPEC.alShimCurrent ;

end
% =========================================================================
function [t0] = estimatekorigintime( Prot ) 
%ESTIMATEKORIGINTIME  Return estimate of time when k-space origin was sampled
% 
% t0 = ESTIMATEKORIGINTIME( Prot )
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
 
nSlices       = Prot.nSlices ;
nEchoes       = length( Prot.echoTime ) ;
nMeasurements = Prot.nMeasurements ;

tAcq = Prot.getacquisitiontime() ;

% Estimate time from excitation to k-space origin:
dt = 0; % [units: ms]

if Prot.Hdr.MrProt.sKSpace.ucTrajectory == 1
    
    if ~strcmp( Prot.Hdr.ScanningSequence, 'EPI' )
    % NOTE: The Siemens DICOM field SliceMeasurementDuration appears to be a misnomer:
    % Rather, it (usually?) refers to the duration of an entire volume. 
        dt = Prot.Hdr.Prot.SliceMeasurementDuration ; % [units: ms]
    else
    % *Exception*: Apparently EPI are handled differently as said DICOM field
    % is oddly large (probably includes dummy volumes?)
        dt = Prot.repetitionTime/nSlices ; % [units: ms]
    end
    
    pf = Prot.partialFourierFactors ;
    
    if all( pf == 1 )
        dt = dt/2 ;
    else
        switch Prot.Hdr.MRAcquisitionType 
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

end
% =========================================================================
% =========================================================================
methods(Sealed=true)
% =========================================================================
% function isSame = iscoincident( Img1, Img2 )
% %ISCOINCIDENT   Check coincidence of 2 images 
% % 
% % isSame = ISCOINCIDENT( Img1, Img2 )
% %
% % Returns TRUE if Img1 and Img2 possess coincident voxel positions
% % and number of measurements/volumes
% %
% % TODO: Check additional properties + add outputs for each?
%
% assert( ( nargin == 2 ) && isa( Img2, 'MrdiUtil' ), 'Missed required input: 2 MrdiUtil-objects' )
%
% isSame = false ;
%
% if Mrdi.compareimggrids( Img1, Img2 )  
%      
%     isSame = true ;
%
%     if ~strcmp( Img1.Hdr.SeriesDescription, Img2.Hdr.SeriesDescription )
%         warning('A computation is being performed based on images acquired in seperate series.')
%     end
% end
%
% end    
% =========================================================================

end
% =========================================================================
% =========================================================================


end
