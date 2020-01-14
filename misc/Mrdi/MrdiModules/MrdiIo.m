classdef MrdiIo < handle
%MrdiIo  MR DICOM Image IO: for instantiating and saving Mrdi objects and images

% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================


% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Img = MrdiIo( imgPath )


end
% =========================================================================
end
% =========================================================================
% =========================================================================    
methods(Static, Hidden=true) 
% =========================================================================
function [] = nii( Img )
%NII - Write MaRdI image to NiFtI file
%
% Wraps to NII( ) (which wraps to the NiFtI toolbox)   
%
%.....
%   Syntax
%
%   nii( Img ) 
%.....
%
% WARNING
%
%   nii() function is convenient for quickly writing a file to throw
%   into an external viewing application (e.g. ImageJ). 
%   The nifti Hdr info (i.e. orientation) is probably all wrong. 
%
%   To save NifTI's properly (takes longer) use Img.write() 

workingDir       = [ pwd '/' ] ;
Params.filename  = [ workingDir Img.Hdr.PatientName.FamilyName '_' num2str( Img.Hdr.SeriesNumber ) '_' Img.Hdr.SeriesDescription  ] ;
Params.voxelSize = Img.getvoxelspacing() ;

nii( squeeze(Img.img), Params )  ;

end
% =========================================================================
function Img = reshapemosaic( Img )
%RESHAPEMOSAIC Reshape Siemens mosaic into volume array and remove padded zeros
%
% Adapted from dicm2nii by
% xiangrui.li@gmail.com 
% http://www.mathworks.com/matlabcentral/fileexchange/42997

assert( ~isempty(strfind( Img.Hdr.ImageType, 'MOSAIC' ) ), 'Corrupt image header?' ) ;       

nImgPerLine = ceil( sqrt( Img.getnumberofslices() ) ); % always nImgPerLine x nImgPerLine tiles

nRows    = size(Img.img, 1) / nImgPerLine; 
nColumns = size(Img.img, 2) / nImgPerLine; 
nEchoes  = size(Img.img, 4) ;  
nVolumes = size(Img.img, 5) ;

img = zeros([nRows nColumns Img.getnumberofslices() nEchoes nVolumes], class(Img.img));

for iImg = 1 : Img.getnumberofslices()

    % 2nd slice is tile(1,2)
    r = floor((iImg-1)/nImgPerLine) * nRows + (1:nRows);     
    c = mod(iImg-1, nImgPerLine) * nColumns + (1:nColumns);

    img(:, :, iImg, :) = Img.img(r, c, :);
end

Img.img = img ;


% -----
% Update header

% -----
% Correct Hdr.ImagePositionPatient 
%   see: http://nipy.org/nibabel/dicom/dicom_mosaic.html

%-------
% Rotation and Scaling matrix: RS  
RS = Img.rotationMatrix * diag(Img.voxelSpacing) ;

Img.Hdr.ImagePositionPatient = Img.Hdr.ImagePositionPatient + ...
    RS(:,1)*(double(Img.Hdr.Rows) - nRows)/2 + ...
    RS(:,2)*(double(Img.Hdr.Columns) - nColumns)/2 ;


Img.Hdr.Rows    = uint16(nRows) ;
Img.Hdr.Columns = uint16(nColumns) ;

tmp = strfind( Img.Hdr.ImageType, '\MOSAIC' ) ;

if ~isempty(tmp) && (tmp > 1)
    Img.Hdr.ImageType = Img.Hdr.ImageType(1:tmp-1) ;
else
    Img.Hdr.ImageType = '' ;    
end

end
% =========================================================================
function [] = write( Img, saveDirectory, imgFormat, isSavingSingleNiis )
%WRITE  Write Mrdi image object to DICOM and/or NifTI
%
%   Usage 
%
%   WRITE( Img )
%   WRITE( Img, saveDirectory )
%   WRITE( Img, saveDirectory, imgFormat ) 
%   WRITE( Img, saveDirectory, imgFormat, isSavingSingleNiis ) 
%
%.....
%
%   Inputs
%
%   default saveDirectory = './tmp'
%   
%   imgFormat can be:
%       'dcm' [default]
%       'nii' (creating temporary DICOMs which are deleted after the system call to dcm2niix)
%       'both' (does not delete the DICOMs)
%
%   isSavingSingleNiis (boolean):
%       false [default] : DICOMs are combined into single NifTI file
%       true : Separate .nii output for each image (passes '-s y' argument to dcm2niix)
%
%.....
%
% Adapted from dicom_write_volume.m (D.Kroon, University of Twente, 2009)
% https://www.mathworks.com/matlabcentral/fileexchange/27941-dicom-toolbox?focused=5189263&tab=function

DEFAULT_SAVEDIRECTORY      = './tmp' ;
DEFAULT_IMGFORMAT          = 'dcm' ;
DEFAULT_ISSAVINGSINGLENIIS = false ;

if nargin < 2 || isempty(saveDirectory)
    saveDirectory = DEFAULT_SAVEDIRECTORY ;
end

if nargin < 3 || isempty(imgFormat)
    imgFormat = DEFAULT_IMGFORMAT ;
end

if nargin < 4 || isempty(isSavingSingleNiis)
    isSavingSingleNiis = DEFAULT_ISSAVINGSINGLENIIS ;
end

fprintf(['\n Writing images to: ' saveDirectory ' ... \n'])

[~,~,~] = mkdir( saveDirectory ) ;

rescaleimg( Img, true ) ;

[X,Y,Z] = Img.getvoxelpositions() ;

%-------
% write Hdr
Hdr.NumberOfSlices          = size( Img.img, 3 ) ;

% Make random series number
SN                          = round(rand(1)*1000);
% Get date of today
today                       = [datestr(now,'yyyy') datestr(now,'mm') datestr(now,'dd')];
Hdr.SeriesNumber            = SN;
Hdr.AcquisitionNumber       = SN;
Hdr.StudyDate               = today;
Hdr.StudyID                 = num2str(SN);
Hdr.PatientID               = num2str(SN);
Hdr.AccessionNumber         = num2str(SN);

% copy from original
Hdr.ImageType               = Img.Hdr.ImageType ; 

Hdr.StudyDescription        = Img.Hdr.StudyDescription ;
Hdr.SeriesDescription       = Img.Hdr.SeriesDescription ;
Hdr.Manufacturer            = Img.Hdr.Manufacturer ;
Hdr.ScanningSequence        = Img.Hdr.ScanningSequence ;
Hdr.SequenceVariant         = Img.Hdr.SequenceVariant ;
Hdr.ScanOptions             = Img.Hdr.ScanOptions ;
Hdr.MRAcquisitionType       = Img.Hdr.MRAcquisitionType ;
Hdr.SliceThickness          = Img.Hdr.SliceThickness ;
Hdr.SpacingBetweenSlices    = Img.Hdr.SpacingBetweenSlices ;
Hdr.PatientPosition         = Img.Hdr.PatientPosition ;
Hdr.PixelSpacing            = Img.Hdr.PixelSpacing ;

Hdr.ImageOrientationPatient = Img.Hdr.ImageOrientationPatient ;
Hdr.SliceLocation           = Img.Hdr.SliceLocation ; 

Hdr.AcquisitionMatrix       = Img.Hdr.AcquisitionMatrix ; 
Hdr.RepetitionTime          = Img.Hdr.RepetitionTime ; 
Hdr.NumberOfAverages        = Img.Hdr.NumberOfAverages ; 
Hdr.PercentSampling         = Img.Hdr.PercentSampling ; 
Hdr.PercentPhaseFieldOfView = Img.Hdr.PercentPhaseFieldOfView ; 
Hdr.InPlanePhaseEncodingDirection  = Img.Hdr.InPlanePhaseEncodingDirection ; 

[rHat, cHat, sHat] = Img.getdirectioncosines( ) ;  

nSlices       = size( Img.img, 3 ) ;
nEchoes       = numel( Img.getechotime() ) ;
nAcquisitions = numel( Img.getacquisitiontime ) ;

nImg = nSlices*nEchoes*nAcquisitions ;

for iSlice = 1 : nSlices 
    for iEcho  = 1 : nEchoes 
        for iAcq = 1 : nAcquisitions  

            iImg = iSlice*iEcho*iAcq ;

            disp( [num2str(100*iImg/nImg) '% ...'] ) ;
    
            %-------
            % filename 
            sliceSuffix = '000000' ;
            sliceSuffix = [ sliceSuffix( length(iSlice) : end ) num2str(iSlice) '-' num2str(iEcho) '-' num2str(iAcq) ] ;
            sliceSuffix = ['-' sliceSuffix '.dcm'] ;
            filename    = [saveDirectory '/' Img.Hdr.PatientName.FamilyName sliceSuffix] ;

            %-------
            % image specific hdr info 
            Hdr.ImageNumber          = iSlice*iEcho*iAcq ;
            Hdr.InstanceNumber       = iSlice*iEcho*iAcq ;
            Hdr.AcquisitionTime      = Img.Hdrs{iSlice,iEcho,iAcq}.AcquisitionTime ;
            
            Hdr.ImagePositionPatient = [(X(1,1,iSlice)) (Y(1,1,iSlice)) (Z(1,1,iSlice))] ;

            Hdr.SliceLocation        = dot( Hdr.ImagePositionPatient, sHat ) ;
           
            dicomwrite( Img.img(:,:,iSlice,iEcho, iAcq) , filename, 'ObjectType', 'MR Image Storage', Hdr ) ;
            
            if( iSlice==1 )
                info                  = dicominfo( filename ) ;
                Hdr.StudyInstanceUID  = info.StudyInstanceUID ;
                Hdr.SeriesInstanceUID = info.SeriesInstanceUID ;
            end

        end
    end
end

%-------
if ~strcmp( imgFormat, 'dcm' )

    if isSavingSingleNiis
        list = findimages( saveDirectory ) ;
        for iImg = 1 : size( list, 1 )
            system(['dcm2niix -s y -f %f_%r ' list{iImg}]) ;
        end
    else
        system(['dcm2niix ' saveDirectory]) ;
    end

    if strcmp( imgFormat, 'nii' )
        delete( [saveDirectory '/*.dcm'] ) ;
    end
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Access=protected)
% =========================================================================
function [] = rescaleimg( Img, isUndoing )
%RESCALEIMG  Rescale image to/from uint16 for saving/loading from DICOM files
%
% []=RESCALEIMG( Img ) 
% []=RESCALEIMG( Img, isUndoing ) 

% Only alters phase images (+ Hdr) by rescaling to radians. Mag. is unaffected.
return ;

end
% =========================================================================

function [img, Hdrs] = loadandsortimages( imgPath, imgTypes, conjunction )
%LOADANDSORTIMAGES Load and sort image files
    arguments
        imgPath (1,:) char = './' ;
        imgTypes = {'*.dcm', '*.IMA'} ;
        conjunction {mustBeMember(conjunction,{'AND', 'OR'})} = 'OR' ;
    end

%% -----
% Validate input 
if isempty(imgPath) || ~isstr( imgPath ) 
    errMsg = ['Function takes a single input: ' ...
        'path string of a single dicom image OR dicom-containing folder' ] ;
    eval( ['help ' mfilename ] ) ; 
    error( errMsg ) ;
end

img  = [] ;
Hdrs = [] ;

%% -----
% Find images 
if isstr( imgTypes )
    
    imgList = findfiles( imgPath, imgTypes ) ;

elseif iscell( imgTypes )

    for iImgType = 1 : numel( imgTypes )
        imgList(iImgType) = findfiles( imgPath, imgTypes{iImgType} ) ;
    end

end
dbstop in loadandsortimages at 37        
%% -----
% Load and sort images 
if isempty( imgList )
    display( ['No images of type ' imgTypes ' found in ' imgPath ] ) ;
    return ;
else
    switch imgTypes
        case strcmpi( imgTypes, {'.dcm', '.IMA'} ) 
            [img, Hdrs] = loadandsortdicoms( imgList ) ;
        otherwise
            error(['Image type ' imgTypes ' not implemented']) ;
    end
end


%% -----
function [] = loadandsortdicoms( imgList ) 
%LOADANDSORTDICOMS  Load and sort dicoms into (ascending?) order

Hdrs = {} ;

for iDir = 1 : numel( imgList )
    for iImg = 1 : numel( imgList{iDir} ) 
        Hdrs{ end+1 } = dicominfo( imgList{iDir, iImg} ) ;
    end
end

nImg = numel( Hdrs ) ; % total number of images

...t.b.c.
% Read protocol info from 1st loaded image (arbitrary):
SpecialHdr = dicominfosiemens( imgList{1} ) ;

nRows      = SpecialHdr.Height ;
nColumns   = SpecialHdr.Width ;

if ~myisfield( SpecialHdr.MrProt, 'lRepetitions' ) 
    nVolumes = 1 ;
else
    nVolumes = SpecialHdr.MrProt.lRepetitions + 1 ;
end

% containers for a few terms varying between dicoms:
sliceLocations = zeros( nImages, 1 ) ; % [units: mm]
echoTimes      = zeros( nImages, 1 ) ; % [units: ms]

% Image headers, not yet organized:
RawHdrs        = cell( nImages, 1 ) ; 

% Read all the image headers and structure the Mrdi object accordingly:
%
% NOTE: It may well be possible to determine all the necessary info from
% the complete Siemens header (e.g. SpecialHdr) of a single image;
% practically, however, the following is easier since the header
% information is abstrusely defined for some sequences.
for iImg = 1 : nImages
    % using dicominfo() here as parse-siemens-shadow() takes much longer
    RawHdrs{ iImg }      = dicominfo( imgList{iImg} ) ;

    sliceLocations(iImg) = RawHdrs{ iImg }.SliceLocation ;
    echoTimes(iImg)      = RawHdrs{ iImg }.EchoTime ;
end

sliceLocations = sort( unique( sliceLocations ), 1, 'ascend' ) ; 
nSlices        = length( sliceLocations ) ;

echoTimes      = sort( unique( echoTimes ), 1, 'ascend' ) ; 
nEchoes        = length( echoTimes ) ; 

Img.img = zeros( nRows, nColumns, nSlices, nEchoes, nVolumes ) ;

% copy of RawHdrs to be reorganized:
Img.Hdrs = cell( nSlices, nEchoes, nVolumes ) ;

for iImg = 1 : nImages

    iHdr    = RawHdrs{ iImg } ;

    iVolume = iHdr.AcquisitionNumber ;
    iSlice  = find( iHdr.SliceLocation == sliceLocations ) ;
    iEcho   = find( iHdr.EchoTime == echoTimes ) ;

    Img.Hdrs{ iSlice, iEcho, iVolume } = iHdr ;

    Img.img( :, :, iSlice, iEcho, iVolume ) = dicomread( imgList{iImg} ) ;

end

% Save the complete 1st Hdr
Img.Hdr = dicominfosiemens( Img.Hdrs{1}.Filename ) ;

% Img.rescaleimg() ;
%
% if ~myisfield( Img.Hdr, 'SpacingBetweenSlices' ) 
%     Img.Hdr.SpacingBetweenSlices = Img.Hdr.SliceThickness ;
% end
%
% if ~isempty( strfind( Img.Hdr.ImageType, 'MOSAIC' ) ) 
%     Img = reshapemosaic( Img ) ;
% end
%
% Img.setslicenormalvector() ;

end %loadandsortdcm()

end
