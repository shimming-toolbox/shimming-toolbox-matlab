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

