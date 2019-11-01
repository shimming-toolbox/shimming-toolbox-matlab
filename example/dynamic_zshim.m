function dynamic_zshim(FM_mag_path, FM_phase_path, MGRE_mag_path, respTrace_path, varargin)

%% ************************************************************************
% function dynaminc_zshim(FM_mag_path, FM_phase_path, MGRE_mag_path, respTrace_path, varargin)
%
% DESCRIPTION: This function will generate static and dynaminc (due to
% respiration) Gz components based on a fieldmap time series (magnitude and
% phase images to be found in 'FM_mag_path' and 'FM_phase_path') and
% respiratory trace information obtained from Siemens bellows (found in
% respTrace_path). An additional multi-gradient echo (MGRE) magnitiude
% image is used (found in MGRE_mag_path) to generate an ROI and resample
% the static and dynaminc Gz component maps to match the MGRE image. 
% Lastly the average Gz values within the ROI are computed for each slice.
%
% INPUTS: 
%
% FM_mag_path : path to DICOM folder containing magnitude magniude images for
%            field mapping timeseries
% 
% FM_phase_path : path to DICOM folder containing magnitude phase images for
%              field mapping timeseries
% 
% MGRE_mag_path : path to DICOM folder containing multi-gradient echo
% (MGRE) magnitude images to be used for segmentation
%
% respTrace_path : path to Siemens respiratory trace recording
% 
% varargin : sort images obtained from a DICOM socket transfer and choose
% whether to move or copy them using 'SocketTransfer_move' or 
% 'SocketTransfer_copy'
%
%
% OUTPUT: text file containing the static and dynamic Gz comnponent values 
% for each slice of the magnitude images found in 'MGRE_mag_path'
% 
% AUTHORS: 
%
% Ryan Topfer, MSc. ryan.topfer@polymtl.ca
% Eva Alonso Ortiz (EAO), PhD. eva.alonso.ortiz@gmail.com 
%
% LAST MODIFIED: Oct. 2019
%
%*************************************************************************


%% ------------------------------------------------------------------------
% Sort DICOM socket transfer images
% EAO todo: enable unsortedDicomDir/sortedDicomDir read-in
%% ------------------------------------------------------------------------
if nargin > 4 
    if strcmp(varargin{1}, 'SocketTransfer_move')
        % isCopying (boolean) is optional (move or copy the .IMA files)
        MaRdI.sortimages( unsortedDicomDir, sortedDicomDir );
    elseif strcmp(varargin{1}, 'SocketTransfer_copy')
        MaRdI.sortimages( unsortedDicomDir, sortedDicomDir, 1 );
    end
end
%% ------------------------------------------------------------------------


%% ------------------------------------------------------------------------
% load MGRE magnitude images for SCT segmentation
% Todo -> Julien
%% ------------------------------------------------------------------------
Mag = MaRdI(MGRE_mag_path);
Params.dataLoadDir = MGRE_mag_path;
Params.centerlineMethod = 'midfov';  % set to 'spinalcord' to create a mask around the spinal cord
shimVoi = Mag.segmentspinalcanal_s(Params);

%% ------------------------------------------------------------------------
% load images and respiratory trace recording as FieldEval + ProbeTracking 
% objects :
%% ------------------------------------------------------------------------

% field map time series
Fields = FieldEval( FM_mag_path, FM_phase_path ); 
% Note: Field.img (static b0) refers to the *mean probe signal (saved in the output as Field.Aux.Data.p) 
% and the respiratory component is a relative deviation from the mean (scaled to RMS PMU)
% 
% That means in the example case, where PMU_mean_value = Field.Aux.Data.p = 1707.13, whenever the PMU 
% reading = 1707.13, there respiratory correction at that moment should be ZERO.
% i.e. correction for iSlice would be: 
% Corrections.static( iSlice ) + ( PMU_current_value - PMU_mean_value ) * Corrections.riro( iSlice )
% 
% in this way, the value of Field.Aux.Data.p needs to be written into the sequence as well
% 
% Alternatively, Corrections.static( iSlice ) could be scaled to refer to the PMU=0 point (then there 
% is no need to keep track of PMU_mean_value) but that would mean tying the static correction at any 
% point in time to the current PMU reading, and that seems less stable to me (e.g. if the belt loosens, 
% or the subject touches it, then both respiratory and static corrections fail, rather than just the former)

% Siemens PMU recording
Pmu   = ProbeTracking(respTrace_path);

%% ------------------------------------------------------------------------
% link the two objects (interpolate PMU to the fieldmap time series)
%% ------------------------------------------------------------------------
Fields.associateaux( Pmu );

% modeled static + respiratory fields (in Field.img and Field.Model.Riro.img respectively)
Field = FieldEval.modelfield( Fields );

%% ------------------------------------------------------------------------
% interpolate static B0 field + Riro images to target slices for shimming:
%% ------------------------------------------------------------------------
[X,Y,Z]  = Mag.getvoxelpositions() ;

% accelerate interpolation by restricting it to the region where signal exists: 
% (NOTE: mask here could instead be the shimVoi output from SCT)
mask = Mag.getreliabilitymask( 0.05 ) ; % returns mask for each echo (4th dim)
mask = logical( sum( mask, 4 ) ) ; % combine echo-specific masks

% 'interp/extrap' (nearest-neighbour substitution for 2d field maps) :
Field.resliceimg( X,Y,Z, mask ) ; % reslice static b0 image 
Field.Model.Riro.resliceimg( X,Y,Z, mask ) ; % reslice RIRO image

%% ------------------------------------------------------------------------
% define the shim system using 'nominal' (ideal) reference maps
%
% alternatively, empitical shim reference maps for the Siemens Prisma at the
% IUGM can be downloaded from :
% https://drive.google.com/open?id=1X3kDizzZeZK2dxs6D_zH8bBbjDaF1veu 
%
% the saved binary location should be referenced in: function shimDir = shimbindir( )
%% ------------------------------------------------------------------------
Shims = ShimOpt_IUGM_Prisma_fit( Field );

% 3rd term along 4th dim of Shims.img is the Z-gradient ref. map [units: Hz/micro-T/m]
Gz = Shims.img(:,:,:,3);

%% ------------------------------------------------------------------------
% Simple (voxelwise) Z-shim optimization:

% 1. static:
% 
% Flip field direction and divide by Gz 
staticTarget   = -Field.img ;
staticGzVoxels = staticTarget ./Gz ;

% 2. riro:
%
% Output image from FieldEval.modelfield() is in units of Hz, scaled to the RMS
% PMU value (this normalizes the output to be independent of the PMU range 
% and it enables a direct comparision between the magnitudes of the static and
% RIRO field components)
% 
% flip RIRO polarity + rescale to units of [Hz/unit-PMU]
riroTarget   = - Field.Model.Riro.img/Field.Model.Riro.Aux.Data.p ;
riroGzVoxels = riroTarget ./Gz ;

%% compute slicewise corrections within shimVoi (spinal cord volume)
nSlices = size( Mag.img, 3 ) ;

% static slicewise Gz correction [units: micro-T]
Corrections.static = zeros( nSlices, 1 ) ; 
% RIRO slicewise Gz correction [units: micro-T/unit-PMU]
Corrections.riro   = zeros( nSlices, 1 ) ; 

for iSlice = 1 : nSlices
    
    sliceVoi               = false( size( shimVoi ) ) ;
    sliceVoi( :,:,iSlice ) = shimVoi(:,:,iSlice ) ;

    Corrections.static( iSlice ) = mean( staticGzVoxels( sliceVoi ) ) ;
    Corrections.riro( iSlice )   = mean( riroGzVoxels( sliceVoi ) ) ;

end

%% ------------------------------------------------------------------------
% write to .txt file readable by sequence
% convert static slicewise Gz correction to units milli-T
% convert RIRO slicewise Gz correction to units milli-T/unit-PMU
%% ------------------------------------------------------------------------

fileID = fopen('Dynamic_Gradients.txt','w');

for iSlice = 1:(nSlices)
    fprintf(fileID,'Vector_Gz[0][%i]= %.6f\n', iSlice-1, 1e-3*Corrections.static(iSlice)); 
    fprintf(fileID,'Vector_Gz[1][%i]= %.6f\n', iSlice-1, 1e-3*Corrections.riro(iSlice)); 
    fprintf(fileID,'Vector_Gz[2][%i]= %.3f\n', iSlice-1, Field.Aux.Data.p); 
end

fclose(fileID);

% Generate Shims_static and Shims_static, the static (B(0)) and respiratory
% (c) components of field: B(t) = c*p(t) + B(0)
% see: Topfer et al. MRM 80:935?946 (2018)

%% ------------------------------------------------------------------------
% Generate z-shim values
%
% Gz(t) = c'*p(t)+Gz(0)
% where;
% Gz(0) = Gz_static (static Gz component) (units?)
% c' = Gz_resp (dynaminc Gz component) (should be unitless ...)
%
%
% Note: MGRE sequence (a_gre_DYNshim) will read in p(t) in realtime, and
% for each slice calculate the expected Gz for that pressure value based on
% the above relationship
%% ------------------------------------------------------------------------
