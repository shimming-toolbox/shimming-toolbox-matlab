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
% define the shim system + target field
% shim reference maps for the Siemens Prisma at the IUGM can be downloaded
% from : https://drive.google.com/open?id=1X3kDizzZeZK2dxs6D_zH8bBbjDaF1veu
% it's location should be referenced in: function shimDir = shimbindir( )
%
% Generate Shims_static and Shims_static, the static (B(0)) and respiratory
% (c) components of field: B(t) = c*p(t) + B(0)
% see: Topfer et al. MRM 80:935?946 (2018)
%% ------------------------------------------------------------------------
Shims_static = ShimOpt_IUGM_Prisma_fit( Field );
Shims_resp = ShimOpt_IUGM_Prisma_fit(Field.Model.Riro);

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

Gz_static = Shims_static.img(:,:,:,3);
Gz_resp = Shims_resp.img(:,:,:,3);

zValues = -Field.img ./Gz_static ;

%% ------------------------------------------------------------------------
% Resample Gz_static and Gz_resp to match images found in MGRE_mag_path and
% compute mean-ROI Gz values for each slice
% Todo -> Ryan/EAO
%% ------------------------------------------------------------------------


