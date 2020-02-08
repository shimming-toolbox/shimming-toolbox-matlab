function realtime_zshim(scan_obj, varargin)

%% ************************************************************************
% function realtime_zshim(scan_obj, varargin)
%
% DESCRIPTION: This function will generate static and dynaminc (due to
% respiration) Gz components based on a fieldmap time series (magnitude and
% phase images to be found in 'FM_mag_path' and 'FM_phase_path') and
% respiratory trace information obtained from Siemens bellows 
% (PMUresp_signal.resp). An additional multi-gradient echo (MGRE) magnitiude
% image is used (found in MGRE_mag_path) to generate an ROI and resample
% the static and dynaminc Gz component maps to match the MGRE image. 
% Lastly the average Gz values within the ROI are computed for each slice.
%
% INPUTS: 
%
% (1)  If working with a DICOM socket transfer (online processing mode):
%      realtime_zshim("scan_obj", 0) or realtime_zshim("scan_obj", 1)
%
% (2)  If working with previously sorted DICOMS (offline processing): 
%      realtime_zshim("scan_obj")
%
% -> "scan_obj" should be either 'phantom' or 'human'
%
% 'phantom' -> a central square-shaped ROI will be generated for averaging
% the Gz field
% 'human' -> The spinal cord toolbox (SCT) will be used to segment the spinal cord

% The boolean 0/1 indicates whether or not the files should be moved or copied.

% The user will then be prompted to input the unsorted
% DICOM directory and the desired sorted DICOM directory. 
%
% There will then be a prompt to input the following paths:
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
%
% OUTPUT: text file containing the static and dynamic Gz comnponent values 
% for each slice of the magnitude images found in 'MGRE_mag_path'
% 
% AUTHORS: 
%
% Ryan Topfer, MSc. ryan.topfer@polymtl.ca
% Eva Alonso Ortiz (EAO), PhD. eva.alonso.ortiz@gmail.com 
%
% LAST MODIFIED: Jan. 2020
%
%*************************************************************************


%% ------------------------------------------------------------------------
% Sort DICOM socket transfer images
% enter the directory path for the unsorted dicom images 
% (ex: '/SYNGO_TRANSFER/SYNGO_TRANSFER/20191101.acdc_95p.201911011532/')
% enter the desired path for the sorted images
%% ------------------------------------------------------------------------
if nargin > 1 
    unsortedDicomDir = input('Enter the directory path for the unsorted dicom images (ADD SLASH AT THE END!): ');
    sortedDicomDir = input('Enter the desired path for the sorted images: ');
    if (varargin{1} == 1)
        % copy the files if optional boolean is 1
        MaRdI.sortimages( unsortedDicomDir, sortedDicomDir, 1 );
    elseif (varargin{1} == 0)
        % move the files is optional boolean is 0
        MaRdI.sortimages( unsortedDicomDir, sortedDicomDir );
    end
    % copy respiratory trace file from mounted drive to local directory
    unix('cp /SYNGO_TRANSFER/SYNGO_TRANSFER/PMUresp_signal.resp .')
end
        

%% ------------------------------------------------------------------------
% Read in paths: FM_mag_path, FM_phase_path, MGRE_mag_path, respTrace_path
% These folders will be generated when using 'MaRdI.sortimages'
%% ------------------------------------------------------------------------
% FM_mag_path = input('(ADD SLASH AT THE END!) Field map mag path: ');
% FM_phase_path = input('Field map phase path: ');
% MGRE_mag_path = input('MGRE mag path: ');
respTrace_path = 'PMUresp_signal.resp';

% Hardcode here if needed:
FM_mag_path = '09_gre_field_mapping_PMUlog/';
FM_phase_path = '10_gre_field_mapping_PMUlog/';
MGRE_mag_path = '11_gre_realtime_zshim_NOSHIM/';



%% ------------------------------------------------------------------------
% load MGRE magnitude images for SCT segmentation
%% ------------------------------------------------------------------------
Mag = MaRdI(MGRE_mag_path);
Params.dataLoadDir = MGRE_mag_path;

if (strcmp('human',scan_obj) == 1)
    Params.centerlineMethod = 'spinalcord';  % 'midfov': middle of FOV; 'spinalcord': to create a mask around the spinal cord
elseif (strcmp('phantom',scan_obj) == 1)
        Params.centerlineMethod = 'midfov';  % 'midfov': middle of FOV; 'spinalcord': to create a mask around the spinal cord
else
    error('\nError: Enter "phantom" or "human" as function argument');
end

shimVoi = Mag.segmentspinalcanal_s(Params);

%% ------------------------------------------------------------------------
% load images and respiratory trace recording as FieldEval + ProbeTracking 
% objects :
%% ------------------------------------------------------------------------

% field map time series

Params.unwrapper = 'FslPrelude'; 
Fields = FieldEval( FM_mag_path, FM_phase_path, Params ); 

% pathToNii = './' ; % where to save the images
% isSavingSingleNiis = true ;
% Fields.write( pathToNii, 'nii', isSavingSingleNiis ) ;

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

% image filtering
% for ind = 1:size(Fields.img,5)
%     %Fields.img(:,:,1,1,ind) = imgaussfilt(squeeze(Fields.img(:,:,1,1,ind)),1);
%     Fields.img(:,:,1,1,ind) = imnlmfilt(squeeze(Fields.img(:,:,1,1,ind)),'DegreeOfSmoothing',20);
% end

% Siemens PMU recording
Pmu   = ProbeTracking(respTrace_path);

%% ------------------------------------------------------------------------
% link the two objects (interpolate PMU to the fieldmap time series)
%% ------------------------------------------------------------------------
Fields.associateaux( Pmu );

%% ------------------------------------------------------------------------
% plot the field map time series
%% ------------------------------------------------------------------------

figure

montage(squeeze(Fields.img));
caxis([0 1000])
colorbar
title('Field map time series (Hz)') ;

print('-djpeg','B0_TimeSeries.jpeg');

%% ------------------------------------------------------------------------
% rescale and compute z-gradients  
%% ------------------------------------------------------------------------

% scaling factor 
g = 1000/(42.576E6) ; % [units: mT/Hz] 

%ImageRes = Fields.getvoxelspacing() ; % [units: mm]
[~,Y0,Z0] = Fields.getvoxelpositions() ; % [units: mm]

for measNo = 1:Fields.getnumberofmeasurements
    [~,Fields.img(:,:,1,1,measNo)] = gradient( ...
     squeeze(g*Fields.img(:,:,1,1,measNo)), Y0(1,:)/1000, Z0(:,1)/1000 ) ; % [units: mT/m]
    %squeeze(g*Fields.img(:,:,1,1,measNo)), ImageRes(1,2)/1000, ImageRes(1,3)/1000 ) ; % [units: mT/m]
end

%% ------------------------------------------------------------------------
% plot the Gz map time series
%% ------------------------------------------------------------------------

figure

montage(1e3*squeeze(Fields.img));
caxis([-200 200])
colorbar
title('Gz map time series (\muT/m)') ;

print('-djpeg','Gz_TimeSeries.jpeg');


%% ------------------------------------------------------------------------
% modeled static + respiratory fields (in Field.img and Field.Model.Riro.img respectively)
%% ------------------------------------------------------------------------
Field = FieldEval.modelfield( Fields );


%% ------------------------------------------------------------------------
% rescale and compute z-gradients  
% Eva todo: Modeling the static+resp. fields then taking the
% gradient is not the same thing as taking the gradient and then modeling the
% static and resp. components. Commented this out and added "rescale and 
% compute z-gradients" above.
%% ------------------------------------------------------------------------

% % scaling factor 
% g = 1000/(42.576E6) ; % [units: mT/Hz] 
% 
% [~,Y0,Z0] = Field.getvoxelpositions() ; % [units: mm]
% 
% [~, Field.img]            = gradient( g*Field.img, Y0/1000, Z0/1000 ) ; % [units: mT/m]
% [~, Field.Model.Riro.img] = gradient( g*Field.Model.Riro.img/Field.Model.Riro.Aux.Data.p, Y0/1000, Z0/1000 ) ; % [units: mT/m]


%% ------------------------------------------------------------------------
% plot some results
%% ------------------------------------------------------------------------

figure

subplot(2,1,1);
imagesc( 1e3*Field.img ) ;
axis equal
caxis([-200 200])
colorbar
title('Static Gz (\muT/m)') ;

subplot(2,1,2);
imagesc( 1e3*Field.Model.Riro.img ) ;
axis equal
caxis([-20 20])
colorbar
title('RMS RIRO') ;

print('-djpeg','Gz_map.jpeg');

% EAO: this does not work for these images. Look into it.
% pathToNii = './' ; % where to save the images
% isSavingSingleNiis = true ;
% Field.write( pathToNii, 'nii', isSavingSingleNiis ) ;

%% ------------------------------------------------------------------------
% interpolate field gradient images to target slices for shimming:
%% ------------------------------------------------------------------------
[X,Y,Z]  = Mag.getvoxelpositions() ;

% accelerate interpolation by restricting it to the region where signal exists: 
% (NOTE: mask here could instead be the shimVoi output from SCT)
mask = Mag.getreliabilitymask( 0.05 ) ; % returns mask for each echo (4th dim)
mask = logical( sum( mask, 4 ) ) ; % combine echo-specific masks

% 'interp/extrap' (nearest-neighbour substitution) :
Field.resliceimg( X,Y,Z, mask ) ; % reslice static b0 image 
Field.Model.Riro.resliceimg( X,Y,Z, mask ) ; % reslice RIRO image


%% ------------------------------------------------------------------------
% plot some results
%% ------------------------------------------------------------------------

figure 

for ind = 1:1:size(Field.img,3)
    subplot_tight(2,size(Field.img,3),ind)
    imagesc( 1e3*Field.img(:,:,ind) ) ;
    caxis([-200 200])
    %colorbar
    set(gca, 'XTickLabel', [],'XTick',[])
    set(gca, 'YTickLabel', [],'YTick',[])
    if ind == 1
        title('resampled static Gz (\muT/m)') ;
        set(get(gca,'title'),'Position',[150 0.3 1])
        cb = colorbar('Location','northoutside');
        pos=get(cb,'Position');
        set(cb,'Position',pos+[0.4,0.08,0.1,0]);
    end
end


for ind = 1:1:size(Field.img,3) 
    subplot_tight(2,size(Field.img,3),size(Field.img,3)+ind)
    imagesc( 1e3*Field.Model.Riro.img(:,:,ind) ) ;
    caxis([-20 20])
    %colorbar
    set(gca, 'XTickLabel', [],'XTick',[])
    set(gca, 'YTickLabel', [],'YTick',[])
    if ind == 1
        title('resampled RMS RIRO (\muT/m)') ;
        set(get(gca,'title'),'Position',[150 0.3 1])
                cb = colorbar('Location','northoutside');
        pos=get(cb,'Position');
        set(cb,'Position',pos+[0.4,0.08,0.1,0]);
    end
end

print('-djpeg','resampled_Gz_map.jpeg');

%% ------------------------------------------------------------------------
% Simple (voxelwise) Z-shim calculation:
%% ------------------------------------------------------------------------

% Flip static gradient field polarity (aim is to cancel it)
staticTarget = -Field.img ;

% Flip RIRO gradient polarity + rescale to units of [mT/m/unit-PMU]
riroTarget   = -Field.Model.Riro.img/Field.Model.Riro.Aux.Data.p ;

% slicewise corrections within shimVoi (spinal cord volume)
nSlices = size( Mag.img, 3 ) ;

% static slicewise Gz correction [units: mT/m]
Corrections.static = zeros( nSlices, 1 ) ; 
% RIRO slicewise Gz correction [units: mT/m/unit-PMU]
Corrections.riro   = zeros( nSlices, 1 ) ; 


for iSlice = 1 : nSlices

    sliceVoi = false( size( shimVoi ) ) ;
    sliceVoi( :,:,iSlice ) = shimVoi(:,:, iSlice ) ;

    Corrections.static( iSlice ) = median( staticTarget( sliceVoi ) ) ;
    Corrections.riro( iSlice )   = median( riroTarget( sliceVoi ) ) ;
    
end

%% ------------------------------------------------------------------------
% write to .txt file readable by sequence
%% ------------------------------------------------------------------------

fileID = fopen('zshim_gradients.txt','w');

for iSlice = 1:(nSlices)
    fprintf(fileID,'Vector_Gz[0][%i]= %.6f\n', iSlice-1, Corrections.static(iSlice)); 
    fprintf(fileID,'Vector_Gz[1][%i]= %.12f\n', iSlice-1, Corrections.riro(iSlice)); 
    fprintf(fileID,'Vector_Gz[2][%i]= %.3f\n', iSlice-1, Field.Aux.Data.p); 
end

fclose(fileID);

% unless in offline processing mode (nargin = 1), copy Dynamic_Gradients.txt to mounted drive
if nargin > 1
   unix('cp zshim_gradients.txt /SYNGO_TRANSFER/SYNGO_TRANSFER/')
end

