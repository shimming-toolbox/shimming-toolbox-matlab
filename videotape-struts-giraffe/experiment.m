disp('Hello');

% hack: install the shimming-toolbox package
addpath(genpath('..'))

data = 'data_testing/'

%% Download data when not already present
if ~isfolder( data )
    url = 'https://osf.io/7d2j5/download?version=3' ;
    fprintf( ['\n Downloading test data...\n URL=' url '\n'] ) ;
    unzip(url, '.') ;
end

tmp = tempname
mkdir(tmp)

niftiPath = fullfile(tmp, 'niftis')
% dicom_to_nifti(fullfile(data, 'dicom_unsorted'), niftiPath)
dicom_to_nifti(fullfile(data, 'ACDC108p'), niftiPath)
acquistionPath = fullfile(niftiPath, 'sub-');
disp(acquistionPath)
ls(acquistionPath)

%% load data
% TODO: Switch to (x,y,z,time,Echo)
manual = false;
if (manual)
    folderMag = input('Choose the magnitude fieldmap data','s')
    folderPhase = input('Choose the phase fieldmap data','s')
else
    folderMag = 'gre_field_mapping_PMUlog_mag';
    folderPhase = 'gre_field_mapping_PMUlog_phase';
%     folderMag = 'a_gre_DYNshim_mag';
%     folderPhase = 'a_gre_DYNshim_phase';
end


listMag = dir(fullfile(acquistionPath, folderMag, '*.nii*'));
nEchos = length(listMag);
if nEchos <= 0 
   error(['No image in acquisition ' folderMag]) 
end

[~ ,tmpInfo, ~] = img.read_nii(fullfile( listMag(1).folder , listMag(1).name ));
% preallocation
% mag  = zeros([tmpInfo.ImageSize nEchos]);
% magInfo =
% magJson = 
if tmpInfo.ImageSize(4) > 1 % If more than one one time
    for iEcho = 1:nEchos
        % Load and make sure it's mag data
        if ~isempty(strfind(listMag(iEcho).name(end-12:end), '_mag'))
            [mag(:,:,:,:,iEcho), magInfo(iEcho), magJson(iEcho)] = img.read_nii( ...
            fullfile( listMag(iEcho).folder , listMag(iEcho).name ) );
        end
    end
else
    for iEcho = 1:nImgs
        % Load and make sure it's mag data
        if ~isempty(strfind(listMag(iEcho).name(end-12:end), '_mag'))
            [mag(:,:,:,1,iEcho), magInfo(iEcho), magJson(iEcho)] = img.read_nii( ...
            fullfile( listMag(iEcho).folder , listMag(iEcho).name ) );
        end
    end
end


% Load phase
listPhase = dir(fullfile(acquistionPath, folderPhase, '*.nii*'));
nEchos = length(listPhase);
if nEchos <= 0
   error(['No image in acquisition ' folderPhase]) 
end

[~ ,tmpInfo, ~] = img.read_nii(fullfile( listPhase(1).folder , listPhase(1).name ));
% preallocation
% phase  = zeros([tmpInfo.ImageSize nEchos]);
% phaseInfo =
% phaseJson = 
if tmpInfo.ImageSize(4) > 1 % If more than one one time
    for iEcho = 1:nEchos
        % Load and make sure it's phase data
        if ~isempty(strfind(listPhase(iEcho).name(end-12:end), '_phase'))
            [phase(:,:,:,:,iEcho), phaseInfo(iEcho), phaseJson(iEcho)] = img.read_nii( ...
            fullfile( listPhase(iEcho).folder , listPhase(iEcho).name ) );
        end
    end
else
    for iEcho = 1:nEchos
        % Load and make sure it's phase data
        if ~isempty(strfind(listPhase(iEcho).name(end-12:end), '_phase'))
            [phase(:,:,:,1,iEcho), phaseInfo(iEcho), phaseJson(iEcho)] = img.read_nii( ...
            fullfile( listPhase(iImg).folder , listPhase(iImg).name ) );
        end
    end
end


%% Unwrap (sunwrap)
disp('Unwrap')
unwrappedPhase = cell(length(phase),1);
for iUnwrap = 1:length(phase)
    magNorm = mat2gray(mag{iUnwrap});
    
    %Assumes there are wraps
    phasePi = mat2gray(phase{iUnwrap})*2*pi - pi;
    
    unwrappedPhase{iUnwrap} = sunwrap(magNorm .* exp( 1i* phasePi ), 0.1);
%     
%     figure(1)
%     subplot(121)
%     imshow(mat2gray(unwrappedPhase{iUnwrap}(:,:,10)))
%     title('unwrapped')
%     subplot(122)
%     imshow(mat2gray(phase{iUnwrap}(:,:,10)))
%     title('wrapped')
    
end    
    
%% Process B0 Field map
% Assumes nEchoes >= 2 
if length(unwrappedPhase)< 2
    error('Less than 2 echo')
end
    
echoTimeDiff = phaseJson{2}.EchoTime - phaseJson{1}.EchoTime
% if using wrapped phase % phasediff = angle( wrappedPhase{1} .* conj(wrappedPhase{2} ) ) ; then unwrap
phasediff = unwrappedPhase{2} - unwrappedPhase{1};
B0Fieldmap = phasediff./(2*pi*echoTimeDiff);
B0Fieldmap = reshape(B0Fieldmap, [size(B0Fieldmap, 1) size(B0Fieldmap, 2) 1 size(B0Fieldmap, 3)]) % montage insists on the format M-by-N-by-1-by-K
figure(2)
montage(B0Fieldmap,'DisplayRange',[0 300])
colorbar
title('B0FieldMap')


disp(['-----'])
% exit;
