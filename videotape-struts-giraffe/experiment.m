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
dicom_to_nifti(fullfile(data, 'dicom_unsorted'), niftiPath)
disp(niftiPath)
ls(niftiPath)

% Load in niftis 
% accounting for the fact that this dataset has separate magnitude/phase channels
% TODO: Switch to (x,y,z,time,Echo)

disp('seperate magnitude and phase')

list = dir(fullfile(niftiPath, 'sub-', 'fmap_mag', '*.nii*'));
% TODO: sort
for i = 1:length(list)
   [mag{i}, magInfo{i}, magJson{i}] = img.read_nii( ...
       fullfile( list(i).folder , list(i).name ) );
end

list = dir(fullfile(niftiPath, 'sub-', 'fmap_phase', '*.nii*'));
% TODO: sort
for i = 1:length(list)
   [phase{i}, phaseInfo{i}, phaseJson{i}] = img.read_nii( ...
       fullfile( list(i).folder , list(i).name ) );
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
