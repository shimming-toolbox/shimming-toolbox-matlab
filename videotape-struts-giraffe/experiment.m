disp('Hello');

% hack: install the shimming-toolbox package
addpath(genpath('..'))

data = 'data_testing/'

%% Download data when not already present
if ~isfolder( data )
    url = 'https://osf.io/7d2j5/?action=download' ;
    fprintf( ['\n Downloading test data...\n URL=' url '\n'] ) ;
    unzip(url) ;
end

tmp = tempname
mkdir(tmp)

nifti_path = fullfile(tmp, 'niftis')
dicom_to_nifti(fullfile(data, 'dicom_unsorted'), nifti_path)
disp(nifti_path)
ls(nifti_path)

% Load in niftis 
% TODO: Switch to (x,y,z,time,Echo)
list = dir(fullfile(nifti_path, ['**' filesep '*.nii*']));
imgs  = cell(length(list),1) ;
infos = cell(length(list),1) ;
jsons = cell(length(list),1) ;
for iImg = 1 : length( list )
   [imgs{iImg}, infos{iImg}, jsons{iImg}] = img.read_nii( ...
       fullfile( list(iImg).folder , list(iImg).name ) );
end

% Seperate in magnitude and phase, potentially doesnt work if other data
% Could load in niftis twice as seperate variables specifying the
% acquisition folder
disp('seperate magnitude and phase')
iMag = 0;
iPhase = 0;
for iList = 1:length(list)
    if ~isempty(strfind(list(iList).name(end-12:end), '_phase'))
        iPhase = iPhase + 1 ;
        phase{iPhase,1} = imgs{iList} ;
        phaseInfo{iPhase,1} = infos{iList} ;
        phaseJson{iPhase,1} = jsons{iList} ;
    else
        iMag = iMag + 1 ;
        mag{iMag,1} = imgs{iList} ;
        magInfo{iMag,1} = infos{iList} ;
        magJson{iMag,1} = jsons{iList} ;
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
    
figure(2)
montage(B0Fieldmap(:,:,:),'DisplayRange',[0 300])
colorbar
title('B0FieldMap')


disp(['-----'])
% exit;
