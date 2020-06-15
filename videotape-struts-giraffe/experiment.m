%% Dependencies:
% dcm2niix version v1.0.20200331 (https://github.com/rordenlab/dcm2niix)
% Note: older versions do not properly convert json files for the phase data
% (one element of ImageType is missing).
% dcm2bids (https://github.com/cbedetti/Dcm2Bids)
% 
% 1)a
% cd to directory of experiment.m and enter the following command:
% 
% matlab -nodisplay -nosplash -r 'experiment;exit'
%
% 1)b
% The experiment.m can also be opened in matlab's editor and run there

%%
clear all;
clc;
close all;
%%

% Add shimming-toolbox to the path
addpath(genpath('..'))

% Relative path for the data
data = 'data_testing/';

% --------------------------------------
%% Download data when not already present
if ~isfolder( data )
    url = 'https://osf.io/7d2j5/download?version=4' ;
    fprintf( ['\n Downloading test data...\n URL=' url '\n'] ) ;
    unzip(url, '.') ;
end

% Create temporary folder for processing
tmp = tempname;
mkdir(tmp)

% ---------
%% Dcm2Bids
% Set output path for niftis
niftiPath = fullfile( tmp, 'niftis' )
% Call dicom_to_nifti which uses dcm2niix and dcm2bids to seperate into
% different nifti acquisitions
dicom_to_nifti(fullfile( data, 'dicom_unsorted' ), niftiPath )

% Set path for patient '' (Patient name is set to nothing for anonymity)
acquisitionPath = fullfile( niftiPath, 'sub-' );

% ----------
%% load data

% Display acquisition path 
disp(acquisitionPath)

% Display the different acquisitions
ls(acquisitionPath)

% Call loadniftis which asks the user which acquisition to load.
disp('Enter magnitude data');
[mag, magInfo, magJson] = loadniftis(acquisitionPath);

% Display the size of the input data
size(mag)

% Call loadniftis which asks the user which acquisition to load.
disp('Enter phase data');
[phase, phaseInfo, phaseJson] = loadniftis(acquisitionPath);

% Display the size of the input data
size(phase)

% -----------------
%% Unwrap (sunwrap)
disp('Unwrap')

% Init Unwrapped Phase
% unwrappedPhase = 
for iAcq = 1:size(phase,4)
    for iEcho = 1:size(phase,5)
        % Get the magnitude for a perticular echo
        magNorm = mat2gray(mag(:,:,:,iAcq,iEcho));

        % Calculate the phase in radians, assumes there are wraps
        phasePi = mat2gray(phase(:,:,:,iAcq,iEcho))*2*pi - pi;
        
        % Unwrap phase using sunwrap
        unwrappedPhase(:,:,:,iAcq,iEcho) = sunwrap(magNorm .* exp( 1i* phasePi ), 0.1);
  
    end
end    
    
% Plot
figure(1)
subplot(121)
imshow(mat2gray(unwrappedPhase(:,:,1,1,1)))
hold on
title('unwrapped')
subplot(122)
imshow(mat2gray(phase(:,:,1,1,1)))
title('wrapped')
hold off

% --------------------
%% Process B0 Field map

% Different process if only 1 echo
if size(unwrappedPhase,5) == 1
    echoTimeDiff = phaseJson(1).EchoTime;
    phasediff    = unwrappedPhase(:,:,:,:);
else
    echoTimeDiff = phaseJson(2).EchoTime - phaseJson(1).EchoTime;
    % if using wrapped phase % phasediff = angle( wrappedPhase(1) .* conj(wrappedPhase(2) ) ) ; then unwrap
    phasediff    = unwrappedPhase(:,:,:,:,2) - unwrappedPhase(:,:,:,:,1);
end
    

B0Fieldmap = phasediff./(2*pi*echoTimeDiff);

% Plot
B0FieldmapPlot = reshape(B0Fieldmap, [size(B0Fieldmap, 1) size(B0Fieldmap, 2) 1 size(B0Fieldmap, 3)]); % montage insists on the format M-by-N-by-1-by-K
figure(2)
montage(B0Fieldmap,'DisplayRange',[min(B0FieldmapPlot,[],'all') max(B0FieldmapPlot,[],'all')])
hold on
colorbar
title('B0FieldMap (Hz)')
hold off

% --------------------
%% TODO: Save images
% RT: NOTE: Would be nice to write the results to file (and display the filepaths upon completion).
% `print()` is probably the easiest function to use to print a figure—
% Not sure if it still works in this headless-commandline mode;
% so maybe `imwrite()` is the better call.
% **Ideal** (necessary TODO at some point soon): Write results *as NIfTI*! 
% NOTE: This means ensuring the header info is correct: 
% Should be easy here since no change has been made to positioning;
% Need to be careful about data type & dimension info.
 
disp(['-----'])
% exit;
