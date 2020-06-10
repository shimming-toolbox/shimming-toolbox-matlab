clear all;
clc;
close all;
%%
disp('Hello');

% hack: install the shimming-toolbox package
addpath(genpath('..'))

data = 'data_testing/'

% --------------------------------------
%% Download data when not already present
if ~isfolder( data )
    url = 'https://osf.io/7d2j5/download?version=4' ;
    fprintf( ['\n Downloading test data...\n URL=' url '\n'] ) ;
    unzip(url, '.') ;
end

tmp = tempname
mkdir(tmp)

% ---------
%% Dcm2Bids
%
niftiPath = fullfile( tmp, 'niftis' )
% TODO : Could add possibility to make own config file : dicom_to_nifti(inputPath, niftiPath, configPath ) 
dicom_to_nifti(fullfile( data, 'acdc_48' ), niftiPath ) 
% dicom_to_nifti(fullfile(data, 'ACDC108p'), niftiPath)
acquisitionPath = fullfile( niftiPath, 'sub-' );
% TODO : Check if there is data

% ----------
%% load data (Could be a function)
% Setting to load data automatically or manually
isManual = true;
disp(acquisitionPath)
ls(acquisitionPath)

if (isManual)
    % List the acquisitions
    acquisitionList = dir(acquisitionPath); % get struct of the folder contents 
    acquisitionList = acquisitionList( ~startsWith( {acquisitionList.name}, '.') ) ; % remove '.' and '..' entries e.g.
    acquisitionList = acquisitionList( [acquisitionList.isdir] ) ;% might also want to remove potential files?
    
    for iDir = 1 : length( acquisitionList )
        fprintf( [ num2str(iDir) ' : ' acquisitionList(iDir).name '\n' ] );
    end
    
    % User input
    isValidInput = false;
    while ~isValidInput
        strMag   = input('Enter the number for the appropriate magnitude fieldmap folder, (type "esc" to quit) : ' , 's');
        if strMag == "esc"
            return
        end
        strPhase = input('Enter the number for the appropriate phase fieldmap folder, (type "esc" to quit) : ' , 's');  
        if strPhase == "esc"
            return
        end
        
        iFolderMag = str2num(strMag) ;
        iFolderPhase = str2num(strPhase) ;
        
        if ~isempty(iFolderMag) && ~isempty(iFolderPhase)
            isValidInput = all( ismember( [iFolderMag iFolderPhase], [1:numel(acquisitionList)] ) );
        end

    end
    folderMag   = acquisitionList(iFolderMag).name;
    folderPhase = acquisitionList(iFolderPhase).name;
        
else
%     folderMag = 'gre_field_mapping_PMUlog_mag'; % ACDC108p
%     folderPhase = 'gre_field_mapping_PMUlog_phase'; % ACDC108p
    folderMag   = 'a_gre_DYNshim_mag'; % dicom_unsorted
    folderPhase = 'a_gre_DYNshim_phase'; % dicom_unsorted
end

% Load mag
listMag = dir(fullfile(acquisitionPath, folderMag, '*.nii*'));
nEchoes = length(listMag);

if nEchoes <= 0 
   error(['No image in acquisition ' folderMag]) 
end

[~ ,tmpInfo, ~] = imutils.read_nii(fullfile( listMag(1).folder , listMag(1).name ));
% preallocation
% mag  = zeros([tmpInfo.ImageSize nEchoes]);
% magInfo =
% magJson = 
if length(tmpInfo.ImageSize) == 3 % If more than one one time
    for iEcho = 1:nEchoes
        % Load and make sure it's mag data
        if ~isempty(strfind(listMag(iEcho).name(end-12:end), '_mag'))
            [mag(:,:,:,:,iEcho), magInfo(iEcho), magJson(iEcho)] = imutils.read_nii( ...
                fullfile( listMag(iEcho).folder , listMag(iEcho).name ) );
        end
    end
else
    for iEcho = 1:nEchoes
        % Load and make sure it's mag data
        if ~isempty(strfind(listMag(iEcho).name(end-12:end), '_mag'))
            [mag(:,:,:,1,iEcho), magInfo(iEcho), magJson(iEcho)] = imutils.read_nii( ...
                fullfile( listMag(iEcho).folder , listMag(iEcho).name ) );
        end
    end
end


% Load phase
listPhase = dir(fullfile(acquisitionPath, folderPhase, '*.nii*'));
nEchoes = length(listPhase);
if nEchoes <= 0
   error(['No image in acquisition ' folderPhase]) 
end

[~ ,tmpInfo, ~] = imutils.read_nii(fullfile( listPhase(1).folder , listPhase(1).name ));
% preallocation
% phase  = zeros([tmpInfo.ImageSize nEchoes]);
% phaseInfo =
% phaseJson = 
if length(tmpInfo.ImageSize) == 3 % If more than one one time
    for iEcho = 1:nEchoes
        % Load and make sure it's phase data
        if ~isempty(strfind(listPhase(iEcho).name(end-12:end), '_phase'))
            [phase(:,:,:,:,iEcho), phaseInfo(iEcho), phaseJson(iEcho)] = imutils.read_nii( ...
            fullfile( listPhase(iEcho).folder , listPhase(iEcho).name ) );
        end
    end
else
    for iEcho = 1:nEchoes
        % Load and make sure it's phase data
        if ~isempty(strfind(listPhase(iEcho).name(end-12:end), '_phase'))
            [phase(:,:,:,1,iEcho), phaseInfo(iEcho), phaseJson(iEcho)] = imutils.read_nii( ...
            fullfile( listPhase(iEcho).folder , listPhase(iEcho).name ) );
        end
    end
end

% -----------------
%% Unwrap (sunwrap)
disp('Unwrap')

% Init Unwrapped Phase
% unwrappedPhase = 
for iAcq = 1:size(phase,4)
    for iEcho = 1:size(phase,5)
        magNorm = mat2gray(mag(:,:,:,iAcq,iEcho));

        %Assumes there are wraps
        phasePi = mat2gray(phase(:,:,:,iAcq,iEcho))*2*pi - pi;

        unwrappedPhase(:,:,:,iAcq,iEcho) = sunwrap(magNorm .* exp( 1i* phasePi ), 0.1);
  
    end
end    
    
% Plot
% figure(1)
% subplot(121)
% imshow(mat2gray(unwrappedPhase(:,:,1,1,1)))
% hold on
% title('unwrapped')
% subplot(122)
% imshow(mat2gray(phase(:,:,1,1,1)))
% title('wrapped')
% hold off

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
% B0FieldmapPlot = reshape(B0Fieldmap, [size(B0Fieldmap, 1) size(B0Fieldmap, 2) 1 size(B0Fieldmap, 3)]); % montage insists on the format M-by-N-by-1-by-K
% figure(2)
% montage(B0Fieldmap,'DisplayRange',[min(min(min(B0FieldmapPlot))) max(max(max(B0FieldmapPlot)))])
% hold on
% colorbar
% title('B0FieldMap (Hz)')
% hold off

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
