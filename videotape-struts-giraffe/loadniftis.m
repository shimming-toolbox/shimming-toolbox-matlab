function [niftis, info, json] = loadniftis(path)
%loadniftis Load a nifti acquisition from dcm2bids 
%
%     niftis = loadniftis(path)
%
% If `path` is a folder containing niftis, directly output niftis. If `path` is
% a folder containing acquisitions, ask the user for which acquisition.

if ~isfolder(path)
    error 'Path does not exist'
end

% Get a list of input path
list = dir(path);
list =  list(~ismember({list.name},{'.','..'}));

% Get if path is an acquisition folder or nifti path
if length(unique([list.isdir])) ~= 1
    error ('Directories and files in input path')
else
    isAcquisitions = unique([list.isdir]);
end

% get path
if isAcquisitions
    % Display all acquisitions
    for iDir = 1 : length( list )
        fprintf( [ num2str(iDir) ' : ' list(iDir).name '\n' ] );
    end
    
    % User input
    isValidInput = false;
    while ~isValidInput
        strNiftis   = input('Enter the number for the appropriate acquisition folder, (type "esc" to quit) : ' , 's');
        if strNiftis == "esc"
            return
        end
        
        iFolder = str2num(strNiftis) ;

        if ~isempty(iFolder)
            isValidInput = all( ismember( iFolder, [1:numel(list)] ) );
        end

    end
    
    niftiPath = fullfile(list(iFolder).folder, list(iFolder).name);
    
% if it is a nifti path
else
    niftiPath = path;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO : looking for the number of files is not a rigorous check, channels
% are also in seperate files. Find a method to assign correctly channels and
% echoes. Also output corresponding json and nifti info.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load
niftiList = dir(fullfile(niftiPath, '*.nii*'));
nEchoes = length(niftiList);

% Check number of echoes
if nEchoes <= 0 
   error(['No image in acquisition ' folderMag]) 
end

% Get info for image size
[~ ,tmpInfo, ~] = imutils.read_nii(fullfile( niftiList(1).folder , niftiList(1).name ));

% Format output data according to (x,y,z,time,echoe,channel)
if length(tmpInfo.ImageSize) == 3 % If more than one time
    for iEcho = 1:nEchoes
        [niftis(:,:,:,:,iEcho), info(iEcho), json(iEcho)] = imutils.read_nii( ...
            fullfile( niftiList(iEcho).folder , niftiList(iEcho).name ) );
    end
else
    for iEcho = 1:nEchoes
        [niftis(:,:,:,1,iEcho), info(iEcho), json(iEcho)] = imutils.read_nii( ...
            fullfile( niftiList(iEcho).folder , niftiList(iEcho).name ) );
    end
end

end




