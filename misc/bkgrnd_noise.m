function sigma = bkgrnd_noise(dataVol)
% Calculate background noise for each slice of an image by computing the 
% average stadard deviation of the signal in four (5 x 5 voxel) corners of 
% an input image (which would presmulably only contain noise) across all 
% echoes
%
% _SYNTAX_
% 
% sigma = bkgrnd_noise(dataVol)
%
% _DESCRIPTION_
%
% _INPUT ARGUMENTS_
%
%    dataVol
%      4D data set dataVol(x,y,z,t)
%
% _OUTPUTS_
%
%   sigma 
%     standard deviation of the noise for each slice in 'dataVol'


% "default" noise level...
sigma(1:size(dataVol,3)) = 1;

% Prepare noise mask
noiseMask = zeros(size(dataVol,1),size(dataVol,2),size(dataVol,3));

% Pick four corners of 5x5
noiseMask(1:5,1:5,:) = 1;
noiseMask(1:5,end-5:end,:) = 1;
noiseMask(end-5:end,1:5,:) = 1;
noiseMask(end-5:end,end-5:end,:) = 1;

% apply mask to the data
    for echo = 1:size(dataVol,4)
        noiseData(:,:,:,echo) = dataVol(:,:,:,echo).*noiseMask(:,:,:); 
    end

noiseData = reshape(noiseData, size(dataVol,1)*size(dataVol,2), size(dataVol,3), size(dataVol,4));

for slice = 1:size(dataVol,3)
    for echo = 1:size(dataVol,4)
        stdNoise(slice,echo) = std(nonzeros(noiseData(:,slice,echo)));
    end
end

for slice = 1:size(dataVol,3)
    sigma(slice) = mean(stdNoise(slice,:));
end

