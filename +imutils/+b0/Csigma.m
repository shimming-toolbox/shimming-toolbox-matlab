function sigma = Csigma(dataVol)
% Calculate background noise for each slice of an image by computing the 
% average standard deviation of the signal in four (5 x 5 voxel) corners of 
% an input image (which would presumably only contain noise) across all 
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
%      4D data set dataVol(x,y,z,nEcho)
%
% _OUTPUTS_
%
%   sigma 
%     standard deviation of the noise for each echo in 'dataVol'


% "default" noise level...


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

%noiseData = reshape(noiseData, size(dataVol,1)*size(dataVol,2), size(dataVol,3), size(dataVol,4));
stdNoise = zeros(1, size(dataVol,4) );

for echo = 1:size(dataVol,4)
        if size(dataVol,3)==1
            stdNoise(echo) = std(squeeze(noiseData(:,:,:,echo)),0,'all');
        else
            stdNoise(echo) = std(noiseData(:,:,:,echo),0,'all');
        end
end
sigma=stdNoise


% SNR = zeros(1, size(dataVol,4) );
% for echo = 1:size(dataVol,4)
%         SNR(echo) = snr(dataVol(:,:,:,echo),noiseData(:,:,:,echo));
% end
% SNR
