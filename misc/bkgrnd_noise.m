function sigma = bkgrnd_noise(data_vol)
% Calculate background noise for each slice of an image by computing the 
% average stadard deviation of the signal in four (5 x 5 voxel) corners of 
% an input image (which would presmulably only contain noise) across all 
% echoes
%
% _SYNTAX_
% 
% sigma = bkgrnd_noise(data_vol)
%
% _DESCRIPTION_
%
% _INPUT ARGUMENTS_
%
%    data_vol
%      4D data set data_vol(x,y,z,t)
%
% _OUTPUTS_
%
%   sigma 
%     standard deviation of the noise for each slice in 'data_vol'


% "default" noise level...
sigma(1:size(data_vol,3)) = 1;

% Prepare noise mask
noise_mask = zeros(size(data_vol,1),size(data_vol,2),size(data_vol,3));

% Pick four corners of 5x5
noise_mask(1:5,1:5,:) = 1;
noise_mask(1:5,end-5:end,:) = 1;
noise_mask(end-5:end,1:5,:) = 1;
noise_mask(end-5:end,end-5:end,:) = 1;

% apply mask to the data
    for echo = 1:size(data_vol,4)
        noise_data(:,:,:,echo) = data_vol(:,:,:,echo).*noise_mask(:,:,:); 
    end

noise_data = reshape(noise_data, size(data_vol,1)*size(data_vol,2), size(data_vol,3), size(data_vol,4));

for slice = 1:size(data_vol,3)
    for echo = 1:size(data_vol,4)
        std_noise(slice,echo) = std(nonzeros(noise_data(:,slice,echo)));
    end
end

for slice = 1:size(data_vol,3)
    sigma(slice) = mean(std_noise(slice,:));
end

