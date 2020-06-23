function mask = threshold_masking(magData, sigma, varargin)
% Create a mask for an input image by thresholding the signal 
%
% _SYNTAX_
% 
% mask = threshold_masking(magData, sigma, varargin)
%
% _DESCRIPTION_
%
% _INPUT ARGUMENTS_
%
%    magData
%      Image from which the mask will be created. Input dimension can be: 
%      (x,y), (x,y,z), (x,y,nEcho) or (x,y,z,nEcho). 
%      In the case of the latter two options, only the first echo time 
%      image is used for thresholding.
%
%    sigma
%      standard deviation of the noise
%
%    optional input: signalThreshold
%      if no value is given signalThreshold = DEFAULT_signalThreshold   
%
% _OUTPUTS_
%
%   mask 
%     Mask with dimension: (x,y,z,nEcho). Note: if input image was 
%     (x,y,nEcho), the function will add a singleton for the "z" dimension 
%     to the output: (x,y,1,nEcho)

DEFAULT_signalThreshold = 30;

if nargin > 2 
    signalThreshold = varargin{1};
else 
    signalThreshold = DEFAULT_signalThreshold;
end

mask = zeros(size(magData));
mask = reshape(mask,size(magData,1)*size(magData,2)*size(magData,3)*size(magData,4),1);
data = reshape(magData,size(magData,1)*size(magData,2)*size(magData,3)*size(magData,4),1);
mask(find(data>signalThreshold*sigma(1))) = 1;
mask = reshape(mask,size(magData,1),size(magData,2),size(magData,3),size(magData,4));
mask = mask(:,:,:,1);
