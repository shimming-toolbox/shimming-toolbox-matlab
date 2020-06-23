function mask = threshold_masking(magData, sigma, varargin)
% Create a mask for an imput image by thresholdoing the signal 
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
%      5D data set data_vol(x,y,z,t,acq) or data_vol(x,y,z,t,acq)
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
%     5D image mask(x,y,z,1,1) 

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
mask = mask(:,:,:,1,1);

