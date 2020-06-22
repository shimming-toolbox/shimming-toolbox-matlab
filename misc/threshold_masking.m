function mask = threshold_masking(mag_data, sigma, varargin)
% Create a mask for an imput image by thresholdoing the signal 
%
% _SYNTAX_
% 
% mask = threshold_masking(mag_data, sigma, varargin)
%
% _DESCRIPTION_
%
% _INPUT ARGUMENTS_
%
%    mag_data
%      4D or 5D data set data_vol(x,y,z,t) or data_vol(x,y,z,t,acq)
%
%    sigma
%      standard deviation of the noise
%
%    optional input: signal_threshold
%      if no value is given signal_threshold = DEFAULT_SIGNAL_THRESHOLD   
%
% _OUTPUTS_
%
%   mask 
%     5D image mask(x,y,z,1,1) 

DEFAULT_SIGNAL_THRESHOLD = 30;

if nargin > 2 
    signal_threshold = varargin{1};
else 
    signal_threshold = DEFAULT_SIGNAL_THRESHOLD;
end

mask = zeros(size(mag_data));
mask = reshape(mask,size(mag_data,1)*size(mag_data,2)*size(mag_data,3)*size(mag_data,4),1);
data = reshape(mag_data,size(mag_data,1)*size(mag_data,2)*size(mag_data,3)*size(mag_data,4),1);
mask(find(data>signal_threshold*sigma(1))) = 1;
mask = reshape(mask,size(mag_data,1),size(mag_data,2),size(mag_data,3),size(mag_data,4));
mask = mask(:,:,:,1,1);

end