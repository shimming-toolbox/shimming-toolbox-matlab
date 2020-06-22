function mask = threshold_masking(mag_data, sigma, varargin)

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