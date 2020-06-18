function mask = threshold_masking(mag_data, sigma)

mask = zeros(size(mag_data));
mask = reshape(mask,size(mag_data,1)*size(mag_data,2)*size(mag_data,3)*size(mag_data,4),1);
data = reshape(mag_data,size(mag_data,1)*size(mag_data,2)*size(mag_data,3)*size(mag_data,4),1);
mask(find(data>100*sigma(1))) = 1;
mask = reshape(mask,size(mag_data,1),size(mag_data,2),size(mag_data,3),size(mag_data,4));
mask = mask(:,:,:,1,1);

end