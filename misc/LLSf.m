function LLSf(mag_path, ph_path)

% open data files
[mag_data, mag_info, mag_json] = imutils.load_niftis(mag_path);
[ph_data, ph_info, ph_json] = imutils.load_niftis(ph_path);

% note the echo index will be 4 once this branch in merged with master
num_echoes = size(mag_data,5);
%echo_times = 

%[mask_img,mask_info,mask_json] = imutils.read_nii(mask_fname);


% create complex data volume
vol_compl(:,:,:,:) = double(mag_data(:,:,:,1,:)).*exp(1i*double(ph_data(:,:,:,1,:)));

% temporally unwrap phase across echoes
delPhaseNet = temporal_phUnwrap;%(vol_compl);

% frequency shift calculation (field map)
[delf, offset, STDX, MSE] = linfitFrequency; %(echo_times, delPhaseNet, num_echoes, mag, mask);

%-------------------------------------------------------------------------%
% Function definitions
%-------------------------------------------------------------------------%

    function delPhaseNet = temporal_phUnwrap %(vol_compl)
    %  Temporally unwrap phase data
    %
    %  Input: vol_compl: complex dataset (mag.*exp(1i*phase))
    %  Output: delPhaseNet: temporally unwrapped phase data
    %
    % Written by Avery Berman (ajberman@mgh.harvard.edu)

    
    %num_echoes = size(vol_compl,4);

    delPhaseNet(:,:,:,1) = angle(vol_compl(:,:,:,1));

    for n=2:num_echoes

        vol_tmp0 = vol_compl(:,:,:,n-1);
        vol_tmp1 = vol_compl(:,:,:,n);

    %               ^ Re
    %               |
    %               |  ^
    %               | /
    %               |/theta
    %               --------> Im                       

        delPhase = angle( vol_tmp1.*conj(vol_tmp0) );
        delPhaseNet(:,:,:,n) = delPhaseNet(:,:,:,n-1) + delPhase;

    end

    end

    function [delf, offset, STDX, MSE] = linfitFrequency %(delt, delPhaseNet, num_echoes, mag, mask)
    %LINFITFREQUENCY3D linearly fits the change in phase vs. echo time
   
    % Ouptput: delf: frequency shift map (field map)'
    %          offset: frequency shift at echo time 0
    %          stdx: standard deviation of the fit
    %          mse: mean square error of the fit
    % Written by Avery Berman (ajberman@mgh.harvard.edu)
    % Modified by Yuhan Ma (yuhanma1@gmail.com)
    % Modification: Use my linearFit instead of Avery's linFit. Accept a 3D
    % volume as an input
    % Modification: use matlab built-in function lscov
    % Last Modified July 2013


    if max(reshape(mask,size(mask,1)*size(mask,2),size(mask,3))) > 1
        mask(find(mask > 0)) = 1;
    end

    delf = zeros(size(delPhaseNet,1),size(delPhaseNet,2),size(delPhaseNet,3));
    offset = zeros(size(delPhaseNet,1),size(delPhaseNet,2),size(delPhaseNet,3));
    STDX = zeros(size(delPhaseNet,1), size(delPhaseNet,2), size(delPhaseNet,3), 2);
    MSE = zeros(size(delPhaseNet,1),size(delPhaseNet,2),size(delPhaseNet,3));

    for i=1:size(delPhaseNet,1)
        for j=1:size(delPhaseNet,2)
            for k = 1:size(delPhaseNet,3)
                if (mask(i,j,k) == 0)
                    delf(i,j,k) = 0;
                else
                    A_design =  [delt(1:num_echoes)',ones(size(delt(1:num_echoes)'))];
                    b = squeeze(delPhaseNet(i,j,k,1:num_echoes));
                    weight = squeeze(mag(i,j,k,1:num_echoes)).^2;
                    [fit_para_tmp, STDX(i,j,k,:), MSE(i,j,k)] = lscov(A_design,b, weight);
                    delf(i,j,k) = fit_para_tmp(1);
                    offset(i,j,k) = fit_para_tmp(2);
                end
            end
        end
    end

    end

    function threshold_mask

    % Calculate background noise 
    sigma = calc_bkgrnd_noise(mag_data);

    mask = zeros(size(mag_data));
    mask = reshape(mask,data_dim(1,1)*data_dim(1,2)*data_dim(1,3)*data_dim(1,4),1);
    data = reshape(data,data_dim(1,1)*data_dim(1,2)*data_dim(1,3)*data_dim(1,4),1);
    mask(find(data>100*sigma(1))) = 1;
    mask = reshape(mask,data_dim);
    mask = mask(:,:,:,1);

    end

    function sigma = calc_bkgrnd_noise(data_vol, data_dim)

    %--------------------------------------------------------------------------
    % Calculate background noise 

    % "default" noise level...
    sigma(1:slices) = 1;

    % Prepare noise mask
    noise_mask = zeros(height,width,slices);

    % Pick four corners of 5x5
    noise_mask(1:5,1:5,:) = 1;
    noise_mask(1:5,end-5:end,:) = 1;
    noise_mask(end-5:end,1:5,:) = 1;
    noise_mask(end-5:end,end-5:end,:) = 1;

    % apply mask to the data
    for slice = 1:slices
        for frame = 1:num_echoes
            noise_data(:,:,slice,frame) = data_vol(:,:,slice,frame).*noise_mask(:,:,slice); 
        end
    end

    noise_data = reshape(noise_data, voxels, slices, num_echoes);

    for slice = 1:slices
        for frame = 1:num_echoes
            std_noise(slice,frame) = std(nonzeros(noise_data(:,slice,frame)));
        end
    end

    for slice = 1:slices
        sigma(slice) = mean(std_noise(slice,:));
    end

    end

end
