function LLSf(mag_path, ph_path)

% open data files
[mag_data, mag_info, mag_json] = imutils.load_niftis(mag_path);
[ph_data, ph_info, ph_json] = imutils.load_niftis(ph_path);


%[mask_img,mask_info,mask_json] = imutils.read_nii(mask_fname);


% create complex data volume
vol_compl(:,:,:,:) = mag_data(:,:,:,:).*exp(1i*ph_data(:,:,:,:));

% temporally unwrap phase across echoes
delPhaseNet = temporal_phUnwrap(vol_compl);

disp('test');


% frequency shift calculation (field map)
%[delf, offset, STDX, MSE] = linfitFrequency(echo_times, delPhaseNet, num_echoes, mag, mask);

%-------------------------------------------------------------------------%
% Function definitions
%-------------------------------------------------------------------------%

    function delPhaseNet = temporal_phUnwrap(vol_compl)
    %  Temporally unwrap phase data
    %
    %  Input: vol_compl: complex dataset (mag.*exp(1i*phase))
    %  Output: delPhaseNet: temporally unwrapped phase data
    %
    % Written by Avery Berman (ajberman@mgh.harvard.edu)

    num_echoes = size(vol_compl,4);

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

    function [delf, offset, STDX, MSE] = linfitFrequency(delt, delPhaseNet, num_echoes, mag, mask)
    %LINFITFREQUENCY3D linearly fits the change in phase vs. echo time
    %   Detailed explanation goes here
    % Input: delt: echo time 
    %        delPhaseNet: temporarily unwrapped phase data
    %        num_echoes: number of echoes used in the fitting
    %        mag: magnitude dataset
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

end
