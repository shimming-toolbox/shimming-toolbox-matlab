function [delf, offset, STDX, MSE] = B0_multiecho_linfit(compl_vol, ph_json)
% function LLSf(compl_vol)
% 
% 
% Inputs: complex data set (compl_vol), phase json (ph_json)
% Outputs: 
%         delf -> field map in Hz 
%         offset -> frequency shift at echo time 0
%         stdx -> standard deviation of the fit
%         mse -> mean square error of the fit

% create magnitude and phase data volumes
mag_data = abs(compl_vol);
ph_data = angle(compl_vol);

% rescale phase - NOTE: this is temporary, scaling should be done bu read_nii but that is not
% working for now
ph_data = rescalePhaseImage(ph_data);

% get number of echoes and echo times
numTE = size(mag_data,4); 
delt = [ph_json.EchoTime];

% pre-allocate memory for variables
delf = zeros(size(ph_data,1),size(ph_data,2),size(ph_data,3));
offset = zeros(size(ph_data,1),size(ph_data,2),size(ph_data,3));
STDX = zeros(size(ph_data,1), size(ph_data,2), size(ph_data,3), 2);
MSE = zeros(size(ph_data,1),size(ph_data,2),size(ph_data,3));

% generate a mask 
sigma = bkgrnd_noise(mag_data);
mask = threshold_masking(mag_data, sigma);

%-------------------------------------------------------------------------%
% temporally unwrap phase across echoes
%-------------------------------------------------------------------------%
delPhaseNet(:,:,:,1) = ph_data(:,:,:,1);

for n=2:numTE

    vol_tmp0 = compl_vol(:,:,:,n-1);
    vol_tmp1 = compl_vol(:,:,:,n);                       

    delPhase = angle( vol_tmp1.*conj(vol_tmp0) );
    delPhaseNet(:,:,:,n) = delPhaseNet(:,:,:,n-1) + delPhase;

end
    

%-------------------------------------------------------------------------%
% frequency shift calculation (field map)
%-------------------------------------------------------------------------%

for i=1:size(ph_data,1)
    for j=1:size(ph_data,2)
        for k = 1:size(ph_data,3)
            if (mask(i,j,k) == 0)
                delf(i,j,k) = 0;
            else
                A =  [delt(1:numTE)',ones(size(delt(1:numTE)'))];
                B = squeeze(delPhaseNet(i,j,k,1:numTE));
                w = squeeze(mag_data(i,j,k,1:numTE)).^2;
                [fit_para_tmp, STDX(i,j,k,:), MSE(i,j,k)] = lscov(A,B,w);
                delf(i,j,k) = fit_para_tmp(1);
                offset(i,j,k) = fit_para_tmp(2);
            end
        end
    end
end


%niftiwrite(delf,'delf_map.nii');

