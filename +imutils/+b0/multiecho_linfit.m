function [delf] = multiecho_linfit(compl_vol, ph_json)
% multiecho_linfit Computes B0 fieldmaps based on a least-squares
% fitting of the phase evolution with respect to time
%
% _SYNTAX_
% 
% [delf] = multiecho_linfit(compl_vol, ph_json)
%
% _DESCRIPTION_
%
% _INPUT ARGUMENTS_
%
%    compl_vol
%      complex 4D data set compl_vol(x,y,z,t)
%    ph_json
%      phase json sidecar struct
%
% _OUTPUTS_
%
%   delf 
%     field map in units of [Hz]
% 
% Code adapted from https://github.com/evaalonsoortiz/LLSf_B0mapping

% create magnitude and phase data volumes
mag_data = abs(compl_vol);
ph_data = angle(compl_vol);

% get number of echoes and echo times
numTE = size(mag_data,4); 
delt = [ph_json.EchoTime];

if ~isrow(delt) 
    delt = delt';
end

% pre-allocate memory for variables
delf = zeros(size(ph_data,1),size(ph_data,2),size(ph_data,3));

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
                [fit_para_tmp] = lscov(A,B,w);
                delf(i,j,k) = fit_para_tmp(1)/(2*pi);
            end
        end
    end
end



