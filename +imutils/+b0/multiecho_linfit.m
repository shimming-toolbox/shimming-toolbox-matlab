function [delf] = multiecho_linfit(complVol, ph_json)
% multiecho_linfit Computes B0 fieldmaps based on a least-squares
% fitting of the phase evolution with respect to time
%
% _SYNTAX_
% 
% [delf] = multiecho_linfit(complVol, ph_json)
%
% _DESCRIPTION_
%
% _INPUT ARGUMENTS_
%
%    complVol
%      complex 4D data set complVol(x,y,z,t)
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
magData = abs(complVol);
phData = angle(complVol);

% get number of echoes and echo times
numTE = size(magData,4); 
delt = [ph_json.EchoTime];

if ~isrow(delt) 
    delt = delt';
end

if numTE ~= size(delt,2)
    error('\nError: The number of echoes in the dataset does not match the number of echo times in the json sidecar');
end

% pre-allocate memory for variables
delf = zeros(size(phData,1),size(phData,2),size(phData,3));

% generate a mask 
sigma = bkgrnd_noise(magData);
mask = threshold_masking(magData, sigma);

%-------------------------------------------------------------------------%
% temporally unwrap phase across echoes
%-------------------------------------------------------------------------%
delPhaseNet(:,:,:,1) = phData(:,:,:,1);

for n=2:numTE

    volTmp0 = complVol(:,:,:,n-1);
    volTmp1 = complVol(:,:,:,n);                       

    delPhase = angle( volTmp1.*conj(volTmp0) );
    delPhaseNet(:,:,:,n) = delPhaseNet(:,:,:,n-1) + delPhase;

end
    

%-------------------------------------------------------------------------%
% frequency shift calculation (field map)
%-------------------------------------------------------------------------%

for i=1:size(phData,1)
    for j=1:size(phData,2)
        for k = 1:size(phData,3)
            if (mask(i,j,k) == 0)
                delf(i,j,k) = 0;
            else
                A =  [delt(1:numTE)',ones(size(delt(1:numTE)'))];
                B = squeeze(delPhaseNet(i,j,k,1:numTE));
                w = squeeze(magData(i,j,k,1:numTE)).^2;
                [fit_para_tmp] = lscov(A,B,w);
                delf(i,j,k) = fit_para_tmp(1)/(2*pi);
            end
        end
    end
end



