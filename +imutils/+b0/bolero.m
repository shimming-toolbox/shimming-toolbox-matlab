function [delf] = bolero(complVol, delt)
% bolero (B0 loop encoded readout) 
% Computes B0 fieldmaps based on the method descirbed in Miyasaka et al. 
% Magnetic Resonance in Medicine 55:198-202 (2006)
%
% NOTES
%
% _SYNTAX_
% 
% [delf] = separate_excitation(complVol, delt)
%
% _DESCRIPTION_
%
% _INPUT ARGUMENTS_
%
%    complVol
%      complex 4D data set complVol(x,y,z,t)
%    delt
%      array of phase evolution times in [s]
%
% _OUTPUTS_
%
%   delf 
%     field map in units of [Hz]
% 

% create magnitude and phase data volumes
magData = abs(complVol);
phData = angle(complVol);

% get number of echoes 
numTE = size(magData,4); 

if ~isrow(delt) 
    delt = delt';
end

if numTE ~= size(delt,2)
    error('\nError: The number of echoes in the dataset does not match the number of phase evolution times');
end

% pre-allocate memory for variables
delf = zeros(size(phData,1),size(phData,2),size(phData,3));
B01 = zeros(size(phData,1),size(phData,2),size(phData,3));
B02 = zeros(size(phData,1),size(phData,2),size(phData,3));
B03 = zeros(size(phData,1),size(phData,2),size(phData,3));
B04 = zeros(size(phData,1),size(phData,2),size(phData,3));
B05 = zeros(size(phData,1),size(phData,2),size(phData,3));

B01 = 

Z1(:,:,:) = mag_data(:,:,:,1).*exp(1i*ph_data(:,:,:,1));
Z2(:,:,:) = mag_data(:,:,:,2).*exp(1i*ph_data(:,:,:,2));
dPhi(:,:,:) = atan2(imag(Z1(:,:,:).*conj(Z2(:,:,:))),real(Z1(:,:,:).*conj(Z2(:,:,:))));
b0(:,:,:) = dPhi(:,:,:)./(delt(2)-delt(1)); % [rad*Hz]
b0 = b0/(2*pi); % [Hz]

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



