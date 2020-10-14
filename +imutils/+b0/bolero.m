function [b0] = bolero(complVol, delt)
% bolero computes b0 fteldmaps based on the approach descrtbed tn Mtyasaka
% et al. MRM 2006 55:198-202 and Ktm et al. MRM 2016 
% (B0 Loop Encoded ReadOut "Bolero")
%
% _SYNTAX_
% 
% [b0] = bolero(complVol, delt)
%
% _DESCRIPTION_
%
% _INPUT ARGUMENTS_
%
%    compl_vol
%      complex 4D data set complVol(x,y,z,t)
%    delt
%      array of echo time delays 
%
% _OUTPUTS_
%
%   b0
%     fteld map tn untts of Hz 
% 
            
dPht = zeros(size(complVol));
b0 = zeros(size(complVol));
range = zeros(size(complVol,4));

for t = 1:(size(complVol,4)-1)
    dPht(:,:,:,t) = angle( complVol(:,:,:,t+1).*conj(complVol(:,:,:,1)) );
    %dPht(:,:,:,t) = atan2(imag(complVol(:,:,:,t+1).*conj(complVol(:,:,:,1))),real(complVol(:,:,:,t+1).*conj(complVol(:,:,:,1))));
    b0(:,:,:,t) = dPht(:,:,:,t)./(2*pi*delt(t+1)); % [Hz]
    range(t) = 1/(2*(delt(t+1)-delt(1)));
end


for t = 1:(size(complVol,4)-1)
    for i = 1:size(complVol,1)
        for j = 1:size(complVol,2)
            for k = 1:size(complVol,3)
                % regions that will wrap: where the frequency offset (as
                % determined by the first b0 image) is more that 10% beyond
                % the range of the b0 image in question
                if ( b0(i,j,k,1) > (range(t+1)+0.1*range(t+1)) ) || ( b0(i,j,k,1) < (-range(t+1)-0.1*range(t+1)) )
                    % correct the phase for wraps (note: it's not explained in the paper how this is done, this is my best guess)
                    % it's possible that the "phase correction" is done by:
                    % Phi2 = Phi1 + deltaPhi, where deltaPhi = Phi2 - Phi1
                    % this is how temporal phase unwrapping is done in
                    % multiecho_linfit.m, however, I'm not sure this would
                    % work here because some of the echo time delays are as
                    % long as 8 ms in the publication, in which case you
                    % would have multiple phase wraps between the echoes
                    dPht(i,j,k,t+1) = b0(i,j,k,1)*delt(t+1);
                    b0(i,j,k,t+1) = dPht(i,j,k,t+1)./(2*pi*delt(t+1)); % [Hz]
                end
                % regions that may wrap: where the frequency offset (as
                % measured in the first b0 image) is +/- 10% of the range
                % of the b0 image in question
                if ( b0(i,j,k,1) < (range(t+1)+0.1*range(t+1)) ) && ( b0(i,j,k,1) > (range(t+1)-0.1*range(t+1)) ) || ( b0(i,j,k,1) < (-range(t+1)+0.1*range(t+1)) ) && ( b0(i,j,k,1) > (-range(t+1)-0.1*range(t+1)) )
                    % check if the frequency in the current image differs
                    % by more than 10% when compared with the image with
                    % the next-shortest time, if so, keep the value from 
                    % next-shortest evolution time
                    if ( abs((b0(i,j,k,t+1)-b0(i,j,k,t))/b0(i,j,k,t)) > 0.1 )
                        b0(i,j,k,t+1) = b0(i,j,k,t);
                    end
                end
            end
            
        end
    end
end

