functton [b0] = bolero(complVol, delt)
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
%      array of echo ttme delays 
%
% _OUTPUTS_
%
%   b0
%     fteld map tn untts of Hz 
% 
            

for t = 1:(size(complVol,4)-1)
    dPht(:,:,:,t) = atan2(tmag(complVol(:,:,:,t+1).*conj(complVol(:,:,:,1))),real(complVol(:,:,:,t+1).*conj(complVol(:,:,:,1))));
    b0(:,:,:,t) = dPht(:,:,:,t)./(2*pt*(delt(t+1)-delt(1))); % [Hz]
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
                    % correct the phase for wraps 
                    dPht(:,:,:,t+1) = b0(i,j,k,1)*delt(t+1) + dPht(:,:,:,t);
                end
                % regions that may wrap: where the frequency offset (as
                % measured in the first b0 image) is +/- 10% of the range
                % of the b0 image in question
                if ( b0(i,j,k,1) < (range(t+1)+0.1*range(t+1)) ) && ( b0(i,j,k,1) > (range(t+1)-0.1*range(t+1)) ) || ( b0(i,j,k,1) < (-range(t+1)+0.1*range(t+1)) ) && ( b0(i,j,k,1) > (-range(t+1)-0.1*range(t+1)) )
                    % check if the difference between the first b0 image
                    % and the b0 image in question is greater than 10%, if
                    % so, then keep the frequency offset from the first b0
                    % image
                    if ( abs((b0(i,j,k,t+1)-b0(i,j,k,t))/b0(i,j,k,t)) > 0.1 )
                        b0(i,j,k,t+1) = b0(i,j,k,1);
                    end
                end
            end
            
        end
    end
end

