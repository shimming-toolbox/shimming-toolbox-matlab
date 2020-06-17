function B0FieldMap = phaseDifferenceMapping(mag, phase, phaseJson)
%phaseDifferenceMapping Compute B0 fieldmaps using 1 or 2 echos
%     
%     B0FieldMap = phaseDifferenceMapping(mag, phase, phaseJson)
% 
% The inputs `mag`, `phase` and `phaseJson` are outputs of the function 
% `load_niftis` (from the `imutils` package). 

% Unwrapping of the phase images (using sunwrap)
disp('Unwrapping...')
% Init Unwrapped Phase
for iAcq = 1:size(phase,5)
    for iEcho = 1:size(phase,4)
        % Get the magnitude for a specific echo
        magNorm = mat2gray(mag(:,:,:,iEcho,iAcq));
        
        % Calculate the phase in radians, assumes there are wraps
        phasePi = mat2gray(phase(:,:,:,iEcho,iAcq))*2*pi - pi;
        
        % Unwrap phase using sunwrap
        unwrappedPhase(:,:,:,iEcho,iAcq) = sunwrap(magNorm .* exp( 1i* phasePi ), 0.1);
        
    end
end

% Plot
figure(1)
subplot(121)
imshow(mat2gray(unwrappedPhase(:,:,1,1,1)))
hold on
title('unwrapped')
subplot(122)
imshow(mat2gray(phase(:,:,1,1,1)))
title('wrapped')
hold off

% Process B0 Field map
disp('Computing B0 maps...')

% Different process if only 1 echo
if size(unwrappedPhase,4) == 1
    echoTimeDiff = phaseJson(1).EchoTime;
    phaseDiff    = unwrappedPhase(:,:,:,:);
else
    echoTimeDiff = phaseJson(2).EchoTime - phaseJson(1).EchoTime;
    % if using wrapped phase % phaseDiff = angle( wrappedPhase(1) .* conj(wrappedPhase(2) ) ) ; then unwrap
    phaseDiff    = unwrappedPhase(:,:,:,2,:) - unwrappedPhase(:,:,:,1,:);
end

B0FieldMap = phaseDiff./(2*pi*echoTimeDiff);

% Plot
B0FieldMapPlot = reshape(B0FieldMap, [size(B0FieldMap, 1) size(B0FieldMap, 2) 1 size(B0FieldMap, 3)]); % montage insists on the format M-by-N-by-1-by-K
figure(2)
montage(B0FieldMap,'DisplayRange',[min(B0FieldMapPlot,[],'all') max(B0FieldMapPlot,[],'all')])
hold on
colorbar
title('B0FieldMap (Hz)')
hold off

disp(['-----'])
end
