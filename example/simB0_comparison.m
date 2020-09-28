% generate a cylindrical susceptibility distribution at 90 degrees using:
% https://github.com/evaalonsoortiz/Fourier-based-field-estimation
% cylindrical_sus_dist = Cylindrical( [128 128 128], [1 1 1], 5, pi/2, [0.36e-6 -8.842e-6]);

% generate susceptibility distribution for my modified Zubal phantom
zubal_sus_dist = Zubal('../zubal_EAO.nii');

% save as nifti
% cylindrical_sus_dist.save('cylindrical90_R5mm_airMineralOil_ChiDist.nii');
zubal_sus_dist.save('zubal_EAO_sus.nii');

% compute deltaB0 for the simulated susceptibility distribution using:
% https://github.com/evaalonsoortiz/Fourier-based-field-estimation
zubal_dBz = FBFest( zubal_sus_dist.volume, zubal_sus_dist.image_res, zubal_sus_dist.matrix, 'Zubal' );
zubal_dBz.save('zubal_EAO_dBz.nii');

% simulate T2* decay for a cylinder of air surrounded by mineral oil with a
% deltaB0 found in an external file 
% cylindrical_vol = NumericalModel('Cylindrical3d',128,1,5,90,'Air', 'SiliconeOil');
% cylindrical_vol.generate_deltaB0('load_external', 'Bdz_cylindrical90_R5mm_airMineralOil_ChiDist.nii');
% cylindrical_vol.simulate_measurement(15, [0.001 0.002 0.003 0.004 0.005 0.006], 100);

% simulate T2* decay for a modified Zubal phantom with a
% deltaB0 found in an external file
zubal_vol = NumericalModel('Zubal','../zubal_EAO.nii');
zubal_vol.generate_deltaB0('load_external', 'zubal_EAO_dBz.nii');
%zubal_vol.simulate_measurement(15, [0.001 0.002 0.003 0.004 0.005 0.006], 100);
zubal_vol.simulate_measurement(15, [0.0005 0.001 0.0015 0.002 0.0025 0.003], 100);

% get magnitude and phase data
magn = zubal_vol.getMagnitude;
phase = zubal_vol.getPhase;
compl_vol = magn.*exp(1i*phase);


% calculate the deltaB0 map from the magnitude and phase data
%[dual_echo_delf] = +imutils.b0.dual_echo(compl_vol(:,:,:,1:2), [0.001 0.002]);
[dual_echo_delf] = +imutils.b0.dual_echo(compl_vol(:,:,:,1:2), [0.0005 0.001]);
%[multi_echo_delf] = +imutils.b0.multiecho_linfit(compl_vol, [0.001 0.002 0.003 0.004 0.005 0.006]); 
[multi_echo_delf] = +imutils.b0.multiecho_linfit(compl_vol, [0.0005 0.001 0.0015 0.002 0.0025 0.003]); 

dual_echo_b0_ppm = 1e6*(dual_echo_delf/3)*(1/42.58e6);
multi_echo_b0_ppm = 1e6*(multi_echo_delf/3)*(1/42.58e6);

% plot results
figure
imagesc(imrotate(squeeze(multi_echo_b0_ppm(128,:,:)),90))
colorbar
title('multi-echo fit: b0 (ppm)')
caxis([-5 5]);

figure
imagesc(imrotate(squeeze(dual_echo_b0_ppm(128,:,:)),90))
colorbar
title('dual-echo fit: b0 (ppm)')
caxis([-5 5]);

figure
imagesc(imrotate(squeeze(1e6.*real(zubal_dBz.volume(128,:,:))),90))
colorbar
title('Fourier-based field estimation for the modified Zubal phantom: b0 (ppm)')
caxis([-5 5]);

% calc diff between dual-echo and multi-echo
diff_dualecho = (dual_echo_b0_ppm-1e6.*real(zubal_dBz.volume));
figure; imagesc(imrotate(squeeze(diff_dualecho(128,:,:)),90)); colorbar; title('dual echo - true dBz');

diff_multiecho = (multi_echo_b0_ppm-1e6.*real(zubal_dBz.volume));
figure; imagesc(imrotate(squeeze(diff_multiecho(128,:,:)),90)); colorbar; title('multi echo - true dBz');

% save b0 maps
nii_vol = make_nii(dual_echo_b0_ppm);
save_nii(nii_vol, ['dualechoB0_ppm_zubal' '.nii']);

nii_vol = make_nii(multi_echo_b0_ppm);
save_nii(nii_vol, ['multiechoB0_ppm_zubal' '.nii']);

% nii_vol = make_nii(imrotate(fliplr(dual_echo_B0_ppm), -90));
% save_nii(nii_vol, ['dualechoB0_ppm_cylindrical90_R5mm_airMineralOil' '.nii']);
% 
% nii_vol = make_nii(imrotate(fliplr(multi_echo_B0_ppm), -90));
% save_nii(nii_vol, ['multiechoB0_ppm_cylindrical90_R5mm_airMineralOil' '.nii']);




B0_hz = 500;
TE = [0.0015 0.0025];
a=NumericalModel('Shepp-Logan2d');
a.generate_deltaB0('2d_linearIP', [B0_hz 0]); 
figure; imagesc(a.deltaB0)

a.simulate_measurement(15, TE, 100);

phaseMeas = a.getPhase();
phaseTE1 = squeeze(phaseMeas(:,:,1,1));
phaseTE2 = squeeze(phaseMeas(:,:,1,2));

B0_meas = (phaseTE2(:, :) - phaseTE1(:, :))/(TE(2) - TE(1));
B0_meas_hz = B0_meas/(2*pi);
figure; imagesc(B0_meas_hz)