% note: adapt ChiDist and FBFest to write json sidecars with info about the
% angle, radius and chi values

% constants
chi_air = 0.36e-6;
chi_mineral_oil = -8.842e-6;
matrix = [128 128 128];
image_res = [1 1 1]; % [mm]
radius = 5; % [mm]
theta = pi/2; 


%%-----------------------------------------------------------------------%%
%% Sphere
%%   - generate susceptibility distribution
%%   - generate deltaB0
%% https://github.com/evaalonsoortiz/Fourier-based-field-estimation
%%-----------------------------------------------------------------------%%

spherical_sus_dist = Spherical( matrix, image_res, radius, [chi_air chi_mineral_oil]);
spherical_sus_dist.save('spherical_sus.nii');

spherical_dBz = FBFest( spherical_sus_dist.volume, spherical_sus_dist.image_res, spherical_sus_dist.matrix, 'Spherical' );
spherical_dBz.save('spherical_dBz.nii');

spherical_vol = NumericalModel('Spherical3d',matrix(1),image_res(1),radius,'air','pure_mineral_oil');
spherical_vol.generate_deltaB0('load_external', 'spherical_dBz.nii');


%%-----------------------------------------------------------------------%%
%% (perpendicular) Cylinder
%%   - generate susceptibility distribution
%%   - generate deltaB0
%% https://github.com/evaalonsoortiz/Fourier-based-field-estimation
%%-----------------------------------------------------------------------%%

cylindrical_sus_dist = Cylindrical( matrix, image_res, radius, theta, [chi_air chi_mineral_oil]);
cylindrical_sus_dist.save('cylindrical_sus.nii');

cylindrical_dBz = FBFest( cylindrical_sus_dist.volume, cylindrical_sus_dist.image_res, cylindrical_sus_dist.matrix, 'Cylindrical' );
cylindrical_dBz.save('cylindrical_dBz.nii');

cylindrical_vol = NumericalModel('Cylindrical3d',matrix(1),image_res(1),radius,(theta*(90/pi)),'air','pure_mineral_oil');
cylindrical_vol.generate_deltaB0('load_external', 'cylindrical_dBz.nii');


%%-----------------------------------------------------------------------%%
%% Modified Zubal Phantom
%%   - generate susceptibility distribution
%%   - generate deltaB0
%% https://github.com/evaalonsoortiz/Fourier-based-field-estimation
%%-----------------------------------------------------------------------%%
zubal_sus_dist = Zubal('../zubal_EAO.nii');
zubal_sus_dist.save('zubal_EAO_sus.nii');

zubal_dBz = FBFest( zubal_sus_dist.volume, zubal_sus_dist.image_res, zubal_sus_dist.matrix, 'Zubal' );
zubal_dBz.save('zubal_EAO_dBz.nii');

% simulate T2* decay for a modified Zubal phantom with a
% deltaB0 found in an external file
zubal_vol = NumericalModel('Zubal','../zubal_EAO.nii');
zubal_vol.generate_deltaB0('load_external', 'zubal_EAO_dBz.nii');


%%-----------------------------------------------------------------------%%
%% Simulate signal decay
%%-----------------------------------------------------------------------%%

% T2* at 3T of pure_mineral_oil is approximated as 0.5*0.063 s
% T1 at 3T of pure_mineral_oil is 181e-3 s

% dual echo
TR = 500e-3;
% Ernst angle for mineral oil = acosd(exp(-TR/T1)) ~= 86
FA = 86; % flip angle [deg]
SNR = 100; % Typical B0 maps SNR values? look it up
% optimal TEs for dual echo are not clear to me. I should play around with
% different timings.
TE = [5e-3 9e-3]; % based on Robinson and Jovicich MRM 2011 66:976-988


% multi-echo: play around with this at the scanner to see what's
% possible
TE = [2e-3 3.5e-3 5e-3 6.5e-3 8e-3 9.5e-3 11e-3 12.5e-3 14e-3];

% bolero
TE = [4e-3 5e-3 6e-3 8e-3 12e-3]; % echo delay times: 0, 1, 2, 4, 8 ms


spherical_vol.simulate_measurement(FA, TE, SNR);

% get magnitude and phase data
magn = spherical_vol.getMagnitude;
phase = spherical_vol.getPhase;
compl_vol = magn.*exp(1i*phase);


%%-----------------------------------------------------------------------%%
%% Compute deltaB0 maps
%%-----------------------------------------------------------------------%%

[dual_echo_delf] = +imutils.b0.dual_echo(compl_vol, TE);
dual_echo_b0_ppm = 1e6*(dual_echo_delf/3)*(1/42.58e6);

% plot results
figure
imagesc(squeeze(dual_echo_b0_ppm(:,:,64)))
colorbar
title('dual-echo fit: b0 (ppm)')

[multi_echo_delf] = +imutils.b0.multiecho_linfit(compl_vol, [0.001 0.002 0.003 0.004 0.005 0.006]); 
multi_echo_b0_ppm = 1e6*(multi_echo_delf/3)*(1/42.58e6);

% plot results
figure
imagesc(squeeze(multi_echo_b0_ppm(:,:,64)))
colorbar
title('multi-echo fit: b0 (ppm)')




figure
imagesc(squeeze(1e6.*real(zubal_dBz.volume(:,:,64))))
colorbar
title('Fourier-based field estimation for the modified Zubal phantom: b0 (ppm)')

% calc diff between dual-echo and multi-echo
diff_dualecho = (dual_echo_b0_ppm-1e6.*real(zubal_dBz.volume));
figure; imagesc(squeeze(diff_dualecho(:,:,64))); colorbar; title('dual echo - true dBz');

diff_multiecho = (multi_echo_b0_ppm-1e6.*real(zubal_dBz.volume));
figure; imagesc(squeeze(diff_multiecho(:,:,64))); colorbar; title('multi echo - true dBz');

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