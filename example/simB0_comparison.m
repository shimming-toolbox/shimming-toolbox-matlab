% Simulate a cylindrical susceptibility distribution at 90 degrees using:
% https://github.com/evaalonsoortiz/Fourier-based-field-estimation
cylindrical_sus_dist = Cylindrical( [128 128 128], [1 1 1], 5, pi/2, [0.36e-6 -8.842e-6]);

% save as nifti
cylindrical_sus_dist.save('cylindrical90_R5mm_airMineralOil_ChiDist.nii');

% compute deltaB0 for the simulated susceptibility distribution using:
% https://github.com/evaalonsoortiz/Fourier-based-field-estimation
Bdz = fourier_based_field_est('cylindrical90_R5mm_airMineralOil_ChiDist.nii');

% simulate T2* decay for a cylinder of air surrounded by mineral oil with a
% deltaB0 found in an external file 
cylindrical_vol = NumericalModel('Cylindrical3d',128,1,5,90,'Air', 'SiliconeOil');
cylindrical_vol.generate_deltaB0('load_external', 'Bdz_cylindrical90_R5mm_airMineralOil_ChiDist.nii');
cylindrical_vol.simulate_measurement(15, [0.001 0.002 0.003 0.004 0.005 0.006], 100);

% save magnitude and phase data
magn = cylindrical_vol.getMagnitude;
phase = cylindrical_vol.getPhase;
compl_vol = magn.*exp(1i*phase);
% cylindrical_vol.save('Magnitude','SimMag_cylindrical90_R5mm_airMineralOil.nii');
% cylindrical_vol.save('Phase','SimPh_cylindrical90_R5mm_airMineralOil.nii');

% calculate the deltaB0 map from the magnitude and phase data
[dual_echo_delf] = B0_dual_echo(compl_vol(:,:,:,1:2), [0.001 0.002]);
[multi_echo_delf] = +imutils.b0.multiecho_linfit(compl_vol, [0.001 0.002 0.003 0.004 0.005 0.006]); 

dual_echo_B0_ppm = (dual_echo_delf/3)*(1/42.58e6);
multi_echo_B0_ppm = (multi_echo_delf/3)*(1/42.58e6);

nii_vol = make_nii(imrotate(fliplr(dual_echo_B0_ppm), -90));
save_nii(nii_vol, ['dualechoB0_ppm_cylindrical90_R5mm_airMineralOil' '.nii']);

nii_vol = make_nii(imrotate(fliplr(multi_echo_B0_ppm), -90));
save_nii(nii_vol, ['multiechoB0_ppm_cylindrical90_R5mm_airMineralOil' '.nii']);



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