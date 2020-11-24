 % generate a "Shepp-Logan" phantom (ie; shape) in 2D
SheppLogan2d_vol = NumericalModel('Shepp-Logan2d',256);
figure; imagesc(SheppLogan2d_vol.starting_volume)
% generate a deltaB0 (B0 inhomogeneity) distribution in 2D (it will be in
% units of Hz)
SheppLogan2d_vol.generate_deltaB0('2d_linearIP', [5 0]); 
% plot it to see how it looks
figure; imagesc(SheppLogan2d_vol.deltaB0)

% simulate MRI magnitude and phase data for your Shepp-Logan phantom in the
% presence of the deltaB0 field you generated
TE = [0.001 0.0015 0.0020 0.0025 0.0030 ]; % echo times
SheppLogan2d_vol.simulate_measurement(15, TE, 100);

% get magnitude and phase data
magn = SheppLogan2d_vol.getMagnitude;
phase = SheppLogan2d_vol.getPhase;

% save the magnitude and phase data as a nifti file
SheppLogan2d_vol.save('Phase', 'SheppLogan2d_simulated_phase.nii');
SheppLogan2d_vol.save('Magnitude', 'SheppLogan2d_simulated_magnitude.nii');

% generate a complex data set from your magnitude and phase data
compl_vol = magn.*exp(1i*phase);

% calculate the deltaB0 map from the magnitude and phase data using the
% SWFM method (it will be in units of Hz)
[dual_echo_delf] = +imutils.b0.dual_echo(compl_vol(:,:,:,:), TE);
% plot it, it should look the same as the deltaB0 you simulated!
figure; imagesc(dual_echo_delf)

% save your deltaB0 map as a nifti file
nii_vol = make_nii(dual_echo_delf);
save_nii(nii_vol, ['dualechoB0_SheppLogan2d' '.nii'])