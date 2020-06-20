classdef (TestTags = {'Simulation', 'Unit'}) NumericalModel_Test < matlab.unittest.TestCase

    properties
        testFileName = 'test'
        testFileNameNii = 'test.nii'
        testFileNameMat = 'test.mat'
    end

    methods (TestClassSetup)
    end

    methods (TestClassTeardown)
        function removeTempFolder(testCase)
            if isfile([testCase.testFileName '.nii'])
                delete([testCase.testFileName '.nii'])
            end
            if isfile([testCase.testFileName '.mat'])
                delete([testCase.testFileName '.mat'])
            end
            
            if isfile([testCase.testFileName '.json'])
                delete([testCase.testFileName '.json'])
            end
        end
    end

    methods (Test)
        %% Class instance tests
        function test_initiate_object_is_expected_class(testCase)

            testObj = NumericalModel();

            expectedClass = 'NumericalModel';
            actualClass = class(testObj);
            testCase.verifyEqual(actualClass, expectedClass)
        end
 
        function test_empty_initialization_returns_expected_starting_volume(testCase)
            testObj = NumericalModel();

            expectedVolume = zeros(128, 128);
            actualVolume = testObj.starting_volume;
            testCase.verifyEqual(actualVolume, expectedVolume)
        end
 
        %% Shepp-Logan type tests instance test
        function test_shepplogan_emtpy_init_returns_expected_starting_volume(testCase)
            testObj = NumericalModel('Shepp-Logan');

            expectedVolume = phantom(128);
            actualVolume = testObj.starting_volume;
            testCase.verifyEqual(actualVolume, expectedVolume)
        end
 
         function test_shepplogan_dims_init_returns_expected_starting_volume(testCase)
            dims = 256;
            testObj = NumericalModel('Shepp-Logan', dims);

            expectedVolume = phantom(dims);
            actualVolume = testObj.starting_volume;
            testCase.verifyEqual(actualVolume, expectedVolume)
        end
        
        function test_shepp_logan_defines_expected_t2star_values(testCase)
            testObj = NumericalModel('Shepp-Logan');

            testCase.assertTrue(all(testObj.volume.T2star(abs(testObj.starting_volume-0.2)<0.001) == testObj.T2star.WM));
            testCase.assertTrue(all(testObj.volume.T2star(abs(testObj.starting_volume-0.3)<0.001) ==  testObj.T2star.GM));
            testCase.assertTrue(all(testObj.volume.T2star(abs(testObj.starting_volume-1)<0.001) == testObj.T2star.CSF));
            testCase.assertTrue(all(testObj.volume.T2star(abs(testObj.starting_volume)<0.001&testObj.starting_volume~=0) == testObj.T2star.WM/2));
            testCase.assertTrue(all(testObj.volume.T2star(abs(testObj.starting_volume-0.1)<0.001) == (testObj.T2star.GM - testObj.T2star.WM/2)/2));
            testCase.assertTrue(all(testObj.volume.T2star(abs(testObj.starting_volume-0.4)<0.001) == testObj.T2star.GM * 1.5));
        end
        
        function test_shepp_logan_defines_expected_magnitude_values(testCase)
            testObj = NumericalModel('Shepp-Logan');

            testCase.assertTrue(all(testObj.volume.protonDensity(abs(testObj.starting_volume-0.2)<0.001) == testObj.protonDensity.WM));
            testCase.assertTrue(all(testObj.volume.protonDensity(abs(testObj.starting_volume-0.3)<0.001) ==  testObj.protonDensity.GM));
            testCase.assertTrue(all(testObj.volume.protonDensity(abs(testObj.starting_volume-1)<0.001) == testObj.protonDensity.CSF));
            testCase.assertTrue(all(testObj.volume.protonDensity(abs(testObj.starting_volume)<0.001&testObj.starting_volume~=0) == testObj.protonDensity.WM/2));
            testCase.assertTrue(all(testObj.volume.protonDensity(abs(testObj.starting_volume-0.1)<0.001) == (testObj.protonDensity.GM - testObj.protonDensity.WM/2)/2));
            testCase.assertTrue(all(testObj.volume.protonDensity(abs(testObj.starting_volume-0.4)<0.001) == testObj.protonDensity.GM * 1.5));
        end

        
        %% simulate_signal method tests
        
        function test_generate_deltaB0_linear_floor_value(testCase)
            testObj = NumericalModel('Shepp-Logan');
            
            m = 0;
            b = 2;
            testObj.generate_deltaB0('linear', [m b]);
            
            deltaB0_map = testObj.deltaB0;
            
            testCase.verifyEqual(mean(deltaB0_map(:)), b/(testObj.gamma/(2*pi)), 'RelTol', 10^-6);
        end
        
        
        function test_generate_deltaB0_linear_slope_value(testCase)
            testObj = NumericalModel('Shepp-Logan');
            
            m = 1;
            b = 0;
            testObj.generate_deltaB0('linear', [m b]);
            
            deltaB0_map = testObj.deltaB0;
            
            dims = size(deltaB0_map);
            [X, Y] = meshgrid(linspace(-dims(1), dims(1), dims(1)), linspace(-dims(2), dims(2), dims(2)));

            testCase.verifyEqual(deltaB0_map(round(dims(1)/2), round(dims(1)/4)), m*X(round(dims(1)/2), round(dims(1)/4))/(testObj.gamma/(2*pi)));
        end
        
        %% simulate_signal method tests
        
        function test_simulate_signal_returns_expected_volume_size(testCase)
            testObj = NumericalModel('Shepp-Logan');
            
            FA = 15;
            TE = [0.003 0.015];
            
            
            testObj.simulate_measurement(FA, TE);
                        
            expectedDims = [128, 128, 1, length(TE)];
            actualDims = size(testObj.measurement);
            testCase.verifyEqual(actualDims, expectedDims);
        end
          
         function test_simulate_signal_get_returns_volume_of_expected_datatype(testCase)
            testObj = NumericalModel('Shepp-Logan');
            
            FA = 15;
            TE = [0.003 0.015];
            
            
            testObj.simulate_measurement(FA, TE);
            
            testCase.assertTrue(isreal(testObj.getMagnitude()));
            testCase.assertTrue(isreal(testObj.getPhase()));
            testCase.assertTrue(isreal(testObj.getReal()));
            testCase.assertTrue(isreal(testObj.getImaginary()));
         end
        
         function test_simulate_signal_SNR_results_in_noisy_backgroun(testCase)
            testObj = NumericalModel('Shepp-Logan');
            
            FA = 15;
            TE = [0.003 0.015];
            SNR = 50;
            
            testObj.simulate_measurement(FA, TE, SNR);
            
            magnitudeData = testObj.getMagnitude;            
            vecMagnitudeROI = magnitudeData(1:10, 1:10, 1, 1);
            vecMagnitudeROI = vecMagnitudeROI(:);

            testCase.assertTrue(std(vecMagnitudeROI)~=0);
            
            phaseData = testObj.getPhase;
            vecPhaseROI = phaseData(1:10, 1:10, 1, 1);
            vecPhaseROI = vecPhaseROI(:);

            testCase.assertTrue(std(vecPhaseROI)~=0);
         end
        
        function test_simulate_signal_dual_echo_calculates_expected_B0(testCase)
            testObj = NumericalModel('Shepp-Logan');
            
            B0_hz = 13;
            testObj.generate_deltaB0('linear', [0.0 B0_hz]); % constant B0 of 10 Hz. deltaB0 is in Tesla.
            TR = 0.025;
            TE = [0.004 0.008];
            testObj.simulate_measurement(TR, TE);
            phaseMeas = testObj.getPhase();
            
            
            phaseTE1 = squeeze(phaseMeas(:,:,1,1));
            phaseTE2 = squeeze(phaseMeas(:,:,1,2));
            
            % Dual echo B0 calculation
            B0_meas = (phaseTE2(64, 64) - phaseTE1(64, 64))/(TE(2) - TE(1));
            B0_meas_hz = B0_meas/(2*pi);
            
            testCase.verifyEqual(B0_hz, B0_meas_hz, 'RelTol', 10^-6);
        end
         
        function test_simulate_signal_dual_echo_righthand_calculates_negative_B0(testCase)
            testObj = NumericalModel('Shepp-Logan');
            testObj.handedness = 'right';

            B0_hz = 13;
            testObj.generate_deltaB0('linear', [0.0 B0_hz]); % constant B0 of 10 Hz. deltaB0 is in Tesla.
            TR = 0.025;
            TE = [0.004 0.008];
            testObj.simulate_measurement(TR, TE);
            phaseMeas = testObj.getPhase();
            
            
            phaseTE1 = squeeze(phaseMeas(:,:,1,1));
            phaseTE2 = squeeze(phaseMeas(:,:,1,2));
            
            % Dual echo B0 calculation
            B0_meas = (phaseTE2(64, 64) - phaseTE1(64, 64))/(TE(2) - TE(1));
            B0_meas_hz = B0_meas/(2*pi);
            
            testCase.verifyEqual(B0_hz, -B0_meas_hz, 'RelTol', 10^-6);
        end
        
        %% save method tests
        function test_save_nii(testCase)
            testObj = NumericalModel('Shepp-Logan');
            
            FA = 15;
            TE = [0.003 0.015];

            testObj.simulate_measurement(FA, TE);
                        
            testObj.save('Phase', testCase.testFileName); % Default option for save is a NIfTI output.
  
            testCase.assertTrue(isfile([testCase.testFileName, '.nii']))
            
            
            % Verify that JSON was written correctly
            testCase.assertTrue(isfile([testCase.testFileName, '.json']))
            
            fname = [testCase.testFileName, '.json'];
            val = jsondecode(fileread(fname));
            testCase.verifyEqual(val.EchoTime', TE)
            testCase.verifyEqual(val.FlipAngle', FA)
            
            
            if isfile([testCase.testFileName '.nii'])
                delete([testCase.testFileName '.nii'])
            end
            if isfile([testCase.testFileName '.json'])
                delete([testCase.testFileName '.json'])
            end
        end

        %% save method tests
        function test_save_nii_with_extension(testCase)
            testObj = NumericalModel('Shepp-Logan');
            
            FA = 15;
            TE = [0.003 0.015];

            testObj.simulate_measurement(FA, TE);
                        
            testObj.save('Phase', [testCase.testFileName, '.nii']); % Default option for save is a NIfTI output.
  
            testCase.assertTrue(isfile([testCase.testFileName, '.nii']))
            
            
            % Verify that JSON was written correctly
            testCase.assertTrue(isfile([testCase.testFileName, '.json']))
            
            fname = [testCase.testFileName, '.json'];
            val = jsondecode(fileread(fname));
            testCase.verifyEqual(val.EchoTime', TE)
            testCase.verifyEqual(val.FlipAngle', FA)
            
            
            if isfile([testCase.testFileName '.nii'])
                delete([testCase.testFileName '.nii'])
            end
            if isfile([testCase.testFileName '.json'])
                delete([testCase.testFileName '.json'])
            end
        end
        
        function test_save_mat(testCase)
            testObj = NumericalModel('Shepp-Logan');
            
            FA = 20;
            TE = [0.003 0.025];

            testObj.simulate_measurement(FA, TE);
                        
            testObj.save('Phase', testCase.testFileName, 'mat');
  
            testCase.assertTrue(isfile([testCase.testFileName, '.mat']))
            
            % Verify that JSON was written correctly
            testCase.assertTrue(isfile([testCase.testFileName, '.json']))
            
            fname = [testCase.testFileName, '.json'];
            val = jsondecode(fileread(fname));
            testCase.verifyEqual(val.EchoTime', TE)
            testCase.verifyEqual(val.FlipAngle', FA)
            
            if isfile([testCase.testFileName '.mat'])
                delete([testCase.testFileName '.mat'])
            end
            if isfile([testCase.testFileName '.json'])
                delete([testCase.testFileName '.json'])
            end
        end
        
        function test_save_mat_with_extension(testCase)
            testObj = NumericalModel('Shepp-Logan');
            
            FA = 20;
            TE = [0.003 0.025];

            testObj.simulate_measurement(FA, TE);
                        
            testObj.save('Phase', [testCase.testFileName, '.mat'], 'mat');
  
            testCase.assertTrue(isfile([testCase.testFileName, '.mat']))
            
            % Verify that JSON was written correctly
            testCase.assertTrue(isfile([testCase.testFileName, '.json']))
            
            fname = [testCase.testFileName, '.json'];
            val = jsondecode(fileread(fname));
            testCase.verifyEqual(val.EchoTime', TE)
            testCase.verifyEqual(val.FlipAngle', FA)
            
            if isfile([testCase.testFileName '.mat'])
                delete([testCase.testFileName '.mat'])
            end
            if isfile([testCase.testFileName '.json'])
                delete([testCase.testFileName '.json'])
            end
        end
        
        %% generate_signal method tests
        function test_generate_signal_case_1(testCase)
            
            protonDensity = 80;
            
            T2star = 100;
            FA = 90;
            TE = 0;
            deltaB0 = 0;
            gamma = 42.58 * 10^6;
            handedness = 'left';
            
            actual_signal = NumericalModel.generate_signal(protonDensity, T2star, FA, TE, deltaB0, gamma, handedness);
            expected_signal = protonDensity;
            
            testCase.verifyEqual(actual_signal, expected_signal);
        end

        function test_generate_signal_case_2(testCase)
            
            FA = 0;
            
            protonDensity = 80;
            T2star = 100;
            TE = 0;
            deltaB0 = 0;
            gamma = 42.58 * 10^6;
            handedness = 'left';

            actual_signal = NumericalModel.generate_signal(protonDensity, T2star, FA, TE, deltaB0, gamma, handedness);
            expected_signal = 0;
            
            testCase.verifyEqual(actual_signal, expected_signal);
        end

        function test_generate_signal_case_3(testCase)
            
            FA = 20;
            protonDensity = 80;
            T2star = 100;
            TE = 0.010;
            deltaB0 = 2;
            gamma = 42.58 * 10^6;
            handedness = 'left';

            switch handedness
                case 'left'
                    sign = -1;
                case 'right'
                    sign = 1;
            end
            
            actual_signal = NumericalModel.generate_signal(protonDensity, T2star, FA, TE, deltaB0, gamma, handedness);
            expected_signal = protonDensity.*sind(FA).*exp(-TE./T2star-sign*1i*gamma*deltaB0.*TE);
            
            testCase.assertTrue(~isreal(expected_signal))
            testCase.verifyEqual(actual_signal, expected_signal);
        end
        
    end

end