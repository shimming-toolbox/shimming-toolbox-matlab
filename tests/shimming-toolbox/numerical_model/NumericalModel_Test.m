classdef (TestTags = {'Simulation', 'Unit'}) NumericalModel_Test < matlab.unittest.TestCase

    properties
    end

    methods (TestClassSetup)
    end

    methods (TestClassTeardown)
    end

    methods (Test)
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
 
        function test_shepplogan_initialization_returns_expected_starting_volume(testCase)
            testObj = NumericalModel('Shepp-Logan');

            expectedVolume = phantom(128, 128);
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

            testCase.assertTrue(all(testObj.volume.magn(abs(testObj.starting_volume-0.2)<0.001) == testObj.protonDensity.WM));
            testCase.assertTrue(all(testObj.volume.magn(abs(testObj.starting_volume-0.3)<0.001) ==  testObj.protonDensity.GM));
            testCase.assertTrue(all(testObj.volume.magn(abs(testObj.starting_volume-1)<0.001) == testObj.protonDensity.CSF));
            testCase.assertTrue(all(testObj.volume.magn(abs(testObj.starting_volume)<0.001&testObj.starting_volume~=0) == testObj.protonDensity.WM/2));
            testCase.assertTrue(all(testObj.volume.magn(abs(testObj.starting_volume-0.1)<0.001) == (testObj.protonDensity.GM - testObj.protonDensity.WM/2)/2));
            testCase.assertTrue(all(testObj.volume.magn(abs(testObj.starting_volume-0.4)<0.001) == testObj.protonDensity.GM * 1.5));
        end
 
        function test_shepp_logan_defines_expected_zero_initial_phase(testCase)
            testObj = NumericalModel('Shepp-Logan');

            testCase.assertTrue(all(all(testObj.volume.phase == 0)));
        end

        function test_generate_signal_case_1(testCase)
            
            protonDensity = 80;
            
            T2star = 100;
            FA = 90;
            TE = 0;
            deltaB0 = 0;
            gamma = 42.58 * 10^6;
            
            actual_signal = NumericalModel.generate_signal(protonDensity, T2star, FA, TE, deltaB0, gamma);
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
            
            actual_signal = NumericalModel.generate_signal(protonDensity, T2star, FA, TE, deltaB0, gamma);
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
            
            actual_signal = NumericalModel.generate_signal(protonDensity, T2star, FA, TE, deltaB0, gamma);
            expected_signal = protonDensity.*sind(FA).*exp(-TE./T2star-1i*gamma*deltaB0.*TE);
            
            testCase.assertTrue(~isreal(expected_signal))
            testCase.verifyEqual(actual_signal, expected_signal);
        end
    end

end