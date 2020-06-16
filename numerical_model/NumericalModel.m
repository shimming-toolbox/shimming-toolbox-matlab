classdef NumericalModel < handle
    %NUMERICALMODEL Numerical Model for B0 data generation.
    
    properties
        starting_volume
        volume
        
        % Default times in seconds @ 3T
        T2star = struct('WM', 0.053, 'GM', 0.066, 'CSF', 0.10) 
    end
   
    methods
    	function obj = NumericalModel(model)
            switch nargin
                case 0
                    obj.starting_volume = zeros(128, 128);
                case 1
                    if strcmp(model, 'Shepp-Logan')
                        obj.shepp_logan([128, 128]);
                        
                    else
                        error('Unknown volume model.')
                    end
                otherwise
                    error('Too many input arguments for class.')
            end
        end
        
        
        function obj = shepp_logan(obj, dims)          
            obj.starting_volume = phantom(dims(1), dims(2));
            
            obj.volume = struct('magn', [], 'phase', [], 'T2star', []);
            obj.volume.magn = obj.starting_volume;
            obj.volume.phase = obj.volume.magn * 0;
            
            obj.volume.T2star = obj.customize_shepp_logan(obj.starting_volume, obj.T2star.WM, obj.T2star.GM, obj.T2star.CSF);

        end
        
    end

    methods (Access = protected)
        function customVolume = customize_shepp_logan(obj, volume, class1, class2, class3)
            customVolume = volume;
            
            % Set regions to T2
            customVolume(abs(volume-0.2)<0.001) = class1;
            customVolume(abs(volume-0.3)<0.001) = class2;
            customVolume(abs(volume-1)<0.001) = class3;

            customVolume((abs(volume)<0.0001)&volume~=0) = class1/2;
            customVolume(abs(volume-0.1)<0.001) = (class2 - class1/2)/2;
            customVolume(abs(volume-0.4)<0.001) = class2*1.5;
        end
    end

end

