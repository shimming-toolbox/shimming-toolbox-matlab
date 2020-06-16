classdef NumericalModel
    %NUMERICALMODEL Numerical Model for B0 data generation.
    
    properties
        starting_volume
    end
   
    methods
    	function obj = NumericalModel(model)
            switch nargin
                case 0
                    obj.starting_volume = zeros(128, 128);
                case 1
                    if strcmp(model, 'Shepp-Logan')
                        obj.starting_volume = phantom(128, 128);
                    else
                        error('Unknown volume model.')
                    end
                otherwise
                    error('Too many input arguments for class.')
            end

        end
        

    end
end

