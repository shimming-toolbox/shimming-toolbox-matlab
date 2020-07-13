classdef NumericalModel < handle
    %NUMERICALMODEL Numerical Model for B0 data generation.
    
    properties
        gamma = 267.52218744 * 10^6; % rad*Hz/Tesla
        fieldStrength = 3.0; % Tesla
        handedness = 'left'; % Siemens & Canon = 'left', GE & Philips = 'right' 
        
        % volume properties
        numVox = 128; % square volume
        pixSize = 1; % default pixel size is 1 mm
        starting_volume
        volume
        measurement
        
        % pulse sequence properties
        TE
        FA
        
        % Default times in seconds @ 3T (Note: Silicone and Mineral oil
        % values are actaully T2, will figure something out to address
        % this)
        T2star = struct('WM', 0.053, 'GM', 0.066, 'CSF', 0.10, 'Air', 0, 'SiliconeOil', 0.566, 'PureMineralOil', 0.063)
        
        % Default proton density in percentage
        protonDensity = struct('WM', 70, 'GM', 82, 'CSF', 100, 'Air', 0, 'SiliconeOil', 71, 'PureMineralOil', 100)
        
        % deltaB0 in Hz
        deltaB0
    end
   
    methods
    	function obj = NumericalModel(model, numVox, varargin)
            % NumericalModel class
            % model: string with the label of the model desired.
            % varargin:  
            %   if 'Spherical3d': iamge_res in [mm], radius in [mm],
            %   material inside the sphere ('Air', 'SiliconeOil' or 
            %   'PureMineralOil'), material outside the sphere ('Air,
            %   'SiliconeOil, or 'PureMineralOil')
            %
            %   if 'Cylindrical3d': iamge_res in [mm], radius in [mm],
            %   theta in [degrees] (angle between principal axis of
            %   cylinder and B0), material inside the sphere ('Air', 
            %   'SiliconeOil' or 'PureMineralOil'), material outside the 
            %   sphere ('Air, 'SiliconeOil, or 'PureMineralOil')
            
            if exist('numVox', 'var')
            	obj.numVox = numVox;
            end
            
            if exist('model', 'var')
                switch model
                    case 'Shepp-Logan2d'
                    	obj.shepp_logan_2d(obj.numVox);
                    case 'Shepp-Logan3d'
                        obj.shepp_logan_3d(obj.numVox);
                    case 'Spherical3d'
                        obj.starting_volume = zeros(obj.numVox, obj.numVox, obj.numVox);
                        obj.spherical_3d(varargin);
                    case 'Cylindrical3d'
                        obj.starting_volume = zeros(obj.numVox, obj.numVox, obj.numVox);
                        obj.cylindrical_3d(varargin);
                    otherwise
                        error('Unknown volume model.')
                end
            else
                obj.starting_volume = zeros(obj.numVox, obj.numVox);
            end
            
            % Define background field
            obj.deltaB0 = obj.starting_volume * 0;
        end
        
        function obj = shepp_logan_2d(obj, numVox)       
            % Create a 2D Shepp_Logan volume for the 
            % dims: [x, y] number of voxels.
            obj.starting_volume = phantom(numVox);
            
            obj.volume = struct('magn', [], 'phase', [], 'T2star', []);

            obj.volume.protonDensity = obj.customize_shepp_logan(obj.starting_volume, obj.protonDensity.WM, obj.protonDensity.GM, obj.protonDensity.CSF);
            obj.volume.T2star = obj.customize_shepp_logan(obj.starting_volume, obj.T2star.WM, obj.T2star.GM, obj.T2star.CSF);
        end
        
        function obj = shepp_logan_3d(obj, numVox)
            % Create a 3D Shepp_Logan volume for the
            % dims: [x, y, z] number of voxels.
            obj.starting_volume = phantom3d(numVox);
            
            obj.volume = struct('magn', [], 'phase', [], 'T2star', []);
            
            obj.volume.protonDensity = obj.customize_shepp_logan(obj.starting_volume, obj.protonDensity.WM, obj.protonDensity.GM, obj.protonDensity.CSF);
            obj.volume.T2star = obj.customize_shepp_logan(obj.starting_volume, obj.T2star.WM, obj.T2star.GM, obj.T2star.CSF);
        end
        
        function obj = spherical_3d(obj, varargin)
            % Create a 3D spherical volume for the
            % dims: [x, y, z] number of voxels.
            
            % define local variables
            matrix = [obj.numVox obj.numVox obj.numVox];
            image_res = [varargin{1,1}{1,1} varargin{1,1}{1,1} varargin{1,1}{1,1}];
            R = varargin{1,1}{1,2}; % [mm]
            
            % define image grid
            [x,y,z] = ndgrid(linspace(-(matrix(1)-1)/2,(matrix(1)-1)/2,matrix(1)),linspace(-(matrix(2)-1)/2,(matrix(2)-1)/2,matrix(2)),linspace(-(matrix(3)-1)/2,(matrix(3)-1)/2,matrix(3)));
            
            % define radial position (in [mm])
            r = sqrt((x.*image_res(1)).^2 + (y.*image_res(2)).^2 + (z.*image_res(3)).^2);
            
            obj.volume = struct('magn', [], 'phase', [], 'T2star', [], 'protonDensity', []);
            
            obj.volume.protonDensity = obj.starting_volume;
            obj.volume.T2star = obj.starting_volume;
            
            switch varargin{1,1}{1,3}
                case 'Air'
                    obj.volume.protonDensity(r <= R ) = obj.protonDensity.Air;
                    obj.volume.T2star(r <= R ) = obj.T2star.Air;
                case 'SiliconeOil'
                    obj.volume.protonDensity(r <= R ) = obj.protonDensity.SiliconeOil;
                    obj.volume.T2star(r <= R ) = obj.T2star.SiliconeOil;
                case 'PureMineralOil'
                    obj.volume.protonDensity(r <= R ) = obj.protonDensity.PureMineralOil;
                    obj.volume.T2star(r <= R ) = obj.T2star.PureMineralOil;     
            end
               
            switch varargin{1,1}{1,4}
                case 'Air'
                    obj.volume.protonDensity(r > R ) = obj.protonDensity.Air;
                    obj.volume.T2star(r > R ) = obj.T2star.Air;
                case 'SiliconeOil'
                    obj.volume.protonDensity(r > R ) = obj.protonDensity.SiliconeOil;
                    obj.volume.T2star(r > R ) = obj.T2star.SiliconeOil;
                case 'PureMineralOil'
                    obj.volume.protonDensity(r > R ) = obj.protonDensity.PureMineralOil;
                    obj.volume.T2star(r > R ) = obj.T2star.PureMineralOil;     
            end
                 
        end
        
        function obj = cylindrical_3d(obj, varargin)
            % Create a 3D cylindrical volume for the
            % dims: [x, y, z] number of voxels.
            
            % define local variables
            matrix = [obj.numVox obj.numVox obj.numVox];
            image_res = [varargin{1,1}{1,1} varargin{1,1}{1,1} varargin{1,1}{1,1}];
            R = varargin{1,1}{1,2}; % [mm]
            theta = varargin{1,1}{1,3}; % angle between main axis of cylinder and z-axis (in degrees)
            
            % define image grid
            [x,y,z] = ndgrid(linspace(-(matrix(1)-1)/2,(matrix(1)-1)/2,matrix(1)),linspace(-(matrix(2)-1)/2,(matrix(2)-1)/2,matrix(2)),linspace(-(matrix(3)-1)/2,(matrix(3)-1)/2,matrix(3)));
            
            % define radial position (in [mm])
            r = sqrt((x.*image_res(1)).^2 + (y.*image_res(2)).^2);
            
            obj.volume = struct('magn', [], 'phase', [], 'T2star', [], 'protonDensity', []);
            
            obj.volume.protonDensity = obj.starting_volume;
            obj.volume.T2star = obj.starting_volume;
            
            switch varargin{1,1}{1,4}
                case 'Air'
                    obj.volume.protonDensity(r <= R ) = obj.protonDensity.Air;
                    obj.volume.T2star(r <= R ) = obj.T2star.Air;
                case 'SiliconeOil'
                    obj.volume.protonDensity(r <= R ) = obj.protonDensity.SiliconeOil;
                    obj.volume.T2star(r <= R ) = obj.T2star.SiliconeOil;
                case 'PureMineralOil'
                    obj.volume.protonDensity(r <= R ) = obj.protonDensity.PureMineralOil;
                    obj.volume.T2star(r <= R ) = obj.T2star.PureMineralOil;     
            end
               
            switch varargin{1,1}{1,5}
                case 'Air'
                    obj.volume.protonDensity(r > R ) = obj.protonDensity.Air;
                    obj.volume.T2star(r > R ) = obj.T2star.Air;
                case 'SiliconeOil'
                    obj.volume.protonDensity(r > R ) = obj.protonDensity.SiliconeOil;
                    obj.volume.T2star(r > R ) = obj.T2star.SiliconeOil;
                case 'PureMineralOil'
                    obj.volume.protonDensity(r > R ) = obj.protonDensity.PureMineralOil;
                    obj.volume.T2star(r > R ) = obj.T2star.PureMineralOil;     
            end
            
            % rotate chi distribution about the y-axis
            t = [cosd(theta)   0      -sind(theta)   0
                 0             1              0      0
                 sind(theta)   0       cosd(theta)   0
                 0             0              0      1];
            tform = affine3d(t);
            obj.volume.T2star = imwarp(obj.volume.T2star,tform);
            obj.volume.protonDensity = imwarp(obj.volume.protonDensity,tform);
            
        end

        function obj = simulate_measurement(obj, FA, TE, SNR)          
            % FA: flip angle value in degrees.
            % TE: Single or array TE value in seconds.
            % SNR: Signal-to-noise ratio
            
            % Set attributes
            obj.FA = FA;
            obj.TE = TE;
            
            % Get dimensions
            numTE = length(TE);
            volDims = size(obj.starting_volume);
               
            % Pre-allocate measurement variables
            if volDims == 2
                obj.measurement = zeros(volDims(1), volDims(2), 1, numTE);
            elseif volDims == 3
                obj.measurement = zeros(volDims(1), volDims(2), volDims(3), numTE);
            end
            
            % Simulate
            for ii = 1:numTE
                obj.measurement(:,:,:,ii) = obj.generate_signal(obj.volume.protonDensity, obj.volume.T2star, FA, TE(ii), obj.deltaB0, obj.gamma, obj.handedness);
            end
            
            if exist('SNR','var')
                obj.measurement = NumericalModel.addNoise(obj.measurement, SNR);
            end
        end
        
        function vol = getMagnitude(obj)
            % Get magnitude data
            vol = abs(obj.measurement);
        end
        
        function vol = getPhase(obj)
            % Get phase data
            vol = angle(obj.measurement);
        end
  
        function vol = getReal(obj)
            % Get real data
            vol = real(obj.measurement);
        end
        
        function vol = getImaginary(obj)
            % Get imaginary data
            vol = imag(obj.measurement);
        end
        
        function vol = save(obj, dataType, fileName, saveFormat)
            % Get magnitude data
            % fileName: String. Prefix of filename (without extension)
            % saveFormat: 'nifti' (default) or 'mat'

            if ~exist('saveFormat', 'var')
                warning('No save format given - saving to NIfTI')
                saveFormat = 'nifti'; 
            end
            
                        
            if strcmp(fileName(end-3:end), '.nii')
                if ~strcmp(saveFormat, 'nifti')
                    warning('File extension and saveFormat do not match - saving to NIfTI format')
                    saveFormat = 'nifti';
                end
                fileName = fileName(1:end-4);
            elseif strcmp(fileName(end-3:end), '.mat')
                if ~strcmp(saveFormat, 'mat')
                    warning('File extension and saveFormat do not match - saving to MAT format')
                    saveFormat = 'mat';
                end
                fileName = fileName(1:end-4);            
            end
                        
            switch dataType
                case 'Magnitude'
                	vol = obj.getMagnitude();
                case 'Phase'
                	vol = obj.getPhase();
                case 'Real'
                	vol = obj.getReal();
                case 'Imaginary'
                	vol = obj.getImaginary();
                otherwise
                    error('Unknown datatype')
            end
            
            switch saveFormat
                case 'nifti'
                    nii_vol = make_nii(imrotate(fliplr(vol), -90));
                    save_nii(nii_vol, [fileName '.nii']);
                    
                    obj.writeJson(fileName)
                case 'mat'
                    save([fileName '.mat'], 'vol')
                    
                    obj.writeJson(fileName)
            end
        end
        
        function obj = writeJson(obj, fileName)
            pulseSeqProperties = struct("EchoTime", obj.TE, "FlipAngle", obj.FA);
            jsonFile = jsonencode(pulseSeqProperties);
            
            fid = fopen([fileName '.json'], 'w');
            if fid == -1, error('Cannot create JSON file'); end
            fwrite(fid, jsonFile, 'char');
            fclose(fid);
        end
 
        function obj = generate_deltaB0(obj, fieldType, params)       
            % fieldType: '2d_linearIP' (2D linear in-plane (IP))
            % params: 
            %    '2d_linearIP': [m, b] where y = mx + b, and the center of the
            %              volume is the origin. m is Hz per voxel. b is
            %              Hz. ang is the angle against the x axis of the
            %              volume
            
            switch fieldType
                case '2d_linearIP'
                    m = params(1);
                    b = params(2);
                    
                    dims = size(obj.starting_volume);

                    % Create coordinates
                    [X, Y] = meshgrid(linspace(-dims(1), dims(1), dims(1)), linspace(-dims(2), dims(2), dims(2)));

                    obj.deltaB0 = m*X+b;
                    
                case '3d_linearTP'
                    m = params(1);
                    b = params(2);
                    
                    dims = size(obj.starting_volume);
                    
                    % Create coordinates
                    [X,Y,Z] = ndgrid(linspace(-(dims(1)-1)/2,(dims(1)-1)/2,dims(1)),linspace(-(dims(2)-1)/2,(dims(2)-1)/2,dims(2)),linspace(-(dims(3)-1)/2,(dims(3)-1)/2,dims(3)));
                    
                    obj.deltaB0 = m*Z+b;
                    
                case 'load_external'
                    % calculate deltaB0 in Hz (external B0 field map should
                    % be in ppm)
                    obj.deltaB0 = imrotate( fliplr( niftiread(params) ) , -90);
                    obj.deltaB0 = (obj.gamma / (2*pi)) * obj.fieldStrength * obj.deltaB0;                    
                    
                otherwise
                    error('Undefined deltaB0 field type')
            end
            
            % Convert field from Hz to T;
            obj.deltaB0 = obj.deltaB0 / (obj.gamma / (2*pi));
        end
    end
    
    methods (Static)
        function signal = generate_signal(protonDensity, T2star, FA, TE, deltaB0, gamma, handedness)
            % FA = flip angle in degrees
            % T2star in seconds
            % TE in seconds
            % B0 in tesla
            % gamma in rad*Hz/Tesla
            
            switch handedness
                case 'left'
                    sign = -1;
                case 'right'
                    sign = 1;
            end
            
            signal = protonDensity.*sind(FA).*exp(-TE./T2star-sign*1i*gamma*deltaB0.*TE);
        end
        function noisyVolume = addNoise(volume, SNR)
            % volume: measurement volume of signals
            % SNR: Signal-to-noise ratio
            
            noiseSTD = max(volume(:))/SNR;
            noisyReal =      real(volume) + randn(size(volume)) * noiseSTD;
            noisyImaginary = imag(volume) + randn(size(volume)) * noiseSTD;

            noisyVolume = noisyReal + 1i*noisyImaginary;
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
            customVolume(abs(volume-0.1)<0.001) = (class2 + class1)/2;
            customVolume(abs(volume-0.4)<0.001) = class2*1.5;
        end
    end

end

