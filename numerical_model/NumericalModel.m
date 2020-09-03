classdef NumericalModel < handle
    %NUMERICALMODEL Numerical Model for B0 data generation.
    
    properties
        gamma = 267.52218744 * 10^6; % rad*Hz/Tesla
        fieldStrength = 3.0; % Tesla
        handedness = 'left'; % Siemens & Canon = 'left', GE & Philips = 'right' 
        
        % volume properties
        type
        numVox = 128; % square volume
        pixSize = 1; % default pixel size is 1 mm
        starting_volume
        volume
        measurement
        
        % pulse sequence properties
        TE
        FA
        
        % T2* values at 3T from
        % https://index.mirasmart.com/ISMRM2019/PDFfiles/4509.html
        % https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.26809
        % Silicone and Mineral oil values are set to 50% of T2
        % T2* CSF is set to 50% of the T2 of CSF (http://mri-q.com/why-is-t1--t2.html)
        % T2* water is set to 50% of the T2 of water (http://mri-q.com/why-is-t1--t2.html)
        % T2* of teeth is set to 50% of the T2 of teeth (https://www.karger.com/Article/Abstract/501901)
        % T2* of fat is set to 50% of the T2 of fat (http://mri-q.com/why-is-t1--t2.html)
        % T2* of muscle is set to 50% of the T2 of muscle (http://mri-q.com/why-is-t1--t2.html)
        % T2* of venous and arterial blood is from:
        % https://onlinelibrary.wiley.com/doi/pdf/10.1002/mrm.21342, table 1 using
        % (1-Ya) = 0.02 and (1-Yv) = 0.39 and a Hct level of 0.44
        % T2* of the "frontal portion eyes" is set to the T2 of the cornea (~ 50
        % ms) from https://onlinelibrary.wiley.com/doi/full/10.1002/jmri.21017
        % -> replaced with T2* of cornea from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4160095/
        % T2* of the "eyeball" is set to the T2 of the Chorioretina (~ 75 ms) from
        % https://onlinelibrary.wiley.com/doi/full/10.1002/jmri.21017
        % -> replaced with T2* of the sclera from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4160095/
        % T2* of the "eye" is set to the 50% of T2 of the vitrous humour (~ 750 ms) from https://onlinelibrary.wiley.com/doi/full/10.1002/jmri.21017
        % T2* of the "lens" is set to the 50% of T2 of the lens nucleus (~ 1130 ms) from https://onlinelibrary.wiley.com/doi/full/10.1002/jmri.21017
        % T2* of the "globus pallidus" is from https://pubs.rsna.org/doi/full/10.1148/radiol.2522081399
        % T2* of cartilage is from https://openmedicinejournal.com/VOLUME/5/PAGE/119/FULLTEXT/
        % T2* of the skin is set to 50% of T2 of the skin at 1.5T from https://www.sciencedirect.com/science/article/pii/S0022202X9190213A
        % T2* of the esophagus is from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5447643/
        % T2* of the parotid gland is set to 50% of T2 of the parotid gland from https://pubmed.ncbi.nlm.nih.gov/23166360/
        % T2* of the pituitary gland is set to 50% of T2 of the pituitary gland from https://www.ajronline.org/doi/full/10.2214/ajr.175.6.1751567
        T2star = struct('sinus',0,'teeth',0.5*150e-3,'bone',3.3e3,'thalamus',58e-3, ...
            'caudate_nucleus',58e-3,'putamen',40e-3,'globus_pallidus',27e-3,'GM',66e-3,'WM',0.053,'CSF',0.5*2,'air',0,'water',0.5*2, ...
            'fat',0.5*70e-3,'muscle',0.5*50e-3,'venous_blood',55e-3,'arterial_blood',20e-3,'frontal_portion_eyes',25e-3,'eyeball',10e-3,...
            'eye',0.5*750e-3,'lens',0.5*1130e-3,'cartilage',25e-3,'skin',0.5*30e-3,'esophagus',17e-3,'parotid_gland',0.5*120e-3,...
            'pituitary_gland',0.5*90e-3,'silicone_oil',0.5*0.566,'pure_mineral_oil',0.5*0.063);
        
        % Default proton density in percentage
        % blood pd from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3023928/
        % pd of fat set to 80, which is a guess
        % pd of muscle set to 65, which is a guess
        % pd of cartilage set to 40, which is a guess
        protonDensity = struct('teeth',20,'bone',20,'GM',82,'WM',70,'CSF',100,...
            'air',0,'water',100,'fat',80,'muscle',65,'blood',83,'cartilage',40,'silicone_oil', 71, 'pure_mineral_oil', 100);
        
        % deltaB0 in Hz
        deltaB0
    end
   
    methods
    	function obj = NumericalModel(model, varargin)
            % NumericalModel class
            % model: string with the label of the model desired ('Cylindrical3d', 'Spherical3d', 'Shepp-Logan3d', 'Shepp-Logan2d', 'Zubal')
            % varargin:  
            %   if 'Spherical3d': number of voxels, iamge_res in [mm], radius in [mm],
            %   material inside the sphere ('air', 'silicone_oil' or 
            %   'pure_mineral_oil'), material outside the sphere ('air,
            %   'silicone_oil, or 'pure_mineral_oil')
            %    ex usage: spherical_vol = NumericalModel('Spherical3d',128,1,5,'Air', 'silicone_oil');
            %
            %   if 'Cylindrical3d': number of voxels, iamge_res in [mm], radius in [mm],
            %   theta in [degrees] (angle between principal axis of
            %   cylinder and B0), material inside the sphere ('air', 
            %   'silicone_oil' or 'pure_mineral_oil'), material outside the 
            %   sphere ('air, 'silicone_oil, or 'pure_mineral_oil')
            %   ex usage: cylindrical_vol = NumericalModel('Cylindrical3d',128,1,5,90,'air', 'silicone_oil');
            %
            %   if 'Zubal': file name containing EAO's modified Zubal phantom
            %   ex usage: zubal_vol = NumericalModel('Zubal',fname_zubalEAO);
            %   
            
            obj.type = model ;
            
            if exist('model', 'var')
                switch model
                    case 'Zubal'
                        obj.numVox = [256 256 128];
                        obj.pixSize = [1.1 1.1 1.4];
                    	obj.Zubal(varargin);
                        
                    case 'Shepp-Logan2d'
                        obj.numVox = varargin{1,1}{1,1};
                    	obj.shepp_logan_2d(obj.numVox);
                        
                    case 'Shepp-Logan3d'
                        obj.numVox = varargin{1,1}{1,1};
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
        
        function obj = Zubal(obj, varargin)
            % load zubal volume
            obj.starting_volume = double(niftiread(varargin{1,1}{1,1}));
            
            % associate voxel values in the zubal phantom with their
            % approproate tissues
            zubal_id = struct('outside_phantom',0,'skin',1,'cerebral_fluid',2,'spinal_cord',3,'skull',4,'spine',5,'dens_of_axis',70, ...
                'jaw_bone',71,'parotid_gland',72,'skeletal_muscle',9,'lacrimal_glands',74,'spinal_canal',75,'hard_palate',76, ...
                'cerebellum',77,'tongue',78,'pharynx',15,'esophagus',16,'horn_of_mandible',81,'nasal_septum',82,'white_matter',83, ...
                'superior_sagittal_sinus',84,'medulla_oblongata',85,'blood_pool',23,'frontal_lobes',89,'bone_marrow',26, ...
                'pons',91,'third_ventricle',92,'trachea',29,'cartilage',30,'occipital_lobes',95,'hippocampus',96,'pituitary_gland',97, ...
                'fat1',98,'fat2',22,'ear_bones',99,'turbinates',100,'caudate_nucleus',101,'zygoma',102,'insula_cortex',103,'sinuses_mouth_cavity',104, ...
                'putamen',105,'optic_nerve',106,'internal_capsule',107,'septum_pellucidium',108,'thalamus',109,'eyeball',110,'corpus_collosum',111, ...
                'special_region_frontal_lobes',112,'cerebral_falx',113,'temporal_lobes',114,'fourth_ventricle',115,'frontal_portion_eyes',116, ...
                'parietal_lobes',117,'amygdala',118,'eye',119,'globus_pallidus',120,'lens',121,'cerebral_aquaduct',122,'lateral_ventricles',123, ...
                'prefrontal_lobes',124,'teeth',125,'sigmoid_sinus',126);
            
            % associate each tissue in the zubal phantom with a proton
            % density
            zubal_pd = struct('outside_phantom',obj.protonDensity.air,'skin',obj.protonDensity.water,'cerebral_fluid',obj.protonDensity.CSF,'spinal_cord',(obj.protonDensity.GM+obj.protonDensity.WM)/2,...
                'skull',obj.protonDensity.bone,'spine',obj.protonDensity.bone,'dens_of_axis',obj.protonDensity.bone,'jaw_bone',obj.protonDensity.bone,'parotid_gland',obj.protonDensity.water,...
                'skeletal_muscle',obj.protonDensity.muscle,'lacrimal_glands',obj.protonDensity.water,'spinal_canal',obj.protonDensity.CSF,'hard_palate',obj.protonDensity.bone,...
                'cerebellum',obj.protonDensity.GM,'tongue',obj.protonDensity.muscle,'pharynx',obj.protonDensity.air,'esophagus',obj.protonDensity.water,'horn_of_mandible',obj.protonDensity.bone,...
                'nasal_septum',obj.protonDensity.bone,'white_matter',obj.protonDensity.WM,'superior_sagittal_sinus',obj.protonDensity.blood,...
                'medulla_oblongata',(obj.protonDensity.GM+obj.protonDensity.WM)/2,'fat1',obj.protonDensity.fat,'fat2',obj.protonDensity.fat,'blood_pool',obj.protonDensity.blood,...
                'frontal_lobes',obj.protonDensity.GM,'bone_marrow',obj.protonDensity.CSF,'pons',obj.protonDensity.WM,'third_ventricle',obj.protonDensity.CSF,'trachea',obj.protonDensity.air,...
                'cartilage',obj.protonDensity.cartilage,'occipital_lobes',obj.protonDensity.GM,'hippocampus',obj.protonDensity.GM,'pituitary_gland',obj.protonDensity.water,'ear_bones',obj.protonDensity.bone,...
                'turbinates',(obj.protonDensity.bone+obj.protonDensity.water)/2,'caudate_nucleus',obj.protonDensity.GM,'zygoma',obj.protonDensity.bone,...
                'insula_cortex',obj.protonDensity.GM,'sinuses_mouth_cavity',obj.protonDensity.air,'putamen',obj.protonDensity.GM,'optic_nerve',obj.protonDensity.WM,...
                'internal_capsule',obj.protonDensity.WM,'septum_pellucidium',(obj.protonDensity.WM+obj.protonDensity.GM)/2,'thalamus',obj.protonDensity.GM,'eyeball',obj.protonDensity.water,...
                'corpus_collosum',obj.protonDensity.WM,'special_region_frontal_lobes',obj.protonDensity.GM,'cerebral_falx',obj.protonDensity.GM,'temporal_lobes',obj.protonDensity.GM,...
                'fourth_ventricle',obj.protonDensity.CSF,'frontal_portion_eyes',obj.protonDensity.water,'parietal_lobes',obj.protonDensity.GM,'amygdala',obj.protonDensity.GM,...
                'eye',obj.protonDensity.water,'globus_pallidus',obj.protonDensity.GM,'lens',obj.protonDensity.water,'cerebral_aquaduct',obj.protonDensity.CSF,...
                'lateral_ventricles',obj.protonDensity.CSF,'prefrontal_lobes',obj.protonDensity.GM,'teeth',obj.protonDensity.teeth,'sigmoid_sinus',obj.protonDensity.blood);
            
            % associate each tissue in the zubal phantom with a T2star
            % value
            % notes:
            % turbinates are made of bone and soft tissue, so the average of bone and water is used
            % T2* of cerebral_falx is set to that of GM
            % T2* of the lacrimal glands is set to that of the parotid gland
            zubal_t2s = struct('outside_phantom',obj.T2star.air,'skin',obj.T2star.skin,'cerebral_fluid',obj.T2star.CSF,'spinal_cord',(obj.T2star.GM+obj.T2star.WM)/2,...
                'skull',obj.T2star.bone,'spine',obj.T2star.bone,'dens_of_axis',obj.T2star.bone,'jaw_bone',obj.T2star.bone,'parotid_gland',obj.T2star.parotid_gland,...
                'skeletal_muscle',obj.T2star.muscle,'lacrimal_glands',obj.T2star.parotid_gland,'spinal_canal',obj.T2star.CSF,'hard_palate',obj.T2star.bone,...
                'cerebellum',obj.T2star.GM,'tongue',obj.T2star.muscle,'pharynx',obj.T2star.sinus,'esophagus',obj.T2star.esophagus,'horn_of_mandible',obj.T2star.bone,...
                'nasal_septum',obj.T2star.bone,'white_matter',obj.T2star.WM,'superior_sagittal_sinus',obj.T2star.venous_blood,...
                'medulla_oblongata',(obj.T2star.GM+obj.T2star.WM)/2,'fat1',obj.T2star.fat,'fat2',obj.T2star.fat,'blood_pool',(obj.T2star.venous_blood+obj.T2star.arterial_blood)/2,...
                'frontal_lobes',obj.T2star.GM,'bone_marrow',obj.T2star.CSF,'pons',obj.T2star.WM,'third_ventricle',obj.T2star.CSF,'trachea',obj.T2star.sinus,...
                'cartilage',obj.T2star.cartilage,'occipital_lobes',obj.T2star.GM,'hippocampus',obj.T2star.GM,'pituitary_gland',obj.T2star.pituitary_gland,'ear_bones',obj.T2star.bone,...
                'turbinates',(obj.T2star.bone+obj.T2star.water)/2,'caudate_nucleus',obj.T2star.caudate_nucleus,'zygoma',obj.T2star.bone,...
                'insula_cortex',obj.T2star.GM,'sinuses_mouth_cavity',obj.T2star.sinus,'putamen',obj.T2star.putamen,'optic_nerve',obj.T2star.WM,...
                'internal_capsule',obj.T2star.WM,'septum_pellucidium',(obj.T2star.WM+obj.T2star.GM)/2,'thalamus',obj.T2star.thalamus,'eyeball',obj.T2star.frontal_portion_eyes,...
                'corpus_collosum',obj.T2star.WM,'special_region_frontal_lobes',obj.T2star.GM,'cerebral_falx',obj.T2star.GM,'temporal_lobes',obj.T2star.GM,...
                'fourth_ventricle',obj.T2star.CSF,'frontal_portion_eyes',obj.T2star.frontal_portion_eyes,'parietal_lobes',obj.T2star.GM,'amygdala',obj.T2star.GM,...
                'eye',obj.T2star.eye,'globus_pallidus',obj.T2star.globus_pallidus,'lens',obj.T2star.lens,'cerebral_aquaduct',obj.T2star.CSF,...
                'lateral_ventricles',obj.T2star.CSF,'prefrontal_lobes',obj.T2star.GM,'teeth',obj.T2star.teeth,'sigmoid_sinus',obj.T2star.venous_blood);
            
            % create a volume 
            obj.volume = struct('magn', [], 'phase', [], 'T2star', [], 'protonDensity', []);
            
            % set the proton density equal to a starting volume which has
            % zero at all voxels
            obj.volume.protonDensity = obj.starting_volume;
            
            % assign the corresponding proton density to each tissue type
            % in the volume
            obj.volume.protonDensity(obj.starting_volume==zubal_id.outside_phantom) = zubal_pd.outside_phantom;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.skin) = zubal_pd.skin;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.cerebral_fluid) = zubal_pd.cerebral_fluid;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.spinal_cord) = zubal_pd.spinal_cord;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.skull) = zubal_pd.skull;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.spine) = zubal_pd.spine;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.dens_of_axis) = zubal_pd.dens_of_axis;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.jaw_bone) = zubal_pd.jaw_bone;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.parotid_gland) = zubal_pd.parotid_gland;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.skeletal_muscle) = zubal_pd.skeletal_muscle;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.lacrimal_glands) = zubal_pd.lacrimal_glands;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.spinal_canal) = zubal_pd.spinal_canal;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.hard_palate) = zubal_pd.hard_palate;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.cerebellum) = zubal_pd.cerebellum;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.tongue) = zubal_pd.tongue;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.pharynx) = zubal_pd.pharynx;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.esophagus) = zubal_pd.esophagus;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.horn_of_mandible) = zubal_pd.horn_of_mandible;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.nasal_septum) = zubal_pd.nasal_septum;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.white_matter) = zubal_pd.white_matter;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.superior_sagittal_sinus) = zubal_pd.superior_sagittal_sinus;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.medulla_oblongata) = zubal_pd.medulla_oblongata;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.blood_pool) = zubal_pd.blood_pool;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.frontal_lobes) = zubal_pd.frontal_lobes;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.bone_marrow) = zubal_pd.bone_marrow;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.pons) = zubal_pd.pons;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.third_ventricle) = zubal_pd.third_ventricle;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.trachea) = zubal_pd.trachea;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.cartilage) = zubal_pd.cartilage;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.occipital_lobes) = zubal_pd.occipital_lobes;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.hippocampus) = zubal_pd.hippocampus;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.pituitary_gland) = zubal_pd.pituitary_gland;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.fat1) = zubal_pd.fat1;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.fat2) = zubal_pd.fat2;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.ear_bones) = zubal_pd.ear_bones;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.turbinates) = zubal_pd.turbinates;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.caudate_nucleus) = zubal_pd.caudate_nucleus;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.zygoma) = zubal_pd.zygoma;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.insula_cortex) = zubal_pd.insula_cortex;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.sinuses_mouth_cavity) = zubal_pd.sinuses_mouth_cavity;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.putamen) = zubal_pd.putamen;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.optic_nerve) = zubal_pd.optic_nerve;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.internal_capsule) = zubal_pd.internal_capsule;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.septum_pellucidium) = zubal_pd.septum_pellucidium;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.thalamus) = zubal_pd.thalamus;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.eyeball) = zubal_pd.eyeball;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.corpus_collosum) = zubal_pd.corpus_collosum;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.special_region_frontal_lobes) = zubal_pd.special_region_frontal_lobes;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.cerebral_falx) = zubal_pd.cerebral_falx;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.temporal_lobes) = zubal_pd.temporal_lobes;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.fourth_ventricle) = zubal_pd.fourth_ventricle;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.frontal_portion_eyes) = zubal_pd.frontal_portion_eyes;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.parietal_lobes) = zubal_pd.parietal_lobes;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.amygdala) = zubal_pd.amygdala;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.eye) = zubal_pd.eye;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.globus_pallidus) = zubal_pd.globus_pallidus;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.lens) = zubal_pd.lens;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.cerebral_aquaduct) = zubal_pd.cerebral_aquaduct;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.lateral_ventricles) = zubal_pd.lateral_ventricles;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.prefrontal_lobes) = zubal_pd.prefrontal_lobes;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.teeth) = zubal_pd.teeth;
            obj.volume.protonDensity(obj.starting_volume==zubal_id.sigmoid_sinus) = zubal_pd.sigmoid_sinus;
            
            % set the T2star equal to a starting volume which has
            % zero at all voxels
            obj.volume.T2star = obj.starting_volume;
            
            % assign the corresponding T2star to each tissue type
            % in the volume
            obj.volume.T2star(obj.starting_volume==zubal_id.outside_phantom) = zubal_t2s.outside_phantom;
            obj.volume.T2star(obj.starting_volume==zubal_id.skin) = zubal_t2s.skin;
            obj.volume.T2star(obj.starting_volume==zubal_id.cerebral_fluid) = zubal_t2s.cerebral_fluid;
            obj.volume.T2star(obj.starting_volume==zubal_id.spinal_cord) = zubal_t2s.spinal_cord;
            obj.volume.T2star(obj.starting_volume==zubal_id.skull) = zubal_t2s.skull;
            obj.volume.T2star(obj.starting_volume==zubal_id.spine) = zubal_t2s.spine;
            obj.volume.T2star(obj.starting_volume==zubal_id.dens_of_axis) = zubal_t2s.dens_of_axis;
            obj.volume.T2star(obj.starting_volume==zubal_id.jaw_bone) = zubal_t2s.jaw_bone;
            obj.volume.T2star(obj.starting_volume==zubal_id.parotid_gland) = zubal_t2s.parotid_gland;
            obj.volume.T2star(obj.starting_volume==zubal_id.skeletal_muscle) = zubal_t2s.skeletal_muscle;
            obj.volume.T2star(obj.starting_volume==zubal_id.lacrimal_glands) = zubal_t2s.lacrimal_glands;
            obj.volume.T2star(obj.starting_volume==zubal_id.spinal_canal) = zubal_t2s.spinal_canal;
            obj.volume.T2star(obj.starting_volume==zubal_id.hard_palate) = zubal_t2s.hard_palate;
            obj.volume.T2star(obj.starting_volume==zubal_id.cerebellum) = zubal_t2s.cerebellum;
            obj.volume.T2star(obj.starting_volume==zubal_id.tongue) = zubal_t2s.tongue;
            obj.volume.T2star(obj.starting_volume==zubal_id.pharynx) = zubal_t2s.pharynx;
            obj.volume.T2star(obj.starting_volume==zubal_id.esophagus) = zubal_t2s.esophagus;
            obj.volume.T2star(obj.starting_volume==zubal_id.horn_of_mandible) = zubal_t2s.horn_of_mandible;
            obj.volume.T2star(obj.starting_volume==zubal_id.nasal_septum) = zubal_t2s.nasal_septum;
            obj.volume.T2star(obj.starting_volume==zubal_id.white_matter) = zubal_t2s.white_matter;
            obj.volume.T2star(obj.starting_volume==zubal_id.superior_sagittal_sinus) = zubal_t2s.superior_sagittal_sinus;
            obj.volume.T2star(obj.starting_volume==zubal_id.medulla_oblongata) = zubal_t2s.medulla_oblongata;
            obj.volume.T2star(obj.starting_volume==zubal_id.blood_pool) = zubal_t2s.blood_pool;
            obj.volume.T2star(obj.starting_volume==zubal_id.frontal_lobes) = zubal_t2s.frontal_lobes;
            obj.volume.T2star(obj.starting_volume==zubal_id.bone_marrow) = zubal_t2s.bone_marrow;
            obj.volume.T2star(obj.starting_volume==zubal_id.pons) = zubal_t2s.pons;
            obj.volume.T2star(obj.starting_volume==zubal_id.third_ventricle) = zubal_t2s.third_ventricle;
            obj.volume.T2star(obj.starting_volume==zubal_id.trachea) = zubal_t2s.trachea;
            obj.volume.T2star(obj.starting_volume==zubal_id.cartilage) = zubal_t2s.cartilage;
            obj.volume.T2star(obj.starting_volume==zubal_id.occipital_lobes) = zubal_t2s.occipital_lobes;
            obj.volume.T2star(obj.starting_volume==zubal_id.hippocampus) = zubal_t2s.hippocampus;
            obj.volume.T2star(obj.starting_volume==zubal_id.pituitary_gland) = zubal_t2s.pituitary_gland;
            obj.volume.T2star(obj.starting_volume==zubal_id.fat1) = zubal_t2s.fat1;
            obj.volume.T2star(obj.starting_volume==zubal_id.fat2) = zubal_t2s.fat2;
            obj.volume.T2star(obj.starting_volume==zubal_id.ear_bones) = zubal_t2s.ear_bones;
            obj.volume.T2star(obj.starting_volume==zubal_id.turbinates) = zubal_t2s.turbinates;
            obj.volume.T2star(obj.starting_volume==zubal_id.caudate_nucleus) = zubal_t2s.caudate_nucleus;
            obj.volume.T2star(obj.starting_volume==zubal_id.zygoma) = zubal_t2s.zygoma;
            obj.volume.T2star(obj.starting_volume==zubal_id.insula_cortex) = zubal_t2s.insula_cortex;
            obj.volume.T2star(obj.starting_volume==zubal_id.sinuses_mouth_cavity) = zubal_t2s.sinuses_mouth_cavity;
            obj.volume.T2star(obj.starting_volume==zubal_id.putamen) = zubal_t2s.putamen;
            obj.volume.T2star(obj.starting_volume==zubal_id.optic_nerve) = zubal_t2s.optic_nerve;
            obj.volume.T2star(obj.starting_volume==zubal_id.internal_capsule) = zubal_t2s.internal_capsule;
            obj.volume.T2star(obj.starting_volume==zubal_id.septum_pellucidium) = zubal_t2s.septum_pellucidium;
            obj.volume.T2star(obj.starting_volume==zubal_id.thalamus) = zubal_t2s.thalamus;
            obj.volume.T2star(obj.starting_volume==zubal_id.eyeball) = zubal_t2s.eyeball;
            obj.volume.T2star(obj.starting_volume==zubal_id.corpus_collosum) = zubal_t2s.corpus_collosum;
            obj.volume.T2star(obj.starting_volume==zubal_id.special_region_frontal_lobes) = zubal_t2s.special_region_frontal_lobes;
            obj.volume.T2star(obj.starting_volume==zubal_id.cerebral_falx) = zubal_t2s.cerebral_falx;
            obj.volume.T2star(obj.starting_volume==zubal_id.temporal_lobes) = zubal_t2s.temporal_lobes;
            obj.volume.T2star(obj.starting_volume==zubal_id.fourth_ventricle) = zubal_t2s.fourth_ventricle;
            obj.volume.T2star(obj.starting_volume==zubal_id.frontal_portion_eyes) = zubal_t2s.frontal_portion_eyes;
            obj.volume.T2star(obj.starting_volume==zubal_id.parietal_lobes) = zubal_t2s.parietal_lobes;
            obj.volume.T2star(obj.starting_volume==zubal_id.amygdala) = zubal_t2s.amygdala;
            obj.volume.T2star(obj.starting_volume==zubal_id.eye) = zubal_t2s.eye;
            obj.volume.T2star(obj.starting_volume==zubal_id.globus_pallidus) = zubal_t2s.globus_pallidus;
            obj.volume.T2star(obj.starting_volume==zubal_id.lens) = zubal_t2s.lens;
            obj.volume.T2star(obj.starting_volume==zubal_id.cerebral_aquaduct) = zubal_t2s.cerebral_aquaduct;
            obj.volume.T2star(obj.starting_volume==zubal_id.lateral_ventricles) = zubal_t2s.lateral_ventricles;
            obj.volume.T2star(obj.starting_volume==zubal_id.prefrontal_lobes) = zubal_t2s.prefrontal_lobes;
            obj.volume.T2star(obj.starting_volume==zubal_id.teeth) = zubal_t2s.teeth;
            obj.volume.T2star(obj.starting_volume==zubal_id.sigmoid_sinus) = zubal_t2s.sigmoid_sinus;

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
            image_res = [varargin{1,1}{1,2} varargin{1,1}{1,2} varargin{1,1}{1,2}];
            R = varargin{1,1}{1,2}; % [mm]
            
            % define image grid
            [x,y,z] = ndgrid(linspace(-(matrix(1)-1)/2,(matrix(1)-1)/2,matrix(1)),linspace(-(matrix(2)-1)/2,(matrix(2)-1)/2,matrix(2)),linspace(-(matrix(3)-1)/2,(matrix(3)-1)/2,matrix(3)));
            
            % define radial position (in [mm])
            r = sqrt((x.*image_res(1)).^2 + (y.*image_res(2)).^2 + (z.*image_res(3)).^2);
            
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
                 
        end
        
        function obj = cylindrical_3d(obj, varargin)
            % Create a 3D cylindrical volume for the
            % dims: [x, y, z] number of voxels.
            
            % define local variables
            matrix = [obj.numVox obj.numVox obj.numVox];
            image_res = [varargin{1,1}{1,2} varargin{1,1}{1,2} varargin{1,1}{1,2}];
            R = varargin{1,1}{1,3}; % [mm]
            theta = varargin{1,1}{1,4}; % angle between main axis of cylinder and z-axis (in degrees)
            
            % define image grid
            [x,y,z] = ndgrid(linspace(-(matrix(1)-1)/2,(matrix(1)-1)/2,matrix(1)),linspace(-(matrix(2)-1)/2,(matrix(2)-1)/2,matrix(2)),linspace(-(matrix(3)-1)/2,(matrix(3)-1)/2,matrix(3)));
            
            % define radial position (in [mm])
            r = sqrt((x.*image_res(1)).^2 + (y.*image_res(2)).^2);
            
            obj.volume = struct('magn', [], 'phase', [], 'T2star', [], 'protonDensity', []);
            
            obj.volume.protonDensity = obj.starting_volume;
            obj.volume.T2star = obj.starting_volume;
            
            switch varargin{1,1}{1,5}
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
               
            switch varargin{1,1}{1,6}
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
                    if strcmp(obj.type,'Zubal')
                        nii_vol = make_nii(vol);
                    else
                        nii_vol = make_nii(imrotate(fliplr(vol), -90));
                    end
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
                    if strcmp(obj.type,'Zubal')
                        obj.deltaB0 = niftiread(params);
                    else
                        obj.deltaB0 = imrotate( fliplr( niftiread(params) ) , -90);                   
                    end
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

