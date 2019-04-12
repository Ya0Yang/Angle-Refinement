classdef Refiner
    properties
        % properties with default values
        ang_step               = 1;
        ang_range              = 5;
        bin_factor             = 1;
        pre_bin_factor         = 3;
        compare_func           = @alignByRFactor;
        forwardProjection_func = @calculate3Dprojection_interp;
        oversampling_ratio     = 2;
        num_refinements        = 1;
        window_half_size       = 5;
        initialProjsNum = 25;
        batchRefinement = 3;
        paddingEdgeSize = 0;
        
        %in case of using real-space forward projector, RealProjection = 1
        RealProjection = 0;
        FullEvolutionRecord = 0; % if this flag is 1, all euler_angles and shift vectors during all iterations will be recorded
        % in AngleEvolution ShiftEvolution
        record_initial = 1;
        
        RefineReferenceAngleInd = 1; % Angle index for reference angle to be fixed
        RefineReferenceAngletoSet = [0 0 0]; % Reference angle will be fixed to this value
        RefineZeroCenterFlag = 0; % 0 for 0 to 180 deg theta convention
        % 1 for -90 to 90 deg theta convention
        
        % mask to be applied to forward projection before comparing with
        % given projections
        Rmethod = 0;
        Rscanres = 1;
        Rscanres_separateR = 0.1;
        RscanRad = 3;
        R_align_redoRec_flag = 0;
        FPmask
        FPind
        
        ang_search_range_cell
        
        % RECONSTRUCTOR is a MATLAB class object for tomographic
        % reconstruction. This can be GENFIRE, EST, SIRT, etc..
        % This will be initialized by its constructor
        RECONSTRUCTOR
        
        % properties which will be initialized in constructor
        ang_search_range
        phi_search_range
        theta_search_range
        psi_search_range
        translation_search_range = -5:5
        x_search_range
        y_search_range
        maximize
        use_parallel
        
        % properties which will be determined and used during the refinement process
        refineModel
        refineProjections  % projections used for refinement, possibly binned
        refineAngles
        refineFullProjections  % full projection, not binned
        measuredProjections
        centers_x
        x_centers
        centers_y
        y_centers
        modelK = [];
        noise_sigma
        bayes_probs
        metrics
        phis
        thetas
        psis
        errK_evolution
        errR_evolution
        ProjsEvolution
        AngleEvolution
        ShiftEvolution
        evaluator_type = 'backprojection';
        saveprefix='';
        
        % save temp reconstruction
        save_temp = 0;
        saveFilename = './results/temp_rec';
        parpool_size = 16
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor method
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % RECONSTRUCTOR should be set with constructor
        function obj = Refiner(val)
            if  nargin > 0
                if strcmp(val,'GENFIRE')
                    obj.RECONSTRUCTOR = GENFIRE_Reconstructor();
                else
                    error('REFINEMENT: unknown reconstructor name!');
                end
            else
                % default RECONSTRUCTOR is GENFIRE
                obj.RECONSTRUCTOR = GENFIRE_Reconstructor();
            end
            
            %set default variables
            obj.ang_search_range    = -obj.ang_range:obj.ang_step:obj.ang_range;%vector of angular displacements to search
            obj.phi_search_range    = obj.ang_search_range;%vector of angular displacements to search phi
            obj.theta_search_range  = obj.ang_search_range;%vector of angular displacements to search theta
            obj.psi_search_range    = obj.ang_search_range;%vector of angular displacements to search psi

            % calculated backprojection. Must be of the form [metric, new_center_x, new_center_y] = function(input_img,calc_img) where
            % metric is the value for R-factor, Xcorr, etc and new_center_x is the
            % optimal location of the center found for input_img based upon comparison
            % with calc_img
            obj.maximize            = false;%determines whether metric from compare_func should be maximized or minimized
            obj.use_parallel        = true;%use parallel (parfor) where applicable
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % declare long methods in external files
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        obj = updateProjections(obj);
        obj = refineOrientation(obj);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % declare short methods in this file
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % method for setting serach range from range and step input
        function obj = set_angle_range_step(obj,ang_range,angle_step)
            obj.ang_search_range    = -ang_range:angle_step:ang_range;%vector of angular displacements to search
            obj.phi_search_range    = obj.ang_search_range;%vector of angular displacements to search phi
            obj.theta_search_range  = obj.ang_search_range;%vector of angular displacements to search theta
            obj.psi_search_range    = obj.ang_search_range;%vector of angular displacements to search psi
        end
        
        function obj = set_angle_range_stepThetaOnly(obj)
            obj.phi_search_range    = 0;%vector of angular displacements to search phi
            obj.theta_search_range  = obj.ang_search_range;%vector of angular displacements to search theta
            obj.psi_search_range    = 0;%vector of angular displacements to search psi
        end
        
        % refineControl method (main workhorse for this REFINE class)
        function obj = refineControl(obj)
            
            obj=obj.CheckParameters();
            
            fprintf('Searching Phi angles between %d and %d with step size %.02f\n', ...
                obj.phi_search_range(1),obj.phi_search_range(end),obj.phi_search_range(2)-obj.phi_search_range(1))
            fprintf('Searching Theta angles between %d and %d with step size %.02f\n', ...
                obj.theta_search_range(1),obj.theta_search_range(end),obj.theta_search_range(2)-obj.theta_search_range(1))
            fprintf('Searching Psi angles between %d and %d with step size %.02f\n', ...
                obj.psi_search_range(1),obj.psi_search_range(end),obj.psi_search_range(2)-obj.psi_search_range(1))
            total_number = numel(obj.phi_search_range) * numel(obj.theta_search_range)...
                * numel(obj.psi_search_range);
            fprintf('Total orientations per projection = %d\n', total_number)
            fprintf('Total number of projections = %d\n', size(obj.refineFullProjections,3))
            
            for refine_num = 1:obj.num_refinements
                fprintf('Refinement iteration #%d\n',refine_num)
                
                if ~isempty(obj.ang_search_range_cell)
                    
                    obj=set_parameters(obj,'ang_search_range',obj.ang_search_range_cell{refine_num});
                    fprintf('Setting angle search range from %f to %f\n',obj.ang_search_range(1),obj.ang_search_range(end));
                end
                
                obj = obj.set_projections_to_RECONSTRUCTOR();
                obj = obj.set_angles_to_RECONSTRUCTOR();
                obj.RECONSTRUCTOR = obj.RECONSTRUCTOR.CheckPrepareData();
                obj.RECONSTRUCTOR = obj.RECONSTRUCTOR.runGridding();
                
                obj.RECONSTRUCTOR = obj.RECONSTRUCTOR.reconstruct();
                
                obj = obj.get_model_from_RECONSTRUCTOR();
                obj.RECONSTRUCTOR = obj.RECONSTRUCTOR.ClearCalcVariables();
                
                obj = obj.apply_binning();
                
                obj.errR_evolution(:,end+1) = obj.RECONSTRUCTOR.errR(end);
                obj.errK_evolution(:,end+1) = min(obj.RECONSTRUCTOR.errK);
                
                obj = obj.refineOrientation();
                % save temp for debugging
                %save(sprintf('tempIter%d.mat',refine_num),'obj');
                % update projection and angles
                obj = obj.updateProjections();
                %save(sprintf('results_%d',refine_num),'obj')
                
            end
        end
        % copy reconstructed model from RECONSTRUCTOR to refineModel property
        function obj = get_model_from_RECONSTRUCTOR(obj)
            obj.refineModel = obj.RECONSTRUCTOR.reconstruction;
        end
        
        % copy projections from RECONSTRUCTOR to refineFullProjections property
        function obj = get_projections_from_RECONSTRUCTOR(obj)
            obj.refineFullProjections = obj.RECONSTRUCTOR.InputProjections;
        end
        
        % copy angles from RECONSTRUCTOR to refineAngles property
        function obj = get_angles_from_RECONSTRUCTOR(obj)
            obj.refineAngles = obj.RECONSTRUCTOR.InputAngles;
        end
        
        % copy current refineFullProjections to RECONSTRUCTOR
        function obj = set_projections_to_RECONSTRUCTOR(obj)
            obj.RECONSTRUCTOR.InputProjections = obj.refineFullProjections;
        end
        
        % copy current refineAngles to RECONSTRUCTOR
        function obj = set_angles_to_RECONSTRUCTOR(obj)
            obj.RECONSTRUCTOR.InputAngles = obj.refineAngles;
        end
        
        % set parameters for REFINE class
        function obj=set_parameters(obj,varargin)
            if mod(length(varargin),2) ~= 0
                error('REFINEMENT: Additional argument list not divisible by 2. Options should be ''key'',''value'' pairs.')
            end
            
            % Apply user-provided options
            par_number = 1;
            while par_number < length(varargin)
                if isprop(obj,varargin{par_number})
                    if strcmp(varargin{par_number},'ang_search_range')
                        obj.(varargin{par_number}) = varargin{par_number+1};
                        obj.phi_search_range = varargin{par_number+1};
                        obj.theta_search_range = varargin{par_number+1};
                        obj.psi_search_range = varargin{par_number+1};
                    elseif strcmp(varargin{par_number},'translation_search_range')
                        obj.x_search_range = varargin{par_number+1};
                        obj.y_search_range = varargin{par_number+1};
                    else
                        obj.(varargin{par_number}) = varargin{par_number+1};
                    end
                    par_number = par_number + 2;
                else
                    error('REFINEMENT: Invalid option %s provided.',varargin{par_number})
                end
            end
        end
        
        function obj = CheckParameters(obj)
            % check size of FPmask
            if isempty(obj.FPmask)
                obj.FPmask = ones(size(obj.refineFullProjections));
                obj.FPmask = My_paddzero(My_stripzero(obj.FPmask,size(obj.refineFullProjections)-[obj.paddingEdgeSize*2 obj.paddingEdgeSize*2 0]),size(obj.refineFullProjections));
                
            else
                if ~isequal(size(obj.FPmask),size(obj.refineFullProjections))
                    error('REFINEMENT: size of FPmask and input projections do not match!')
                end
            end
            
            if obj.FullEvolutionRecord==1 && obj.record_initial == 1
                obj.AngleEvolution = obj.refineAngles;
                obj.ShiftEvolution = zeros(size(obj.refineFullProjections,3),2);
            end
            
            if ~isempty(obj.ang_search_range_cell)
                if length(obj.ang_search_range_cell)~=obj.num_refinements
                    error('REFINEMENT: number of elements in ang_search_range_cell should be the same with num_iterations!')
                end
            end
        end
        
        % apply binning to obj.refineModel and obj.refineProjections
        function obj = apply_binning(obj)
            obj.refineModel = bin(obj.refineModel,obj.bin_factor,3);
            obj.refineProjections = bin(obj.refineFullProjections,obj.bin_factor,2);
        end
        
        % applybinning to projections before reconstruction
        function obj = pre_binning(obj)
            obj.refineFullProjections = TianBins(obj.refineFullProjections,obj.pre_bin_factor,2);
        end
    end
end

