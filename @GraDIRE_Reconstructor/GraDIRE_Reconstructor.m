classdef GraDIRE_Reconstructor
    
    properties
        % internal variables to be cleared after run
        InputProjections
        InputAngles
        
        % internal variables to be kept after run
        reconstruction
        errR
        Dim1  % projection size1
        Dim2  % projection size2
        n1_oversampled  % projection size1 after oversampling
        n2_oversampled  % projection size2 after oversampling
        NumProjs
        Rarr_record
        Rarr2_record
        
        % R factor monitoring
        monitor_R = 0
        monitorR_loopLength = 10
        
        % filenames
        filename_Projections = ''
        filename_Angles = '';
        filename_Results = './results/GraDIRE_rec.mat';
        
        % reconstruction parameters
        numIterations = 50;
        oversamplingRatio = 3;
        griddingMethod = 1;
        
        % axes vectors for phi, theta, psi
        vector1 = [0 0 1];
        vector2 = [0 1 0];
        vector3 = [1 0 0];
        
        % save temp reconstruction
        save_temp = 0;
        save_loopLength = 50;
        saveFilename = './results/temp_rec';
        
        step_size = 2.;         % normalized step_size is chosen in range (1,3)
        sigma  = 000;           % smoothing constraint, smoother with smaller sigma, 0 if no smooth
        dtype  = 'single';      % data type, single for memory efficiency
        kernel;
        sum_rot_pjs;
        use_parallel = 0
        xj;
        yj;
        zj;
        Rot_x;
        Rot_y;
        Support;
        initial_model = 0;
        shape_constrain;
        shape_limit = 0;
    end
    
    methods
        
        % declare long methods in external files
        obj = reconstruct(obj)
        
        % declare short methods in this file
        function obj = readFiles(obj)
            if FileExist(obj.filename_Projections)
                obj.InputProjections = importdata(obj.filename_Projections);
            else
                error('GraDIRE: Projections file does not exist!')
            end
            if FileExist(obj.filename_Angles)
                obj.InputAngles = importdata(obj.filename_Angles);
            else
                error('GraDIRE: Angles file does not exist!')
            end
        end
        
        function obj = CheckPrepareData(obj)
            
            % set number of projections
            obj.NumProjs = size(obj.InputProjections,3);
            obj.Dim1 = size(obj.InputProjections,1);
            obj.Dim2 = size(obj.InputProjections,2);

            % input angle and projection size check
            if obj.NumProjs~=size(obj.InputAngles,1)
                error('GraDIRE: Number of projections and Angles does not match!')
            end
            
            % input angle check
            if size(obj.InputAngles,2) >3
                error('GraDIRE: Input Angle 2nd dimenstions larger than three!')
            end
            
            % if only one angle set is given, make them three Euler angles
            if size(obj.InputAngles,2) == 1
                obj.InputAngles = [zeros(obj.NumProjs ,1); obj.InputAngles; zeros(obj.NumProjs ,1)];
            end
            
            % check if interp method is legitimate
            if obj.griddingMethod ~= 1
                error('GraDIRE: Unrecognized gridding method.')
            end
            
            obj.n1_oversampled = round(obj.Dim1 * obj.oversamplingRatio);
            obj.n2_oversampled = round(obj.Dim2 * obj.oversamplingRatio);
            
        end
        function obj = runGridding(obj)
            fprintf('GraDIRE: Interpolate real space projections...\n\n');
            switch obj.griddingMethod
                case 1
                    obj = obj.interp_pj_realspace();
            end
        end

        % clear big temporary arrays
        function obj = ClearCalcVariables(obj)
            obj.InputProjections = [];
            obj.InputAngles      = [];
            obj.sum_rot_pjs      = [];
            obj.xj               = [];
            obj.yj               = [];
            obj.zj               = [];
            obj.Rot_x            = [];
            obj.Rot_y            = [];
        end
        function SaveResults(obj)
            save(obj.filename_Results, 'obj')
        end
        % set parameters for GraDIRE class
        function obj=set_parameters(obj,varargin)
            if mod(length(varargin),2) ~= 0
                error('GraDIRE: Additional argument list not divisible by 2. Options should be ''key'',''value'' pairs.')
            end
            % Apply user-provided options
            par_number = 1;
            while par_number < length(varargin)
                if isprop(obj,varargin{par_number})
                    obj.(varargin{par_number}) = varargin{par_number+1};
                else
                    error('GraDIRE: Invalid option %s provided.',varargin{par_number})
                end
                par_number = par_number + 2;
            end
        end
    end
end