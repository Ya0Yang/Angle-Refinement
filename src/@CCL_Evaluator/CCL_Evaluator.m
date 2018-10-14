classdef CCL_Evaluator
    % evaluate projection orientation by calculating cross common lines
    % between each projection and all other projections
    properties
        REFINER % a Refiner object
        %%%%
        %% more stuff goes here
        %%%%
        oversampling_ratio       = 2;
        compare_func             = @mse; %% currently this is ignored and mse is inlined in evaluateOrientation
        pj_k                     = [];
        translation_phase_shifts = [];
        angle_start              = 0;
        angle_end                = 360;
        angle_step               = 0.1;
        lines                    = []; % contains the 1d Fourier lines
        line_angles              = []; % contains angles for each line determined by angle_start,angle_end, and angle_step
        line_angles_allpos       = [];
         %         maximize           = false;
    end
    methods
        function obj = CCL_Evaluator(refiner)
            obj.REFINER = refiner;
        end
        function obj = prepData(obj)
            
            % compute FFTs of projections and precompute 1d Fourier lines
            % and translational phase shifts
            
            OS                   = obj.REFINER.oversampling_ratio;
            [dimx,dimy,num_proj] = size(obj.REFINER.refineProjections);
            pad_dimx = dimx*OS;
            pad_dimy = dimy*OS;
            if (mod(pad_dimx,2)==0)
                pad_dimx = pad_dimx - 1;
            end
            if (mod(pad_dimy,2)==0)
                pad_dimy = pad_dimy - 1;
            end
            obj.pj_k             = zeros(pad_dimx,pad_dimy,num_proj);
            
            obj.line_angles = obj.angle_start : obj.angle_step  : obj.angle_end;
            num_lines       = length(obj.line_angles);
            tmp_lines       = zeros(pad_dimx,num_lines,num_proj); % layout so each 1d line is continuous in column major ordering
            
            for i = 1:num_proj
                fprintf('Precomputing projections for projection #%d/%d\n',i,num_proj)
%                 obj.pj_k(:,:,i) = my_fft(padarray(,[floor(( (OS-1) * dimx)/2), floor(( (OS-1) * dimy)/2)]));
                obj.pj_k(:,:,i) = my_fft(My_paddzero(obj.REFINER.refineProjections(:,:,i),[pad_dimx, pad_dimy]));
                if obj.REFINER.use_parallel
                    parfor line_num = 1:num_lines
                       tmp_lines(:,line_num,i) = calculate2Dprojection_interp(obj.pj_k(:,:,i),obj.line_angles(line_num))';
                    end
                else
                    for line_num = 1:num_lines
                       tmp_lines(:,line_num,i) = calculate2Dprojection_interp(obj.pj_k(:,:,i),obj.line_angles(line_num))';
                    end
                end
            end 
            obj.lines = tmp_lines;
           
                                    
            % calculate phase shifts for an x/y shift for each line
            n = size(tmp_lines,1);
            obj.translation_phase_shifts = zeros(n, 2, num_proj);
            
            for ang_num = 1:length(obj.line_angles)
                phi = obj.line_angles(ang_num);    
                obj.translation_phase_shifts(:, 1, ang_num) = translationalFourierPhaseShift(n, 1, 0, phi);     
                obj.translation_phase_shifts(:, 2, ang_num) = translationalFourierPhaseShift(n, 0, 1, phi);       
            end
            % take log now because the shifts are applied using exp
            obj.translation_phase_shifts = log(obj.translation_phase_shifts); 
            obj.line_angles_allpos = obj.line_angles;
            obj.line_angles(obj.line_angles > 180) = obj.line_angles(obj.line_angles > 180) - 360;
        end
        
        function idx = getNearestLine(obj, phi, pj_num)
           %% get Fourier line closest to phi 
           [~,idx] = min(abs(obj.line_angles-phi));
%            if phi < 0;
%                phi = phi + 180;
%            end
%            num_lines = numel(obj.line_angles_allpos);
%            spacing = 360/(num_lines-1);
%            idx1 = floor(phi / spacing) + 1;
%            idx2 = idx1 + 1;
%            d1 = phi - obj.line_angles_allpos(idx1);
%            d2 = obj.line_angles(idx2) - phi;
%            if d1 > d2
%                idx = idx2;
%            else
%                idx = idx1;
%            end
%            closest_line = obj.lines(:,idx,pj_num);
        end
        
        function [metric, suggested_center_x, suggested_center_y] = evaluateOrientation(obj,pj_num,phi,theta,psi) % compare projections
                                 
            [dimx,dimy,num_proj] = size(obj.pj_k);
%             ncx = round((dimx+1)/2);
%             ncy = round((dimy+1)/2);
            
            OR = obj.REFINER.oversampling_ratio;
            dimx_small = (dimx-1)/OR;
            dimy_small = (dimx-1)/OR;
            ncx_small = round((dimx_small+1)/2);
            ncy_small = round((dimy_small+1)/2);
            
            % get the rotation matrix for the orientation we are testing
            R_test = get_Euler_matrix([phi,theta,psi]); % rotation matrix of test orientation
            best_total_error = 1e30;
            
            % cache some variables for performance
            error_function = obj.compare_func; 
            lines = obj.lines;
            for xshift = obj.REFINER.x_search_range
                for yshift = obj.REFINER.y_search_range
                    % setup error metric
                    total_err = 0;
                    for compare_pj_num = 1:num_proj % compare with all projections (excluding self)
                        if compare_pj_num==pj_num
                            continue
                        end
                        R_ref = get_Euler_matrix(obj.REFINER.refineAngles(compare_pj_num,:)); % rotation matrix of test orientation
                        [v_ref, v_test] = accel_get_commonline_vectors(R_ref, R_test); % implied common lines in each projection
%                         [v_ref, v_test] = get_commonline_vectors(R_ref, R_test, [0;0;1]); % implied common lines in each projection
                        idx_ref      = obj.getNearestLine(atan2d(v_ref(2),v_ref(1)), compare_pj_num);
                        l_ref        = lines(:,idx_ref,compare_pj_num);
                        idx_test     = obj.getNearestLine(atan2d(v_test(2),v_test(1)),pj_num);
                        l_test       = lines(:,idx_test,pj_num);
                        
%                         l_test = l_test .* obj.translation_phase_shifts(:, 1, idx_test) .^ xshift ...
%                                         .* obj.translation_phase_shifts(:, 2, idx_test) .^ yshift;
                        l_test = l_test .* exp(xshift*(obj.translation_phase_shifts(:, 1, idx_test)) + ...
                                                yshift.*(obj.translation_phase_shifts(:, 2, idx_test)));

                        % accumulate error -- it's faster to inline this
                        % function but more flexible to use a function
                        % handle
                        total_err = total_err + error_function(l_ref,l_test);
%                         total_err = total_err + mean((abs(l_test-l_ref)).^2);
                    end
                    if (total_err < best_total_error)
                        best_total_error = total_err;
                        suggested_center_x = ncx_small + xshift;
                        suggested_center_y = ncy_small + yshift;
                    end
                end
            end
            metric = best_total_error ./ num_proj; % take average
        end
        
      function obj=set_parameters(obj,varargin)
        if mod(length(varargin),2) ~= 0
            error('Additional argument list not divisible by 2. Options should be ''key'',''value'' pairs.')
        end

        % Apply user-provided options
        par_number = 1;
        while par_number < length(varargin)
            if isprop(obj,varargin{par_number})
                  obj.(varargin{par_number}) = varargin{par_number+1};
                par_number = par_number + 2;
            else
                error('REFINEMENT: Invalid option %s provided.',varargin{par_number})
            end
        end
      end
        
    end
end