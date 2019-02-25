      % refineControl method (main workhorse for this REFINE class)
      function obj = refineControl_XT(obj)
        
        obj=obj.CheckParameters();
       
        fprintf('Searching Phi angles between %d and %d with step size %d\n', ...
            obj.phi_search_range(1),obj.phi_search_range(end),obj.phi_search_range(2)-obj.phi_search_range(1))
        fprintf('Searching Theta angles between %d and %d with step size %d\n', ...
            obj.theta_search_range(1),obj.theta_search_range(end),obj.theta_search_range(2)-obj.theta_search_range(1))
        fprintf('Searching Psi angles between %d and %d with step size %d\n', ...
            obj.psi_search_range(1),obj.psi_search_range(end),obj.psi_search_range(2)-obj.psi_search_range(1))
        fprintf('Searching X translations between %d and %d with step size %d\n', ...
            obj.x_search_range(1),obj.x_search_range(end),obj.x_search_range(2) -obj.x_search_range(1)) 
        fprintf('Searching Y translations between %d and %d with step size %d\n', ...
            obj.y_search_range(1),obj.y_search_range(end),obj.y_search_range(2) - obj.y_search_range(1))
        total_number = numel(obj.phi_search_range) * numel(obj.theta_search_range)...
            * numel(obj.psi_search_range) * numel(obj.x_search_range) ...
            * numel(obj.y_search_range);
        fprintf('Total orientations per projection = %d\n', total_number)
        fprintf('Total number of projections = %d\n', size(obj.refineFullProjections,3))
        
        for refine_num = 1:obj.num_refinements
            fprintf('Refinement iteration #%d\n',refine_num)
            
            allProjs = obj.refineFullProjections;
            allAngles = obj.refineAngles;
            
            num_projs = size(allProjs,3);
            
            for proj_num = 1:num_projs
              
              fprintf('Refineing proj # %d / %d\n',proj_num,num_projs)
            
              if ~isempty(obj.ang_search_range_cell)

                  obj=set_parameters(obj,'ang_search_range',obj.ang_search_range_cell{refine_num});
                  fprintf('Setting angle search range from %f to %f\n',obj.ang_search_range(1),obj.ang_search_range(end));
              end
              
              currProjs = allProjs;
              refFullProj = currProjs(:,:,proj_num);
              refProj = bin(refFullProj,obj.bin_factor,2);
              currProjs(:,:,proj_num) = [];
              
              currAngles = allAngles;
              refAngle = currAngles(proj_num,:);
              currAngles(proj_num,:) = [];
              
              obj.refineFullProjections = currProjs;
              obj.refineAngles = currAngles;
                            

              % set initial (or updated) projection and angles in RECONSTRUCTOR for next
              % reconstruction
              obj = obj.set_projections_to_RECONSTRUCTOR();      
              obj = obj.set_angles_to_RECONSTRUCTOR();              

              % run RECONSTRUCTOR reconstruction. Note that the
              % RECONSTRUCTOR here can be any tomographic reconstruction
              % object, including EST, GENFIRE, SIRT, etc.
              obj.RECONSTRUCTOR = obj.RECONSTRUCTOR.CheckPrepareData();
              obj.RECONSTRUCTOR = obj.RECONSTRUCTOR.runGridding();
              obj.RECONSTRUCTOR = obj.RECONSTRUCTOR.reconstruct();


              % get reconstructed model from RECONSTRUCTOR
              obj = obj.get_model_from_RECONSTRUCTOR();            

              obj.RECONSTRUCTOR = obj.RECONSTRUCTOR.ClearCalcVariables();

              % apply binning for refinement
              obj = obj.apply_binning();

              % run orientation refinement
              obj = obj.refineOrientation_XT(refProj,refAngle); 
              % save temp for debugging
              %save(sprintf('tempIter%d.mat',refine_num),'obj');

              % update projection and angles 
              [obj,updated_projection,updated_angle] = obj.updateProjections_XT(refFullProj,refProj, proj_num);    
              %save(sprintf('results_%d',refine_num),'obj')
              
              allProjs(:,:,proj_num) = updated_projection;
              allAngles(proj_num,:) = updated_angle;

            end            
            
            tempSaveShift = obj.ShiftEvolution;
            save(sprintf('%s_refineControlTempSaveIter%d.mat',obj.saveprefix,refine_num),'tempSaveShift','allProjs','allAngles');
            
            obj.refineFullProjections = allProjs;
            obj.refineAngles = allAngles;

        end
      end