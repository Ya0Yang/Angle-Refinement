      % refineControl method (main workhorse for this REFINE class)
      function obj = refineControl_progressive(obj)
        
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
        
%         originalFullProjections = obj.refineFullProjections;
        projectionSeq = obj.projectionSeq;
        initialProjsNum = obj.initialProjsNum;
        batchRefinement = obj.batchRefinement; 
        
        obj = obj.pre_binning();  %so refineFullProjection is binned here, to speed up
        
        allProjs = obj.refineFullProjections; %temporary store the full dataset (binned)
        allAngles = obj.refineAngles;

        currProjs = allProjs(:,:,projectionSeq(1:initialProjsNum));
        currAngles = allAngles(projectionSeq(1:initialProjsNum),:);
        
        obj.refineFullProjections = currProjs;
        obj.refineAngles = currAngles;
        
        %refine the first 25 projections
        for refine_num = 1:obj.num_refinements
            fprintf('Initial 25 projections refinement iteration #%d\n',refine_num)
            
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
            obj = obj.refineOrientation_initial(); 

            % update projection and angles
            obj = obj.updateProjections_initial();    
        end
        allProjs(:,:,projectionSeq(1:initialProjsNum)) = obj.refineFullProjections; 
        allAngles(projectionSeq(1:initialProjsNum),:) = obj.refineAngles;
        
%         currProjsNum = initialProjsNum;
        for currProjsNum=initialProjsNum:batchRefinement:numel(projectionSeq)
            fprintf('Progressive refining Current Projections: # %d\n',currProjsNum);
            if currProjsNum==initialProjsNum
                
                obj.refineFullProjections = allProjs(:,:,projectionSeq(1:currProjsNum+batchRefinement));
                obj.refineAngles = allAngles(projectionSeq(1:currProjsNum+batchRefinement),:);
                
                obj = obj.get_model_from_RECONSTRUCTOR();    
                obj.RECONSTRUCTOR = obj.RECONSTRUCTOR.ClearCalcVariables();
                obj = obj.apply_binning();
                
                %refine the projs in batchRefinemetn with earlier model
                obj = obj.refineOrientationProgressive(currProjsNum,batchRefinement);
                obj = obj.updateProjectionsProgressive(currProjsNum,batchRefinement);
                
                allProjs(:,:,projectionSeq(1:currProjsNum+batchRefinement)) = obj.refineFullProjections;
                allAngles(projectionSeq(1:currProjsNum+batchRefinement),:) = obj.refineAngles;
                
            elseif currProjsNum+batchRefinement < numel(projectionSeq)
                %refine the projs in batchRefinemetn with new
                %reconstruction
                currProjs = allProjs(:,:,projectionSeq(1:currProjsNum));
                currAngles = allAngles(projectionSeq(1:currProjsNum),:);
                
                obj.refineFullProjections = currProjs;
                obj.refineAngles = currAngles;

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
                
                % refine and update
                obj.refineFullProjections = allProjs(:,:,projectionSeq(1:currProjsNum+batchRefinement));
                obj.refineAngles = allAngles(projectionSeq(1:currProjsNum+batchRefinement),:);
                
                obj = obj.apply_binning();
                obj = obj.refineOrientationProgressive(currProjsNum,batchRefinement);
                obj = obj.updateProjectionsProgressive(currProjsNum,batchRefinement);
                
                allProjs(:,:,projectionSeq(1:currProjsNum+batchRefinement)) = obj.refineFullProjections;
                allAngles(projectionSeq(1:currProjsNum+batchRefinement),:) = obj.refineAngles;
                
            elseif currProjsNum+batchRefinement >= numel(projectionSeq)
                currProjs = allProjs(:,:,projectionSeq(1:currProjsNum));
                currAngles = allAngles(projectionSeq(1:currProjsNum),:);
                obj.refineFullProjections = currProjs;
                obj.refineAngles = currAngles;

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
                
                obj.refineFullProjections = allProjs(:,:,projectionSeq(1:end));
                obj.refineAngles = allAngles(projectionSeq(1:end),:);
                obj = obj.apply_binning();
                obj = obj.refineOrientationProgressive(currProjsNum,numel(projectionSeq)-currProjsNum);
                obj = obj.updateProjectionsProgressive(currProjsNum,numel(projectionSeq)-currProjsNum);
                
%                 allProjs(:,:,projectionSeq(1:end)) = obj.refineFullProjections;
%                 allAngles(projectionSeq(1:end),:) = obj.refineAngles;
                
                fprintf('FINISHING ALL REFINEMENT PROGRESSIVELY.\n')
                break;
            end
            
            
        end
        
      end
            

