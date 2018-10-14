classdef Backprojection_Evaluator
    % evaluate projection orientation by calculating backprojection and
    % comparing to reference via compare_func
    properties
        REFINER % a Refiner object
        forwardProjection_func = @calculate3Dprojection_interp_wrapper; % function for calculating forward projection
        compare_func           = @alignByNormXCorr; % function for comparing two projections
        modelK                 = [];
        oversampling_ratio     = 2;
    end
    methods
        function obj = Backprojection_Evaluator(refiner)
            obj.REFINER = refiner;
        end
        function [metric, suggested_center_x, suggested_center_y] = evaluateOrientation(obj,pj_num,phi,theta,psi) % compare projections
            ref_pj  = obj.REFINER.refineProjections(:,:,pj_num); % todo: change refineProjections name to something like projections
            calc_pj = obj.forwardProjection_func(obj,phi,theta,psi);
            [metric, suggested_center_x, suggested_center_y] = obj.compare_func(ref_pj,calc_pj);
        end
        function obj = prepData(obj)
           if isempty(obj.modelK) % if empty, need to populate it. This way you only have to do this once
                obj.modelK = my_fft(My_paddzero(obj.REFINER.refineModel,round(size(obj.REFINER.refineModel)*obj.oversampling_ratio)));
            end
        end
    end
end