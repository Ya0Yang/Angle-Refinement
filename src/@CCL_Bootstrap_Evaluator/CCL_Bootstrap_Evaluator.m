classdef CCL_Bootstrap_Evaluator
    % evaluate projection orientation by calculating cross common lines
    % between a set of projections that is sampled with replacement from
    % the full set to prevent overfitting
    properties
        REFINER % a Refiner object
        %%%%
        %% more stuff goes here
        %%%%
        oversampling_ratio     = 2;
    end
    methods
        function obj = Backprojection_Evaluator(refiner)
            obj.REFINER = refiner;
        end
        function [metric, suggested_center_x, suggested_center_y] = evaluateOrientation(obj,pj_num,phi,theta,psi) % compare projections
        %%%%
        %% more stuff goes here
        %%%%
        end
    end
end