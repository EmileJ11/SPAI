classdef Source < handle
    %SOURCE Summary of this class goes here
    %   Detailed explanation goes here
    % Created by: Aaron Geldert & Nils Meyer-Kahlen
    % Last modified: 5 Oct 2022
    
    properties
        position
        index
        previousReflector = -1
        reflectionHistory
        order = 0
        parent
        reflector
        reflectors
        valid = 1
        path
    end
    
    methods
        function obj = Source(position)
            %SOURCE Constructor
            obj.position = position;
            obj.valid = 1; % assume true
%             global sourceIndex;
%             obj.index = sourceIndex;
%             sourceIndex = sourceIndex + 1;
        end
    end
end

