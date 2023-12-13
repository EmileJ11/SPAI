classdef SDMPC < handle
    %SDMPC Spatial Decomposition Method Point Cloud
    %   Processes SDM and prepares point cloud subsets
    % Created by: Aaron Geldert
    % Last modified: 16 Aug 2022
    
    properties
        % source data
        fs      % sample rate
        srir    % B-format audio (n x 4)
        nOrig   % original number of samples
        tVec    % time vector for all samples
        
        % preprocessed data
        pos     % cartesian coords
        doa     % spherical coords
        dst     % distance (via c*tVec)
        pMass  % energy
        dbMass  % decibels
        
        % selected data
        nSel    % selected number of samples
        selInds    % indices of selected points
        
        % point cloud
        PC
       
        % constants
        c = 343; % speed of sound
        
    end
    
    methods
        function obj = SDMPC(fs, arg1, arg2, distOffset)
            %SDMPC Construct an instance of this class
            % Either done from: 
            %   - fs, srir data (arg2 is a string), distOffset
            %   - fs, rir and doa data, distOffset
            % note: distOffset should be in meters (default 0)
            
            if nargin<4
                distOffset = 0;
            end
            
            assert(distOffset>=0, 'distOffset must be a non-negative value in meters');
            assert(isscalar(fs), 'First argument must be a scalar sample rate.');
            
            obj.fs = fs;
            obj.nOrig = length(arg1);
            obj.tVec = ((0:(obj.nOrig-1)).')./fs;
            obj.dst = obj.c.*obj.tVec + distOffset;
            
            if ischar(arg2)
                assert(size(arg1, 2)>=4, 'SRIR data must contain 4 channels of FOA data.');
                
                % Load from B-format data
                [B,A] = butter(4, 5200/(0.5*obj.fs)); % lpf at 5.2k for Eigenmike
                srirFilt = filtfilt(B,A,arg1);
                pos = srirFilt(:, 1) .* [srirFilt(:, 4), srirFilt(:, 2), srirFilt(:, 3)];
                pos = pos./vecnorm(pos,2,2);
                
                obj.pos = pos.*obj.dst;
                obj.doa = C2S(obj.pos);
                
                obj.srir = arg1;
                obj.pMass = abs(arg1(:,1));
                obj.dbMass = db(arg1(:,1));
                
            else
                if isrow(arg1)
                    arg1 = arg1.';
                end
                assert(size(arg2, 2)==3, 'DOA data must be Nx3 in spherical coordinates.')
                
                % Load from SDM data
                obj.srir = arg1; % just omni channel
%                 obj.doa = arg2;
%                 obj.pos = S2C(obj.doa);
                obj.pos = arg2;

                obj.pMass = abs(arg1); % bipolar values?
                obj.dbMass = db(abs(arg1));
            end
            
        end
        
        function obj = selectPoints(obj, threshold, num)
            % Selects up to num points that exceed db threshold
            [vals, ids] = sort(obj.dbMass, 'descend');
            firstInd = find(vals<threshold, 1);
            vecEnd = min(num, firstInd-1);
            ids = ids(1:vecEnd);
            ids = sort(ids);
            obj.selInds = ids;
            obj.nSel = numel(ids);
            obj.getPC;
        end
        
        function PC = getPC(obj)
            % Returns the point cloud of data
            PC.n = obj.nSel;
            PC.pos = obj.pos(obj.selInds, :);
            PC.mass = obj.pMass(obj.selInds); % mass is based on energy!
            obj.PC = PC;
        end
        
        function DF = getDF(obj)
            % Returns the diffuse field (opposite of point cloud)
            % stays full - no need for sparsity
%             dfInds = setdiff(1:obj.nOrig, obj.selInds);
            DF.n = obj.nOrig;
            DF.pos = obj.pos;
            DF.pos(obj.selInds,:) = 0;
            DF.mass = obj.pMass;
            DF.mass(obj.selInds) = 0;
        end
    end
    
    methods (Static)
        function [sdmpc1, sdmpc2] = selectBalancedPCs(sdmpc1, sdmpc2, nTarget)
            % This function takes in two SDMPC objects and updates their
            % point clouds to best balance their energies, aiming for about
            % nTarget selected points
            
            assert(isa(sdmpc1, 'SDMPC'), 'First argument must be an SDMPC object.');
            assert(isa(sdmpc2, 'SDMPC'), 'Second argument must be an SDMPC object.');
            assert(isscalar(nTarget), 'nTarget must be a scalar.');
            
            n1 = 1; n2 = 1; % number of points used
            
            % sort the masses
            [m1, i1] = sort(abs(sdmpc1.pMass), 'descend');
            [m2, i2] = sort(abs(sdmpc2.pMass), 'descend');
            
            while ((n1*n2 < nTarget*nTarget) &&...
                    (n1 < sdmpc1.nOrig) && (n2 < sdmpc2.nOrig))
                if sum(m1(1:n1)) < sum(m2(1:n2)) % dist 1 is smaller 
                    n1 = n1+1;
                else % dist 2 is smaller
                    n2 = n2+1;
                end
            end
            
            sdmpc1.nSel = n1;
            sdmpc1.selInds = i1(1:n1);
            sdmpc1.getPC;
            
            sdmpc2.nSel = n2;
            sdmpc2.selInds = i2(1:n2);
            sdmpc2.getPC;
        end
        
    end
end

