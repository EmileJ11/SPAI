function [DOA, M] = interpSDM(k, SDM1, SDM2, T, virtuals)
%INTERPSDM Summary of this function goes here
%   Interpolate at k, using transport plans T that may include virtual mass
%   
%   T: transport plan (masses placed places, zero otherwise)
% Created by: Aaron Geldert
% Last Modified: 24 May 2022

    assert(k>=0 && k<=1, 'k must be between 0 and 1.');

    if virtuals == 0
        [DOA, M] = transportPlanInterp3D(k, SDM1.pos, SDM2.pos, T);
    else
        %% There are 3 types of mass to account for:
        % 1. Transported Mass (TM)
        % - The interior of the transport matrix defines the masses of all
        % couplings of real points to other real points. Their corresponding
        % positions are simply linearly interpolated
        %
        % 2. Created Mass (CM)
        % - When a real point is coupled with a virtual point and one or 
        % more other real points, the mass of the coupling is linearly 
        % added to the other transported masses coupled with the real point
        %
        % 3. Spontaneous Mass (SM)
        % - When a real point is coupled only with a virtual point, the 
        % mass of the coupling spontaneously appears at position the real 
        % point posisiton, mass proportional to the interpolation parameter
        
        % 1. Transported Mass (TM)
        Ttm = T(1:end-1, 1:end-1); % ignore virtual mass columns
        
        % 2. Created Mass (CM)
        % From X->Y 
        xRat = Ttm./vecnorm(Ttm, 1, 2); % ratio of CM distribute, rowsum=1
        Tcm1 = xRat .* T(1:end-1, end); % calculate CM values
        Tcm1(isnan(Tcm1)) = 0; % remove NaN elements
        
        % From Y->X
        yRat = Ttm./vecnorm(Ttm, 1, 1); % colsum = 1
        Tcm2 = T(end, 1:end-1) .* yRat;
        Tcm2(isnan(Tcm2)) = 0;
        
        % Calculate transported & created masses at interp value k
        Ttmk = Ttm + (1-k)*Tcm1 + k*Tcm2;
        [I, J, TM] = find(Ttmk);
        numTm = numel(TM);
        
        % Calculate transported mass positions at interp value k
        DOAtm = zeros(numTm,3);
        for ii = 1:numTm
            DOAtm(ii,:) = (1-k).*SDM1.pos(I(ii),:) + k.*SDM2.pos(J(ii),:);
        end
               
        % 3. Spontaneous Mass (SM)
        % Which virtual masses are coupled only to Y?
        whereY = isnan(sum(yRat, 1));
        Sm1 = T(end,whereY).'; % ordered vector
        DOAsm1 = SDM2.pos(whereY,:);

        whereX = isnan(sum(xRat, 2)); % virtual masses coupled only to X
        Sm2 = T(whereX,end);
        DOAsm2 = SDM1.pos(whereX,:);
        
        % Concatenate all points
        DOA = [DOAtm; DOAsm1; DOAsm2];
        M = [TM; Sm1; Sm2];
        
    end
    
end
