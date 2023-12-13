function pcNN = nearestNeighborInterp(pc1, pc2, k)
% Created by: Aaron Geldert
% Last modified: 18 Oct 2022
    
    % sort based on mass
    [~, sortInds1] = sort(pc1.mass, 'descend');
    [~, sortInds2] = sort(pc2.mass, 'descend');
    
    pc1.mass = pc1.mass(sortInds1);
    pc1.pos = pc1.pos(sortInds1, :);
    pc2.mass = pc2.mass(sortInds2);
    pc2.pos = pc2.pos(sortInds2, :);
    
    % determine number of mappings
    nMap = min(pc1.n, pc2.n);
    
    % Make iterative deterministic mappings:
    % each iteration, largest unmapped of pc1 to nearest in pc2
    pcNN.n = 0;
    pcNN.mass = NaN(nMap, 1);
    pcNN.pos = NaN(nMap, 3);
    
    for ii = 1:nMap
        % find nearest mass
        [~, jj] = min(vecnorm(pc2.pos - pc1.pos(ii,:), 2, 2));
        
        % interpolate and add to results
        pcNN.mass(ii) = (1-k)*pc1.mass(ii) + k*pc2.mass(jj);
        pcNN.pos(ii,:) = (1-k)*pc1.pos(ii,:) + k*pc2.pos(jj,:);
        pcNN.n = pcNN.n + 1;
        
        % remove this mass from pc2!
        remainingInds = setdiff(1:pc2.n, jj);
        pc2.mass = pc2.mass(remainingInds);
        pc2.pos = pc2.pos(remainingInds,:);
        pc2.n = pc2.n-1;
    end
end

