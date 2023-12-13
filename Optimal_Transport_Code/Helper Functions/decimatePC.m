function PCdec = decimatePC(PCfull, nTarget, nRange)
% Created by: Aaron Geldert
% Last modified: 17 Oct 2022

    n = round(nTarget + randn * nRange);
    n = max(n, 1);
    n = min(PCfull.n, n);
    
    % truncate based on mass
    [~, sortInds] = sort(PCfull.mass, 'descend');
    sortInds = sortInds(1:n);
    
    PCdec.mass = PCfull.mass(sortInds, :);
    PCdec.pos = PCfull.pos(sortInds, :);
    PCdec.n = n;
end