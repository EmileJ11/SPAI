function rankByMethod = determineRankings(ERR)

% counter for each method
% 1=POT, 2=NN, 3=ALN, 4=LIN
rankByMethod = zeros(4,4); % ranks are rows, methods are cols

validK = 2:20;

% iterate through the 100 trials
for trialInd = 1:100

    errs = [ERR.otc(validK,trialInd),...
        ERR.nnc(validK,trialInd),...
        ERR.alc(validK,trialInd),...
        ERR.lic(validK,trialInd)];

    [~, sortInd] = sort(errs, 2);
    count = histc(sortInd, 1:4); % count the number of times 
    rankByMethod = rankByMethod + count;
end

end