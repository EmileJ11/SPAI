function numDuplicates = countDuplicatesInPC(PC)
% Created by: Aaron Geldert
% Last modified: 16 Oct 2022
    uniquePos = unique(PC.pos,'rows');
    numDuplicates = PC.n - size(uniquePos,1);
end