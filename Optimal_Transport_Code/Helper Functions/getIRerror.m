function E = getIRerror(PCtest, PCtrue, fs, pl)
%GETIRERROR calculates the error by applying a spreading function and
%sum of absolute difference, and more :)
% Created by: Aaron Geldert
% Last modified: 17 Oct 2022

if nargin<3 
    fs = 48e3; % default
end
if nargin < 4
    pl = 0;
end

t_er = 0.060; % 60 ms default early reflection assumption
intTime = 0.004; % 2 ms spreading on each side
convLen = round(intTime*fs)+1; % odd number for peak symmetry

% get MONO IRs from the PCs, always 200ms
hTest = monoIRfromPC(PCtest, 0.2, fs);
hTrue = monoIRfromPC(PCtrue, 0.2, fs);

% apply spreading function
% sf = spreadingFunction(convLen, -10);
sf = hann(convLen);
hsTest = conv(hTest, sf);
hsTrue = conv(hTrue, sf);
% yes, the convolution shifts things in time, but it matches

% truncate the IRs to the early reflection portion
n_ds = find(hsTrue,1); % find earliest non-zero sample using TRUE IR
n_er = round(n_ds + fs*t_er); % end of ER portion
hsTest = hsTest(1:n_er); % yes, consider from t=0!
hsTrue = hsTrue(1:n_er);

% sum of squared absolute difference, normalized by reference
E = sum((hsTest-hsTrue).^2) ./ sum(hsTrue.^2);

if pl ~=0
    figure();
    subplot(121);
    plot(1:n_er, hsTest); hold on; 
    plot(1:n_er, hsTrue);
    legend('hs1','hs2');
    subplot(122);
    plot(1:n_er, (hsTest-hsTrue).^2); title('Squared Difference');
    xlabel('Time (samples)');
    sgtitle(['Error: ' num2str(E)]);
end

end