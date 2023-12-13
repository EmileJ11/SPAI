function PCnoisy = addNoiseToPC(PCclean, posNoise, magNoise)
% Created by: Aaron Geldert
% Last modified: 17 Oct 2022

PCnoisy = PCclean;
PCnoisy.pos = PCclean.pos + posNoise*randn(size(PCclean.pos));
PCnoisy.mass = PCclean.mass .* (1+ magNoise*randn(size(PCclean.mass)));

end