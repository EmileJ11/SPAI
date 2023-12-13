function reflectedSource = reflectSource(source, surface)
%REFECTSOURCE Summary of this function goes here
%   Checks for validity of reflection as well
% Created by: Aaron Geldert
% Last modified: 29 Oct 2022

% if source is not on reflective side, then invalid
if dot(surface.normal, (source.position-surface.center)) <= 0
    reflectedSource = NaN;
    return;
end

pos = (source.position-surface.center) - 2*dot(surface.normal, (source.position-surface.center)) * surface.normal + surface.center;

reflectedSource = Source(pos);

end

