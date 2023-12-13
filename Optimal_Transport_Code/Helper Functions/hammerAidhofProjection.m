function [xy] = hammerAidhofProjection(aziEle)
    % Hammer-Aidhof projection
    azi = aziEle(:, 1);
    ele = aziEle(:, 2);
    hap = @(azi, ele) [-cos(ele).*sin(azi/2) 0.5 * sin(ele)] ./ ...
        (sqrt(1+cos(ele).*cos(azi/2)));
    xy = hap (mod (azi+ pi, 2 * pi) - pi, ele);
end