function [surfaces, furthestPt] = loadGeometry(filename)
%LOADGEOMETRY Loads a set of quadrilateral surfaces from file
%   Quadrilaterals are defined as 3D pts in lines of a .txt file
% Created by: Aaron Geldert & Nils Meyer-Kahlen
% Last modified: 15 Oct 2022

warning('Using deprecated standalone function - loadGeometry() in class IsmSimulator is preferred.')
iPoint = 0;
numSurfaces = 0;
pointsForSurface = zeros(4,3);
surfaces = [];
furthestPt = zeros(1,3);

lines = readlines(filename);
for iLine = 1:length(lines)
    % try to read coordinates, NaN if invalid
    coords = str2double(split(lines(iLine),' '));
    if any(isnan(coords))
        continue;
    else
        % update furthest point (yes, per x y z dimension!)
        furthestPt = max(coords.', furthestPt);
        
        iPoint = iPoint + 1;
        pointsForSurface(iPoint,:) = coords.';
        if iPoint == 4
            % add a surface
            newSurface = Surface(pointsForSurface);
            surfaces = [surfaces; newSurface];
            numSurfaces = numSurfaces + 1;
            iPoint = 0;
        end
    end
end
if isempty(surfaces)
    warning(['No geometry data could be loaded from ' filename]);
end
disp(['Loaded ' num2str(numSurfaces) ' surfaces from ' filename]);

end

