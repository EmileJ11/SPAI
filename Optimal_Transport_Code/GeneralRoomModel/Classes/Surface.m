classdef Surface < handle
    %SURFACE A 3D quadrilateral planar surface
    %   The normal vector should point towards the interior of a geometry
    % Created by: Aaron Geldert & Nils Meyer-Kahlen
    % Last modified: 15 Oct 2022
    
    properties
        points % 4 points in 3D space [4x3]
        normal % normal vector [1x3]
        center % mean of points [1x3]
        area % surface area
        index
    end
    
    methods
        function obj = Surface(inPoints)
            %SURFACE Constructor
            obj.points = inPoints;
            A = inPoints(1,:);
            B = inPoints(2,:);
            C = inPoints(3,:);
            n = cross(B-A, C-A);
            obj.normal = n ./ norm(n);
            obj.center = mean(obj.points);
            
            % surface area is a sum of 2 triangles (convex assumption)
            obj.area = areaOfTriangle3D(inPoints(1,:), inPoints(2,:), inPoints(3,:))...
                + areaOfTriangle3D(inPoints(3,:), inPoints(4,:), inPoints(1,:));
            
            global globalSurfaceCount;
            obj.index = globalSurfaceCount;
            globalSurfaceCount = globalSurfaceCount+1;
            
        end
        
        function ret = inFront(obj, otherSurface)
            % check if any vertex of a given polygon is in front of the surface
            ret = 0;
            for ii = 1:length(otherSurface.points)
                compare = otherSurface.points(ii,:) - obj.center;
                if dot(compare, obj.normal) > 0
                    ret = 1; return;
                end
            end
        end
        
        function ret = doesLineIntersect(obj, pointA, pointB)
            % returns 0 if not, returns point otherwise
            % check if the line from pointA to pointB intersects plane
            l = pointB - pointA;
            if dot(l, obj.normal) == 0
                ret = 0; return;
            end
            d = dot((obj.center - pointA), obj.normal) ./ dot(l, obj.normal);
            if d == 0
                ret = 0; return;
            end
            % intersection point with plane
            p = pointA + l * d;
            % now intersect the line through intersection pt and center pt
            % with all edge, check if the number is even
            intersectionCount = 0;
            numPoints = 4;
            for iEdge = 1:numPoints
                A = obj.points(iEdge,:);
                B = obj.points(mod(iEdge,numPoints)+1,:);
                if lineToLineSegmentIntersection(A, B, obj.center, p) ~= 0 
                    intersectionCount = intersectionCount+1;
                end
            end
            
            if mod(intersectionCount, 2) == 0
                ret = p; % return intersection point
            else
                ret = 0; 
            end
        end
    end
end

