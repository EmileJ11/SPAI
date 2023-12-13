function ret = lineToLineSegmentIntersection(A, B, C, D)
%LINETOLINESEGMENTINTERSECTION check if line AB intersects line segment CD
% Created by: Aaron Geldert & Nils Meyer-Kahlen
% Last modified: 5 Oct 2022


coPlanarThreshold = 1e-3;
lengthErrorThreshold = 1e-3;

da = D-C;
db = B-A;
dc = A-C;

if abs(dot(dc, cross(da, db))) >= coPlanarThreshold
    ret = 0; return;
end

s = 0;
if norm(cross(da,db))^2 ~= 0
    s = dot(cross(dc,db), cross(da,db)) ./ norm(cross(da,db)).^2;
end

if s >= 0 && s <= 1
    intersection = C + s * da;
    
    if ((norm(intersection - A)^2 + norm(intersection - B)^2)<= (norm(A-B)^2 + lengthErrorThreshold))
        ret = 1; return;
    end
end

ret = 0;

end

