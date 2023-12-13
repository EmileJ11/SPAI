function area = areaOfTriangle3D(p1, p2, p3)
%AREAOFTRIANGLE3D Calculates the area of a triangle with 3 cartesian points
% Created by: Aaron Geldert
% Last modified: 15 Oct 2022

    assert(isequal(numel(p1), numel(p2), numel(p3), 3),'Coordinate vectors must consist of 3 elements');
    x = [p1(1), p2(1), p3(1)];
    y = [p1(2), p2(2), p3(2)];
    z = [p1(3), p2(3), p3(3)];
    oneVec = ones(size(x));
    area = 0.5*sqrt(det([x;y;oneVec])^2 + det([y;z;oneVec])^2 + det([z;x;oneVec])^2);
end