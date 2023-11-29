function [distance] = euclidean_distance(position1,position2)
%UNTITLED Summary of this function goes here
%   Calculate the Euclidean distance between 2 points in space
%   Postion1 and Position2 should have same length
distance = 0;
dimensions = length(position1);

for i=1:dimensions
    distance = sqrt((position1(1) - position2(1))^2 + (position1(2) - position2(2))^2 + (position1(3) - position2(3))^2);
end

end

