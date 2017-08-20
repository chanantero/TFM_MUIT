function [N1, N2] = multilateration(refPoints, distances)
% refPoints. numPoints x 3 matrix.
% distances. numPoints-element column vector

numPoints = numel(distances);

A = [ones(numPoints, 1), -2*refPoints];
b = distances.^2 - sum(refPoints.^2, 2);

if  (numPoints > 3)
    Xpdw=pinv(A)*b;
 
    N1 = Xpdw; %Xpdw(2:end);
    N2 = N1;
else
    error('The number of reference points must be 5 or bigger')
end

end