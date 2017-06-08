function greenValues = GreenFunction(x, k)
% Green's function
% dim. Dimension with size 3 that corresponds to the x, y, and z
% coordinates of the points in x

if nargin < 3
    dim = 2;
end

dist = sqrt(sum(x.^2, dim));
greenValues = exp(-1i*k*dist)./dist;

end