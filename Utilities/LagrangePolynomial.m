function polynom = LagrangePolynomial( x, y )
% Lagrange polynomial.
% Explanation at: http://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html

x = x(:);
y = y(:);

numDataPoints = length(x);
polynom = zeros(1, numDataPoints);
for i=1:numDataPoints
    p=1;
    for j=1:numDataPoints
        if j~=i
            c = poly(x(j))/(x(i)-x(j));
            p = conv(p,c);
        end
    end
    term = p*y(i);
    polynom = polynom + term;
end

end

