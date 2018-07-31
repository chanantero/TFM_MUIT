function [ y ] = complex2polar( x )

y = [abs(x), rad2deg(angle(x))];

end

