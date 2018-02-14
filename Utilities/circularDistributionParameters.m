function [meanAlpha, varianceAlpha, stdAlpha] =  circularDistributionParameters( alpha, dim )
% https://en.wikipedia.org/wiki/Directional_statistics
% https://en.wikipedia.org/wiki/Mean_of_circular_quantities
% http://www.fiserlab.org/manuals/procheck/manual/man_cv.html

if exist('dim', 'var') == 0
    dim = 1;
end

x = exp(1i*alpha);

meanX = mean(x, dim);
meanAlpha = angle(meanX);
R = abs(meanX);
varianceAlpha = 1 - R;
stdAlpha = sqrt(-2*log(R));

end