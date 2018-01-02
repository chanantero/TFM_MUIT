coef = [(0:0.1:1)', randn(11, 1)];
coefOnlyNoise = ones(11,1);
f = visualizeSignalCoefficients(coefOnlyNoise, coef);