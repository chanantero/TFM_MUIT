function f = visualizeSignalCoefficients(recOnlyNoiseCoef, recCoef, WFScoef, NScoef)
% recOnlyNoiseCoef. P-element vector. Coefficients when there is only the effect of the noise source.
% coef. (P x N) matrix. Coefficients after cancellation

[P, N] = size(recCoef);

powOnlyNoise = (abs(recOnlyNoiseCoef).^2)/2;
pow = (abs(recCoef).^2)/2;
powWFS = (abs(WFScoef).^2)/2;
powNS = (abs(NScoef).^2)/2;
cancel = pow./repmat(powOnlyNoise, 1, N);
cancelGlobal = sum(pow, 1)./sum(powOnlyNoise);

% Terminar
powOnlyNoiseLog = 10*log10(powOnlyNoise);
powWFSLog = 10*log10(powWFS);
powNSLog = 10*log10(powNS);
powLog = 10*log10(pow);
cancelLog = 10*log10(cancel);

strRecOnlyNoise = getStatsInfoStr( powOnlyNoise );
strWFS = getStatsInfoStr( powWFS );
strRec = getStatsInfoStr( pow );
strCancel = getStatsInfoStr( cancel );
for n = 1:N
    strCancel{n} = [strCancel{n}, ', $C_{\mathit{global}} = ', num2str(cancelGlobal(n), 2), '$'];
end


f = figure;
axCoefHist = axes(f);
axCoefHist.NextPlot = 'Add';
axCoefHist.OuterPosition = [0 0 0.5 1];
axCancelHist = axes(f);
axCancelHist.NextPlot = 'Add';
axCancelHist.OuterPosition = [0.5 0 0.5 1];

% B) Logaritmic Histogram
edgesCoef = -110:10:40; % Edges for the coefficients histogram
edgesCancel = -110:10:20; % Edges for the cancellation histogram

powOnlyNoiseLog(powOnlyNoiseLog < edgesCancel(1)) = edgesCancel(1);
powWFSLog(powWFSLog < edgesCoef(1)) = edgesCoef(1);
powLog(powLog < edgesCoef(1)) = edgesCoef(1);
cancelLog(cancelLog < edgesCancel(1)) = edgesCancel(1);

powOnlyNoiseLog(powOnlyNoiseLog > edgesCancel(end)) = edgesCancel(end);
powWFSLog(powWFSLog > edgesCoef(end)) = edgesCoef(end);
powLog(powLog > edgesCoef(end)) = edgesCoef(end);
cancelLog(cancelLog > edgesCancel(end)) = edgesCancel(end);

map = colormap('lines');

histogram(axCoefHist, powOnlyNoiseLog, edgesCoef, 'Normalization', 'Probability', 'DisplayStyle', 'stairs', 'EdgeColor', 'black');
for n = 1:N
    histogram(axCoefHist, powWFSLog(:, n), edgesCoef, 'Normalization', 'Probability', 'DisplayStyle', 'stairs', 'LineStyle', '--', 'EdgeColor', map(n, :));
    histogram(axCoefHist, powLog(:, n), edgesCoef, 'Normalization', 'Probability', 'DisplayStyle', 'stairs', 'EdgeColor', map(n, :));
    histogram(axCancelHist, cancelLog(:, n), edgesCancel, 'Normalization', 'Probability', 'DisplayStyle', 'stairs', 'EdgeColor', map(n, :));
end

axCoefHist.XTick = edgesCoef;
axCoefHist.XTickLabels{1} = '-Inf';
axCoefHist.XTickLabels{end} = 'Inf';
axCoefHist.XLim = [edgesCoef(1), edgesCoef(end)];
axCoefHist.XLabel.String = 'Normalized Power (dB)';
legend(axCoefHist, [strRecOnlyNoise; strWFS; strRec], 'Interpreter', 'Latex');

axCancelHist.XTick = edgesCancel;
axCancelHist.XTickLabels{1} = '-Inf';
axCancelHist.XTickLabels{end} = 'Inf';
axCancelHist.XLim = [edgesCancel(1), edgesCancel(end)];
axCancelHist.XLabel.String = 'Cancellation (dB)';
legend(axCancelHist, strCancel, 'Interpreter', 'Latex');

% 2D map
axMapOnlyNoise = 
visualizeSignalCoefficients2Dmap(ax, recCoef, WFScoef, NScoef);

end


function strs = getStatsInfoStr( coefMat )

N = size(coefMat, 2);

coefMean = mean(coefMat, 1);
coefMax = max(coefMat, [], 1);
coefMin = min(coefMat, [], 1);
coefStd = std(coefMat, 0, 1);

strs = cell(N, 1);
for n = 1:N
    strs{n} = ['$\bar{p} = ', num2str(coefMean(n), 2), '$, $\max{p} = ', num2str(coefMax(n), 2), '$, $\min{p} = ', num2str(coefMin(n), 2), '$, $\sigma_p = ', num2str(coefStd(n), 2), '$'];
end

end