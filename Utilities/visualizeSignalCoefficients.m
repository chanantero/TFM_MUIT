function f = visualizeSignalCoefficients(coefRef, coef)
% coefRef. P-element vector. Coefficients when there is only the effect of the noise source.
% coef. (P x N) matrix. Coefficients after cancellation

[P, N] = size(coef);

powRef = (abs(coefRef).^2)/2;
pow = (abs(coef).^2)/2;
cancel = pow./repmat(powRef, 1, N);

coefMean = mean(pow, 1);
coefMax = max(pow, [], 1);
coefMin = min(pow, [], 1);
coefStd = std(pow, [], 1);

cancelMean = mean(cancel, 1);
cancelMax = max(cancel, 1);
cancelMin = min(cancel, 1);
cancelStd = std(cancel, 1);

cancelGlobal = sum(pow, 1)./sum(powRef);

powLog = 10*log10(pow);
cancelLog = 10*log10(cancel);

strParamCoef = cell(1, N);
strParamCancel = cell(1, N);
for n = 1:N
    strParamCoef{n} = ['$\bar{p} = ', num2str(coefMean(n), 2), '$, $\max{p} = ', num2str(coefMax(n), 2), '$, $\min{p} = ', num2str(coefMin(n), 2), '$, $\sigma_p = ', num2str(coefStd(n), 2), '$'];
    strParamCancel{n} = ['$\bar{C} = ', num2str(cancelMean(n), 2), '$, $\max{C} = ', num2str(cancelMax(n), 2), '$, $\min{C} = ', num2str(cancelMin(n), 2), '$, $\sigma_C = ', num2str(cancelStd(n), 2), '$, $C_{\mathit{global}} = ', num2str(cancelGlobal(n), 2), '$'];
end


f = figure;
axCoefHist = axes(f);
axCoefHist.NextPlot = 'Add';
axCoefHist.OuterPosition = [0 0 0.5 1];
axCancelHist = axes(f);
axCancelHist.NextPlot = 'Add';
axCancelHist.OuterPosition = [0.5 0 0.5 1];

% B) Logaritmic Histogram
edgesCoef = -110:10:0; % Edges for the coefficients histogram
edgesCancel = -110:10:20; % Edges for the cancellation histogram

powLog(powLog < edgesCoef(1)) = edgesCoef(1);
powLog(powLog > edgesCoef(end)) = edgesCoef(end);
cancelLog(cancelLog < edgesCancel(1)) = edgesCancel(1);
cancelLog(cancelLog > edgesCancel(end)) = edgesCancel(end);

for n = 1:N
    histogram(axCoefHist, powLog(:, n), edgesCoef, 'Normalization', 'Probability', 'DisplayStyle', 'stairs');
    histogram(axCancelHist, cancelLog(:, n), edgesCancel, 'Normalization', 'Probability', 'DisplayStyle', 'stairs');
end

axCoefHist.XTick = edgesCoef;
axCoefHist.XTickLabels{1} = '-Inf';
axCoefHist.XLim = [edgesCoef(1), edgesCoef(end)];
axCoefHist.XLabel.String = 'Normalized Power (dB)';
legend(axCoefHist, strParamCoef, 'Interpreter', 'Latex');

axCancelHist.XTick = edgesCancel;
axCancelHist.XTickLabels{1} = '-Inf';
axCancelHist.XLim = [edgesCancel(1), edgesCancel(end)];
axCancelHist.XLabel.String = 'Cancellation (dB)';
legend(axCancelHist, strParamCancel, 'Interpreter', 'Latex');

end