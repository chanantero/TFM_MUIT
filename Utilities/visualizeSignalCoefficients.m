function f = visualizeSignalCoefficients(coefRef, coef, plotType)
% coefRef. P-element vector. Coefficients when there is only the effect of the noise source.
% coef. (P x N) matrix. Coefficients after cancellation

if nargin == 1
    plotType = 'linear';
end

[P, N] = size(coef);

powRef = (abs(coefRef).^2)/2;
pow = (abs(coef).^2)/2;
cancel = pow./repmat(powRef, 1, N);

coefMean = mean(pow, 1);
coefMax = max(pow, 1);
coefMin = min(pow, 1);
coefStd = std(pow, 1);

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
    strParamCoef{n} = ['$\bar{c} = ', num2str(coefMean(n), 2), '$, $\max{c} = ', num2str(coefMax(n), 2), '$, $\min{c} = ', num2str(coefMin(n), 2), '$, $\sigma_c = ', num2str(coefStd(n), 2), '$'];
    strParamCancel{n} = ['$\bar{c} = ', num2str(cancelMean(n), 2), '$, $\max{c} = ', num2str(cancelMax(n), 2), '$, $\min{c} = ', num2str(cancelMin(n), 2), '$, $\sigma_c = ', num2str(cancelStd(n), 2), '$'];
end

f = figure;
axCoefHist = axes(f);
axCoefHist.NextPlot = 'Add';

coef = mat2cell(coef, P, ones(1, N));

axHist = axes(f);
axHist.NextPlot = 'Add';

switch plotType
    case 'linear'
        
        % A) Histogram
        edges = 0:0.1:1; % Edges for the histogram: divide the 0 to 1 interval into equal intervals.
        
        strParam = cell(N, 1);
        
        for n = 1:N
            powRef = (abs(coefRef{n}).^2)/2;
            
            histogram(axHist, powRef, edges, 'DisplayStyle', 'stairs');
            
            meanVal = mean(powRef);
            maxVal = max(powRef);
            minVal = min(powRef);
            stdVal = std(powRef);
            strParam{n} = ['$\bar{c} = ', num2str(meanVal), '$, $\max{c} = ', num2str(maxVal), '$, $\min{c} = ', num2str(minVal), '$, $\sigma_c = ', num2str(stdVal), '$'];
        end
        
        legend(strParam, 'Interpreter', 'Latex');
        
    case 'logarithmic'
        % B) Logaritmic Histogram
        edges = -110:10:0; % Edges for the histogram
        
        strParam = cell(N, 1);
        
        powRef = (abs(coefRef).^2)/2;
        
        for n = 1:N           
            pow = (abs(coef{n}).^2)/2;
            cancel = pow./powRef;
            
            powLog = 10*log10(pow);
            cancelLog = 10*log10(cancel);
            
            powLog(powLog < edges(1)) = edges(1);
            
            histogram(axHist, powLog, edges, 'Normalization', 'Probability', 'DisplayStyle', 'stairs');
            
            meanVal = mean(powRef);
            maxVal = max(powRef);
            minVal = min(powRef);
            stdVal = std(powRef);
            cancelGlobal = sum(pow)/sum(powRef);
            strParam{n} = ['$\bar{c} = ', num2str(meanVal), '$, $\max{c} = ', num2str(maxVal), '$, $\min{c} = ', num2str(minVal), '$, $\sigma_c = ', num2str(stdVal), '$'];
        end
                
        axHist.XTick = edges;
        axHist.XTickLabels{1} = '-Inf';
        axHist.XLim = [edges(1), edges(end)];
        axHist.XLabel.String = 'Normalized Power (dB)';
        legend(strParam, 'Interpreter', 'Latex');
end

end