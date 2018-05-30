function [transfFunct, transfFunctMinimized, fsel, axAbsInd, axPhaseInd, axAbsMin, axPhaseMin] = calculateTransferFunction(inputSignals, outputSignals, fs, varargin)
% Fuction for the analysis of correction factors based on time signals
% Input arguments:
% - outputSignals. Numeric N-dimensional array.
% - inputSignals. Numeric N-dimensional array. Same size as outputSignals
% - fs. Numeric scalar. Sampling frequency.
% - dimNames. Cellstring vector with N elements. The i-th element
% contains the name of the i-th dimension. The dimension corresponding to
% time should be called 'time'. The dimension where minimization of the
% quadratic error should be applied must be called 'minimizDim'.
% - fLim. Limit of the frequencies of interest.

numDims = ndims(outputSignals);
dimSizes = size(outputSignals);

p = inputParser;

addOptional(p, 'dimNames', [{'minimizDim'}, {'time'}, cell(1, numDims - 2)])
addParameter(p, 'freqEdges', [])
addParameter(p, 'fLim', [0 1000])

parse(p, varargin{:})

dimNames = p.Results.dimNames;
fLim = p.Results.fLim;

if ismember('freqEdges', p.UsingDefaults)
    predef_freqEdges = false;
else
    predef_freqEdges = true;
end

[~, indTimeDim] = ismember('time', dimNames);

numSamp = size(outputSignals, indTimeDim);
% t = (0:numSamp-1)/fs;
f = (0:numSamp-1)*fs/numSamp;
fIndSel = f >= fLim(1) & f <= fLim(2);
fsel = f(fIndSel);

% Perform FFT and calculate the correction factor for each frequency of
% interest

newDimOrder = [indTimeDim, 1:indTimeDim-1, indTimeDim+1:numDims]; % Move time frequency to be the 1

output_freq = fft(permute(outputSignals, newDimOrder)); %/fs;
output_freq = output_freq(fIndSel, :, :, :);

input_freq = fft(permute(inputSignals, newDimOrder)); %/fs;
input_freq = input_freq(fIndSel, :, :, :);

transfFunct = output_freq./input_freq;

% Histograms
if ~predef_freqEdges
    numBins = 100;
    binWidth = (fLim(2) - fLim(1))/numBins;
    freqEdges = 10:binWidth:1000;
else
    freqEdges = p.Results.freqEdges;
end
    
absCorrFactEdges = 0:0.1:4;
phaseCorrFactEdges = -10:90;

freqMat_ind = pointWiseExtend(fsel', transfFunct);

axAbsInd = histogram2D( abs(transfFunct), [], freqMat_ind, freqEdges, absCorrFactEdges );
axAbsInd.Title.String = '|H| for individual cases';
axAbsInd.XLabel.String = 'Frequency (Hz)';
axAbsInd.YLabel.String = '|H|';
colorbar

axPhaseInd = histogram2D( rad2deg(angle(transfFunct)), [], freqMat_ind, freqEdges, phaseCorrFactEdges );
axPhaseInd.Title.String = '\angle{H} for individual cases';
axPhaseInd.XLabel.String = 'Frequency (Hz)';
axPhaseInd.YLabel.String = '\angle{H} (º)';
colorbar

% Global cancellation
% Find the value of the transfer function that minimizes the quadratic
% error for all cases in the dimension minimizDim
[~, minimizDim] = ismember('minimizDim', dimNames);
minimizDim_adapted = find(newDimOrder == minimizDim);
newDimOrder_Minimiz = [minimizDim_adapted, 1, 2:minimizDim_adapted - 1, minimizDim_adapted + 1:numDims];
output_freq = permute(output_freq, newDimOrder_Minimiz);
input_freq = permute(input_freq, newDimOrder_Minimiz);
dimSizesRed = dimSizes; dimSizesRed(indTimeDim) = length(fsel);
globSize = dimSizesRed;
globSize(minimizDim) = 1;
transfFunctMinimized = zeros(globSize);
numCases = prod(globSize);

for k = 1:numCases
    outputVec = output_freq(:, k);
    inputVec = input_freq(:, k);
    transfFunctMinimized(k) = inputVec\outputVec;    
end

[~, aux] = ismember(1:numDims, newDimOrder);
freqMat_Min = pointWiseExtend(permute(fsel', aux), transfFunctMinimized);

axAbsMin = histogram2D( abs(transfFunctMinimized), [], freqMat_Min, freqEdges, absCorrFactEdges );
axAbsMin.Title.String = '|H| for minimized case';
axAbsMin.XLabel.String = 'Frequency (Hz)';
axAbsMin.YLabel.String = '|H|';
colorbar

axPhaseMin = histogram2D( rad2deg(angle(transfFunctMinimized)), [], freqMat_Min, freqEdges, phaseCorrFactEdges );
axPhaseMin.Title.String = '\angle{H} for minimized case';
axPhaseMin.XLabel.String = 'Frequency (Hz)';
axPhaseMin.YLabel.String = '\angle{H} (º)';
colorbar


end
