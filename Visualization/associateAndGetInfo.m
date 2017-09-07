function [ varargout ] = associateAndGetInfo( refValues, calcValues, selectedInfo )
% Input arguments:
% - refValues
% - calcValues

refValues = refValues(:);
calcValues = calcValues(:);

% Correspondence
[indRef, indCalc] = peakAssociation(refValues, calcValues);

varargout = cell(1, numel(selectedInfo));
[varargout{:}] = getAssociationInfo( refValues, calcValues, indRef, indCalc, selectedInfo );

end