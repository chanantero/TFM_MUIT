function varargout = getAssociationInfo( refValues, calcValues, indRef, indCalc, selectedOutputs )

calcValues = calcValues(:);
refValues = refValues(:);
[indRef, order] = sort(indRef);
indCalc = indCalc(order);
numRefVal = numel(refValues);
numCalcVal = numel(calcValues);
numAssociations = numel(indRef); % numel(indCalc) is the same

numOutputs = numel(selectedOutputs);
results = cell(1, numOutputs);
for k = 1:numOutputs
    switch selectedOutputs{k}
        case 'indRef'
            results{k} = indRef;
        case 'indCalc'
            results{k} = indCalc;
        case 'orderedCalc'
            orderedCalc = NaN(size(refValues));
            orderedCalc(indRef) = calcValues(indCalc);
            results{k} = orderedCalc;
        case 'numAssociations'
            % Number of associated reference values
            results{k} = numAssociations;
        case 'numPhantomsRef'
            % Number of not associated reference values
            numNotAssocRef = numRefVal - numAssociations;
            results{k} = numNotAssocRef;
        case 'numPhantomsCalc'
            % Number of not associated calculated values (phantom values)
            numNotAssocCalc = numCalcVal - numAssociations;
            results{k} = numNotAssocCalc;
        case 'existAssocRef'
            % Existence of a correspondence for each reference value
%             existAssocRef = ismember(1:numRefVal, indRef);
existAssocRef = false(size(refValues));
existAssocRef(indRef) = true;
            results{k} = existAssocRef;
        case 'existAssocCalc'
            % Existence of a correspondence for each calculated value
%             existAssocCalc = ismember(1:numCalcVal, indCalc);
existAssocCalc = false(size(calcValues));
existAssocCalc(indCalc) = true;
            results{k} = existAssocCalc;
        case 'deviation'
            % Deviation of associated values. For the non associated, the result is NaN
            deviation = NaN(1, numRefVal);
            deviation(indRef) = calcValues(indCalc) - refValues(indRef);
            results{k} = deviation;
        otherwise
            error('getAssociationInfo:notRecognizedOutput', 'The output selection is not recognized')
    end
end

varargout = results;

end