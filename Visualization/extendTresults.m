function T = extendTresults( Tresults, realDoAname, calcDoAname, selectedInfo, newVariableNames )

tab = rowfun(@(realDoA, calcDoA) associateAndGetInfo( realDoA, calcDoA, selectedInfo ), Tresults, 'InputVariables', {realDoAname, calcDoAname}, 'NumOutputs', numel(selectedInfo), 'OutputFormat', 'cell', 'ExtractCellContents', true);

arrayVariables = {'numAssociations', 'numPhantomsRef', 'numPhantomsCalc'};
cellVariables = {'indRef', 'indCalc', 'existAssocRef', 'existAssocCalc', 'deviation'};

indArray = ismember(selectedInfo, arrayVariables);
indCell = ismember(selectedInfo, cellVariables);

% Suena anti-intuitivo, pero es así
tabArray = cell2table(tab(:, indArray), 'VariableNames', newVariableNames(indArray)); % Convert to array if possible
tabCell = array2table(tab(:, indCell), 'VariableNames', newVariableNames(indCell));

tab = [tabArray, tabCell];
[~, order] = sort([find(indArray) find(indCell)]);
tab = tab(:, order);

tab(:, indArray) = tabArray;
tab(:, indCell) = tabCell;

% Add columns to Tresults
T = [Tresults tab];

end