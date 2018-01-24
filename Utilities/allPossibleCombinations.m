params = {'SourceFilter', 'AcousticPath', 'Constraint', 'Grouping', 'ZerosFixed'};

sourceFilter = {'None', 'Loudspeaker'};
acPath = {'Current', 'Theoric'};
constraint = {'None', 'Magnitude'};
grouping = {'Independent', 'Alltogether'};
zerosFixed = {'No', 'Yes'};

paramValues = {sourceFilter, acPath, constraint, grouping, zerosFixed};
lengths = cellfun(@(p) numel(p) , paramValues);

numCombs = prod(lengths);
numCol = numel(lengths);

% Create all possible combinations matrix
% % Old way
% combs = zeros(numCombs, numCol);
% for k = 1:numCol
%     numCombsAfter = prod(lengths(k+1:end));
%     a = repmat(1:lengths(k), numCombsAfter, 1);
%     a = a(:);
%     numCombsBefore = prod(lengths(1:k-1));
%     col = repmat(a, numCombsBefore, 1);
%     combs(:, k) = col;
% end

% Better way
combs = cell(numCol, 1);
[combs{:}] = ind2sub(lengths, 1:numCombs);
combs = cell2mat(combs)';

combsString = cell(numCombs, numCol);
for k = 1:numCombs
    comb = combs(k, :);
    for c = 1:numCol
        combsString(k, c) = paramValues{c}(comb(c));
    end
end

combStringCell = mat2cell(combsString, numCombs, ones(numCol, 1));
T = table(combStringCell{:}, 'VariableNames', params)
