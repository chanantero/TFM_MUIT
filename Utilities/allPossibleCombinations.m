sourceFilter = {'None', 'Loudspeaker'};
acPath = {'Current', 'Theoric'};
constraint = {'None', 'Magnitude'};
grouping = {'Independent', 'Alltogether'};
zerosFixed = {'No', 'Yes'};

params = {sourceFilter, acPath, constraint, grouping, zerosFixed};
lengths = cellfun(@(p) numel(p) , params);

% Create all possible combinations matrix
numCombs = prod(lengths);
numCol = numel(lengths);
combs = zeros(numCombs, numCol);
for k = 1:numCol
    numCombsAfter = prod(lengths(k+1:end));
    a = repmat(1:lengths(k), numCombsAfter, 1);
    a = a(:);
    numCombsBefore = prod(lengths(1:k-1));
    col = repmat(a, numCombsBefore, 1);
    combs(:, k) = col;
end

combsString = cell(numCombs, numCol);
for k = 1:numCombs
    comb = combs(k, :);
    for c = 1:numCol
        combsString(k, c) = params{c}(comb(c));
    end
end
