function [ insertionCost, deletionCost, substitutionCost ] = genCostMatrices( refSeq, newSeq )
% New way
lengthRefSequence = numel(refSeq);
lengthNewSequence = numel(newSeq);

if(lengthRefSequence > lengthNewSequence)
    % Insertions are prohibited
    insertionCost = -Inf(lengthRefSequence + 1, lengthNewSequence);
    deletionCost = zeros(lengthRefSequence, lengthNewSequence + 1);
elseif(lengthRefSequence < lengthNewSequence)
    % Deletions are prohibited
    insertionCost = zeros(lengthRefSequence + 1, lengthNewSequence);
    deletionCost = -Inf(lengthRefSequence, lengthNewSequence + 1);
else
    % Deletions and insertions are prohibited
    insertionCost = -Inf(lengthRefSequence + 1, lengthNewSequence);
    deletionCost = -Inf(lengthRefSequence, lengthNewSequence + 1);
end

[refMesh, newMesh] = ndgrid(refSeq, newSeq);
substitutionCost = -(refMesh - newMesh).^2;

end

