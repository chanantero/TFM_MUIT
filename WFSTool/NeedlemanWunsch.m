function alignment = NeedlemanWunsch(insertionCost, deletionCost, substitutionCost)

lengthRefSequence = size(substitutionCost, 1);
lengthNewSequence = size(substitutionCost, 2);

acumulatedCost = zeros(lengthRefSequence + 1, lengthNewSequence + 1);
sequenceAction = categorical(zeros(size(acumulatedCost)));

% Initialize variables
acumulatedCost(2:end, 1) = cumsum(deletionCost(:, 1));
acumulatedCost(1, 2:end) = cumsum(insertionCost(1, :));
sequenceAction(:) = 'undefined';
sequenceAction(1, 2:end) = 'insert';
sequenceAction(2:end, 1) = 'delete';

for row = 2:lengthRefSequence + 1
    for column = 2:lengthNewSequence + 1
        
        insertAcumCost = acumulatedCost(row, column - 1) + insertionCost(row, column - 1);
        deletAcumCost = acumulatedCost(row - 1, column) + deletionCost(row - 1, column);
        subsAcumCost = acumulatedCost(row - 1, column - 1) + substitutionCost(row - 1, column - 1);
        [acumulatedCost(row, column), ind] = max([insertAcumCost, deletAcumCost, subsAcumCost]);
        
        if(ind == 1)
            sequenceAction(row, column) = 'insert';
        elseif(ind == 2)
            sequenceAction(row, column) = 'delete';
        elseif(ind == 3)
            sequenceAction(row, column) = 'substitute';
        end
        
    end
end

% Trace back
alignment = categorical(zeros(1, lengthRefSequence + lengthNewSequence));
step = 0;
row = lengthRefSequence + 1;
column = lengthNewSequence + 1;
while ~(column == 1 && row == 1)
    step = step + 1;
    switch sequenceAction(row, column)
        case 'insert'
            alignment(step) = 'insert';
            column = column - 1;
        case 'delete'
            alignment(step) = 'delete';
            row = row - 1;
        case 'substitute'
            alignment(step) = 'substitute';
            row = row - 1;
            column = column - 1;
        case 'undefined'
            if ~(column==1 && row ==1)
                error('NeedlemanWunsch:unknownError', 'Unknown error')
            end
    end
end

alignment = flip(alignment(1:step));

end