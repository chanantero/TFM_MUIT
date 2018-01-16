% Apply the inverse operation (third argument flag to true)
simulFieldFormatted = mergeAndPermute(data, modVec, true, sizeObjective);

% Visualize
visualObj = animation({xVec, yVec, zVec, amplitudeVec, phaseVec, 1:obj.numReceivers},...
    {simulFieldFormatted}, labels, {'Cancellation (dB)'}, [], []);

