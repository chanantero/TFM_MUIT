function array = stretchArray( array, factors )
% It takes an array and modifies it en the next way:
% size(arrayOutput) = size(array).*factors
% Example: stretchArray([1:3, 4:6], [2 2])
% One factor for each dimension

numDims = numel(factors);

for stretchingDim = 1:numDims
    sizeArray = ones(1, numDims);
    sizeArrayAux = size(array);
    sizeArray(1:numel(sizeArrayAux)) = sizeArrayAux;
    
    repeatVector = ones(1, numDims);
    repeatVector(stretchingDim) = factors(stretchingDim);
    repeatedArray = repmat(array, repeatVector);
    subindices = cell(1, numDims);
    for d = 1:numDims
        if(d ~= stretchingDim)
            shiftedVector = shiftdim((1:sizeArray(d))', -(d - 1));
            repeatVector2 = sizeArray;
            repeatVector2(d) = 1;
            repeatVector2(stretchingDim) = factors(stretchingDim)*sizeArray(stretchingDim);
            subindices{d} = repmat(shiftedVector, repeatVector2);
        else
            aux = repmat((1:sizeArray(stretchingDim):sizeArray(stretchingDim)*(factors(stretchingDim) - 1) + 1)', 1, sizeArray(d));
            aux = aux + repmat((0:size(aux, 2)-1), size(aux, 1), 1);
            aux = aux(:);
            shiftedVector = shiftdim(aux, -(d - 1));
            repeatVector2 = sizeArray;
            repeatVector2(d) = 1;
            subindices{d} = repmat(shiftedVector, repeatVector2);
        end
    end 
    linearIndices = sub2ind(size(repeatedArray), subindices{:});
    array = repeatedArray(linearIndices);
end

end

