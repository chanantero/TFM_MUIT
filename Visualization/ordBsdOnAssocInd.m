function [ orderedVector ] = ordBsdOnAssocInd( sizeRef, values, refInd, assocInd  )

orderedVector = NaN(sizeRef);
orderedVector(refInd) = values(assocInd);
    
end