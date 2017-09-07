function [ dataVector, newParameters ] = cellContents2Vector( cellMat, parameters )

numElem = cellfun(@(thres) numel(thres), cellMat); % Number of elements in each cell
acum = [0; cumsum(numElem(:))]; % Acumulated sum
totalPoints = sum(numElem(:));
numParam = numel(parameters);
newParameters = mat2cell(zeros(totalPoints, numParam), totalPoints, ones(1, numParam));
dataVector = zeros(totalPoints, 1);

for k = 1:numel(cellMat)
    dataVector(acum(k)+1:acum(k+1)) = cellMat{k};
    for p = 1:numParam
        newParameters{p}(acum(k)+1:acum(k+1)) = parameters{p}(k); 
    end
end


end

