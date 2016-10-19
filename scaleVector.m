function vectors = scaleVector(vectors, scalingFactors)

vectors = vectors.*repmat(scalingFactors, 1, 3);

end