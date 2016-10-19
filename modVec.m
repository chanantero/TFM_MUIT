function [ modules ] = modVec( vectors )

modules = sqrt(sum(vectors.^2, 2));

end

