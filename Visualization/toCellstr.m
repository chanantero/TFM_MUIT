function [ output ] = toCellstr( x )

if isnumeric(x)
    output = numToCellstr(x);
elseif iscategorical(x)
    output = cellstr(x);
end

end

function output = numToCellstr( x )

output = cell(size(x));

for k = 1:numel(x)
    output{k} = num2str(x(k));
end

end