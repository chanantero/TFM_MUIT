function [ params, gofs ] = fitInterface( data, models )
% Simplification of the use of fit, so it is more comfortable
% Input arguments:
% - data. N element structure. There are two fields: x and y.  
% - models. Cell string array with numModels elements.

N = numel(data);
numModels = numel(models);

[~, gof] = fit([0; 1], [0; 1], 'poly1');
gofs = repmat(gof, N, numModels);
fs = cell(N, numModels);
len = 0;
for n = 1:N
    fprintf(repmat('\b', 1, len));
    msg = sprintf('%d/%d', n, N);
    fprintf(msg);
    len = numel(msg);
    
    x = data(n).x;
    y = data(n).y;
    for m = 1:numModels
        [f, gof] = fit(x, y, models{m});        
        fs{n, m} = f;
        gofs(n, m) = gof;        
    end
end
fprintf('\n');

% Generate structure arrays of parameters for the fitted data
params = cell(1, numModels);
for m = 1:numModels
    names = coeffnames(fs{1, m});
    numParams = numel(names);
    aux = [names'; cell(1, numParams)];
    paramStruct = repmat(struct(aux{:}), N, 1);
    for n = 1:N
        values = coeffvalues(fs{n, m});
        for p = 1:numParams
            paramStruct(n).(names{p}) = values(p);
        end
    end
    params{m} = paramStruct;
end

end

