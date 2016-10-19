function fig = singleFreqRepr( r, U, t, f, varargin )

p = inputParser;

addOptional(p, 'clim', []);
addOptional(p, 'axs', []);

parse(p, varargin{:});

clim = p.Results.clim;
axs = p.Results.axs;

if ~iscell(U)
    U = {U};
end
N = numel(U);

if isempty(axs)
    fig = figure;
    axs = gobjects(N, 1);
    for k = 1:N
        axs(k) = axes(fig, 'Units', 'normalized', 'OuterPosition', [(k-1)/N, 0, 1/N, 1]);
    end
end

s = gobjects(N,1);
for k = 1:N
    s(k) = scatter3(axs(k), r(:,1), r(:,2), r(:,3), [], real(U{k}), '.');
end

if ~isempty(clim)
    for k = 1:N
        axs(k).CLim = clim;
    end
else
    for k = 1:N
        axs(k).CLim = [-max(abs(U{k})), max(abs(U{k}))];
    end
end

incrT = diff(t);
incrT = [incrT, incrT(end)];
for l = 1:numel(t)
    for k = 1:N
        U_repr = U{k}*exp(1i*f*2*pi*t(l));
        s(k).CData = real(U_repr);
    end
    pause(incrT(l))
end

end

