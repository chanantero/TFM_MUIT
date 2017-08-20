function [maxtab, mintab] = peakdet(v, delta, x)
% Modified version of the original function by Eli Billauer

%PEAKDET Detect peaks in a vector
%        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
%        maxima and minima ("peaks") in the vector V.
%        MAXTAB and MINTAB consists of two columns. Column 1
%        contains indices in V, and column 2 the found values.
%
%        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
%        in MAXTAB and MINTAB are replaced with the corresponding
%        X-values.
%
%        A point is considered a maximum peak if it has the maximal
%        value, and was preceded (to the left) by a value lower by
%        DELTA.

% Eli Billauer, 3.4.05 (Explicitly not copyrighted).
% This function is released to the public domain; Any use is allowed.

v = v(:); % Just in case this wasn't a proper vector

if nargin < 3
    x = (1:length(v))';
else
    x = x(:);
    if length(v)~= length(x)
        error('Input vectors v and x must have same length');
    end
end

if (length(delta(:)))>1
    error('Input argument DELTA must be a scalar');
end

if delta <= 0
    error('Input argument DELTA must be positive');
end

mn = Inf; mx = -Inf;
mnind = NaN; mxind = NaN;
numpicos=0;
numvalles=0;

lookformax = false;
orderStablished = false;
maxtab=zeros(length(v),2);
mintab=zeros(length(v),2);

for i=1:length(v)
    this = v(i);
    if this > mx, mx = this; mxind = i; end
    if this < mn, mn = this; mnind = i; end
    
    if(orderStablished)
        if lookformax
            if this < mx - delta
                maxtab(numpicos+1,:) = [mxind mx];
                mn = this; mnind = i;
                lookformax = false;
                numpicos=numpicos+1;
            end
        else
            if this > mn + delta
                mintab(numvalles+1,:) = [mnind mn];
                mx = this; mxind = i;
                lookformax = true;
                numvalles=numvalles+1;
            end
        end
    else
        if this < mx - delta
            maxtab(numpicos+1,:) = [mxind mx];
            mn = this; mnind = i;
            lookformax = false;
            numpicos = numpicos+1;
            orderStablished = true;
        end
        if this > mn + delta
            mintab(numvalles+1,:) = [mnind mn];
            mx = this; mxind = i;
            lookformax = true;
            numvalles=numvalles+1;
            orderStablished = true;
        end
    end
end

if(orderStablished)
    if(lookformax)
        maxtab = [maxtab(1:numpicos, :); [mxind mx]];
        mintab = mintab(1:numvalles, :);
    else
        maxtab = maxtab(1:numpicos,:);
        mintab = [mintab(1:numvalles,:); [mnind mn]];
    end
else
    maxtab = maxtab(1:numpicos,:);
    mintab = mintab(1:numvalles, :);
end

onlyFullPeaks = false;
if(onlyFullPeaks) % Only maximum surrounded by 2 minimum. Only minimum surrounded by 2 maximum.
    if(orderStablished)
        % Cut the first peak
        [~, typeFirstPeak] = min([maxtab(1,1), mintab(1,1)]);
        if(typeFirstPeak == 1) % Primero hubo un máximo
            maxtab = maxtab(2:numpicos, :);
        else % Primero hubo un mínimo
            mintab = mintab(2:numpicos, :);
        end
        
        % Cut the last peak
        %     [~, typeLastPeak] = min([maxtab(end,1), mintab(end,1)]);
        if(lookformax) % Last peak was a maximum
            maxtab = maxtab(1:end-1, :);
        else % Last peak was a minimum
            mintab = mintab(1:end-1, :);
        end
    end
end

maxtab(:,1) = x(maxtab(:,1));
mintab(:,1) = x(mintab(:,1));

end
