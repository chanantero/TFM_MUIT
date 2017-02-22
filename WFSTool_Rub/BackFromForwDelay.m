function backDelay = BackFromForwDelay(prevForwDelay, forwDelay, storedBackDelay, Fs)
% We assume that the delay follows the theory in...

numPrevDelay = numel(prevForwDelay);
numBackDelay = numel(storedBackDelay);

firstSample = numBackDelay + 1;
lastSample = numel(forwDelay) + floor(forwDelay(end)*Fs); % With the floor(...) we make sure each interpolated point is calculated with one point at each side

backDelay = forward2BackwardDelay( [prevForwDelay; forwDelay], [prevForwDelay; forwDelay], Fs, firstSample + numPrevDelay, lastSample + numPrevDelay );

backDelay = [storedBackDelay; backDelay'];

end
