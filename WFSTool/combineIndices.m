function [combKey, mapping] = combineIndices(baseKey, priorityKey)
% We have two vectors of integers: baseKey and priorityKey. Each integer is
% the identification of a destination (mapping). For example, if baseKey is
% [2, 4, 10], it could mean that is the mapping of three signals over three
% audio channels (2nd, 4th and 10th audio channels) of an audio device.
% From those two vectors, we extract a third one: combKey. Each element of
% this third vector also is a key onto the same destination (e.g., audio
% channels). But, we want this third vector to be a combination of the
% original ones.
% For example, we have two set of signals s1 and s2. The first one has
% associated a mapping baseKey = [2, 4, 10]. The second one has a mapping
% priorityKey = [10, 5]. If we wanted to generate a signal set s3, it would
% have a mapping vector combKey = [2, 4, 5, 10]. The signal set s3 would be
% a combination of the sets s1 and s2. Especifically, s3 = [s1(1), s1(2),
% s2(2), s2(1)]. Ase we see, the signals from s1 are the base of signals of
% s3, but in channels that are occupied by both s1 and s2, s2 has priority.
% The output of the function is, not only the combination key vector
% combKey, but also the mapping of the original vectors baseKey and
% priorityKey over combKey.
% Following the previous example, the mapping for the base vector would
% consist on two vectors. The first one is the combination indices and the
% second one is the corresponding base indices. This is, combInd_base = [1,
% 2] and baseInd = [1, 2].
% The mapping for the priority vector would be: combInd_priority = [3, 4]
% and priorityInd = [2, 1];
% The mathematical relation is:
% combKeys(combInd_base) = baseKey(baseInd);
% combKeys(combInd_priority) = priorityKey(priorityInd);
% combInd_base and combInd_priority don't have common elements, i.e., they
% are exclusive.
% If there are two repeated keys in any of the vectors, the first one (the 
% one with the lowest index) has priority.
% Output arguments:
% - mapping. Structure vector of 2 elements. The first one contains the
% mapping for the base vector and the second one for the priority one. It
% has two fields: originInd and destinationInd.

baseKey = baseKey(:);
priorityKey = priorityKey(:);

combKey = unique([baseKey; priorityKey], 'sorted');

[flagComb, priorityInd] = ismember(combKey, priorityKey);
combInd_priority = find(flagComb);
priorityInd = priorityInd(flagComb);

% Do the same but only with the combination key elements that are not
% occupied by the priority assignation
[flagComb, baseInd] = ismember(combKey(~flagComb), baseKey);
combInd_base = find(flagComb);
baseInd = baseInd(flagComb);

mapping = struct('destinationInd', combInd_base, 'originInd', baseInd);
mapping = repmat(mapping, [2, 1]);
mapping(2).destinationInd = combInd_priority;
mapping(2).originInd = priorityInd;

end