function [x,y] = find_middle()
% Computes the middle point of the array of speakers. The x coordinate is
% taken as the average of the maximum and minimum x coordinates of the
% array. Ídem for the y coordinate.
alt = generate_array();
x = (min(alt(1,:)) + max(alt(1,:))) / 2;
y = (min(alt(2,:)) + max(alt(2,:))) / 2;
end

