function [x,y] = find_middle()
alt = generate_array();
x = (min(alt(1,:)) + max(alt(1,:))) / 2;
y = (min(alt(2,:)) + max(alt(2,:))) / 2;
end

