x = 0:10;
y = 11:15;

[X, Y] = ndgrid(x, y);
data = X+Y;

obj = animation({x, y}, {data}, {'x', 'y'}, {'data'}, [], []);