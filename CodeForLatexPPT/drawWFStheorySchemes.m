%% Draw WFS theory schemes

ax = axes(figure);
x = [0 1 1 0 0];
y = [0 0 1 1 0];
z = [0 1 2 1 0];
plot3(ax, x, y, z)

printfig(ax.Parent, '', 'prueba', {'eps', 'svg', 'pdf'})