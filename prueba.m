% No funciona con auriculares!!! ¿Por qué??
obj.setNumWFSarraySources(2);

obj.reproduceAndRecord('main', 'soundTime', 2);

x = signalFunc(0, 44100);
ax = axes(figure);
plot(ax, x(:,1))
