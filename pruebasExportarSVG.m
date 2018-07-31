%% Pruebas para exportar figuras a SVG, seleccionando qu� se exporta como texto y qu� como path

ax = axes(figure);

% Cuando el XLabel.Interpreter es latex, se exporta como path. Cuando es
% none o es tex, depende del contenido del texto se exporta como path o no.
% Si el texto contiene alg�n s�mbolo de latex o tex, entonces se exporta
% como path. Si contiene texto plano, como texto

ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = 'hola';
printfig(ax.Parent, imagesPath, 'PruebaA', 'svg');

ax.XLabel.Interpreter = 'none';
printfig(ax.Parent, imagesPath, 'PruebaB', 'svg');

% La idea es, por tanto, exportar como path (interpreter = 'latex') todo
% aquello que no se quiera como texto en Latex. Para ello, hay que eliminar
% todo aqu�l texto que no se quiera como path.

% Despu�s, se exporta una imagen con todo el texto como texto plano (hay
% que sustituir los comandos de latex por etiquetas identificables), y se
% busca en el archivo svg por ese fragmento de texto para inclu�rlo en el
% primer svg.

h = figure;
ax = axes(h);
plot(ax, cos(0:0.1:2*pi));
ax.Title.String = 't�tulo';
ax.XLabel.String = 'Eje x';
ax.YLabel.String = 'Eje y';
legend(ax, 'hopoa')
colorbar

TexObj = findall(h,'Type','Text'); % normal text, titels, x y z labels
LegObj = findall(h,'Type','Legend'); % legend objects
axPath = findall(h,'Type','Axes');  % axes containing x y z ticklabel
ColObj = findall(h,'Type','Colorbar'); % containg color bar tick

for k = 1:numel(TexObj)
    TexObj(k).String = [];
end

for k = 1:numel(LegObj)
    pos = LegObj(k).Position;
    LegObj(k).String = {''};
    drawnow
    LegObj(k).Position = pos;
end
