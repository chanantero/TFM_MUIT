%% Pruebas para exportar figuras a SVG, seleccionando qué se exporta como texto y qué como path

ax = axes(figure);

% Cuando el XLabel.Interpreter es latex, se exporta como path. Cuando es
% none o es tex, depende del contenido del texto se exporta como path o no.
% Si el texto contiene algún símbolo de latex o tex, entonces se exporta
% como path. Si contiene texto plano, como texto

ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = 'hola';
printfig(ax.Parent, imagesPath, 'PruebaA', 'svg');

ax.XLabel.Interpreter = 'none';
printfig(ax.Parent, imagesPath, 'PruebaB', 'svg');

% La idea es, por tanto, exportar como path (interpreter = 'latex') todo
% aquello que no se quiera como texto en Latex. Para ello, hay que eliminar
% todo aquél texto que no se quiera como path.

% Después, se exporta una imagen con todo el texto como texto plano (hay
% que sustituir los comandos de latex por etiquetas identificables), y se
% busca en el archivo svg por ese fragmento de texto para incluírlo en el
% primer svg.

h = figure;
ax = axes(h);
plot(ax, cos(0:0.1:2*pi));
ax.Title.String = 'titulo';
ax.Title.Interpreter = 'none';
ax.XLabel.String = 'Eje x';
ax.YLabel.String = 'Eje y';
legend(ax, 'hopoa')
colorbar

options.TickLabels2Latex = false;
Plot2LaTeX(h, [imagesPath, 'pruebaPlot2LaTeXfalse'], options)

options.TickLabels2Latex = true;
Plot2LaTeX(h, [imagesPath, 'pruebaPlot2LaTeX'], options)


h2 = copyobj(h, 0);
hBackup = copyobj(h, 0);
delete(h2.Children)
copyobj(hBackup.Children, h2)
