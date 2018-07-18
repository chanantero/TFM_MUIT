ax = axes(figure);

currentFolder = pwd;
cd(imagesPath); % Needed for inkscape to link svg files properly
options.TickLabels2Latex = false;
Plot2LaTeX(ax.Parent, 'PruebaPlot2LaTeX', options)
cd(currentFolder)

ax.XTickLabel{1} = '$0.3131$';
ax.XTickLabel{2} = '$\alpha$';
ax.TickLabelInterpreter= 'latex';
printfig(ax.Parent, imagesPath, 'PruebaPrintFigLatex', 'svg')
ax.TickLabelInterpreter= 'tex';
printfig(ax.Parent, imagesPath, 'PruebaPrintFigTex', 'svg')
ax.TickLabelInterpreter= 'none';
printfig(ax.Parent, imagesPath, 'PruebaPrintFigNone', 'svg')

printfig(ax.Parent, imagesPath, 'PruebaPrintFig', 'pdf')

Figure2PDFplusLatex( ax.Parent, [imagesPath, 'PruebaFigure2PDFplusLatex']);
options.TickLabels2Latex = false;
Plot2LaTeX(ax.Parent, [imagesPath, 'PruebaPlot2LaTeXfalse'], options)
options.TickLabels2Latex = true;
Plot2LaTeX(ax.Parent, [imagesPath, 'PruebaPlot2LaTeXtrue'], options)
