function Figure2PDFplusLatex( h, filename )

%% Config and checks
%   Specify location of your inkscape installation
DIR_INKSC = 'c:\Program Files\Inkscape\Inkscape.exe'; 

%test if installation is correct
if ~exist(DIR_INKSC,'file')
    error([DIR_INKSC, ' cannot be found, check installation location'])
end

if verLessThan('matlab', '8.4.0.')
	error('Older versions than Matlab 2014b are not supported')
end

if ~strcmp(h.Type,'figure')
    error('h object is not a figure')
end

%% Save to fig and SVG
saveas(h, filename, 'svg'); % export to svg

%% Invoke Inkscape to generate PDF + LaTeX
[pathFolder, name, ~] = fileparts(filename);
currentFolder = pwd;

if isempty(pathFolder)
    DIR_FIG = [currentFolder,'\'];
else
    DIR_FIG = [pathFolder, '\'];
end

cd(DIR_FIG)
[status,cmdout] = system(['"', DIR_INKSC, '"',...
                ' "', DIR_FIG, name,'.svg"', ...
                ' ','--export-pdf',...
                ' "', DIR_FIG, name,'.pdf"',...
                ' ','--export-latex',...
                ' ','-export-area-drawing']);
cd(currentFolder)
         
% test if a .pdf and .pdf_tex file exist
if exist([filename,'.pdf'],'file')~= 2 || exist([filename,'.pdf_tex'],'file')~= 2
    cmdout
    warning('No .pdf or .pdf_tex file produced, please check your Inkscape installation and specify installation directory correctly.')
end

