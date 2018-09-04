function Plot2LaTeX( h, filename, options )
% Modified by Rub�n Chuli� (chanantero@hotmail.com)
%PLOT2LATEX saves matlab figure as a pdf file in vector format for
%inclusion into LaTeX. Requires free and open-source vector graphics
%editor Inkscape.
%
%   PLOT2LATEX(h,filename) saves figure with handle h to a file specified by
%   filename, without extention. Filename can contain a relative location
%   (e.g. 'images\title') to save the figure to different location.
%
%   PLOT2LATEX(h,filename, options) saves figure with specified options.
%   The y-offset of all text can be modified using options.yCorrFactor.
%   The default is options.yCorrFactor = 0.8. The units are px. With
%   options.Renderer the renderer of the figure can be specified:
%   ('opengl', 'painters').
%
%   PLOT2LATEX requires a installation of Inkscape. The program's
%   location has to be 'hard coded' into this matlab file if it differs
%   from 'c:\Program Files (x86)\Inkscape\Inkscape.exe'. Please specify
%   your inscape file location by modifying DIR_INKSC variable on the
%   first line of the actual code.
%
%   PLOT2LATEX saves the figures to .svg format. It invokes Inkscape to
%   save the svg to a .pdf and .pdf_tex file to be incorporated into LaTeX
%   document using \begin{figure} \input{image.pdf_tex} \end{figure}.
%   More information on the svg to pdf conversion can be found here:
%   ftp://ftp.fu-berlin.de/tex/CTAN/info/svg-inkscape/InkscapePDFLaTeX.pdf
%
%   PLOT2LATEX produces three files: .svg, .pdf, .pfd_tex. The .svg-file
%   contains vector image. The .pdf-file contains figure without the text.
%   The .pdf_tex-file contains the text including locations and other
%   type setting.
%
%   The produced .svg file can be manually modified in Inkscape and
%   included into the .tex file using the using the built-in "save to pdf"
%   functionality of Inkscape.
%
%   PLOT2LATEX saves the figure to a svg and pdf file with the
%   approximately the same width and height. Specify the Font size and size
%   within Matlab for correct conversion.
%
%   Workflow
%   - Matlab renames all strings of the figure to labels. The strings are
%   stored to be used later. To prevent a change in texbox size, labels are
%   padded to match the size of the texbox.
%   - Matlab saves the figure with labels to a svg file.
%   - Matlab opens the svg file and restores the labels  wiht the original
%   string
%   - Matlab invokes Inkscape to save the svg file to a pdf + pdf_tex file.
%   - The pdf_tex is to be included into LaTeX.
%
%   Features:
%   - Complex figures such as plotyy, logarithmic scales.
%   - It parses LaTeX code, even if it is not supported by Matlab LaTeX.
%   - Supports real transparency.
%   - SVG is a better supported, maintained and editable format than eps
%   - SVG allows simple manual modification into Inkscape.
%
%   Limitation:
%   - Text resize is still done in PLOT2LATEX. The LaTeX fonts in matlab do
%   not correspond completely with the LaTeX font size.
%   - Legend size is not always correct, use \hspace or \vspace in matlab
%   legend to achieve a nicer fit. Requires some iterations.
%   - Rotating 3D images using toolbar does not work, using view([]) works.
%   - Text boxes wiht LaTeX code which is not interpretable by matlab
%   results in too long text boxes.
%   - Very large figures sometimes result in very large waiting times.
%   - Older versions than matlab 2014b are not supported.
%   - PLOT2LATEX currently does not work with titles consisting of multiple
%   lines.
%   - PLOT2LATEX does not work with annotation textbox objects.
%   - PLOT2LATEX does not suport colored text.
%
%   Trouble shooting
%   - For Unix users: use the installation folder such as:
%   '/Applications/Inkscape.app/Contents/Resources/script ' as location.
%   - For Unix users: For some users the bash profiles do not allow to call
%   Inkscape in Matlab via bash. Therefore change the bash profile in Matlab
%   to something similar as setenv('DYLD_LIBRARY_PATH','/usr/local/bin/').
%   The bash profile location can be found by using '/usr/bin/env bash'

%   To do:
%   - Restore Interpreter instead of putting it to LaTeX
%   - Annotation textbox objects
%   - Allow multiple line text
%   - Use findall(h,'-property','String')
%   - Speed up code by smarter string replacent of SVG file
%   - Resize of legend box using: [h,icons,plots,str] = legend(); (not so simple)
%   - PLOT2LATEX does not suport colored text. (Matlab limitation in saving to sgv)
%   - Size difference .svg and .fig if specifying units other than px.
%       (Matlab limitation?)

%   Version:  1.2
%   Autor:    J.J. de Jong, K.G.P. Folkersma
%   Date:     20/04/2016
%   Contact:  j.j.dejong@utwente.nl

%   Change log
%   v 1.1 - 02/09/2015 (not released)
%   - Made compatible for Unix systems
%   - Added a waitbar
%   - Corrected the help file
%   v 1.2 - 20/04/2016
%   - Fixed file names with spaces in the name. (Not adviced to use in latex though)
%   - Escape special characters in XML (<,>,',",&) -> (&lt;,&gt;,&apos;,&quot;,&amp;)

%% Config and checks
%   Specify location of your inkscape installation
DIR_INKSC = 'c:\Program Files\Inkscape\Inkscape.exe';

% initize waitbar
nStep = 5; Step = 0;
hWaitBar = waitbar(Step/nStep,'Initializing');

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

%% process options, first set default
yCorrFactor = 0.8; % default
% do not set default renderer

if nargin > 2
    if isfield(options,'yCorrFactor')
        yCorrFactor = options.yCorrFactor; % y offset of the text in px
    end
    
    if isfield(options,'Renderer') %WARNING: large size figures can become very large
        h.Renderer = options.Renderer; % set render
    end
    
    if isfield(options, 'TickLabels2Latex')
        TickLabelToLatex = options.TickLabels2Latex;
    else
        TickLabelToLatex = true;
    end
    
    if isfield(options, 'Legend2Latex')
        LegendToLatex = options.Legend2Latex;
    else
        LegendToLatex = true;
    end
else
    TickLabelToLatex = true;
    LegendToLatex = true;
end

% TickLabelToLatex should always be true, as I checked. The reason for this
% is that Plot2LaTeX edits the alignment of text so when inserted in a
% LaTeX document, the text appears in the right position. You can set TickLabelToLatex
% to false, but the result will always be worse. I leave it here just to
% gain insight into how the function works in future reviews. Before
% executing the function ChangeInterpreter, the saved svg makes number as
% paths. But after executing it, it exports all the text as text. It is a
% line of future research.
h = copyobj(h, 0); % Don't touch original

h2 = copyobj(h, 0);

%% Configure the figure that will be exported as path to svg
TexObj = findall(h2,'Type','Text'); % normal text, titels, x y z labels
LegObj = findall(h2,'Type','Legend'); % legend objects
axPath = findall(h2,'Type','Axes');  % axes containing x y z ticklabel
ColObj = findall(h2,'Type','Colorbar'); % containg color bar tick

for k = 1:numel(TexObj)
    TexObj(k).String = [];
end

if LegendToLatex
for k = 1:numel(LegObj)
    pos = LegObj(k).Position;
    N = length(LegObj(k).String);

    for n = 1:N
        LegObj(k).String{n} = '';
    end
    
    while pos(3) >= LegObj(k).Position(3) % first label of legend should match box size
        LegObj(k).String{1} = [LegObj(k).String{1},'.'];
    end
    LegObj(k).String{1} = LegObj(k).String{1}(1:end-1);
    
    LegObj(k).TextColor = [1 1 1];
    
    drawnow
    LegObj(k).Position = pos;
end    
end

for k = 1:numel(ColObj)
    ColObj(k).TickLabels = [];
end

if TickLabelToLatex
    for k = 1:numel(axPath)
        axPath(k).XTickLabel = [];
        axPath(k).YTickLabel = [];
    end
end

% Export
ChangeInterpreter(h2, 'Latex') % In order to be exported as path
saveas(h2, [filename, '_path'], 'svg')
% saveas(h2, [filename, '_pathCopy'], 'svg')
close(h2)

%% Find all objects with text
TexObj = findall(h,'Type','Text'); % normal text, titels, x y z labels
LegObj = findall(h,'Type','Legend'); % legend objects
AxeObj = findall(h,'Type','Axes');  % axes containing x y z ticklabel
ColObj = findall(h,'Type','Colorbar'); % containg color bar tick

PosAnchSVG      = {'start','middle','end'};
PosAligmentSVG  = {'start','center','end'};
PosAligmentMAT  = {'left','center','right'};

ChangeInterpreter(h,'Latex')
h.PaperPositionMode = 'auto'; % Keep current size

n_Axe = length(LegObj);
for i = 1:n_Axe % scale text omit in next version
    LegPos(i,:) = LegObj(i).Position;
end

%% Replace text with a label
Step = Step + 1;
waitbar(Step/nStep,hWaitBar,'Replacing text with labels');

iLabel = 0; % generate label iterator

n_TexObj = length(TexObj);
for i = 1:n_TexObj % do for text, titles and axes labels
    iLabel = iLabel + 1;
    
    % find text string
    Labels(iLabel).TrueText = TexObj(i).String;
    
    % find text aligment
    Labels(iLabel).Alignment = PosAligmentSVG(...
        find(ismember(...
        PosAligmentMAT,...
        TexObj(i).HorizontalAlignment)));
    % find achor aligment svg uses this
    Labels(iLabel).Anchor = PosAnchSVG(...
        find(ismember(...
        PosAligmentMAT,...
        TexObj(i).HorizontalAlignment)));
    % generate label
    Labels(iLabel).LabelText = LabelText(iLabel);
    
    %find text posiont
    Labels(iLabel).Position = TexObj(i).Position;
    
    % replace string with label
    TexObj(i).String = LabelText(iLabel);
end

% do similar for legend objects
if LegendToLatex
n_LegObj = length(LegObj);
iLegEntry = 0;
for i = 1:n_LegObj
    n_Str = length(LegObj(i).String);
    
    iLegEntry = iLegEntry + 1;
    iLabel = iLabel + 1;
    
    Labels(iLabel).TrueText = LegObj(i).String{1};
    Labels(iLabel).Alignment = PosAligmentSVG(1); % legends are always left aligned
    Labels(iLabel).Anchor = PosAnchSVG(1);
    
    % generate legend label padded with dots to fill text box
    LegObj(i).String{1} = LegText(iLegEntry);
    while LegPos(i,3) >= LegObj(i).Position(3) % first label of legend should match box size
        LegObj(i).String{1} = [LegObj(i).String{1},'.'];
    end
    if length(LegObj(i).String{1}) > 1
    LegObj(i).String{1} = LegObj(i).String{1}(1:end-1);
    end
    Labels(iLabel).LabelText = LegObj(i).String{1}; % write as label
    
    for j = 2:n_Str % do short as possible label for other entries
        iLegEntry = iLegEntry + 1;
        iLabel = iLabel + 1;
        Labels(iLabel).TrueText = LegObj(i).String{j};
        Labels(iLabel).Alignment = PosAligmentSVG(1);
        Labels(iLabel).Anchor = PosAnchSVG(1);
        Labels(iLabel).LabelText = LegText(iLegEntry);
        LegObj(i).String{j} = LegText(iLegEntry);
    end
end
end

% do similar for axes objects, XTick, YTick, ZTick
if TickLabelToLatex
    n_AxeObj = length(AxeObj);
    for i = 1:n_AxeObj
        n_Str = length(AxeObj(i).XTickLabel);
        for j = 1:n_Str
            iLabel = iLabel + 1;
            Labels(iLabel).TrueText = AxeObj(i).XTickLabel{j};
            Labels(iLabel).Alignment = PosAligmentSVG(2);
            Labels(iLabel).Anchor = PosAnchSVG(2);
            Labels(iLabel).LabelText = LabelText(iLabel);
            AxeObj(i).XTickLabel{j} = LabelText(iLabel);
        end
        
        isRightAx = strcmp(AxeObj(i).YAxisLocation,'right'); % exeption for yy-plot
        n_Str = length(AxeObj(i).YTickLabel);
        for j = 1:n_Str
            iLabel = iLabel + 1;
            Labels(iLabel).TrueText = AxeObj(i).YTickLabel{j};
            if isRightAx % exeption for yy-plot, aligment is left for the right axis
                Labels(iLabel).Alignment = PosAligmentSVG(1);
                Labels(iLabel).Anchor = PosAnchSVG(1);
            else % normal y labels are right aligned
                Labels(iLabel).Alignment = PosAligmentSVG(3);
                Labels(iLabel).Anchor = PosAnchSVG(3);
            end
            Labels(iLabel).LabelText = LabelText(iLabel);
            AxeObj(i).YTickLabel{j} = LabelText(iLabel);
        end
        
        n_Str = length(AxeObj(i).ZTickLabel);
        for j = 1:n_Str
            iLabel = iLabel + 1;
            Labels(iLabel).TrueText = AxeObj(i).ZTickLabel{j};
            Labels(iLabel).Alignment = PosAligmentSVG(3);
            Labels(iLabel).Anchor = PosAnchSVG(3);
            Labels(iLabel).LabelText = LabelText(iLabel);
            AxeObj(i).ZTickLabel{j} = LabelText(iLabel);
        end
    end
else
    for k = 1:length(AxeObj)
        XLabelPos = AxeObj(k).XLabel.Position;
        YLabelPos = AxeObj(k).YLabel.Position;
        AxeObj(k).XTickLabel = [];
        AxeObj(k).YTickLabel = [];
        AxeObj(k).XLabel.Position = XLabelPos;
        AxeObj(k).YLabel.Position = YLabelPos;
    end
end

% do similar for color bar objects
n_ColObj = length(ColObj);
for i = 1:n_ColObj
    isAxIn = strcmp(ColObj(i).AxisLocation,'in'); % find internal external text location
    isAxEast = strcmp(ColObj(i).Location,'east'); % find location
    isRightAx = isAxIn ~= isAxEast;
    
    n_Str = length(ColObj(i).TickLabels);
    for j = 1:n_Str
        iLabel = iLabel + 1;
        Labels(iLabel).TrueText = ColObj(i).TickLabels{j};
        if isRightAx % if text is right aligned
            Labels(iLabel).Alignment = PosAligmentSVG(1);
            Labels(iLabel).Anchor = PosAnchSVG(1);
        else % if text is left aligned
            Labels(iLabel).Alignment = PosAligmentSVG(3);
            Labels(iLabel).Anchor = PosAnchSVG(3);
        end
        Labels(iLabel).LabelText = LabelText(iLabel);
        ColObj(i).TickLabels{j} = LabelText(iLabel);
    end
end
nLabel = iLabel;

% set text interpreter to plain text
ChangeInterpreter(h,'none');

%% Save to fig and SVG
Step = Step + 1;
waitbar(Step/nStep,hWaitBar,'Saving figure to .svg file');

% savefig(h,[filename,'_temp']); % to see the intermediate situation
saveas(h,filename,'svg'); % export to svg
% saveas(h,[filename, 'Copy'],'svg'); % export to svg
%% Modify SVG file to replace labels with original text
Step = Step + 1;
waitbar(Step/nStep,hWaitBar,'Restoring text in .svg file');

for iLabel = 1:nLabel
    Labels(iLabel).XMLText = EscapeXML(Labels(iLabel).TrueText);
end

try
    fin = fopen([filename,'.svg'], 'r'); % open svg file
    strNormal = char(fread(fin)');
    fclose(fin);
    
    [startInd, endInd, tokens] = regexp(strNormal, '<g[^>]*?><text[^>]*>(\S*?)</text\s*>\s*</g\s*>', 'start', 'end', 'tokens', 'all');
    
    % For each label, find the substring than contains it
    numTexts = length(startInd);
    textFragments = cell(numTexts, 1);
    textContent = cell(numTexts, 1);
    for t = 1:numTexts
        textFragments{t} = strNormal(startInd(t):endInd(t));
        textContent{t} = tokens{t}{1};
    end
    labels = {Labels.LabelText};
    [isFound, index] = ismember(labels, textContent);
    LabelsFound = Labels(isFound);
    index = index(isFound);
    
    % Substitute the label by the original text
    numLabels = length(LabelsFound);
    newTextFragments = cell(numLabels, 1);
    for l = 1:numLabels
        text = textFragments{index(l)};
        % Substitute the label by the original text
        text = strrep(text, LabelsFound(l).LabelText, LabelsFound(l).XMLText);
        
        % Add some parameters to the style attribute
        en = regexp(text, '<g[^>]*style=".*?"', 'end');
        text = [text(1:en-1), ' text-align:', LabelsFound(l).Alignment{1},...
            '; text-anchor:', LabelsFound(l).Anchor{1}, text(en:end)];
        
        % Change the value of x and y coordinates
        text = regexprep(text, 'x="\S*"', 'x="0"');
        
        [token, tokenExtents] = regexp(text,'<text[^>]*?y="(\S*)"', 'tokens', 'tokenExtents');
        yOffsetStr = token{1}{1};
        yOffsetStrStartInd = tokenExtents{1}(1);
        yOffsetStrEndInd = tokenExtents{1}(2);
        yOffset = str2double(yOffsetStr);
        text = [text(1:yOffsetStrStartInd-1), num2str(yOffset*yCorrFactor), text(yOffsetStrEndInd+1:end)];
        
        newTextFragments{l} = text;
    end
    newTextFragment = strjoin(newTextFragments, '\n');
    
    % Write newTextFragments in fout
    fout = fopen([filename, '_path.svg'], 'r');
    strPath = char(fread(fout)');
    fclose(fout);
    
    start = regexp(strPath, '</g\s*>\s*</svg\s*>', 'start');
    strPath = [strPath(1:start-1), newTextFragment, strPath(start:end)];
    
    fout = fopen([filename, '_path.svg'], 'w');
    fwrite(fout, strPath);
    fclose(fout);
    
catch ME
    fclose(fin);
    fclose(fout);
    rethrow(ME)
end

movefile([filename,'_path.svg'],[filename,'.svg'])
%% Invoke Inkscape to generate PDF + LaTeX
Step = Step + 1;
waitbar(Step/nStep,hWaitBar,'Saving .svg to .pdf file');

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

%% Restore figure in matlab
oldWay = false;
if oldWay
% for nicety replace labels with the original text
Step = Step + 1;
waitbar(Step/nStep,hWaitBar,'Restoring Matlab figure');

if nLabel > 0
    iLabel = 0;
    n_TexObj = length(TexObj);
    for i = 1:n_TexObj %
        iLabel = iLabel + 1;
        TexObj(i).String = Labels(iLabel).TrueText;
    end
    
    n_LegObj = length(LegObj);
    for i = 1:n_LegObj
        n_Str = length(LegObj(i).String);
        for j = 1:n_Str
            iLabel = iLabel + 1;
            LegObj(i).String{j} = Labels(iLabel).TrueText;
        end
    end
    
    if TickLabelToLatex
        n_AxeObj = length(AxeObj);
        for i = 1:n_AxeObj
            n_Str = length(AxeObj(i).XTickLabel);
            for j = 1:n_Str
                iLabel = iLabel + 1;
                AxeObj(i).XTickLabel{j} = Labels(iLabel).TrueText;
            end
            
            n_Str = length(AxeObj(i).YTickLabel);
            for j = 1:n_Str
                iLabel = iLabel + 1;
                AxeObj(i).YTickLabel{j} = Labels(iLabel).TrueText;
            end
            
            n_Str = length(AxeObj(i).ZTickLabel);
            for j = 1:n_Str
                iLabel = iLabel + 1;
                AxeObj(i).ZTickLabel{j} = Labels(iLabel).TrueText;
            end
        end
    end
    
    n_AxeObj = length(ColObj);
    for i = 1:n_AxeObj
        n_Str = length(ColObj(i).TickLabels);
        for j = 1:n_Str
            iLabel = iLabel + 1;
            ColObj(i).TickLabels{j} = Labels(iLabel).TrueText;
        end
    end
end

% restore interpreter
ChangeInterpreter(gcf,'Latex')
else
    delete(h)
end

close(hWaitBar);
end

function Str = LabelText(iLabel)
% LABELTEXT generates labels based on label number
Str = 'X000';
idStr = num2str(iLabel);
nStr = length(idStr);
Str(end - nStr + 1 : end ) = idStr;
end

function Str = LegText(iLedEntry)
% LEGTEXT generates legend labels based on legend entry number
Str = num2str(iLedEntry);
end

function ChangeInterpreter(h,Interpreter)
% CHANGEINTERPRETER puts interpeters in figure h to Interpreter

if strcmp('Interpreter', 'none')
    saveas(h, 'prueba_a','svg');
end
TexObj = findall(h,'Type','Text');
LegObj = findall(h,'Type','Legend');
AxeObj = findall(h,'Type','Axes');
ColObj = findall(h,'Type','Colorbar');

Obj = [TexObj;LegObj]; % Tex and Legend opbjects can be treated similar



n_Obj = length(Obj);
for i = 1:n_Obj
    Obj(i).Interpreter = Interpreter;
end

Obj = [AxeObj;ColObj]; % Axes and colorbar opbjects can be treated similar

n_Obj = length(Obj);
for i = 1:n_Obj
    Obj(i).TickLabelInterpreter = Interpreter;
end
if strcmp('Interpreter', 'none')
    saveas(h,'prueba_b','svg');
end
end

function strXML = EscapeXML(str)
% ESCAPEXML repaces special characters(<,>,',",&) -> (&lt;,&gt;,&apos;,&quot;,&amp;)
escChar = {'&','<','>','''','"'};
repChar = {'&amp;','&lt;','&gt;','&apos;','&quot;'};
strXML = regexprep(str,escChar,repChar);
end