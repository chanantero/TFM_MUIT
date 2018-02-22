function drawWFSarrayFun( fileName, varargin )

svgText = WFSarraySVG(varargin{:});

% Write file
destFile = fopen(fileName, 'w', 'n', 'UTF-8');
fwrite(destFile, svgText, 'char');
fclose(destFile);

end

