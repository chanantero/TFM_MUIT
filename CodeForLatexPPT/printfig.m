function printfig(fig, path, name, formats)
    if ~iscellstr(formats)
        formats = {formats};
    end

    numFormats = numel(formats);
    for f = 1:numFormats
        format = formats{f};
        ext = ['.', format];
        formatOption = getFormatOption(format);
        print(fig, [path, name, ext], formatOption);
    end
end

function formatOption = getFormatOption(format)

formats = {'pdf', 'eps', 'svg', 'emf'};
formatOptions = {'-dpdf', '-depsc', '-dsvg', '-dmeta'};

formatOption = formatOptions{ismember(formats, format)};

end