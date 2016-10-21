function new_pageSize = check_pageSize(pageSize)
% Check if pageSize is correct
default_pageSize = 4096;
if pageSize < 2048 || pageSize > 16384
    fprintf('Invalid page size input,will be set to default page size %d',default_pageSize);
    new_pageSize = default_pageSize;
else
    fprintf('Page Size set to %d.\n',pageSize);
    new_pageSize = pageSize;
end

end