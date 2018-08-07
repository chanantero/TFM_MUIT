workspaceInfo = whos;

bytes = [workspaceInfo.bytes];
names = {workspaceInfo.name};
total = sum(bytes);

[bytesSorted, ind] = sort(bytes, 'descend');
namesSorted = names(ind);

bytesSortedRel = bytesSorted/total;

acumRel = cumsum(bytesSortedRel);

ind = 1:find(acumRel > 0.9, 1, 'first');

ax = axes(figure);
p = pie(ax, bytesSortedRel(ind), namesSorted(ind));
for k = 1:numel(p)
    if isa(p(k), 'matlab.graphics.primitive.Text')
        p(k).Interpreter = 'none';
    end
end
ax.Title.String = ['Total = ', num2str(total)];