function plotlog_pn(filename)

file = fopen(filename);
sizeX = fread(file, 1, 'int64')
sizeY = fread(file, 1, 'int64')
lengthX = fread(file, 1, 'double')
lengthY = fread(file, 1, 'double')
numMoments = fread(file, 1, 'int64')
data = fread(file, sizeX * sizeY * numMoments, 'double');
fclose(file);

M = zeros(sizeX, sizeY);
for i = 1:sizeX
    for j = 1:sizeY
        k = (i-1) * sizeY * numMoments + (j-1) * numMoments + 1;
        M(i,j) = data(k) / (2.0 * sqrt(pi));
        if M(i,j) < 0
            M(i,j) = 0;
        end
    end
end

dx = lengthX / sizeX;
dy = lengthY / sizeY;
x = dx * (0:(sizeX - 1))' + dx / 2 - lengthX / 2;
y = dy * (0:(sizeY - 1))' + dy / 2 - lengthY / 2;

figure;
set(gca, 'FontName', 'Helvetica');
imagesc(y, x, log10(M), [max(max(log10(M))) - 7, max(max(log10(M)))]);
%surf(y, x, M);
%shading interp;
axis square;
h = colorbar;
set(h, 'FontName', 'Helvetica');

%title('Solution', 'FontSize', 20, 'FontName', 'Helvetica');
print('-depsc2', 'solution.eps');
system('epstopdf solution.eps');
