clear all;
fd=fopen('rst.bin','r');
data = fread(fd,[1 200], 'double');
%imagesc(data);
%colorbar;
fd=fclose(fd);
