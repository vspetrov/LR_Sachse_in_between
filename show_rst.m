clear all;
fd = fopen('rst.bin','r');
Size = 15;
Count = 100;
data=fread(fd,[Size Count],'double');
imagesc(data');
colorbar;
