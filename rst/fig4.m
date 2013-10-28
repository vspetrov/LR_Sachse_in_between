#!/usr/bin/octave
clear all;

fd=fopen('chain_30_100_30_time5000ms.bin','r');
x=fread(fd, [160 2500], 'double');
imagesc(x, [-80 20]);
colorbar;

set(gca(),'xtick', [1 1250 2500]);
set(gca(),'xticklabel', [0 2500 5000]);
xlabel('time, ms')
set(gca(),'ytick', [1 31 130 160]);
set(gca(),'yticklabel', [1 31  131 160 ]);
ylabel('i')


print -dpng fig4.png');
% print('fig4.pdf');
