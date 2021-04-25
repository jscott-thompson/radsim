function plot_rd_map(range_swath,doppler_swath,rd_map)
%PLOT_RD_MAP Plots a range Doppler map

imagesc(doppler_swath,range_swath,10*log10(rd_map'));
title('Range Doppler Map');
xlabel('Doppler [Hz]');
ylabel('Range [m]');
