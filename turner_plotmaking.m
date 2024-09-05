%turner_plotmaking.m

% using dan's plotting functions.
% author TEJ
% date 21 August 2024
clear;clc;clf;
%% loading processed glider data (each 'bindata')
glider33 = load('23503301_bin.mat').bindata;
glider37 = load('23503701_bin.mat').bindata;
glider45 = load('23504501_bin.mat').bindata;
glider56 = load('23505601_bin.mat').bindata;

sglider141 = load('sg141.mat').bindata;
sglider526 = load('sg526.mat').bindata;
sglider528 = load('sg528.mat').bindata;
sglider687 = load('sg687.mat').bindata;

slocum684 = load('osu684.mat').bindata;
slocum685 = load('osu685.mat').bindata;
slocum686 = load('osu686.mat').bindata;


% load mapped data
load('map2023.mat')

%% looking at time

dt_array = datetime(map.time, 'ConvertFrom', 'posixtime');
disp(dt_array(1)); disp(dt_array(end));
%31-May-2023 to 31-Jul-2023

%% psection
% Plot section with x=distance, y=depth, color=var1, contour=var2.

figure(1);
subplot(4,1,1); 
psection(glider33, 't', 'sigma', ':', [], 20:0.5:30, 20:30, 'dist')

subplot(4,1,2)
psection(glider37, 't', 'sigma', ':', [], 20:0.5:30, 20:30, 'dist')

subplot(4,1,3)
psection(glider45, 't', 'sigma', ':', [], 20:0.5:30, 20:30, 'dist')

subplot(4,1,4)
psection(glider56, 't', 'sigma', ':', [], 20:0.5:30, 20:30, 'dist')

%% psection3d
% Plot 3-d section with x=lon, y=lat, z=z, color=var1, contour=var2.

figure(2)
subplot(4,1,1); 
psection3d(glider33, 't', 'sigma',':',[])
subplot(4,1,2)
%psection3d(glider37, 't', 'sigma',':',[]) %throwing error.
subplot(4,1,3)
psection3d(glider45, 't', 'sigma',':',[])
subplot(4,1,4)
psection3d(glider56, 't', 'sigma',':',[])

%% pxz_xy
% Plot sections where x=dist, y=depth, color=var1, contour=var2.
% map vs mapem
% mapw vs mapwem

figure(3)
subplot(4,1,1)
pxz_xy(map,'t','s',3,50,[],[],[])

subplot(4,1,2)
pxz_xy(mapem,'t','s',3,50,[],[],[])

subplot(4,1,3)
pxz_xy(mapw,'w','rogeo',3,50,[],[],[])

subplot(4,1,4)
pxz_xy(mapwem,'w','rogeo',3,50,[],[],[])

%% psectime
% Plot section where x=time, y=depth, color=var1, contour=var2.

figure(4)
subplot(4,1,1);
psectime(glider33, 't', 'sigma',':',[],20:0.5:30,20:30)

subplot(4,1,2);
%psectime(glider37, 't', 'sigma',':',[],20:0.5:30,20:30)
%throwing error. ask dan.

subplot(4,1,3);
psectime(glider45, 't', 'sigma',':',[],20:0.5:30,20:30)

subplot(4,1,4);
psectime(glider56, 't', 'sigma',':',[],20:0.5:30,20:30)

%%
figure
pxyll_xyuv_turner(map,'pvvert','sigma','obs',5,1,[],[],[],4)
pxyll_xyuv_turner(map,'pvvert','sigma','obs',5,2,[],[],[],1)
pxyll_xyuv_turner(map,'pvvert','sigma','obs',5,3,[],[],[],1)

%% movie making.
movietimell1cuv(map,'t',ualong,23,[],[],ctd,[])
movietimell1cuv(map,'t','ualong',23,[],[],ctd,[])
movietimell1cuv(mapem,'t','obs',11,[],1,ctd,6)


%% looking at lat/lon, color by time.

%sux
figure(5);
subplot(2,1,1)
scatter(ctd.lat, ctd.time)
subplot(2,1,2)
scatter(ctd.lon, ctd.time)
%%
figure(6);
scatter( ctd.lon, ctd.lat,4,ctd.time, 'filled');
cb = colorbar;
ylabel('latitude')
xlabel('longitude')
cb.Label.String = 'time (dn)';

%%
figure(7)
scatter(glider37.lon, glider37.lat, [], glider37.time, 'filled')


%%
figure(8);
osu = true;%true;
uw = true;%true;
spray = true;

% light green, spring, kelly, sea green
greens = ["#B0FF93" "#1BE77D" "#00AD00" "#0D824B"];
% maize, naples yellow, mikado, metallic gold
yellows = ["#FFEC70", "#F7DB4D", "#FFC800", "#CCA702"];
% sandy brown, orange peel, web dark orange, pumpkin
oranges = ["#FFB264", "#FFA200", "#FF8800", "#FF6A00"];

figure;
hold on;

plot_handles = []; legend_labels = {};

if osu
    % OSU SLOCUM GLIDERS: GREEN
    scatter(slocum684.lon, slocum684.lat, [], slocum684.time, 'filled', 'HandleVisibility', 'off');
    p1 = plot(slocum684.lon, slocum684.lat, 'Color', greens(1), 'LineWidth', 2, 'MarkerIndices', 1:15:length(slocum684.lat));
    scatter(slocum685.lon, slocum685.lat, [], slocum685.time, 'filled', 'HandleVisibility', 'off');
    p2 = plot(slocum685.lon, slocum685.lat, 'Color', greens(2), 'LineWidth', 2, 'MarkerIndices', 1:15:length(slocum685.lat));
    scatter(slocum686.lon, slocum686.lat, [], slocum686.time, 'filled', 'HandleVisibility', 'off');
    p3 = plot(slocum686.lon, slocum686.lat, 'Color', greens(3), 'LineWidth', 2, 'MarkerIndices', 1:15:length(slocum686.lat));
    plot_handles = [plot_handles, p1, p2, p3];
    legend_labels = [legend_labels, {'osu1', 'osu2', 'osu3'}];
end

if uw
    % UW SEAGLIDERS: YELLOW
    scatter(sglider141.lon, sglider141.lat, [], sglider141.time, 'filled', 'HandleVisibility', 'off');
    p4 = plot(sglider141.lon, sglider141.lat, 'Color', yellows(1), 'LineWidth', 2, 'MarkerIndices', 1:15:length(sglider141.lat));
    scatter(sglider526.lon, sglider526.lat, [], sglider526.time, 'filled', 'HandleVisibility', 'off');
    p5 = plot(sglider526.lon, sglider526.lat, 'Color', yellows(2), 'LineWidth', 2, 'MarkerIndices', 1:15:length(sglider526.lat));
    scatter(sglider528.lon, sglider528.lat, [], sglider528.time, 'filled', 'HandleVisibility', 'off');
    p6 = plot(sglider528.lon, sglider528.lat, 'Color', yellows(3), 'LineWidth', 2, 'MarkerIndices', 1:15:length(sglider528.lat));
    scatter(sglider687.lon, sglider687.lat, [], sglider687.time, 'filled', 'HandleVisibility', 'off');
    p7 = plot(sglider687.lon, sglider687.lat, 'Color', yellows(4), 'LineWidth', 2, 'MarkerIndices', 1:15:length(sglider687.lat));
    plot_handles = [plot_handles, p4, p5, p6, p7];
    legend_labels = [legend_labels, {'sg1', 'sg2', 'sg3', 'sg4'}];
end

if spray
    % SPRAYS: ORANGE
    scatter(glider33.lon, glider33.lat, [], glider33.time, 'filled', 'HandleVisibility', 'off');
    p8 = plot(glider33.lon, glider33.lat, 'Color', oranges(1), 'LineWidth', 2, 'MarkerIndices', 1:15:length(glider33.lat));
    scatter(glider37.lon, glider37.lat, [], glider37.time, 'filled', 'HandleVisibility', 'off');
    p9 = plot(glider37.lon, glider37.lat, 'Color', oranges(2), 'LineWidth', 2, 'MarkerIndices', 1:15:length(glider37.lat));
    scatter(glider45.lon, glider45.lat, [], glider45.time, 'filled', 'HandleVisibility', 'off');
    p10 = plot(glider45.lon, glider45.lat, 'Color', oranges(3), 'LineWidth', 2, 'MarkerIndices', 1:15:length(glider45.lat));
    scatter(glider56.lon, glider56.lat, [], glider56.time, 'filled', 'HandleVisibility', 'off');
    p11 = plot(glider56.lon, glider56.lat, 'Color', oranges(4), 'LineWidth', 2, 'MarkerIndices', 1:15:length(glider56.lat));
    plot_handles = [plot_handles, p8, p9, p10, p11];
    legend_labels = [legend_labels, {'spray1', 'spray2', 'spray3', 'spray4'}];
end

legend(plot_handles, legend_labels);

hold off; 
xlabel('Longitude');
ylabel('Latitude');
title('Sprays, Seagliders, and Slocums');

cb = colorbar;
ylabel(cb, 'Date and Time');
datetick(cb, 'y', 'dd-mmm', 'keeplimits');
numTicks = 10;  % adjust if desired
set(cb, 'YTick', linspace(min(glider37.time), max(glider37.time), numTicks));
cbTicks = get(cb, 'YTick');
cbTickLabels = cellstr(datestr(datetime(cbTicks, 'ConvertFrom', 'posixtime'), 'dd-mmm'));
set(cb, 'YTickLabel', cbTickLabels);


%% 3d scatter

figure(10);
scatter3(glider33.lon, glider33.lat, glider33.depth); hold on; % doesnt like depth being only 100x1
plot3(glider33.lon, glider33.lat, glider33.depth, 'Color', oranges(1), 'LineWidth', 2);

view(45, 30);  %hmmm

%% more psection
% Plot section with x=distance, y=depth, color=var1, contour=var2.

figure(11);
subplot(4,1,1); 
psection(glider33, 's', 'sigma', ':', [], 20:0.5:30, 20:30, 'dist')

subplot(4,1,2)
psection(glider37, 's', 'sigma', ':', [], 20:0.5:30, 20:30, 'dist')

subplot(4,1,3)
psection(glider45, 's', 'sigma', ':', [], 20:0.5:30, 20:30, 'dist')

subplot(4,1,4)
psection(glider56, 's', 'sigma', ':', [], 20:0.5:30, 20:30, 'dist')



% Plot section with x=distance, y=depth, color=var1, contour=var2.

figure(12);
subplot(4,1,1); 
psection(glider33, 'sigma', 't', ':', [],3:4:32, [], 'dist')

subplot(4,1,2)
psection(glider37, 'sigma', 't',  ':', [], 20:0.5:30, 20:30, 'dist')

subplot(4,1,3)
psection(glider45, 'sigma', 't',  ':', [], 20:0.5:30, 20:30, 'dist')

subplot(4,1,4)
psection(glider56, 'sigma', 't', ':', [], 20:0.5:30, 20:30, 'dist')



