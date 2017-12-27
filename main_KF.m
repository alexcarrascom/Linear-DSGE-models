%% Main: Kalman Filter
% this driver is used to test some of the functions that I am creating.
% By Alex Carrasco, november 2017

clc; clear all; close all;
addpath(genpath([pwd '\Estimation']));

path_g = [ pwd '\Graphs'];
% imprime = @(x) print( gcf, '-depsc2', [path_g filesep x]);
% imprpdf = @(x) eps2pdf( [path_g filesep x '.eps']);
formataxis  = @(x) set(x, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 15, 'Box', 'Off', 'PlotBoxAspectRatio', [1 0.75 1]);
formatlegend  = @(x) set(legend, 'Location', x, 'Orientation', 'Vertical', 'Box', 'Off', 'Fontsize', 18, 'Fontangle', 'normal');

%% [I] Kalman filter
% loading data
load KFexample.mat;
lambda = 1600;

% cycle I: using GDP 
Psi = [1 0]; Su = lambda/16;
A1  = [2 -1;1 0]; Ae = [1;0]; Se = Su*lambda^(-1);

z0=[thp(1) thp(1)]; % initial
P0=diag([1e+10,1e+10]);       % initial  

[KF1,LF1,Lt1]=KalmanFilter(y,A1,Ae,Psi,Se,Su,'init_P',P0,'init_z',z0);
tkf1 = KF1.smoothZ.z(1,:)';
tkf1_ns = KF1.filterZ.z(1,2:end)';
ckf1 = y'-tkf1;
P1_t = vec(KF1.smoothZ.P(1,1,:));

% cycle II: using annualized change
Psi = [1 4 -4]; Su = 0;
A1 = [1 0 0; 0 0 0; 0 1 0];
Ae = [1 0; 0 1; 0 0]; Se = diag([1 lambda/16]);


z0 = [Dthp(2) chp(2) chp(1)]';
P0 = diag([1e+10 lambda/16 lambda/16]);

[KF2,LF2,Lt2]=KalmanFilter(Dy',A1,Ae,Psi,Se,Su,'init_P',P0,'init_z',z0);
tkf2    = KF2.smoothZ.z(1,:)';
tkf2_ns = KF2.filterZ.z(1,2:end)';
ckf2 = KF2.smoothZ.z(2,:)';
ckf2_ns = KF2.filterZ.z(2,2:end)';
P2_t = vec(KF2.smoothZ.P(1,1,:));
P2_c = vec(KF2.smoothZ.P(2,2,:));

%% [II] Plotting
figure(1)
h=plot(ran,[chp(2:end),ckf1(2:end),ckf2]);
set(h(1),'color',[.6 .6 .6],'linestyle','-','linewidth',6);
set(h(2),'color',[.2 .2 .9],'linestyle','-','linewidth',2.5);
set(h(3),'color','k','linestyle',':','linewidth',1.25);
aux=gca; 
set(aux,'xtick',ranx,'xticklabel',ranx_lab);
formataxis(aux);
legend('IrisToolBox','KF - ln(Y_t)','KF - 4\Delta ln(Y_t)');
formatlegend('NorthEast');
xlim( [min(ran)-1 max(ran)+1]), ylim([-18 18]);
% imprime('HP_compar');
% imprpdf('HP_compar');

figure(2)
hold on;
a = area( ran, [(ckf2-P2_c) 2*P2_c ] );
h = plot( ran, [ckf2 ckf2_ns]);
hold off;
set(a(1), 'EdgeColor', 'none', 'FaceColor', 'none');
set(a(2), 'EdgeColor', 'none', 'FaceColor', [1 0.65 0.65]);
set(h(1), 'Linewidth', 2.5, 'Color', [0.75 0.1 0.1]);
set(h(2), 'Linewidth', 1.5,  'Linestyle', '--', 'Color', 'k');
aux=gca; 
set(aux,'xtick',ranx,'xticklabel',ranx_lab);
formataxis(aux);
legend([h;a(2)],'smoothed cycle','filtered cycle','bands (\pm \sigma)');
xlim( [min(ran)-1 max(ran)+1]), ylim([-25 25]);
formatlegend('best');
% imprime('cycle_HPKF');
% imprpdf('cycle_HPKF');

figure(3)
hold on;
a = area( ran, [(tkf1(2:end)-P1_t(2:end)) 2*P1_t(2:end) ] );
h = plot( ran, [tkf1(2:end) tkf1_ns(2:end)]);
hold off;
set(a(1), 'EdgeColor', 'none', 'FaceColor', 'none');
set(a(2), 'EdgeColor', 'none', 'FaceColor', [0.65 0.65 1]);
set(h(1), 'Linewidth', 2.5, 'Color', [0.1 0.1 0.75]);
set(h(2), 'Linewidth', 1.5,  'Linestyle', '--', 'Color', 'k');
aux=gca; 
set(aux,'xtick',ranx,'xticklabel',ranx_lab);
formataxis(aux);
legend([h;a(2)],'smoothed trend','filtered trend','bands (\pm \sigma)'); 
xlim( [min(ran)-1 max(ran)+1]), ylim([1050 1200]);
formatlegend('best');
% imprime('trend_HPKF');
% imprpdf('trend_HPKF');

figure(4)
hold on;
a = area( ran, [(tkf2-P2_t) 2*P2_t] );
h = plot( ran, [tkf2 tkf2_ns]);
hold off;
set(a(1), 'EdgeColor', 'none', 'FaceColor', 'none');
set(a(2), 'EdgeColor', 'none', 'FaceColor', [0.65 0.65 1]);
set(h(1), 'Linewidth', 2.5, 'Color', [0.1 0.1 0.75]);
set(h(2), 'Linewidth', 1.5,  'Linestyle', '--', 'Color', 'k');
aux=gca; 
set(aux,'xtick',ranx,'xticklabel',ranx_lab);
formataxis(aux);
legend([h;a(2)],'smoothed \Delta y_t','filtered trend \Delta y_t','bands (\pm \sigma)'); 
xlim( [min(ran)-1 max(ran)+1]), ylim([-10 15]);
formatlegend('NorthWest');
% imprime('DY_HPKF');
% imprpdf('DY_HPKF');

figure(5)
plot(ran,Lt1(2:end),'-b','linewidth',2);
aux=gca; 
set(aux,'xtick',ranx,'xticklabel',ranx_lab);
formataxis(aux); xlim( [min(ran)-1 max(ran)+1]);
% imprime('Like_KF1');
% imprpdf('Like_KF1');

figure(6)
plot(ran,Lt2,'-b','linewidth',2);
aux=gca; 
set(aux,'xtick',ranx,'xticklabel',ranx_lab);
formataxis(aux); xlim( [min(ran)-1 max(ran)+1]);
% imprime('Like_KF2');
% imprpdf('Like_KF2');

% End
% close all;