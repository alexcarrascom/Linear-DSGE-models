%% main.m
% Estimates a SOE models: the objective for this driver is to test the
% created functions.
% This file performs RWMH estimation using own codes.
% Some changes must be done, basically changing paths
% **************************************
% Alex Carrasco, december 2017
% **************************************
clear; clc; close all;

% This part needs to be changed updating folders
addpath(genpath([pwd '\Solution']));
addpath(genpath([pwd '\Estimation']));

load data.mat;

path_g = [pwd '\Graphs']; 
% imprime = @(x) print( gcf, '-depsc2', [path_g filesep x]);
% imprpdf = @(x) eps2pdf( [path_g filesep x '.eps']); % function must be add
formataxis  = @(x) set(x, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 20, 'Box', 'Off', 'PlotBoxAspectRatio', [1 0.75 1]);
formatlegend  = @(x) set(legend, 'Location', x, 'Orientation', 'Vertical', 'Box', 'Off', 'Fontsize', 20, 'Fontangle', 'normal');
label_x   =@(x) xlabel(x,'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 20);
label_y   =@(x) ylabel(x,'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 20);
titul   =@(x) title(x,'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 20);

est = 0; % estimate or load past results?

global Y X
%% [I] Data
data=rmfield(data,'i_f1');
obs_sort={'c','pi','i','Ds','pif','i_f2'};
vars = fieldnames(data);

Y=[];
s = 1:numel(obs_sort);
for ii=1:numel(obs_sort)
    ind=strcmp(obs_sort{ii},vars);
    Y(ii,:)=data.(vars{s(ind)});    
end
T=numel(Y(1,:));
X = zeros(1,T);

%% [II] Baseline calibration
clear theta;
% for initialization of csminwel I use Dynare mode computation
mod_pastcomputation = [pwd '\Mod\model_basic_mode.mat'];
load(mod_pastcomputation);
theta0 = xparam1;
Sigma0 = inv(hh);

%% [III] IRF's
[Gamma0,Gamma1,Const,Psie,Pi,Se]=sims_basic_model(theta0);
[A1,~,Ae,~,~,~,~,eu] = gensys(Gamma0,Gamma1,Const,Psie,Pi);
IR2=IRF(A1,Ae,zeros(7),'noprint','nograph','periods',20,'transform',@(x) 100*x);
[Psi2,Su,Psi1]=measure_eq_basic(theta0);
[~,ll,like]=KalmanFilter(Y,A1,Ae,Psi2,Se,Su,'nosmooth','exog',X,Psi1,'init_P',10*eye(size(A1,1)));

%% [II] Plotting data
cc=1;

for ii=1:numel(obs_sort)
    figure(cc); hold on;
    h=plot(rngQ,Y(ii,1:end)'); 
    plot([rngQ(1) rngQ(end)],[0 0],':k');
    set(h,'linewidth',2,'color',[0.1 0.1 0.6]); axis tight;
    set(gca,'xtick',rng_plot,'xticklabel',rng_plotlab);
    formataxis(gca);
%     imprime(['fig_' num2str(cc) ]);
%     imprpdf(['fig_' num2str(cc) ]);
    cc=cc+1;
end

%% [III] Estimation
if est    
tic
    c = .125; % jump scale for sampling
    N = 5000;   % 100000
    [THETA,thetamode,Sigma,~,aratio,ll] = RWMH_dsge(Y,X,c,theta0,Sigma0, ...
                        @priors, @measure_eq_basic,@sims_basic_model,'n_replic',N, ...
                     'n_drop',.75,'init_scal',0,'mode_file','mode_comp'); 
toc
    save BE_result THETA aratio thetamode
else
    load BE_result.mat
end

%% [IV] Plots  posterior distribution
load([pwd '\Mod\model_basic_results.mat']);
density = oo_.posterior_density.parameters;
par_names=fieldnames(density);

for ii=1:26  % 26 estimated parameters
    figure(cc); hold on;
    [temp,~,x]=mykernel(THETA(:,ii),'kernel',2,'ngrids',300,'eachseries');
    h1=plot(x,temp);
    set(h1,'color',[0.7 0.1 0.1],'linestyle','none','marker','x','markersize',7);
    temp=density.(par_names{ii});
    h2=plot(temp(:,1),temp(:,2));
    set(h2,'color',[0.1 0.1 0.7],'linestyle','-','linewidth',2);
    axis tight;
    legend('Own code','Dynare');
    formatlegend('best');
    formataxis(gca);    
%     imprime(['fig_' par_names{ii}]);
%     imprpdf(['fig_' par_names{ii}]);
    cc=cc+1;
end




