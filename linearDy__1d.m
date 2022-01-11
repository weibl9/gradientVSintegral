% Copyright 2020, All Rights Reserved
% Code by Baolei Wei
% For paper, "Parameter estimation for grey system models: 
%                 Gradient matching versus integral matching"
% by Baolei Wei

%%
clc
clear
close all

%% true model: dx/dt=ax x(t0)=xi
dyn = 1;
switch dyn 
    case 1
        a = 0.75; xi = 0.35;    % growth
    otherwise 
        a = -0.75; xi = 3.50;   % decay 
end
t0 = 0; tn = 5;                 % time interval 

syms x(t)
eqn = diff(x,t) == a*x;
cond = x(t0) == xi;
xt = dsolve(eqn,cond);

npar = 1;                    % number of parameters

%% main process
sobs = [21,51,101,251,501];  % set of sample sizes
snvr = [0.01,0.04,0.09];     % set of noise-variance ratios
nrep = 5000;                 % repetitions in each (obs,nvr) scenario

par_est = cell(length(snvr),length(sobs));

for invr = 1:length(snvr) 
    nvr = snvr(invr);

    for iobs = 1:length(sobs) 
        nobs = sobs(iobs);                  % sample size
        tobs = linspace(t0,tn,nobs)';       % observed time instants 
        xtru = eval(subs(xt,tobs));         % true time-series observations
        
        est_gm = nan(nrep,npar);
        est_im = nan(nrep,npar);
        ini_cnd = nan(nrep,1);
        
        parfor irep=1:nrep
            % data generation
            rng(irep)
            nois = sqrt(nvr)*std(xtru)*randn(size(xtru));
            xobs = xtru + nois;              % noisy observations 

            % gradient matching
            dt = diff(tobs);
            dx = diff(xobs)./dt;
            Xd = 0.5*xobs(1:end-1)+0.5*xobs(2:end);
            est_gm(irep,:) = regress(dx, Xd)';      % convert to row vector

            % integral mathcing
            cx = cumsum(xobs(1:end-1)+xobs(2:end)).*dt/2;
            Xc = [cx ones(nobs-1,1)];
            beta = regress(xobs(2:end), Xc)';       % convert to row vector
            est_im(irep,:) = beta(1:end-1);
            ini_cnd(irep,:) = beta(end);            
        end
        
        par_est{invr,iobs} = {nvr; nobs; est_gm; est_im; ini_cnd};        
    end
    
end

%% group boxplot for gradient and integral matching-based structure parameters
fnts = 13;

h = figure;
[nrow, ncol] = size(par_est);
for irow=1:nrow
    gm = []; im = []; xlab = [];
    for icol=1:ncol
        gm = horzcat(gm, par_est{irow,icol}{3});
        im = horzcat(im, par_est{irow,icol}{4});
        xlab = horzcat(xlab,par_est{irow,icol}{2}); 
    end
    colr = [0, 255, 0, 250; 
            255, 255, 0, 250]'/256;
    if irow == 1
        lgd = {'gradient matching','integral matching'}; 
    else 
        lgd = 0;
    end
    
    subplot(1,nrow,irow)
    gboxplot({gm,im}, xlab, colr, lgd, 'northwest')
    grid on; grid minor
    ylim([0.3 1.3])
    yline(a,':',"Color",'b','linewidth',1,"HandleVisibility","off");
    
    xlabel('$n$','interpreter','latex')
    ylabel('$\hat{\beta}$','interpreter','latex')
    title(['nvr = ',num2str(sqrt(par_est{irow,1}{1})*100),'%'])
    set(gca,'fontsize',fnts)
end
set(gcf,'Position',[50 200 1000 250])
set(h,'PaperSize',[15 10])
print(['../LaTexSourceFiles DeIn/figs/1d-linear-',num2str(dyn)],'-bestfit','-dpdf')

%% group boxplot for integral matching-based initial cinditions
h = figure;
subplot(1,nrow,1)
ini_nvrs = cell(1,ncol);
for icol=1:ncol
    ini_nvrs{icol} = horzcat(par_est{1,icol}{5}, par_est{2,icol}{5}, par_est{3,icol}{5});
end

colr = [255, 255, 0, 255;  
        255, 0, 255, 255; 
        0, 255, 255, 255;
        0, 255, 0, 255;                
        255, 163, 0, 255]'/256;
    
xlab = {'10%', '20%', '30%'}; % snvr
lgd = {'n=  21', 'n=  51', 'n=101', 'n=251', 'n=501'}; % cellstr(split(string(sobs))); 
    
gboxplot(ini_nvrs, xlab, colr, lgd, 'northwest')

grid on;
grid minor
ylim([-0.4 1.5])
yline(xi,':',"Color",'b','linewidth',1,"HandleVisibility","off");

set(gca,'fontsize',fnts)
xlabel('$nvr$','interpreter','latex')
ylabel('$\hat{\eta}$','interpreter','latex')
    
set(gcf,'Position',[50 200 1000 250])
set(h,'PaperSize',[15 10])
print(['../LaTexSourceFiles DeIn/figs/1d-linear-init',num2str(dyn)],'-bestfit','-dpdf')

%% tables of the estimates
[nrow, ncol] = size(par_est);
par_tab = cell(nrow,ncol);

for irow=1:nrow
    for icol=1:ncol
        irc_par = par_est{irow,icol};
        par_tab{irow,icol} = {irc_par{1}; ...
            irc_par{2}; ...
            [mean(irc_par{3}); std(irc_par{3})]; ...
            [mean(irc_par{4}); std(irc_par{4})]; ...
            [mean(irc_par{5}); std(irc_par{5})] };
    end
end

T = [];
for irow=1:nrow
    gm_im_init = [];
    for icol=1:ncol
        gm_im_init = vertcat(gm_im_init, ...
            [par_tab{irow,icol}{3}, par_tab{irow,icol}{4}, par_tab{irow,icol}{5}]);
    end
    T = vertcat(T, gm_im_init);
end
