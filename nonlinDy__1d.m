% Copyright 2020, All Rights Reserved
% Code by Baolei Wei
% For paper, "Parameter estimation for grey system models: 
%                 Gradient matching versus integral matching"
% by Baolei Wei

%%
clc
clear
close all

%% true model: dx/dt=ax+bx^2 x(t0)=xi
dyn = 1;
switch dyn 
    case 1
        r = 2.0;
        a = 1.0; b = -1.1; xi = 0.1;   % growth
    otherwise 
        r = 2.0;
        a = 1.0; b = -0.5; xi = 0.2;   % decay 
end
t0 = 0; tn = 10; % time interval 

syms x(t)
eqn = diff(x,t) == a*x+b*x^r;
cond = x(t0) == xi;
xt = dsolve(eqn,cond);

npar = 2;                     % number of parameters

%% main process
sobs = [21,51,101,251,501];   % set of sample sizes
snvr = [0.01,0.04,0.09];      % set of noise-variance ratios
nrep = 5000;                  % repetitions in each (obs,nvr) scenario

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
            Xd = [ xobs(1:end-1)+xobs(2:end),...
                   xobs(1:end-1).^r+xobs(2:end).^r ]/2;
            est_gm(irep,:) = regress(dx, Xd)';      % convert to row vector

            % integral mathcing
            cx = [ cumsum(xobs(1:end-1)+xobs(2:end)).*dt,...
                   cumsum(xobs(1:end-1).^r+xobs(2:end).^r).*dt ]/2;
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

for ipar = 1:npar
    h = figure(ipar);
    [nrow, ncol] = size(par_est);
    for irow=1:nrow
        gm = []; im = []; xlab = [];
        for icol=1:ncol
            gm = horzcat(gm, par_est{irow,icol}{3}(:,ipar));
            im = horzcat(im, par_est{irow,icol}{4}(:,ipar));
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

        set(gca,'fontsize',fnts)
        xlabel('$n$','interpreter','latex')
        switch ipar
            case 1
                yline(a,':',"Color",'b','linewidth',1,"HandleVisibility","off");
                ylabel('$\hat{\beta}_1$','interpreter','latex')
                ylim([0 2])
            case 2
                yline(b,':',"Color",'b','linewidth',1,"HandleVisibility","off");
                ylabel('$\hat{\beta}_2$','interpreter','latex')
                ylim([-2.2 0])
        end
        title(['nvr = ',num2str(sqrt(par_est{irow,1}{1})*100),'%'])
    end
    set(gcf,'Position',[50 200 1050 250])
    set(h,'PaperSize',[15 10])
    print(['../LaTexSourceFiles DeIn/figs/1d-nonlinear-',num2str(ipar)],'-bestfit','-dpdf')
end

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
yline(xi,':',"Color",'b','linewidth',1,"HandleVisibility","off");
grid on; grid minor

set(gca,'fontsize',fnts)
xlabel('$nvr$','interpreter','latex')
ylabel('$\hat{\eta}$','interpreter','latex')
    
set(gcf,'Position',[50 200 1000 250])
set(h,'PaperSize',[15 10])
print('../LaTexSourceFiles DeIn/figs/1d-nonlinear-init','-bestfit','-dpdf')


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



