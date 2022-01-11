% Copyright 2020, All Rights Reserved
% Code by Baolei Wei
% For paper, "Parameter estimation for grey system models: 
%                 Gradient matching versus integral matching"
% by Baolei Wei

%%
clc
clear
close all

%% true model
t0 = 0; tn = 10;                    % time interval
syms x1(t) x2(t)                    % sate variabels
x = [x1; x2];                       % state vector
A = [ -0.5  -1; 
       1.2  -0.6];                  % structural parameter matrix
xi = [8.0; 7.5];
cond = x(t0) == xi;         % initial condition vector 

odes = diff(x) == A*x;              % ODEs in Matrix Form
[x1Sol(t), x2Sol(t)] = dsolve(odes, cond); % xSol(t) = simplify(xSol(t))  

npar = 4;                    % number of parameters

%% main process
sobs = [21,51,101,251,501];  % set of sample sizes
snvr = [0.01,0.04,0.09];     % set of noise-variance ratios
nrep = 500;                 % repetitions in each (obs,nvr) scenario

par_est = cell(length(snvr),length(sobs));

for invr = 1:length(snvr) 
    nvr = snvr(invr);

    for iobs = 1:length(sobs) 
        nobs = sobs(iobs);                  % sample size
        tobs = linspace(t0,tn,nobs)';       % observed time instants
        xtru = [eval(subs(x1Sol,tobs)), eval(subs(x2Sol,tobs))];	% true time-series observations 
        
        est_gm = nan(nrep,npar);
        est_im = nan(nrep,npar);
        ini_cnd = nan(nrep,2);
        
        parfor irep=1:nrep
            % data generation
            rng(irep)
            nois = sqrt(nvr)*std(xtru).*randn(size(xtru));
            xobs = xtru + nois;              % noisy observations 
            Xall = xobs;                     % all features in vector field: x1 x2
            
            % gradient matching
            dt = diff(tobs);
            dx = diff(xobs)./dt;            
            Xdll = Xall(1:end-1,:)/2+Xall(2:end,:)/2;
            est_gm(irep,:) = [ regress(dx(:,1), Xdll)',...     % A(1,1) A(1,2)
                               regress(dx(:,2), Xdll)' ];      % A(2,1) A(2,2)

            % integral mathcing
            cXall = cumsum(Xall(1:end-1,:)+Xall(2:end,:)).*dt/2;        % x1 x2 x1*x2
            Xcll = [cXall ones(nobs-1,1)];
            beta = [ regress(xobs(2:end,1), Xcll)', ...      % A(1,1) A(1,2) \eta_1
                     regress(xobs(2:end,2), Xcll)' ];        % A(2,1) A(2,2) \eta_2
            est_im(irep,:) = beta([1,2,4,5]);       % A(1,1) A(1,2) A(2,1) A(2,2)
            ini_cnd(irep,:) = beta([3,6]);          % \eta_1 \eta_2          
        end
        
        par_est{invr,iobs} = {nvr; nobs; est_gm; est_im; ini_cnd};        
    end
    
end

%% group boxplot for gradient and integral matching-based structure parameters
fnts = 15;
A2r = reshape(A',1,[]);

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
        yline(A2r(ipar),':',"Color",'b','linewidth',1,"HandleVisibility","off");

        set(gca,'fontsize',fnts)
        xlabel('$n$','interpreter','latex')
        switch ipar
            case 1
                ylabel('$\hat{\beta}_{1,1}$','interpreter','latex')
                ylim([-1.2 0.3])
            case 2
                ylabel('$\hat{\beta}_{1,2}$','interpreter','latex')
                ylim([-1.5 0.0])
            case 3
                ylabel('$\hat{\beta}_{2,1}$','interpreter','latex')
                ylim([0.3 2.7])
            case 4
                ylabel('$\hat{\beta}_{2,2}$','interpreter','latex')
                ylim([-1.3 0.3])
        end
        title(['nvr = ',num2str(sqrt(par_est{irow,1}{1})*100),'%'])
    end
    set(gcf,'Position',[50 200 1050 450])
    set(h,'PaperSize',[15 10])
end 
%% group boxplot for integral matching-based initial cinditions
for icmp = 1:2
    h = figure;
    subplot(1,2,1)
    ini_nvrs = cell(1,ncol);
    for icol=1:ncol
        ini_nvrs{icol} = horzcat(par_est{1,icol}{5}(:,icmp), par_est{2,icol}{5}(:,icmp), par_est{3,icol}{5}(:,icmp));
    end

    colr = [255, 255, 0, 255;  
            255, 0, 255, 255; 
            0, 255, 255, 255;
            0, 255, 0, 255;                
            255, 163, 0, 255]'/256;

    xlab = {'10%', '20%', '30%'}; % snvr
    lgd = {'n=  21', 'n=  51', 'n=101', 'n=251', 'n=501'}; % cellstr(split(string(sobs)));

    gboxplot(ini_nvrs, xlab, colr, lgd, 'southwest')
    yline(xi(icmp),':',"Color",'b','linewidth',1,"HandleVisibility","off");

    set(gca,'fontsize',fnts)
    xlabel('$nvr$','interpreter','latex')
    ylabel('$\hat{\eta}$','interpreter','latex')    
    switch icmp
        case 1
            ylabel('$\hat{\eta}_{1}$','interpreter','latex')
        case 2
            ylabel('$\hat{\eta}_{2}$','interpreter','latex')
    end
    set(gcf,'Position',[50 200 1050 450])
    set(h,'PaperSize',[15 10])
end


%% tables of the estimates
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



