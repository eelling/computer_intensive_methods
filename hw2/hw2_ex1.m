%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%           Comp. intens. HW2 EX1         %%%%%
%%%%%               6 may 2018                %%%%%
%%%%%       Eva Elling Linn Öhström           %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hybrid MCMC (Gibbs with MH-steps)algorithm sampeling 
% from the given posterior
clc; close all;

% load data
load coal-mine.csv

% given parameters
t1 = 1851; % fixed end points of the dataset in years
tD = 1963; % D = d+1
v = 0.1; % hyper parameter, to be changed
bp = 4; % # of breakpoints,, to be changed
d = bp+1; % time points
rho = 0.05; % tuning parameter of the proposal distributions
N = 100000; % %steps
burn_in = 10000; % how much should be discarded?
M = 10000 + burn_in;
t = zeros(d+1,M); % the breakpoints, d+1 to have d intervals
t(:,1) = t1;
t(:,d+1) = tD;
t(:,1) = linspace(t1,tD,d+1); % end- and breakpoints
lam = zeros(d,M); % intensities lambda for each interval of t (thus d-dim.)
n = zeros(d,M); % # of disasters for each interval of t
theta = zeros(d,M);

% implementation Gibbs
for j = 1:M-1
    % # number of disasters in the sub-interval
    for i = d
        n(i,j) = sum(coal_mine >= t(i,j) & coal_mine < t(i+1,j));
    end
    % draw theta
    for i = d
        theta(:,j) = gamrnd(2*d+2, 1./(v + sum(lam(:,j)))); %why theta first? 
    end
    % draw lambda
    for i = 1:d
       lam(:,j) = gamrnd(n(:,j)+2, 1./(diff(t(:,j))+theta(:,j)));
    end
    
    % rwmh for t 
    t = RWMH(t,lam,j,n,d,rho,coal_mine);
end


figure(1)
for i = 1:d
    %subplot(d,1,i)
    %plot(theta(i,burn_in:end));
    histfit(theta(i,burn_in:end),20,'gamma'); % ska man få olika theta?
    title('\theta distribution')
end

figure(2)
 for i = 1:bp
    subplot(bp,1,i)
    hold on
    %histogram(t(i,burn_in:end),20)
    %plot(lam(i,burn_in:end))
    ksdensity(lam(i,burn_in:end))
    title('\lambda distribution')
 end 

figure(3)
 for i = 1:d
    subplot(d,1,i)
    hold on
    %histogram(t(i,burn_in:end),50)
    %plot(t(i,burn_in:end))
    ksdensity(t(i,burn_in:end))
    title('t distribution')
 end 

% rwmh implementation
function t_prop = RWMH(t,lam,j,n,d,rho,coal_mine)
    for i = 2:d
        R = rho*(t(i+1,j)-t(i-1,j)); 
        t_prop =  t(i,j) + floor((2*R)*rand) - R;
        t_star = t(:,j);
        t_star(i) = t_prop;
        if t_prop > t(i-1,j) && t_prop < t(i+1,j)
            for i = 1:d
                n_star(i,1) = sum(coal_mine >= t_star(i) & ...
                    coal_mine <= t_star(i+1));
            end
            % # of disasters year i are assumed to be Po
            a = prod(lam(:,j).^n_star + exp(-lam(:,j).*diff(t_star)).*diff(t_star)) / ...
                prod(lam(:,j).^n(:,j) + exp(-lam(:,j).*diff(t(:,j))).*diff(t(:,j)));
            if rand <= a
                t(i,j+1) = t_prop;
            else
                t(i,j+1) = t(i,j);
            end            
        else
            t(i,j+1) = t(i,j);
        end
    end
t_prop = t; 
end
