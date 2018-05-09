
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%           Comp. intens. HW1 EX3         %%%%%
%%%%%             29 april 2018               %%%%%
%%%%%       Eva Elling Linn Öhström           %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; format short; %rng default

% data
load stations.mat % basis stations position
load RSSI-measurements.mat % signal

% given
v = 90;
s = 6;
eta = 3;
std = 1.5;
sigma = 0.5;
dt = 0.5;
d_t=dt;
alpha = 0.6;

phi_tilde = [1,dt,dt^2/2;0,1,dt;0,0,alpha];
Phi = [phi_tilde,zeros(3,3);zeros(3,3),phi_tilde];
psi_z = [dt^2/2,dt,0]';
psi_w = [dt^2/2,dt,1]';
Psi_z = [psi_z,zeros(3,1);zeros(3,1),psi_z];
Psi_w = [psi_w,zeros(3,1);zeros(3,1),psi_w];

% simulation
m = 500;
N = 10000;
tau = zeros(2,m);

W = @() mvnrnd(zeros(2,1),sigma*eye(2),N);
X0 =  mvnrnd(zeros(6,1),diag([500,5,5,200,5,5]),N)';

exp = @(X,l) v - 10*eta*log10(sqrt((X(1,:)-pos_vec(1,l)).^2+ (X(4,:)-pos_vec(2,l)).^2))';
p = @(Y,X,l) normpdf(Y,exp(X,l),std); % observation density, for weights
p_yx = @(Y,X) prod([p(Y(1),X,1),p(Y(2),X,2),p(Y(3),X,3),...
    p(Y(4),X,4),p(Y(5),X,5),p(Y(6),X,6)],2);

w0 = p_yx(Y(:,1),X0);
W_gen = @(r,c) mvnrnd(zeros(r,1),sigma*eye(r,r),c)';

tau(:,1) = [sum(X0(1,:)*w0)/sum(w0) sum(X0(4,:)*w0)/sum(w0)]';

X = X0;
w = w0;

state_vector = transpose(datasample([1 2 3 4 5]', N))';

for n = 2:m
    
    [Zn,state_vector] = Z_N(state_vector,N);
    X = Phi*X + Psi_z*Zn + Psi_z*W()';
    w = p_yx(Y(:,n),X).*w;
    tau(:,n) = [sum(X(1,:)*w)/sum(w) sum(X(4,:)*w)/sum(w)]';

end
X_state = tau;
figure, plot(X_state(1,:),X_state(2,:),'-bl'), hold on, plot(X_state(1,1),X_state(2,1),'xr'), hold on, grid on
plot(pos_vec(1,:),pos_vec(2,:),'bl*'), legend('Trajectory','Initial position (X_0^1,X_0^2)', 'Location of the basis stations')

function [Zn,state_vector] = Z_N(state_vector,N)
    P = (1/20)*(ones(5) + 15*eye(5)); %transition  probability matrix
    Z = [[0,0]', [3.5,0]', [0,3.5]', [0,-3.5]', [-3.5,0]'];
    pi = rand(1,N);
    pi_states = P(:,state_vector);
    trans = cumsum(pi_states);

    for i = 1:N
         state_vector(i) = find(trans(:,i)'>pi(i),1); 
    end 

    Zn = Z(:,state_vector);
end
