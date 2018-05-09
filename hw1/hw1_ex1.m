disp('           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('           %%%%%           Comp. intens. HW1 EX1         %%%%%')
disp('           %%%%%             16 nov 2017                 %%%%%')
disp('           %%%%%       Eva Elling Linn Öhström           %%%%%')
disp('           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')
clc; clear all; format short; %rng default

% given
sigma = 0.5;
dt = 0.5;
alpha = 0.6;
m = 800;

P = (1/20)*(ones(5) + 15*eye(5)); %transition  probability matrix 
W = @() mvnrnd(zeros(2,1),sigma^2*eye(2));
Z = [[0,0]', [3.5,0]', [0,3.5]', [0,-3.5]', [-3.5,0]'];

phi_tilde = [1,dt,dt^2/2;0,1,dt;0,0,alpha];
Phi = [phi_tilde,zeros(3,3);zeros(3,3),phi_tilde];
psi_z = [dt^2/2,dt,0]';
psi_w = [dt^2/2,dt,1]';
Psi_z = [psi_z,zeros(3,1);zeros(3,1),psi_z];
Psi_w = [psi_w,zeros(3,1);zeros(3,1),psi_w];

X0 =  mvnrnd(zeros(6,1),diag([500,5,5,200,5,5]))';
X_state = zeros(2,m);
X = X0;

state0 = randi([1,5]);
STATE = state0
state = num2str(state0);

for n = 1:m
    switch state
        case '1' 
            pi = rand;
            pi_state = P(STATE,:);
            trans = cumsum(pi_state);
            new_state = find(trans > pi,1);
            X = Phi*X + Psi_z*Z(new_state) + Psi_z*W()';
            X_state(:,n) = [X(1) X(4)];
            STATE = new_state
            state = num2str(new_state);

        case '2'
            pi = rand;
            pi_state = P(STATE,:);
            trans = cumsum(pi_state);
            new_state = find(trans > pi,1);
            X = Phi*X + Psi_z*Z(new_state) + Psi_z*W()';
            X_state(:,n) = [X(1) X(4)]; 
            STATE = new_state
            state = num2str(new_state);
            
        case '3'
            pi = rand;
            pi_state = P(STATE,:);
            trans = cumsum(pi_state);
            new_state = find(trans > pi,1);
            X = Phi*X + Psi_z*Z(new_state) + Psi_z*W()';
            X_state(:,n) = [X(1) X(4)]; 
            STATE = new_state
            state = num2str(new_state);
            
        case '4'
            pi = rand;
            pi_state = P(STATE,:);
            trans = cumsum(pi_state);
            new_state = find(trans > pi,1);
            X = Phi*X + Psi_z*Z(new_state) + Psi_z*W()';
            X_state(:,n) = [X(1) X(4)];
            STATE = new_state
            state = num2str(new_state);
            
        case '5'
            pi = rand;
            pi_state = P(STATE,:);
            trans = cumsum(pi_state);
            new_state = find(trans > pi,1);
            X = Phi*X + Psi_z*Z(new_state) + Psi_z*W()';
            X_state(:,n) = [X(1) X(4)];
            STATE = new_state
            state = num2str(new_state);
    end
        
end
X_state = X_state';   
space = 8;
plot(X_state(1:space:end,1),X_state(1:space:end,2),'--bl'), xlabel('X1'), ...
    ylabel('X2'), hold on, plot(X_state(1,1),X_state(1,2),'r*'), grid on, ...
    legend('Trajectory path','Initial position (X_0^1,X_0^2)')
%%

