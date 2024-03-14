% kiet vu
%contain:
%differential equations 1 2 3 4 5 6
%Laplace transform 1 2
%Probability 1 -> 18

%% differential equations
% ex1
clear
close all
R=1;
L=3;
U=5;
%R*q'(t)+1/C*q(t)=U, q(0)=0
%first degree equation
%a*y'(t)+b*y(t)=A,y(0)=y0
%solution_formulas.pdf, 1b i) p.2
a=L;
b=R;
A=U;
y0=0;
la=b/a;
%solution
i=@(t) A/b+(y0-A/b)*exp(-la*t);
%i(t)=q'(t)
d=@(t) (y0-A/b)*exp(-la*t)*(-la);

tau=1/la; %time constant
tmax=5*tau;
t=0:tau/100:tmax;
plot(t,R*i(t),'linewidth',1.5)
hold on
plot(t,L*d(t),'linewidth',1.5)
hold off
grid
xlabel('time t')
xlim([0,tmax])
legend({'UR','UL)'},'fontsize',12)

%% ex2
clear
close all

a=1;
k = 50;
V = 1000;
b=k/V;
C = 2;
D = 2;
w = 0.2;
A=k*C;
B=k*D;
y0=0;
la=b/a;
K=1/sqrt(b^2+(a*w)^2);
phi=atan2(-a*w,b);
E=y0-K*A*sin(phi)-B/b;

q=@(t) K*A*sin(w*t+phi)+B/b+E*exp(-la*t);
q0 = K*A*sin(w*0+phi)+B/b+E*exp(-la*0);
%dq=@(t) K*A*cos(w*t+phi)*w+E*exp(-la*t)*(-la) %q'(t)=i(t)

%u=@(t) A*sin(w*t)+B

p=@(t) C*sin(w*t)+D;

tau=1/la; %time constant
tmax=10*tau;
T=2*pi/w; %period of the sine-curves
t=0:tau/100:tmax;

figure(1)
plot(t,q(t)/V,'r','linewidth',1.5)
hold on
plot(t,p(t),'b','linewidth',1.5)
hold off
grid
   xlabel('time t')
legend({'x(t)/V','p(t)'},'fontsize',12)
xlim([0,200])

% amplification and phase shift
w=0:0.1:10;
L=1./sqrt(b^2+(a*w).^2);
phi=atan2(-a*w,b);

figure(2)
subplot(2,1,1)
plot(w,L/20,'linewidth',1.5)
grid
title('amplification K')
subplot(2,1,2)
plot(w,phi,'linewidth',1.5)
grid
title('phase shift \phi')
xlabel('angular frequency \omega')
%% ex3
clear
close all

R = 1; 
L = 2; 
C = 0.5; 
U = 1; 
T = 25; 

t = linspace(0, T, 1000);

sys = @(t, x) [x(2); (U - R*x(2) - x(1)/C)/L];

% Set the initial conditions
x0 = [0; 0];

[t_sol, x_sol] = ode45(sys, t, x0);


q = x_sol(:, 1);
q_prime = x_sol(:, 2);
UR = R * q_prime;
UL = L * diff(q_prime)./diff(t_sol);
UC = q / C;

figure;
plot(t_sol, UR, 'r', 'LineWidth', 2); hold on;
plot(t_sol(1:end-1), UL, 'g', 'LineWidth', 2);
plot(t_sol, UC, 'b', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('U_R', 'U_L', 'U_C');
title('Voltages in the RLC Circuit');
grid on;
hold off;
%% ex4
clear
close all

M = 100; 
m = 10; 
k = 500000; 
b = 100; 
R = 0.1; 
omega = 100; % Angular speed (rad/s)

% Define the time vector
T = 1; % Time interval (Seconds)
t = linspace(0, T, 1000);

sys = @(t, x) [x(2); (m * R * omega^2 * sin(omega * t) - b * x(2) - k * x(1)) / M];

x0 = [0; 0];

[t_sol, x_sol] = ode45(sys, t, x0);

figure;
plot(t_sol, x_sol(:, 1), 'b', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('y(t)');
grid on;

% Calculate KA for the given range of omega
omega_range = linspace(0, 10 * sqrt(k / M), 1000);
KA = (m * R * omega_range.^2) ./ sqrt((k - M * omega_range.^2).^2 + (b * omega_range).^2);
% Plot KA(omega)
figure;
plot(omega_range, KA, 'r', 'LineWidth', 2);
xlabel('Angular Speed \omega (rad/s)');
ylabel('Amplitude KA');
title('Amplitude KA vs Angular Speed');
grid on;
%% ex5
clear
close all
R = 1;
L = 2;
C = 0.1;
U = 10;
T = 5; 
q0 = 0; 
i0 = 0; 

dt = 0.01; % time step

% Increasing part
t1 = 0:dt:T;
u1 = U / T * t1;

% Decreasing part
t2 = T + dt:dt:T - dt;
u2 = -U / T * (t2 - T) + U;

% One sawtooth
t = [t1, t2];
u = [u1, u2];

% Five sawtooths
t = [t, t + T, t + 2 * T, t + 3 * T, t + 4 * T];
u = [u, u, u, u, u];


N = length(t);
q = zeros(1, N);
i = zeros(1, N);
di = zeros(1, N - 1);
q(1) = q0;
i(1) = i0;

for n = 1:N - 1
    di(n) = (1 / L) * (u(n) - R * i(n) - (1 / C) * q(n)); % i'(t) = di(t)
    q(n + 1) = q(n) + i(n) * dt + 1 / 2 * di(n) * dt^2;
    i(n + 1) = i(n) + di(n) * dt; 
end


UR = R * i;
UL = L * [di, di(end)]; 
UC = q / C;


figure;
subplot(4, 1, 1);
plot(t, u, 'linewidth', 1.5);
grid;
title('Input Voltage u(t)');
grid on;

subplot(4, 1, 2);
plot(t, UC, 'linewidth', 1.5);
grid;
title('Charge UC');
grid on;

subplot(4, 1, 3);
plot(t, UR, 'linewidth', 1.5);
grid;
title('Voltage UR = R * q''(t)');
grid on;

subplot(4, 1, 4);
plot(t, UL, 'linewidth', 1.5);
grid;
title('Voltage UL = L * q''''(t)');
grid on;

%% ex6
clear
close all
J = 0.1; 
L = 0.5; 
R = 0.3;
b = 0.08; 
K = 0.1; 
Kp = 5;
Ki = 2;
Kd = 0.001;


num = K;
den = [J*L, R*J + b*L, b*R + K^2];
sys = tf(num, den);


pid_controller = pid(Kp, Ki, Kd);


cl_sys = feedback(pid_controller*sys, 1);


t = linspace(0, 10, 1000);

% Setpoint
r = ones(size(t));


[y, t_out] = lsim(cl_sys, r, t);
r_out = ones(size(t_out));

e = r_out - y;

u = Kp*e + Ki*cumsum(e)*(t_out(2)-t_out(1)) + Kd*[diff(e); 0]/(t_out(2)-t_out(1));

figure;
subplot(3, 1, 1);
plot(t_out, y);
xlabel('Time (s)');
ylabel('\omega(t)');
title('Angular Speed \omega(t)');
grid on;

subplot(3, 1, 2);
plot(t_out, e);
xlabel('Time (s)');
ylabel('e(t)');
title('Error e(t) = r(t) - \omega(t)');
grid on;

subplot(3, 1, 3);
plot(t_out, u);
xlabel('Time (s)');
ylabel('u(t)');
title('Output u(t) of the PID-controller');
grid on;

%% Probability
%%ex1
clear
close all
p=zeros(1,11); %p(k+1) = probability to have k=0,1,...,7 winning numbers
for k=0:10
   p(k+1)=nck(20,k)*nck(50,10-k)/nck(70,10); 
end

figure(1)
bar(0:10,p,1)
grid
xlabel('k')
ylabel('P(k winning numbers)')
% simulation 

N=100000; %number of rounds
%winning numbers
oikea_rivi=randperm(70,20);%random permutation, 
%7 numbers from 1,2,...,39, all different
k_oikein=zeros(1,11); %k_oikein(k+1) = number of rounds with  k=0,1,2,...,7 winning numbers 
for n=1:N
    oma_rivi=randperm(70,10); %chosen numbers
    k=length(intersect(oikea_rivi,oma_rivi)); %number of winning numbers
    %intersect(u,v)=common numbers in vectors u and v 
    k_oikein(k+1)=k_oikein(k+1)+1;
end


figure(2)
bar(0:10,k_oikein/N,1)
hold on
plot(0:10,p,'r.','markersize',20)
hold off
grid
xlabel('k')
ylabel('P(k winning numbers)')
legend('simu','exact')

%k oikein osuudet
k_oikein/N; %simulointi
p; %tarkka
%% ex2
clear
close all
N=100000; %number of simulations
case1=0; %4 times at lease 1 6 occurs
case2=0; %2 dice 24 times at lease 1 double 6 occurs


for n=1:N
   result1=randi(6,1,4); %results of the six throws
   if max(result1)== 6 %unique(result) = number of different numbers in vector result
      case1=case1+1;
   end
   result2=randi(6,2,24);
   if max(sum(result2))== 12
      case2=case2+1; 
   end
end

P1=(6^4-5^4)/6^4;
P2=(36^24-35^24)/((6^24)^2);

bar([1,2],[case1,case2]/N)
hold
plot([1,2],[P1,P2],'r.','markersize',12)
hold off
grid
legend('simu','exact','location','northwest')
%% ex3
clear
close all
N = 100000;
counter = 0;
    for i = 1:N
        numbers = sort(randperm(39, 7));
        consecutive = false;
        for j = 2:7
            if numbers(j) - numbers(j-1) == 1
                consecutive = true;
                break;
            end
        end
        if ~consecutive
            counter = counter + 1;
        end
    end
 
 counter/N;

%% ex4
clear 
close all
% Parameters
days = 365;
target_prob = 0.5;
num_simulations = 10000;

% Initialize variables
n = 0;
prob = 0;

% Find the minimum value of n for which the probability is greater than 0.5
while prob <= target_prob
    n = n + 1;
    count = 0;
    
    for i = 1:num_simulations
        % Run a single simulation
        birthdays = randi(days, [1, n]);
        unique_birthdays = unique(birthdays);
        found_all = (length(unique_birthdays) == days);

        if found_all
            count = count + 1;
        end
    end

    prob = count / num_simulations; 
end

fprintf('Minimum value of n for probability > 0.5: %d\n', n);

%% ex 5
clear
close all
n = 100;
N = 10000;
count = 0;
for i = 1:N
    birthdays = randi(365, 1, n);
    unique_birthdays = unique(birthdays);
    for j = 1:length(unique_birthdays)
        if sum(birthdays == unique_birthdays(j)) >= 3
            count = count + 1;
            break;
        end
    end
end
prob = count / N;

%% ex6a
clear
close all
N = 100000;
count = 0; 
M=[(1:13)'; (1:13)'; (1:13)' ;(1:13)'];

for i = 1:N 
    deck = randperm(52);
    SHm=M(deck);
    if (sum(SHm(1:13) == 1) == 4) 
        count = count + 1;
    end
    if (sum(SHm(14:26) == 1) == 4) 
        count = count + 1;
    end    
    if (sum(SHm(27:39) == 1) == 4) 
        count = count + 1;
    end    
    if (sum(SHm(40:52) == 1) == 4) 
        count = count + 1;       
    end
end

prob = count / N; 
%% ex6b
clear
close all
N = 100000;
count = 0; 
M=[(1:13)'; (1:13)'; (1:13)' ;(1:13)'];

for i = 1:N 
    deck = randperm(52);
    SHm=M(deck);
    if (sum(SHm(1:13) == 1) == 1) && (sum(SHm(14:26) == 1) == 1) && (sum(SHm(27:39) == 1) == 1) && (sum(SHm(40:52) == 1) == 1) 
        count = count + 1;       
    end
end
prob = count / N; 
%% EX7
clear
close all
n = 20; % number of elements
N = 10000; 

counts = zeros(n+1, 1); 

for i = 1:N
    order = randperm(n); 
    count = sum(order == (1:n)); 
    counts(count+1) = counts(count+1) + 1; 
end
probs = counts / N;
bar(0:n, probs, 1)
%% EX 8
clear
close all
n = 13;
k = 0:n; 
p = 1/3; % probability one game right
q = 2/3; % probability one game wrong

probs = zeros(1, n+1);
for i = 0:n
    probs(i+1) = nchoosek(n, i) * p^i * q^(n-i);
end

bar(0:n, probs, 1)

%% EX 9
clear
close all
n = 20; 
p = 0.5; 
N=100000;
heads = zeros(1, n+1);
k = 4;
ok=0;
for i = 1:N
    heads_row = 0;
    for j = 1:n
        if rand() < p
            heads_row = heads_row + 1;
            if heads_row==k
                ok=ok+1;
                break
            end
        else
            heads_row = 0;
        end
    end
end
ok/N;

%% ex10
clear
close all
n = 10; % coupons
t = 30; %boxes to be bought

trials = 10000;
count = 0; 

for i = 1:trials
    collected = zeros(1, n);
    for j = 1:t
        coupon = randi(n);
        collected(coupon) = 1;
        if sum(collected) == n
            count = count + 1;
            break
        end
    end
end


probability = count / trials;
%% ex11
clear
close all
p = 0.5;
N = 10;
m = 1;
n = 5;
num_simulations = 10000;
num_A_wins = 0;
    for i = 1:num_simulations
        A_wins = N - m;
        B_wins = N - n;
        while A_wins < N && B_wins < N
            if rand() < p
                A_wins = A_wins + 1;
            else
                B_wins = B_wins + 1;
            end
        end
        if A_wins == N
            num_A_wins = num_A_wins + 1;
        end
    end
probability = num_A_wins / num_simulations;

%% ex12
clear
close all
% Define variables
b1 = 5;
w1 = 7;
b2 = 3;
w2 = 4;
b3 = 4;
w3 = 6;
num_simulations = 100000;

num_black = 0;
for i = 1:num_simulations
    % Choose a box randomly
    box = randi(3);
    % Draw a ball randomly from the chosen box
    if box == 1
        if rand() < b1 / (b1 + w1)
            num_black = num_black + 1;
        end
    elseif box == 2
        if rand() < b2 / (b2 + w2)
            num_black = num_black + 1;
        end
    else
        if rand() < b3 / (b3 + w3)
            num_black = num_black + 1;
        end
    end
end
P_black = num_black / num_simulations;

fprintf("The probability of drawing a black ball is %.4f\n", P_black);
%% ex13
clear
close all
num_simulations = 10000;
num_wins = 0;
for i = 1:num_simulations
    % Roll two dice
    roll1 = randi(6);
    roll2 = randi(6);
    sum = roll1 + roll2;
    if sum == 7 || sum == 11
        % Player wins
        num_wins = num_wins + 1;
    elseif sum == 2 || sum == 3 || sum == 12
        % Player loses
    else
        first_roll = sum;
        while true
            roll1 = randi(6);
            roll2 = randi(6);
            sum = roll1 + roll2;
            if sum == first_roll
                % Player wins
                num_wins = num_wins + 1;
                break;
            elseif sum == 7
                % Player loses
                break;
            end
        end
    end
end
P = num_wins / num_simulations;

fprintf("The probability of winning is %.4f\n", P);
%% ex14
clear
close all
% Define variables
num_simulations = 10000;

A = 0;
for i = 1:num_simulations
    % Generate random arrival times for A and B
    t1 = rand() * 60;
    t2 = rand() * 60;
    % Calculate the absolute difference between the arrival times
    diff = abs(t1 - t2);
    if diff <= 10
        % A and B are in the cafe simultaneously
        A = A + 1;
    end
end
P_simulated = A / num_simulations;

fprintf("The probability of being in the cafe simultaneously is %.4f\n", P_simulated);
%% ex15
clear
close all
m = 2; 
M = 5; 
n = 10; 
num_samples = 100000; 


X = m + (M-m)*rand(num_samples*n,1);
X_bar = reshape(X,[n,num_samples]);
X_bar = mean(X_bar);


mu = (m+M)/2;
O = sqrt(1/12)*(M-m)/sqrt(n);

% Define the normal distribution
x = mu-3*O:0.001:mu+3*O;
y = normpdf(x,mu,O);

histogram(X_bar,100,'Normalization','pdf');
hold on;
plot(x,y,'r','LineWidth',2);

xlabel('x');
ylabel('pdf');
%% ex16
clear
close all
% Define variables
mu1 = 2; 
sigma1 = 1;                       
mu2 = 5;    
sigma2 = 1.414; 
a1 = 2/3; 
a2 = 1/3; 
num_samples = 100000; 


x1 = normrnd(mu1, sigma1, [1,num_samples]);
x2 = normrnd(mu2, sigma2, [1,num_samples]);


x = a1*x1 + a2*x2;


mu = a1*mu1 + a2*mu2;
sigma = sqrt(a1^2*sigma1^2 + a2^2*sigma2^2);


x1_dist = @(x) normpdf(x,mu1,sigma1);
x2_dist = @(x) normpdf(x,mu2,sigma2);
x_dist = @(x) normpdf(x,mu,sigma);

figure;
fplot(x_dist,'LineWidth',2);
hold on;
fplot(x1_dist,'LineWidth',2);
hold on;
fplot(x2_dist,'LineWidth',2);
hold off;
grid
xlim([-2,12])

figure;
histogram(x1,100,'Normalization','pdf');
hold on;
fplot(x1_dist,'LineWidth',2);
xlabel('x');
ylabel('pdf');
title(sprintf('Normal distribution N(%.2f, %.2f^2)',mu1,sigma1));

figure;
histogram(x2,100,'Normalization','pdf');
hold on; 
fplot(x2_dist,'LineWidth',2);
xlabel('x');
ylabel('pdf');
title(sprintf('Normal distribution N(%.2f, %.2f^2)',mu2,sigma2));

figure;
histogram(x,100,'Normalization','pdf');
hold on;
fplot(x_dist,'LineWidth',2);
xlabel('x');
ylabel('pdf');
title(sprintf('Normal distribution N(%.2f, %.2f^2)',mu,sigma));
%% ex17
clear
close all

% Define variables
mu = 20; % mean of the normal distribution
sigma = 10; % standard deviation of the normal distribution
n = 30; % sample size
num_samples = 100000; % number of samples

% Generate 100,000 samples of size n from a normal distribution
X = normrnd(mu, sigma, [n,num_samples]);

% Calculate the sample standard deviations
S = std(X);

% Calculate the chi-squared values
chi_sq = ((n-1)*S.^2)/sigma^2;

% Define the theoretical chi-squared distribution
df = n-1; % degrees of freedom
x = 0:0.1:df*4;
y = chi2pdf(x,df);

alpha = 0.05;
chi_sq_a = chi2inv(alpha/2,df);
chi_sq_1a = chi2inv(1-alpha/2,df);
lower_bound = sqrt((n-1)/(chi_sq_1a))*S;
upper_bound = sqrt((n-1)/(chi_sq_a))*S;
num_within_ci = sum(lower_bound <= sigma & sigma <= upper_bound);
percent_within_ci = num_within_ci/num_samples*100;
% Plot the histogram and theoretical chi-squared distribution

figure;
histogram(X,100,'Normalization','pdf');
hold on;

figure;
histogram(chi_sq,100,'Normalization','pdf');
hold on;
plot(x,y,'r','LineWidth',2);
plot(chi_sq_a ,0,'g.','markersize',20)
plot(chi_sq_1a,0, 'r.','markersize',20)
xlabel('\chi^2');
ylabel('pdf');
title(sprintf('Chi-squared distribution with %d degrees of freedom',df));
legend('Sample distribution','Χn-1^2');
%% ex18
clear
close all
% Parameters
lambda = 3; % Average number of events per unit time
T = 1000; % Time interval
num_simulations = 100000;
k_max = 12; % Maximum value of k


% Run simulations
counts = zeros(1, k_max + 1);
for i = 1:num_simulations
    N = lambda * T;
    occurrences = rand(1, N) * T;
    num_occurrences_in_interval = sum(occurrences >= 0 & occurrences <= 1);
    
    if num_occurrences_in_interval <= k_max
        counts(num_occurrences_in_interval + 1) = counts(num_occurrences_in_interval + 1) + 1;
    end
end

% Calculate simulated probabilities
simulated_probs = counts / num_simulations;

% Plot theoretical and simulated probabilities
figure;
bar(0:k_max, simulated_probs', 'grouped');
xlabel('k');
ylabel('Probability');
title('Poisson Distribution Probabilities (Theoretical and Simulated)');
legend('Simulated');
grid on;



