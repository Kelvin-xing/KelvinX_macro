%**************************************************************************
%                   LSE Macroeconomics Summer Program
%                   Part II: Heterogeneous Agents
%                   Instructor: Wouter J. Den Haan
%
%                   use of this program in any fee-based program requires
%                   explicit permission (wjdenhaan@gmail.com)
%**************************************************************************

%Cleaning
clear all; close all; clc

%% 1. Simulate data with true law of motion
%--------------------------------------------------------------------------

%Settings
T = 10000;

%Stochastics
%randn('state',1712)                                              %#ok<RAND>
rng(1712,'v5normal'); %this is the same as randn('state',1712) but not recommended by Matlab
                     %Matlab prefers rng(1712)
z = 0.01*randn(T,1);

%Simulation
k = zeros(T,1);
k(1) = 1;
for t = 2:T,
    k(t) = motion(k(t-1),z(t));
end

%% 2. Linear approximation of law of motion
%--------------------------------------------------------------------------

%Preparing left- and right-hand
%Y = ...;
%x = ...;

%Regression
b = X\y;

%% 3. Check R-squared: close to one?
%--------------------------------------------------------------------------

%Sum of squared
sse = norm(y-X*b)^2;                        %Error sum of squares
tss = norm(y-mean(y))^2;                    %Total sum of squares

%R-squared
R2 = 1-sse/tss;
fprintf('\n R-squared = %f \n \n',R2)

%% 4. Check simulation: similar?
%--------------------------------------------------------------------------

%Simulation
k2 = zeros(T,1);
k2(1) = 1;          %k2 should have the same starting value as k and the shocks should also be same.
for t = 2:T,
    k2(t) = [1 k2(t-1) z(t)] * b;
end

%Plot
figure(1)
plot(1:T,[k k2])
title('Simulation: similar?'), box off
legend('k','k2')
xlabel('\it{t}'), ylabel('\it{k_{t}}','Rotation',0)
axis([9000 10000 -inf inf])