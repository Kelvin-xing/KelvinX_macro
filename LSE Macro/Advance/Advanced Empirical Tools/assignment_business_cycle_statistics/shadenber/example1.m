%
% program to plot series with NBER bands
%
% the dates are stored in dateX.mat (they can be created with
% create_dates_xlabels.m
%

clear all; close all; clc
load dateX
TT = length(dateX);
r1 = randn(TT,1);
r2 = randn(TT,1);
figure('Name','example graph')
hold all
title('\bf{random series with NBER recessions}')
plot(dateX,r1,'LineWidth',1.5,'Color','k')
hold all
plot(dateX,r2,'LineWidth',1.5,'Color','b')
xlabel('year')
datetick('x','YY')
%set(gcf,'Renderer','zbuffer');
shadenber
maximize