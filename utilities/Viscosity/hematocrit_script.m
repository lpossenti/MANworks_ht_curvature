%script for evaluation of phase separation effect

clear all
close all
clc

Hd=0.45;
Dalfa=0.0317*1E-4*2*1E6;
Dbeta=0.0317*1E-4*2*1E6;
DF=0.04*1E-4*2*1E6;
QB=0.6423353*DF^2*3.14/4;
Qalfa=0.52566659*Dalfa^2*3.14/4;
Qbeta=0.52549672*Dbeta^2*3.14/4;

[A,B,x0,Halfa,Hbeta]=phase_separation(DF, Dalfa, Dbeta, QB, Qalfa, Qbeta, Hd);