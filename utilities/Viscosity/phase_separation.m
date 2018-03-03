function [A,B,x0,Halfa,Hbeta]=phase_separation(DF, Dalfa, Dbeta, QB, Qalfa, Qbeta, Hd)

%%  phase_separation computes the phase separation effect using the empirical formula found by Pries et al. 2005
%     A,B,x0 coefficients of empirical formula
%       Halfa hematocrit in branch alfa [-]
%       Hbeta hematocrit in branch beta [-]
%       DF diameter of mother vessel [µm]
%       Dalfa diameter of branch alfa [µm]
%       Dbeta diameter of branch beta [µm]
%       QB blood flow in mother vessel [same unit of measure of other blood flows]
%       Qalfa blood flow in mother vessel [same unit of measure of other blood flows]
%       Qbeta blood flow in mother vessel [same unit of measure of other blood flows]
%       Hd hematocrit of mother vessel
% 
%   Author: Simone Di Gregorio
%   Politecnico di Milano, 07/07/2017
%   Contact: simone.digre@gmail.com  
%%


A=-13.29*((Dalfa/Dbeta)^2-1)/((Dalfa/Dbeta)^2+1)*(1-Hd)/DF;
B=1+6.98*(1-Hd)/DF;
x0=0.964*(1-Hd)/DF;
FQB=Qalfa/QB;

if FQB <= x0
    FQE=0;
elseif FQB < 1 - x0 && FQB > x0
    X=(FQB-x0)/(1-2*x0);

    logitFQE=A+B*log(X/(1-X));

    FQE=exp(logitFQE)/(1+exp(logitFQE));
elseif FQB >= 1-x0
    FQE=1;
end

Halfa=FQE*Hd*QB/Qalfa;

Hbeta=(Hd*QB-Halfa*Qalfa)/Qbeta;