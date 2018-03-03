function visco=Kiani(raggio, ematocrito, mu_plasma)


mu_c=exp(0.48+2.35*ematocrito);
D=raggio*2;
delta=2.03-2*ematocrito;
d_m=2.7;
visco=(1-(1-(mu_plasma/mu_c))*(1-2*delta/D)^4)^(-1)*(1-(d_m/D)^4)^(-1);