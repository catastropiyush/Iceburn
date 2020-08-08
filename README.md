# Iceburn
# Stars burn
G = 6.67*10^(-11);       //Gravitational constant
   sig = 5.670373*10^(-8);  //Stefan-Boltzmann constant
   con= 2.287*10^(-28);     //constant used in Luminosity equation
   k = 1.3806488*10^(-23);  //Boltzmann constant
   miu = 0.6;               //mean molecular weight
   pi = 3.142596;           
   mp = 1.672*10^(-27);     //mass of proton
   n = 10000;               //number of iterations
   ro = 695500000;          //radius
   h=ro/n
//Initialize $t$ and $y$ as length n+1 column vectors.
    r = zeros(n+1,1);
    P = zeros(n+1,1);
    M = zeros(n+1,1);
    rho= zeros(n+1,1);
    L = zeros(n+1,1);
    T = zeros(n+1,1);
    opa = zeros(n+1,1);
//Declare the initial values of the variable on the surface of the star.
r(1) = 0;
P(1) = 2.477*10^(16);
M(1) = 0;
L(1) = 0;
T(1) = 1.57*10^7;
rho(1) = P(1)* miu *mp/(k*T(1));
opam = 0.01*0.02;
opaH(1) = (1.1*10^(-25)*0.02^(0.5)*rho(1)^0.5*T(1)^0.5);
opae(1) = (0.2*(1+0.71)*(1+2.7*10^11*rho(1)/(T(1)^2))^-1*(1+ (T(1) / 4.5*10^8)^0.86)^(-1));
opak(1) = (4*10^25*(1+0.71)*(0.02+0.001)*rho(1)/(T(1)^3.5));
opa(1) = opam(1) + (opaH(1)^(-1) + ( opae(1) + opak(1) )^(-1))^(-1);
r(2) = r(1) + h;
P(2) = P(1) - ((h^2)/2)*G*(4/3)*pi*rho(1)^2;
T(2) = T(1) - h^2*opa(1)*rho(1)^2*T(1)*con/(32*sig);
M(2) = h^3*(4*pi/3)*rho(1);
L(2) = (4*pi/3)*h^3*rho(1)*con*T(1)^4;
for i = 2:n
 rho(i) = P(i)* miu *mp/(k*T(i));
 opaH(i) = (1.1*10^(-25)*0.02^(0.5)*rho(i)^0.5*T(i)^0.5);
 opae(i) = (0.2*(1+0.71)*(1+2.7*10^11*rho(i)/(T(i)^2))^-1*(1+ (T(i) / 4.5*10^8)^0.86)^(-1));
 opak(i) = (4*10^25*(1+0.71)*(0.02+0.001)*rho(i)/(T(i)^3.5));
 opa(i) = opam + (opaH(i)^(-1) + ( opae(i) + opak(i) )^(-1))^(-1);
 r(i+1) = r(i) + h;
 P(i+1) = P(i) - h*G*M(i)*rho(i)/(r(i)^2);
 M(i+1) = M(i) + h*4*pi*r(i)^2*rho(i);
 T(i+1) = T(i) - h*3*opa(i)*rho(i)*L(i)/(64*pi*sig*r(i)^2*T(i)^3);
 L(i+1) = L(i) + h*4*pi*r(i)^2*rho(i)*con*T(i)^4;
end
    rho(n+1) = P(n+1)* miu *mp/(k*T(n+1));
subplot(2,3,1);
plot(r/ro,log(T),'r','LineWidth',2)
set(gca(),"grid",[2 2]);
title('Temperature');
('square')

subplot(2,3,2)
plot(r/ro,log(L),'b','LineWidth',2)
set(gca(),"grid",[2 2])
title('Luminosity');
('square')

subplot(2,3,3)
plot(r/ro,log(rho),'g','LineWidth',2);
set(gca(),"grid",[2 2])
title('Density');
('square')
