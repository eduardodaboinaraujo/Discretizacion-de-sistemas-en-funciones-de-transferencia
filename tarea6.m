clear 
clc 
close all

syms s t
J=0.00048115; %Kg*m^2
Kt=0.22076;   %N*m/A  
Ka=0.22076;   %V/rad*s
Bo=0.0026829;  %N*m*s
Rf=4.08;      %Ohms
Lf=0.011307;  %H

num=[Kt];
den=[J*Lf (Bo*Lf+J*Rf) Bo*Rf Ka*Kt];
g=tf(num,den);

[A B C D]=tf2ss(num,den);

% ---- Tiempo de muestreo adecuado ---- %
Frmx=bandwidth(g);
T=5*Frmx;
To=1/T;
%----------------------------------------%
% ----- Funcion discreta metodo ZOH -----%
numz=[0.113];
denz=[1 -1.808 0.8328 -5.163e-06];
gzoh=tf(numz,denz,To);
%-----------------------------------------%
% ----- Comando c2d metodo ZOH -----------%
Gzoh=c2d(g,To,'zoh');
%-----------------------------------------%
% ----- Funcion discreta metodo trapezoidal -----%
numt=[0.22 0.66 0.66 0.22];
dent=[8.19 -8.75 -3.95 4.91];
gtus=tf(numt,dent,To);
%------------------------------------------------%
% -------- Comando c2d metodo Tustin -----------%
Gtus=c2d(g,To,'tustin');
% ----- Funcion discreta metodo transformación de polos y ceros -----%
numc=[0.2829e-1 0.5658e-1 0.2829e-1];
denp=[1 -1.808 0.8328 -0.5163e-5];
gpc=tf(numc,denp,To);
%--------------------------------------------------------------------%
% ------- Ubicacion de polos -------%
% -- ZOH --%
 [p,z]=pzmap(numz,denz);
%-- Tustin --%
 [pt,zt]=pzmap(numt,dent);
% -- Polos y ceros --%
 [pp,zz]=pzmap(numc,denp);
%------- discretizacion exacta -------%
i=A*inv(A);
I=s*i;
Ai1=I-A;
Aj=inv(Ai1);
Phi=ilaplace(Aj);
Bi1=Phi*B;
gamma=int(Bi1,0,t);

Ad1=double(subs(Phi,t,To));
Bd1=double(subs(gamma,t,To));

sys1=ss(Ad1,Bd1,C,D,To);
Gz1=tf(sys1);

figure(1)
step(g)
hold on 
step(gzoh)
hold on 
step(Gzoh,'k')
title('Representacion del sistema continua vs discretizado metodo ZOH');
legend('continua','discreta ZOH','discreta C2d');
grid on

figure(2)
step(g)
hold on 
step(gtus)
hold on 
step(Gtus,'k')
title('Representacion del sistema continua vs discretizado metodo Tustin');
legend('continua','discreta Tustin','discreta C2d');
grid on

figure(3)
step(g)
hold on  
step(gpc,'k')
title('Representacion del sistema continua vs discretizado por transformación de polos y ceros');
legend('continua','discreta P-Z');
grid on

figure(4)
axis([-1 1 -1 1]);
zplane(z,p);
title('Ubicacion de Ceros-Polos metodo Zoh');
grid on
 
figure(5)
axis([-1 1 -1 1]);
zplane(zt,pt);
title('Ubicacion de Ceros-Polos metodo Tustin');
grid on

figure(6)
axis([-1 1 -1 1]);
zplane(zz,pp);
title('Ubicacion de Ceros-Polos Metodo polos y ceros');
grid on

figure(7)
step(g)
hold on 
step(gzoh,'r')
hold on
step(gtus,'g')
hold on
step(gpc,'k')
title('Analisis temporal ');
legend('continua','discreta ZOH','discreta tustin','discreta P-Z');
grid on

figure(8)
step(g)
hold on 
step(Gz1,'y')
hold on 
step(gzoh,'r')
hold on
step(gtus,'g')
hold on
step(gpc,'k')
title('Representacion del sistema continua vs discretizado Exacta,ZOH,Tustin,P-Z');
legend('continua','discreta Exacta','discreta ZOH','discreta Tustin','discreta P-Z');
grid on