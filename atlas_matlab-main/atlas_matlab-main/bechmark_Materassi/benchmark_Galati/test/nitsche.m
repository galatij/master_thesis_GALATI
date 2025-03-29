close all
clear all
clc
sym=0; % 1 sym version
consistency=1; % 0 no consistency term
if sym == 1
    disp('symmetric version');
else
    disp('unsymmetric version');
end
if consistency == 0
    disp('Just penalization');
end
f=1;
mu=1;
L=2;
Nhp=20;
Nhm=20;
Ndof=Nhp+Nhm;
Ndoftot=Ndof+2; % + b points
ap=[1,Nhp];
am=[Nhp+1,Ndof];
h=L/Ndof;
A=sparse(Ndof,Ndof);
As=sparse(Ndof,Ndof);
gamma=10;
g=1;
%% stiffness part
A(ap(1),ap(1))=1; 
for i = ap(1):ap(2)-1  % plus part
    A(i,i    )=A(i,i  )+1;
    A(i,i+1  )=A(i,i+1)-1;
    A(i+1,i  )=A(i+1,i)-1;
    A(i+1,i+1)=A(i+1,i+1)+1;
end
% minus part
A(am(2),am(2))=1; 
for i = am(1):am(2)-1
    A(i,i    )=A(i,i  )+1;
    A(i,i+1  )=A(i,i+1)-1;
    A(i+1,i  )=A(i+1,i)-1;
    A(i+1,i+1)=A(i+1,i+1)+1;
end
%% Save just stiffness
As=A;
As=(mu/h)*As; % scale it 
%% Consistency term
p2=ap(2)-1;
p1=ap(2);
m1=am(1);
m2=am(1)+1;
if consistency ~=0
    % first row of plus domain
    %-(p2-p1+m1-m2) 1/2
    A(p1,p2)=A(p1,p2)+1/2;
    A(p1,p1)=A(p1,p1)-1/2;
    A(p1,m1)=A(p1,m1)+1/2;
    A(p1,m2)=A(p1,m2)-1/2;
    % last row of minus domain
    A(m1,p2)=A(m1,p2)-1/2;
    A(m1,p1)=A(m1,p1)+1/2;
    A(m1,m1)=A(m1,m1)-1/2;
    A(m1,m2)=A(m1,m2)+1/2;
end
%% Penalization
% first row of plus domain
A(p1,p1)=A(p1,p1) + gamma;
A(p1,m1)=A(p1,m1) - gamma;
% last row of minus domain
A(m1,p1)=A(m1,p1) - gamma;
A(m1,m1)=A(m1,m1) + gamma;
%%
% symmetry
if sym == 1 & consistency ~=0
    % column of plus (p1)
    % avg(v')*salto(u)
    % - -1/h 1/2 up1 - um1
    A(p2,p1)=A(p2,p1)+1/2;
    A(p1,p1)=A(p1,p1)-1/2;
    A(m1,p1)=A(m1,p1)+1/2;
    A(m2,p1)=A(m2,p1)-1/2; 
    % column of minus
    A(p2,m1)=A(p2,m1)-1/2;
    A(p1,m1)=A(p1,m1)+1/2;
    A(m1,m1)=A(m1,m1)-1/2;
    A(m2,m1)=A(m2,m1)+1/2; 
end
A=(mu/h)*A; % scale it
%% rhs
b=zeros(Ndof,1);
for i = 1:Ndof
    b(i)=f*h;
end
% correct
b(m1)=b(m1)/2;
b(p1)=b(p1)/2;
bs=b; %save rhs without penalty & symmetry
%% penalty
b(m1)=b(m1)-gamma*(mu/h)*g;
b(p1)=b(p1)+gamma*(mu/h)*g;
%% symmetry
if sym == 1  & consistency ~=0
    %-0.5 1/h mu
    b(p2)=b(p2)+0.5*(mu/h)*g;
    b(p1)=b(p1)-0.5*(mu/h)*g;
    b(m1)=b(m1)+0.5*(mu/h)*g;
    b(m2)=b(m2)-0.5*(mu/h)*g;
end
%%            
bp=(L^2*f)/(mu*(L + 2)) + (2*g)/(L + 2);
bm=bp;
cp=-(2*L*g)/(L + 2) - (L^2*f*(L - 2))/(2*mu*(L + 2));
 
a =-f/(2*mu);
uexp = @(x) a.*x.*x + bp.*x +cp;
uexm = @(x) a.*x.*x + bm.*x;
xval=[linspace(0,L/2,Nhp+1),linspace(L/2,L,Nhm+1)];
uex =[uexm([xval(1:Nhp+1)]),uexp([xval(Nhp+2:end)])];
plot(xval,uex);
grid on;
uh = A\b;
uh_ext=[0;uh;0];
hold on;
plot(xval,uh_ext);
legend('Exact','computed');

%%calcolo errore
xx=[h/2:h:L-h/2];
uexx =uexm(xx);
uexxp =uexp(xx);
Nm=length(xx)/2;

diffp=uexx(1:Nm) - 0.5*(uh_ext(1:Ndof/2)+uh_ext(2:Ndof/2+1))'
diffm=uexxp(Nm+1:end) - 0.5*(uh_ext(Ndof/2+2:end-1)+uh_ext(Ndof/2+3:end))'

eep=h*(diffp).^2;
eem=h*(diffm).^2;

err=(sum(eep)+sum(eem))^0.5;

disp('L2 error')
err
%%
%%
%%Derivative mu du/dx(L/2)
der_exact=-f*L/2+mu*bm;
disp(['Exact tension=',num2str(der_exact)]);
jump_uh=uh(p1)-uh(m1);
disp(['Jump = ',num2str(jump_uh), ' Exact jump= ',num2str(g), ' Error=',num2str(g-jump_uh)]);
nitschecontrib=(jump_uh-g)*gamma*mu/h;
disp(['Nitsche contribution to tension=',num2str(nitschecontrib)]);
avg_der=0.5*mu*(uh(p1)-uh(p2)+uh(m2)-uh(m1))/h;
disp(['tension computed via derivative=',num2str(avg_der)]);
jump_der=mu*(uh(p1)-uh(p2)-uh(m2)+uh(m1))/h;
disp(['error in tension continuity=',num2str(jump_der)]);

%% Now go vaiational
% Extract the rows we are interested
R=sparse(Ndof,Ndof);
R(Nhm,Nhm)=1;
R(m1,m1)=1;
R(p1,p1)=1;
Ar=R*As;
br=R*bs;
t=Ar*uh-br;
tp=t(p1);
tm=t(m1);
disp(['left and right tension computed variationally= ',num2str(tm),', ',num2str(tp)]);

disp(['avcerage tension computed variationally= ',num2str((-tm+tp)/2)]);





