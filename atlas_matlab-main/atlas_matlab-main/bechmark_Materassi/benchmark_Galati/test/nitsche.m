close all
clear all
clc
%% Params
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
gamma=10;
g=1;

%% Mesh params
L=2;
Nhp=20;  % Number of points in the plus part
Nhm=20;  % Number of nodes in the minus part
Ndof=Nhp+Nhm;  % number of DOFs
Ndoftot=Ndof+2;  % + b points
ap=[1,Nhp]; % Nodes in the plus part
am=[Nhp+1,Ndof];  % Nodes in the minus part
h=L/Ndof;  % Mesh spacing
A=sparse(Ndof,Ndof);
As=sparse(Ndof,Ndof);

%% Stiffness part
% Local stiffness: 1D, linear FE, grad-grad matrix local matrix
%          |  1   -1 |
% Ke = 1/h*|         |  --> add contributions to global elements
%          | -1    1 |

A(ap(1),ap(1))=1;  % Enforcing Dirichlet
for i = ap(1):ap(2)-1  % Plus part
    A(i,i    )=A(i,i  )+1;
    A(i,i+1  )=A(i,i+1)-1;
    A(i+1,i  )=A(i+1,i)-1;
    A(i+1,i+1)=A(i+1,i+1)+1;
end

% minus part
A(am(2),am(2))=1;  % Enforcing Dirichlet
for i = am(1):am(2)-1
    A(i,i    )=A(i,i  )+1;
    A(i,i+1  )=A(i,i+1)-1;
    A(i+1,i  )=A(i+1,i)-1;
    A(i+1,i+1)=A(i+1,i+1)+1;
end


%% Save just stiffness
As=A;  % I will modify this adding consistency and other terms, see below...
As=(mu/h)*As;  % Scale it (Poisson problem)

%% Consistency term
p2 = ap(2)-1;  % one before the last node in plus domain
p1 = ap(2);    % last node in plus domain (interface node)
m1 = am(1);    % first node in minus domain (interface node)
m2 = am(1)+1;  % second node in minus domain

if consistency ~=0
    % last row of plus domain
    %-(p2-p1+m1-m2) 1/2
    A(p1,p2)=A(p1,p2)+1/2;  % interaction with neigh in plus domain
    A(p1,p1)=A(p1,p1)-1/2;  % diagonal correction
    A(p1,m1)=A(p1,m1)+1/2;  % coupling with minus domain interface
    A(p1,m2)=A(p1,m2)-1/2;  % interaction with neigh in minus domain

    % first row of minus domain
    A(m1,p2)=A(m1,p2)-1/2;  % interaction with neigh in plus domain
    A(m1,p1)=A(m1,p1)+1/2;  % coupling with plus domain interface
    A(m1,m1)=A(m1,m1)-1/2;  % diagonal correction
    A(m1,m2)=A(m1,m2)+1/2;  % interaction with neigh in minus domain
end


%% Penalization
% Penalty term: 1D, linear FE, jump-jump local matrix
%                 |  1   -1 |
% Pe = gamma*mu/h*|         |  --> add contributions to global elements
%                 | -1    1 |

% first row of plus domain
A(p1,p1)=A(p1,p1) + gamma;
A(p1,m1)=A(p1,m1) - gamma;
% last row of minus domain
A(m1,p1)=A(m1,p1) - gamma;
A(m1,m1)=A(m1,m1) + gamma;


%% Symmetry
% Simply the transpose of the consistency term
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


%% Rhs
b=zeros(Ndof,1);
for i = 1:Ndof
    b(i)=f*h;  % Piecewise constant quadrature rule
end

%% Correct
% To avoid double-counting the shared point
b(m1)=b(m1)/2;
b(p1)=b(p1)/2;
bs=b;  % Save rhs without penalty & symmetry

% Penalty: impose weakly the jump
b(m1)=b(m1)-gamma*(mu/h)*g;
b(p1)=b(p1)+gamma*(mu/h)*g;

% Symmetry: from IBP for a nonhomogeneous jump
if sym == 1  & consistency ~=0
    % -0.5 1/h mu
    b(p2)=b(p2)+0.5*(mu/h)*g;
    b(p1)=b(p1)-0.5*(mu/h)*g;
    b(m1)=b(m1)+0.5*(mu/h)*g;
    b(m2)=b(m2)-0.5*(mu/h)*g;
end


%% Exact solution
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

%% Numerical solution
uh = A\b;
uh_ext=[0;uh;0];
hold on;
plot(xval,uh_ext);
legend('Exact','Computed');

%% Error analysis
xx=[h/2:h:L-h/2];
uexx =uexm(xx);
uexxp =uexp(xx);
Nm=length(xx)/2;

diffp = uexx(1:Nm) - 0.5*(uh_ext(1:Ndof/2)+uh_ext(2:Ndof/2+1))';
diffm = uexxp(Nm+1:end) - 0.5*(uh_ext(Ndof/2+2:end-1)+uh_ext(Ndof/2+3:end))';

eep=h*(diffp).^2;
eem=h*(diffm).^2;

err=(sum(eep)+sum(eem))^0.5;

disp('L2 error')
disp(err)


%% Derivative mu du/dx(L/2)
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





