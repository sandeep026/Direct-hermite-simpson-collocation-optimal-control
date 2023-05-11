clear
close all

import casadi.*

T=2;
a=1;
b=-2.694528;
c=-1.155356;
N=50;
opti=Opti();
X=opti.variable(2,N+1);
U=opti.variable(1,N+1);

X0=[0;0];

opti.subject_to(X(:,1)-X0==0)

dt = T/N; % length of a control interval

for k=1:N % loop over control intervals
% collocation point
Xmid=0.5*(X(:,k)+X(:,k+1))+dt/8*(ode_fun(X(:,k),U(:,k))-ode_fun(X(:,k+1),U(:,k+1)));
Umid=0.5*(U(:,k)+U(:,k+1));
%simpson quadrature
opti.subject_to( X(:,k+1)-X(:,k)-dt/6*(ode_fun(X(:,k),U(:,k))+4*ode_fun(Xmid,Umid)+ode_fun(X(:,k+1),U(:,k+1)))==0)
end

s=0;
for k=1:N
Xmid=0.5*(X(:,k)+X(:,k+1))+dt/8*(ode_fun(X(:,k),U(:,k))-ode_fun(X(:,k+1),U(:,k+1)));
Umid=0.5*(U(:,k)+U(:,k+1));
s=s+dt/6*(obj(X(:,k),U(k))+obj(X(k+1),U(k+1))+4*obj(Xmid,Umid));
end
obj=s;


opti.minimize(obj)
opti.subject_to(a*X(1,end)+b*X(2,end)-c==0)

opti.solver('ipopt')
sol=opti.solve();

% Compare the numerical and analytical solution
t=linspace(0,2,N+1);
X=sol.value(X);
U=sol.value(U);
[Xa,Ua]=analytical_solution(t);
plot(Xa(1,:),Xa(2,:))
hold on
scatter(X(1,:),X(2,:))
grid
xlabel('$x_1$')
ylabel('$x_2$')
legend('analytical','single shooting','Location','northwest')


 
 
figure
plot(t,Ua)
hold on
plot(t,U)
xlabel('time [s]')
ylabel('f')
legend('analytical','single shooting','Location','northwest')

figure
Lag=opti.f+ opti.lam_g'*opti.g;
H=hessian(Lag,opti.x);
spy(H) 
legend('Hessian sparsity') 
figure 
Jac=jacobian(opti.g,opti.x);
spy(Jac)
legend('Jacobian sparsity') 
