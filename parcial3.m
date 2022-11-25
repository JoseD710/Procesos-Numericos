clc
clear
close all
%datos de entrada
r2=3;
r3=8;
r4=5;
phi1=0.0;
phi4=pi/2;
q=[9;
   0.9;
   0.3];
%variables de actuación
phi2=0;
r10=1;
rp10=1;
rpp10=-1;
%grado de libertad
counter=0;
for t=0:1.25:1.25
    counter=counter+1;
    tol=100;
    iter=0;
    while tol>1e-8 && iter<100
        iter=iter+1;
        fq=[r2*cos(q(2))+r3*cos(q(3))-r4*cos(phi4)-q(1)*cos(phi1);
            r2*sin(q(2))+r3*sin(q(3))-r4*sin(phi4)-q(1)*sin(phi1);
            q(2)-r10-rp10*t-0.5*rpp10*t^2];

        J=[-cos(phi1), -r2*sin(q(2)), -r3*sin(q(3));
           -sin(phi1),  r2*cos(q(2)), r3*cos(q(3));
            0,          1,             0;];
        qi=-inv(J)*fq+q;
        q=qi;
        tol=norm(fq);
    end
    P(:,counter)=q;
    %calculo de la velocidad
    v=-inv(J)*[0; 0; -rp10-rpp10*t];
    V(:,counter)=v;
    %calculo de la aceleracion
    Jp=[0, -r2*cos(q(2))*v(2), -r3*cos(q(3))*v(3);
        0, -r2*sin(q(2))*v(2), -r3*sin(q(3))*v(3);
        0,  0,                 0;];
    phi_t_pp=[0; 0; -rpp10];
    a=-inv(J)*(Jp*v+phi_t_pp);
    A(:,counter)=a;
end
t=0:1.25:1.25;
%plot de la posicion
figure
plot(t,P(1,:),'r')
hold on
plot(t,P(2,:),'g')
plot(t,P(3,:),'b')
legend('r1','phi2','phi3')
xlabel('t')
%plot de la velocidad
figure
plot(t,V(1,:),'r')
hold on
plot(t,V(2,:),'g')
plot(t,V(3,:),'b')
%plot de la aceleracion
figure
plot(t,A(1,:),'r')
hold on
plot(t,A(2,:),'g')
plot(t,A(3,:),'b')

% figure
% for i=1:1:size(P,2)
%     O2=[0;0];
%     A=[r2*cos(P(2,i));  r2*sin(P(2,i))];


%     B=A-[r3*cos(P(3,i)); r3*sin(P(3,i))];
%     
%     clf
%     line([O2(1) A(1)],[O2(2) A(2)],'color','blue');
%     line([A(1) B(1)],[A(2) B(2)],'color','red');    
%     hold off
%     axis([-8 8 -8 8])
%     %axis equal
%     pause(0.1)
% end