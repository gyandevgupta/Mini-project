% step size
h=0.1;  

t_final=15;
x = 0:h:t_final;

m =0.5;
r =0.05;
d =0.08;
I=[(3/5)*m*(r^2/4+d^2), 0, 0;0, (3/5)*m*(r^2/4+d^2), 0;0, 0, (3/10)*m*r^2];

% Calculates upto y(t_final)
% creates a matrix of length(x)*3
y = zeros(length(x),3);                            

% initial condition          
w = [4,2,6];
y(1,:) = w;
r = [1,0,0;0,1,0;0,0,1];
wcross=[0,-w(1,3),w(1,2);w(1,3),0,-w(1,1);-w(1,2),w(1,1),0];
%
Q=zeros(4,length(x));
R=r*expm(wcross.*h);
Q(:,1)=rotm2quat(R);
Lin=zeros(length(x),3);
Lin(1,:)=transpose(inv(R)*(I*transpose(w)));

% plot(x,y(:,1),'g',x,y(:,2),'b',x,y(:,3),'r')
% xlabel('time')
% ylabel('W_x,W_y,W_z')
% legend({'x,y(:,1)','x,y(:,2)','x,y(:,3)'},'Location','southwest')
% legend({'W_x','W_y','W_z'},'Location','southwest')

for i=1:(length(x)-1) 
    
    k_1 = RHS_euler_eq(y(i,:));
    k_2 = RHS_euler_eq(y(i,:)+0.5*h*k_1);
    k_3 = RHS_euler_eq(y(i,:)+0.5*h*k_2);
    k_4 = RHS_euler_eq(y(i,:)+k_3*h);
    y(i+1, :) = y(i, :) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
    R=R*expm([0,-y(i+1,3),y(i+1,2);y(i+1,3),0,-y(i+1,1);-y(i+1,2),y(i+1,1),0].*h);
    Lin(i+1,:)=inv(R)*(I*transpose(y(i+1,:)));
    Q(:,i+1)=rotm2quat(R);
end
L=y*I;

function [result] = RHS_euler_eq(w)
m =0.5;
r =0.05;
d =0.08;
I=[(3/5)*m*(r^2/4+d^2), 0, 0;0, (3/5)*m*(r^2/4+d^2), 0;0, 0, (3/10)*m*r^2];
% Making column vector
w = transpose(w);
result = inv(I)*cross(I*w,w);

% Making a row vector
result = transpose(result);

end