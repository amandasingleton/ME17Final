function outvar=Amanda_Singleton_Advection_Diffusion2D(D,Source_Term,BC,Initial_Data,a,b,c,d,I,J,tfinal)
x=linspace(a,b,I); dx=x(2)-x(1);
y=linspace(c,d,J); dy=y(2)-y(1);

un=zeros(I*J);
unp1=zeros(I*J);
A=zeros(I*J,I*J);
RHS=zeros(I*J,1);
exact=zeros(I*J,1);
rho=zeros(I,J);
dt=0.25*dx; t=0;

u = zeros(I*J,1);
v = zeros(I*J,1);

for i=1:I
    for j=1:J
        p=(j-1)*I+i;
%         un(p)=BC(x(i),y(j),t+dt);
          un(p)=0;
        u(p)=0;
        v(p)=y(j)*(1-y(j));
    end
end

C = 1+D*2*dt*(1/dx/dx+1/dy/dy);
L= -D*dt/dx/dx;
R= L;
T= -D*dt/dy/dy;
B= T;

% u = @(x,y) 0;
% v = @(x,y) y*(1-y);
while t<tfinal
   if t+dt>tfinal
       dt=tfinal-t;
       C = 1 + 2*D*dt/dx/dx + 2*D*dt/dy/dy;
        L =     -D*dt/dx/dx;
        R =     -D*dt/dx/dx;
        B =                   -D*dt/dy/dy;
        T =                   -D*dt/dy/dy;
   end
    
    %Define A and RHS for interior points
    for i=2:I-1
        for j=2:J-1
            p=(j-1)*I+i;
            A(p,p)=C;
            A(p,p+1)=R;
            A(p,p-1)=L;
            A(p,p+I)=T;
            A(p,p-I)=B;
           % RHS(p)=un(p)*(1-u(x(i),y(j))*dt/dx-v(x(i),y(j))*dt/dy)+u(x(i),y(j))*dt/dx*un(p-1)+v(x(i),y(j))*dt/dy*un(p-I)+dt*Source_Term(x(i),y(i),t,u(x(i),y(j)),v(x(i),y(j)));
           %if v is positive
           if v(p)>0
           RHS(p)=un(p)-dt/dx*u(p)*(un(p)-un(p-1))-dt/dy*v(p)*(un(p)-un(p-I))+dt*Source_Term(x(i),y(j),t,u(p),v(p));
           else
           RHS(p)=un(p)-dt/dx*u(p)*(un(p)-un(p-1))-dt/dy*v(p)*(un(p+I)-un(p))+dt*Source_Term(x(i),y(j),t,u(p),v(p));    
           end
         end
    end
    
    %Define A and RHS for the botom & top walls
    for i=1:I
        %top wall
        j=J;
        p=(j-1)*I+i;
        A(p,p)=1; RHS(p)=BC(x(i),y(j),t+dt);
        
        %bottom wall
        j=1;
        p=(j-1)*I+i;
        A(p,p)=1; RHS(p)=BC(x(i),y(j),t+dt);
    end
    
    %Define A and RHS for the left & right walls
    for j=1:J
        %left wall
        i=1;
        p=(j-1)*I+i;
        A(p,p)=1; RHS(p)=BC(x(i),y(j),t+dt);
        
        %right wall
        i=I;
        p=(j-1)*I+i;
        A(p,p)=1; RHS(p)=BC(x(i),y(j),t+dt);
    end  
    A;
    RHS;
    unp1=A\RHS;
    for i=1:I
        for j=1:J
            p=(j-1)*I+i;
            BCtest(i,j)=BC(x(i),y(j),t+dt);
            rho(i,j)=unp1(p);
            %rho(i,j)=BC(x(i),y(j),t+dt);
            exact(p)=BC(x(i),y(j),t+dt);
        end
    end
    
    
    figure(1)
    subplot(2,1,1);
    mesh(x,y,BCtest);
    axis([0 1 0 1 -2 0]);
    title('exact');
    subplot(2,1,2);
    mesh(x,y,rho);
    axis([0 1 0 1 -2 0]);
    title('scheme')
    pause(dt);
    t=t+dt;
    un=unp1;
end
% Error analysis:
maxError = max(abs(unp1 - exact));
fprintf('The maximum error on a %dx%d grid is %2.2e. \n',I, J, maxError);
end