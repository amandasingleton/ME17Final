function outvar=Heat2D(D,Source_Term,BC,Initial_Data,a,b,c,d,I,J,dt,tfinal)
x=linspace(a,b,I); dx=x(2)-x(1);
y=linspace(c,d,J); dy=y(2)-y(1);

un=zeros(I,J);
unp1=zeros(I,J);
A=sparse(I*J,I*J);
RHS=zeros(I*J,1);
exact=zeros(I*J,1);
rho=zeros(I,J);


for i=1:I
    for j=1:J
    un(i,j)=Initial_Data(x(i),y(j));
    end
end

C = 1+D*(2*dt/dx/dx+2*dt/dy/dy);
L=-D*dt/dx/dx;
R=L;
T=-D*dt/dy/dy;
B=T;
t=0; 
while t<tfinal
    if t+dt>tfinal
        dt=tfinal-t;
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
            RHS(p)=un(p)+dt*Source_Term(x(i),y(i),t);
        end
    end
    
    %Define A and RHS for the botom & top walls
    for i=1:I
        %top wall
        j=J;
        p=(j-1)*I+i;
        A(p,p)=1; RHS(p)=BC(x(i),y(i),t);
        
        %bottom wall
        j=1;
        p=(j-1)*I+i;
        A(p,p)=1; RHS(p)=BC(x(i),y(i),t);
    end
    
    %Define A and RHS for the left & right walls
    for j=1:J
        %left wall
        i=1;
        p=(j-1)*I+i;
        A(p,p)=1; RHS(p)=BC(x(i),y(i),t);
        
        %right wall
        i=I;
        p=(j-1)*I+i;
        A(p,p)=1; RHS(p)=BC(x(i),y(i),t);
    end  
    
    unp1=A\RHS;
    for i=1:I
        for j=1:J
            BCtest(i,j)=BC(x(i),y(j),t);
            rho(i,j)=unp1((j-1)*I+i);
        end
    end
    
    t=t+dt;
    un=unp1;
    
    subplot(2,1,1);
    mesh(x,y,BCtest);
    axis([-1 1 -1 1 -1 1]);
    title('exact');
    subplot(2,1,2);
    mesh(x,y,rho);
    axis([-1 1 -1 1 -1 1]);
    title('scheme')
    pause(dt);
end

end