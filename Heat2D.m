function outvar=Heat2D(D,Source_Term,BC,Initial_Data,a,b,c,d,I,J,dt,tfinal)
X=linspace(a,b,I); dx=X(2)-X(1);
Y=linspace(c,d,J); dy=Y(2)-Y(1);

un=zeros(I*J);
unp1=zeros(I*J);
A=sparse(I*J,I*J);
RHS=zeros(I*J,1);
exact=zeros(I*J,1);
rho=zeros(I,J);


for i=1:I
    for j=1:J
        p=((j-1)*I+i);
    un(p)=Initial_Data(X(i),Y(j));
    end
end
x=1/dx/dx;
y=1/dy/dy;

C = 1+D*dt*(x+y);
L = -D*dt/2*x;
R = L;
T=-D*dt/2*y;
B=T;

nC = 1-D*dt*(x+y);
nL = D*dt/2*x;
nR = nL;
nT = D*dt/2*y;
nB = nT;
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
            RHS(p)= nC*un(p)+nL*un(p-1)+nR*un(p+1)+nT*un(p+I)+nB*un(p-I)+0.5*dt*Source_Term(X(i),Y(i),t)+0.5*dt*Source_Term(X(i),Y(i),t+dt);
            
        end
    end
    
    %Define A and RHS for the botom & top walls
    for i=1:I
        %top wall
        j=J;
        p=(j-1)*I+i;
        A(p,p)=1; RHS(p)=BC(X(i),Y(i),t+dt);
        
        %bottom wall
        j=1;
        p=(j-1)*I+i;
        A(p,p)=1; RHS(p)=BC(X(i),Y(i),t+dt);
    end
    
    %Define A and RHS for the left & right walls
    for j=1:J
        %left wall
        i=1;
        p=(j-1)*I+i;
        A(p,p)=1; RHS(p)=BC(X(i),Y(i),t+dt);
        
        %right wall
        i=I;
        p=(j-1)*I+i;
        A(p,p)=1; RHS(p)=BC(X(i),Y(i),t+dt);
    end  
    
    unp1=A\RHS;
    for i=1:I
        for j=1:J
            p=(j-1)*I+i;
            BCtest(i,j)=BC(X(i),Y(j),t+dt);
            exact(p)=BC(X(i),Y(j),t+dt);
            rho(i,j)=unp1(p);
        end
    end
    
    t=t+dt;
    un=unp1;
    
    subplot(2,1,1);
    mesh(X,Y,BCtest);
    axis([-1 1 -1 1 -1 1]);
    title('exact');
    subplot(2,1,2);
    mesh(X,Y,rho);
    axis([-1 1 -1 1 -1 1]);
    title('scheme')
    pause(dt);
end
% Error analysis:
maxError = max(abs(unp1 - exact));
fprintf('The maximum error on a %dx%d grid is %2.2e. \n',I, J, maxError);
end
