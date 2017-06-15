%clf;
I=40; J=40;
a=0;b=1;
c=0;d=1;

t=0; tfinal=1; 
D=1;

% Source_Term=@(x,y,t,u,v) -(exp(-t)-1)*(sin(pi*x)+sin(pi*y))+u*(exp(-t)-1)*pi*cos(pi*x)...
%     +v*(exp(-t)-1)*pi*cos(pi*y)+D*(exp(-t)-1)*pi^2*(sin(pi*x)+sin(pi*y));
Source_Term=@(x,y,t,u,v) -(exp(-t))*(sin(pi*x)+sin(pi*y))+u*(exp(-t)-1)*pi*cos(pi*x)...
    +v*(exp(-t)-1)*pi*cos(pi*y)+D*(exp(-t)-1)*pi^2*(sin(pi*x)+sin(pi*y));

BC=@(x,y,t) (exp(-t)-1)*(sin(pi*x)+sin(pi*y));

Initial_Data=@(x,y) 0;
%Source_Term=@(x,y,t,u,v) 0;
Amanda_Singleton_Advection_Diffusion2D(D,Source_Term,BC,Initial_Data,a,b,c,d,I,J,tfinal);