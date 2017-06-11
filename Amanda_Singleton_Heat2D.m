clf;
I=80; J=80;
a=-1;b=1;
c=-1;d=1;

t=0; tfinal=4; dt=0.1;
D=0.5;

Source_Term=@(x,y,t) sin(x)*(0.1+D*0.03+D*0.1*t)+0.1*cos(y)*(1+D*t);
BC=@(x,y,t) 0.03*sin(x)+0.1*(sin(x)+cos(y))*t;
Initial_Data=@(x,y) 0.03*sin(x);

Heat2D(D,Source_Term,BC,Initial_Data,a,b,c,d,I,J,dt,tfinal);
