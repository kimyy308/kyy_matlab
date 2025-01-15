function [Y,C,W,PM,F]=PRO(w,g,gamp,ffreq,n,namp,cubic,f,dt,L,IC)

% w=2*pi*w;
% ffreq=2*pi*ffreq;

%initialize

Y=zeros(2,length(0:dt:L));
C=Y; F=Y;
Y(:,1)=IC;
W=zeros(length(0:dt:L),1);

PM=zeros(2,2,length(0:dt:L));

C(1,1)=-cubic*IC(1)^3;

F(1,1)=parfun(f,ffreq,0);

N(1,1)=n*(1+parfun(namp,ffreq,0));
dW=sqrt(dt)*randn;
W(1)=dW;

%initialize propagator matrix

M=[g+parfun(gamp,ffreq,0) w; -w 0];
PM(:,:,1)=M;

% use Euler method for first step

Y(:,2)=Y(:,1)+dt*(M*Y(:,1)+C(:,1))+N(:,1)*dW+0.5*N(:,1)*((dW)^2-dt);
W(2)=W(1)+sqrt(dt)*randn;

%integrate

for t=dt:dt:L;
    
    dW=sqrt(dt)*randn;
    
    m=round(t/dt+1);
    
    M=[g+parfun(gamp,ffreq,t) w; -w 0];
    C(1,m)=-cubic*Y(1,m)^3;
    N(1,m)=n*(1+parfun(namp,ffreq,t));
    F(1,m)=parfun(f,ffreq,t);
    
    PM(:,:,m)=M;
    W(m+1)=W(m)+dW;
 
    Y(:,m+1)=Y(:,m)+dt*(M*Y(:,m)+C(:,m)+F(:,m))+N(:,m)*dW;
    
end

Y=Y(:,2:end-1);
C=C(:,2:end);
PM=PM(:,:,2:end);
W=W(2:end-1);
F=F(:,2:end-1);

end

function y=parfun(A,f,x)
y=A*cos(f*x);
end
