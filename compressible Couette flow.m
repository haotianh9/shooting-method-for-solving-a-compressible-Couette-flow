%using a shooting method to solve a compressible Couette flow
%cold wall
clc;
clear all;
miu0=1.7161E-5;
T0=273.16;
k0=0.02415;
Ts=111;
ue=500;%the speed of the upper wall
D=0.1;%the height of the duct
Te=300;%temperature of the upper wall
Tw=300;%temperature of the lower wall
divN=501;%precision
y=linspace(0,D,divN);
dy=(D-0)/(divN-1);
u=y*ue/D;
tau0=miu0*ue/D;
T=linspace(Tw,Tw,divN);
eTey0=0;%partialT/partialy,y=0
tau=tau0;
tau1=tau0;
for i=1:10000;
    eTey(1)=eTey0;
    for j=1:10000;
        T(1)=Tw;
        for m=2:divN;
            k=k0*(T(m-1)/T0)^(3/2)*(T0+Ts)/(T(m-1)+Ts);
            T(m)=T(m-1)+(eTey+tau/k*(u(1)-u(m)))*dy;
        end
        if abs(T(divN)-Te)<1E-10
            break
        end
        eTey=eTey+(Te-T(divN))/D*0.01;       
    end
    for m=1:divN
        miu(m)=miu0*(T(m)/T0)^(3/2)*(T0+Ts)/(T(m)+Ts);
        k(m)=k0*(T(m)/T0)^(3/2)*(T0+Ts)/(T(m)+Ts);
    end
      for j=1:10000
          u(1)=0;
        for m=2:divN
           u(m)=u(m-1)+dy*(tau/miu(m));
        end
        if abs(u(divN)-ue)<1E-10
            break
        end
        tau=tau-(u(divN)-ue)/D*miu(floor(divN/2)); 
    end
      if abs(tau1-tau) <1E-10
          break
      end
      tau1=tau;
end
plot(u/ue-y/D,y/D),xlabel('u(y)/ue-y/D'),ylabel('y/D');
plot(T/Tw,y/D),xlabel('T(y)/Tw'),ylabel('y/D');


