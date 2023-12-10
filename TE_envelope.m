%% --- Basic gearbox data

z1=30;
z2=47;
ns=38*60/2.5;
n1=2279;
n2=n1*z1/z2;
fs=ns/60;
f1=n1/60;
f2=n2/60;
fz=f1*z1;
fk=fz/lcm(z1,z2);
fmax=12000;
om1=n1*2*pi/60;
om2=n2*2*pi/60;
counter=ceil(om2/(2*pi/z2));

load TE_results.mat

fpodz2=[-0.75,1.71,-0.8,-0.8,-5.89,5.8,-1.47,-2.19,...
    -5.45,4.51,-4.55,-1.56,-5,1.88,-2.75,...
    -2.75,-2.75,-2.75,3.68,1.53,-5.97,-6.51,...
    5.48,-1.48,-1.48,-3.64,4.89,-0.49,0.27,-0.81,9.07,...
    -3.1,1.8,-0.72,-2.11,9.25,0.9,0.9,-2.25,8.85,...
    0,3.86,-4.13,7.72,-1.66,1.8,-1.08]/1000;
Fpodz2=cumsum(fpodz2);

fpodz1=[-1.63,1.77,-0.53,2.84,-4.04,-0.43,2.31,...
    -0.82,-1.38,3.51,4.86,-1.53,-0.82,1.45,...
    3.98,-0.68,-0.71,2.56,-1.96,-1,-0.68,-0.89,...
    2.52,-5.83,-1.99,3.48,-4.62,-0.46,0.96,-1.32]/1000;
Fpodz1=cumsum(fpodz1);

j1=1;
j2=1;
ODLfi2p=FI2;
ODLerrp=odl;
fikorp=fi_kor;

for i=1:counter
    if i==1
        ODLfi2=ODLfi2p;
        ODLerr=ODLerrp;
        FIkor=fikorp;
    else 
        
        ODLfi2=[ODLfi2,ODLfi2p+(i-1)*(2*pi/z2)];
        ODLerr=[ODLerr,ODLerrp];

        if j2==47
            j2=1;
        else
            j2=j2+1;
        end
        if j1==30
            j1=1;
        else
            j1=j1+1;
        end
    end
end

%% --- Envelope of TE signal

[FIsort,Isort] = sort(ODLfi2);
ODLsort=ODLerr(Isort);

[ODLupper,ODLlower] = envelope(ODLsort,1,'peak');

figure
hold on
plot(FIsort,ODLsort,'-r','Linewidth',2)
plot(ODLfi2,ODLerr,'-k','Linewidth',1)
hold off
grid on
ylabel('dist');
xlabel('fi2');

x=ODLlower/1000;
dt=gradient(FIsort)/om2;
t=cumsum(dt);
v=gradient(x)./dt;
a=gradient(v)./dt;

TE=[t;x;v;a];
csvwrite('TE_min.csv',TE)
