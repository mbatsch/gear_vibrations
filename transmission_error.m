load data.mat

fi1p=fi1pt;
fi1k=fi1kt;
N=200;
dfi1=(fi1k-fi1p)/(N-1);
figure
k=1;
for fi1=fi1p:dfi1:fi1k
    fi2=fi1/i12;
    Mf1=rot(0,0,-fi1);
    Mfh=[cos(kappay) sin(kappax)*sin(kappay) cos(kappax)*sin(kappay) 0; 0 cos(kappax) -sin(kappax) 0; sin(kappay) sin(kappax)*cos(kappay) cos(kappax)*cos(kappay) 0; 0 0 0 1];
    Mf2=Mfh*trans(ar+ax,ay,az)*rot(0,0,fi2);

    %% --- Surface of pinion tooth in f
    [imax, jmax]=size(x11m);
    n=1;
    for j=1:jmax
        for i=1:imax
            LI1(i,j)=i;
            LJ1(i,j)=j;
            r11=[x11m(i,j);y11m(i,j);z11m(i,j);1];
            r1f=Mf1*r11;
            x1f(i,j)=r1f(1,:);
            y1f(i,j)=r1f(2,:);
            z1f(i,j)=r1f(3,:);
        end
    end

    %% --- Surface of gear tooth in f
    [imax, jmax]=size(x22m);
    n=1;
    for j=1:jmax
        for i=1:imax
            LI2(i,j)=i;
            LJ2(i,j)=j;
            r22=[x22m(i,j);y22m(i,j);z22m(i,j);1];
            r2f=Mf2*r22;
            x2f(i,j)=r2f(1,:);
            y2f(i,j)=r2f(2,:);
            z2f(i,j)=r2f(3,:);
        end
    end
    
    %% --- Interpolant
    
    [th,ro,z]=cart2pol(x2f,y2f,z2f);
    F = scatteredInterpolant(ro(:),z(:),th(:));
    
    %% --- EaseOff plot
    
    [imax, jmax]=size(x11m);

    for j=1:jmax
        for i=1:imax

            z02=z1f(i,j);
            ro02=sqrt(x1f(i,j)^2+y1f(i,j)^2);
            th02=F(ro02,z02);

            [x2f0,y2f0,z2f0]=pol2cart(th02,ro02,z02);

            k0=sqrt((x1f(i,j)-x2f0)^2+(y1f(i,j)-y2f0)^2);
            K0(i,j)=k0;
            FIkor(i,j)=2*asin(k0/2/ro02);
            FIKOR(i,j,k)=FIkor(i,j);
            Kkor(i,j,k)=ro02*FIkor(i,j);

        end
    end
    odl(k)=min(K0,[],'all');
    fi_kor(k)=min(FIkor,[],'all');
    FI1(k)=fi1+fi_kor(k);
    FI2(k)=fi2;

    surf(z1f,sqrt(x1f.^2+y1f.^2),K0,'Facecolor',[1 (N-k)/N 0]);%,'Meshstyle','none')
    hold on
    xlabel 'z1f'
    ylabel 'ry1'
    zlabel 'fikor'
    grid on

    k=k+1;
end

figure
plot(FI2,fi_kor)
hold on
plot(FI2+1*2*pi/z2,fi_kor)
plot(FI2+2*2*pi/z2,fi_kor)
plot(FI2+3*2*pi/z2,fi_kor)

save('TE_results.mat','FI1','FI2','fi_kor','odl')
