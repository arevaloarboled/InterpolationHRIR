d=[30,40,50,75,100,130];
a=[30,90,150,210,270,330];
e=[-30,0,30,60];
%d=[40,50,75,100,130];
%a=[210,270];
%e=[-30,0,30,60];
dtam=size(d);
atam=size(a);
etam=size(e);
r=zeros(dtam(2)+atam(2)+etam(2),2+3);
r2=zeros(dtam(2)+atam(2)+etam(2),2+3);
c=1;
for i=1:dtam(2)
    for j=1:atam(2)
        for k=1:etam(2)
            display([a(j) e(k) d(i)])
            [el1,er1,sdl1,sdr1,y1,soHs1]=comparingInter([a(j) e(k) d(i)],'r');
            [el2,er2,sdl2,sdr2,y2,soHs2]=comparingInter([a(j) e(k) d(i)],'l');
            r(c,:)=[sdl1 sdr1 a(j) e(k) d(i)];
            r2(c,:)=[sdl2 sdr2 a(j) e(k) d(i)];
            c=c+1;
        end
    end
end

