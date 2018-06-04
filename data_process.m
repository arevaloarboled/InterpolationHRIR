[x,y,z]=sph2cart(deg2rad(r(:,3)),deg2rad(r(:,4)),deg2rad(r(:,5)));
[x2,y2,z2]=sph2cart(deg2rad(r(:,3)),deg2rad(r(:,4)),deg2rad(r(:,5)));
scatter3(x,y,z,50,r(:,1)./max(r(:,1)),'filled');
catter3(x2,y2,z2,50,r(:,2)./max(r(:,2)),'filled');

%%%plot mesh
d=[20 30 40 50 75 100 130 160];
%e=[-40:10:90];
e=[-40:10:90];
coor=zeros(8*4,3);
c=1;
for i=1:8
    for j=1:14
        if j==11
            A=[0:10:359];
        elseif j==12
            A=[0:15:359];
        elseif j==13
            A=[0:30:359];
        elseif j==14
            A=0;
            coor(c,:)=[A e(j) d(i)];
            c=c+1;
            break;
        else
            A=[0:5:359];
        end
        a=size(A);
        for k=1:a(2)
            coor(c,:)=[A(k) e(j) d(i)];
            c=c+1;
        end
    end
end
[x,y,z]=sph2cart(deg2rad(coor(:,1)),deg2rad(coor(:,2)),deg2rad(coor(:,3)));

