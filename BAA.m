%[晶相]=BAA[坐标数据]
function [Data,Y]=BAA(data)
Point=data;
N=length(Point);
%----------------------------------------------------------------------
%[x,y,z]
x=Point(:,1);
y=Point(:,2);
z=Point(:,3);
clear Point
%----------------------------------------------------------------------
%N0
for i=1:N
    x0=x(i);y0=y(i);z0=z(i);
    d=sqrt((x-x0).^2+(y-y0).^2+(z-z0).^2);
    t=sort(d);
    p(i,1)=(sqrt(sum(t(2:7).^2)/6));
    max(i)=length(find((d<=sqrt(1.45*p(i,1)^2))&(d>0)));
end
max_line=fliplr(sort(max));
clear i x0 y0 z0 d t p
for i=1:N
    x0=x(i);y0=y(i);z0=z(i);
    d=sqrt((x-x0).^2+(y-y0).^2+(z-z0).^2);
    t=sort(d);
    p(i,1)=(sqrt(sum(t(2:7).^2)/6));
    s=find(d<=sqrt(1.45*p(i,1)^2)&d>0);
    if length(s)==max_line(1)
        N_0(i,:)=s';
    else
       t=(max_line(1)-length(s));
       for j=1:t
           s=[s;0];
       end
        N_0(i,:)=s';
    end
end
clear x0 y0 z0 i  t max s p j d max_line
%----------------------------------------------------------------------
%N1
for i=1:N
    x0=x(i);y0=y(i);z0=z(i);
    d=sqrt((x-x0).^2+(y-y0).^2+(z-z0).^2);
    t=sort(d);
    p(i,1)=(sqrt(sum(t(2:7).^2)/6));
    max(i)=length(find((d<=sqrt(1.55*p(i,1)^2))&(d>0)));
end
max_line=fliplr(sort(max));
clear i x0 y0 z0 d t p
for i=1:N
    x0=x(i);y0=y(i);z0=z(i);
    d=sqrt((x-x0).^2+(y-y0).^2+(z-z0).^2);
    t=sort(d);
    p(i,1)=(sqrt(sum(t(2:7).^2)/6));
    s=find(d<=sqrt(1.55*p(i,1)^2)&d>0);
    if length(s)==max_line(1)
        N_1(i,:)=s';
    else
       t=(max_line(1)-length(s));
       for j=1:t
           s=[s;0];
       end
        N_1(i,:)=s';
    end
end
clear x0 y0 z0 i  t max s p j d max_line
%----------------------------------------------------------------------
%angle_function
Y=zeros(N,8);
for i=1:N
    for j=1:length(nonzeros(N_0(i,:)))
        for k=j+1:length(nonzeros(N_0(i,:)))
            a=[x(N_0(i,j))-x(i),y(N_0(i,j))-y(i),z(N_0(i,j))-z(i)];
            b=[x(N_0(i,k))-x(i),y(N_0(i,k))-y(i),z(N_0(i,k))-z(i)];
            cos_kij=dot(a,b)/norm(a)/norm(b);
        if ((cos_kij>=-1.0)&&(cos_kij<-0.945))
                    Y(i,1)=Y(i,1)+1;
            elseif ((cos_kij>=-0.945)&&(cos_kij<-0.915))
                    Y(i,2)=Y(i,2)+1;
            elseif ((cos_kij>=-0.915)&&(cos_kij<-0.755))
                    Y(i,3)=Y(i,3)+1;
            elseif ((cos_kij>=-0.755)&&(cos_kij<-0.195))
                    Y(i,4)=Y(i,4)+1;
            elseif ((cos_kij>=-0.195)&&(cos_kij<0.195))
                    Y(i,5)=Y(i,5)+1;
            elseif ((cos_kij>=0.195)&&(cos_kij<0.245))
                    Y(i,6)=Y(i,6)+1;
            elseif ((cos_kij>=0.245)&&(cos_kij<0.795))
                    Y(i,7)=Y(i,7)+1;
            elseif ((cos_kij>=0.795)&&(cos_kij<1.0))
                    Y(i,8)=Y(i,8)+1;
        end
        end
    end
end
clear i j k cos_kij a b
%----------------------------------------------------------------------
%delta
delta(:,1)=0.35*Y(:,5)./(Y(:,6)+Y(:,7)-Y(:,5));                     %bcc
delta(:,2)=0.61*(abs(Y(:,1)+Y(:,2)-6)+Y(:,3))/6;                    %fcc
delta(:,3)=(abs(Y(:,1)-3)+abs(Y(:,1)+Y(:,2)+Y(:,3)+Y(:,4)-9))/12;   %hcp
delta(:,4)=abs(1-Y(:,7)./24);                                       %cp
%----------------------------------------------------------------------
%judge
str=zeros(N,1);
for i=1:N
    if length(nonzeros(N_0(i,:)))<6
        str(i:1)=0;
    end
    if Y(i,1)==7
         delta(i,1)=0;
         str(i,1)=1;
    elseif Y(i,1)==6
         delta(i,2)=0;
         str(i,1)=2;
    elseif Y(i,1)<=3
         delta(i,3)=0;
         str(i,1)=3;
     end
end
clear i
for i=1:N
        if Y(i,8)>0
             str(i,:)=0;
        elseif Y(i,5)<3
            if (length(nonzeros(N_1(i,:)))>13||length(nonzeros(N_1(i,:)))<11)
                str(i,:)=0;
            else
                str(i,:)=4;
            end
        elseif (delta(i,1)<delta(i,4))
            if (length(nonzeros(N_1(i,:)))<11)
                str(i,:)=0;
            else
                str(i,:)=1;
            end
        elseif ((length(nonzeros(N_1(i,:)))>12)||(length(nonzeros(N_1(i,:)))<11))
                str(i,:)=0;
        elseif (delta(i,2)<delta(i,3))
                str(i,:)=2;
        else
                str(i,:)=3;
        end
end
clear i

Data=data;
Data(:,4)=str;
return 
end