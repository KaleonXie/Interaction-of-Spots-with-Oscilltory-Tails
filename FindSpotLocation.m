function [x_c,y_c]=FindSpotLocation(xq,yq,uq,spots_num)

x_c=[];
c=[];
uq_n=uq;
x_n=xq;
y_n=yq;

for i=1:spots_num
        if i>1
            u_o=uq_n;
            x_o=x_n;
            y_o=y_n;
            uq_n=u_o(sqrt((x_o-x_c(i-1)).^2+(y_o-y_c(i-1)).^2)>0.1);
            x_n=x_n(sqrt((x_o-x_c(i-1)).^2+(y_o-y_c(i-1)).^2)>0.1);
            y_n=y_n(sqrt((x_o-x_c(i-1)).^2+(y_o-y_c(i-1)).^2)>0.1);
        end
    c(i)=max( max(uq_n) );
    x_cc=x_n( uq_n==c(i) );
    y_cc=y_n( uq_n==c(i) );
    x_c(i)=x_cc(1);
    y_c(i)=y_cc(1);
end
end