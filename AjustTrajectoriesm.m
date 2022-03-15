newxc=zeros(size(xcs));
newyc=zeros(size(ycs));
spot_num=size(xcs,2); 
tn=size(xcs,1);
tt=10*[1:tn]';
newxc(1,:)=xcs(1,:);
newyc(1,:)=ycs(1,:);
for i=2:tn
    for j=1:spot_num
        o=sqrt(abs(newxc(i-1,j)-xcs(i,:)).^2+abs(newyc(i-1,j)-ycs(i,:)).^2);
        index=find(o==min(o));
        newxc(i,j)=xcs(i,index);
        newyc(i,j)=ycs(i,index);
    end
end


% tau=1/k3+0.01; N=6; r=0.3125; w= 0.0019;
% tau=1/k3+0.01; N=7; r=0.3593; w= 0.0010;