function imret = blendImagePoisson3(im1, im2, roi, targetPosition)

%% Preprocessing
persistent A1 L1 U1 roi0
[m,n,~] = size(im1);
[m1,~,~] = size(im2);

roi(:,[1,2]) = roi(:,[2,1]);
roi = [roi;roi(1,:)];
targetPosition(:,[1,2]) = targetPosition(:,[2,1]);
targetPosition = [targetPosition;targetPosition(1,:)];

%dx dy denote the displacement from roi to target
dx = targetPosition(1,1)-roi(1,1);
dy = targetPosition(1,2)-roi(1,2);
%Find the point in targetPosition
A = repmat([1:m]',[1,n]);
B = repmat([1:n],[m,1]);
in = inpolygon(A,B,targetPosition(:,1),targetPosition(:,2));
%get these point position
[x,y] = find(in == 1);
L0= length(x);
disp(L0);
[judge,location] = ismember([x-1,y;x,y-1;x+1,y;x,y+1],[x,y],'row');
Loc=[[1:L0,1:L0,1:L0,1:L0]',location];
%% Construct Sparse Matrix A1
if isempty(A1)||norm(roi0-roi,'fro') > 1e-6|| numel(roi0)~=numel(roi)
index=1:4*L0; 
M =Loc(index(all(Loc')),:);
A11 = sparse(M(:,1),M(:,2),-1,L0,L0);
A12 = sparse((1:L0)',(1:L0)',4,L0,L0);
A1 = A11 +A12; 
roi0 = roi;
%% LU decompostion
[L1,U1] = lu(A1);
end
%% Construct b
[M1,M2,M3] = deal(im2(:,:,1),im2(:,:,2),im2(:,:,3));
[i0,j0]= deal(round(x-dx),round(y-dy));
% vpq
s1 = 4*double(M1(i0+(j0-1)*m1))-double(M1(i0-1+(j0-1)*m1))-double(M1(i0+1+(j0-1)*m1))-double(M1(i0+(j0-2)*m1))-double(M1(i0+j0*m1));
s2 = 4*double(M2(i0+(j0-1)*m1))-double(M2(i0-1+(j0-1)*m1))-double(M2(i0+1+(j0-1)*m1))-double(M2(i0+(j0-2)*m1))-double(M2(i0+j0*m1));
s3 = 4*double(M3(i0+(j0-1)*m1))-double(M3(i0-1+(j0-1)*m1))-double(M3(i0+1+(j0-1)*m1))-double(M3(i0+(j0-2)*m1))-double(M3(i0+j0*m1));
b1 = [s1,s2,s3];
%f*q
[N1,N2,N3] = deal(im1(:,:,1),im1(:,:,2),im1(:,:,3));
%outside element 1,inside 0
judge =1-reshape(judge,[length(x),4]);
t1 = [double(N1(x-1+(y-1)*m)),double(N1(x+1+(y-1)*m)),double(N1(x+(y-2)*m)),double(N1(x+y*m))];
t2 = [double(N2(x-1+(y-1)*m)),double(N2(x+1+(y-1)*m)),double(N2(x+(y-2)*m)),double(N2(x+y*m))];
t3 = [double(N3(x-1+(y-1)*m)),double(N3(x+1+(y-1)*m)),double(N3(x+(y-2)*m)),double(N3(x+y*m))];
[b21,b22,b23] = deal(sum(t1.*judge,2),sum(t2.*judge,2),sum(t3.*judge,2));
b2 = [b21,b22,b23];
b = b1+b2;
%% LU QR Chol decompostion using '\'
X = U1\(L1\b);
 %Comparing time comsuming, LU perform better than QR(4s to 23s)
%% TODO: compute blended image
imret = im1;
for i =1: length(x)
    imret(x(i),y(i),:)= X(i,:);
end

