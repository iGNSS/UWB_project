function theta = chan2d(BS,range,noise)
%BS为基站的坐标；
%range为距离；
%noise为噪声误差；
[BSN,~] = size(BS);%基站的个数；
Q = diag(ones(BSN,1)*noise);
k = sum(BS.^2,2);
h1 = range.^2 - k;
G1 = [-2*BS,ones(BSN,1)];
theta0 = (G1'/Q*G1)\G1'/Q*h1;%LS；

dis_s = (sum((BS - ones(BSN,1)*theta0(1:2,:)').^2,2)).^(1/2);%计算真实距离；
B1 = diag(dis_s);
cov1 = B1*Q*B1;%测量误差和测量距离的协方差矩阵；
theta1 = (G1'/cov1*G1)\G1'/cov1*h1;%第一次WLS估计值Za；

h2 = [theta1(1,1)^2; theta1(2,1)^2; theta1(3,1)];
G2 = [1 0; 0 1; 1 1];
B2 = diag([theta1(1,1); theta1(2,1); 1/2]);
cov_theta1 = G1'/cov1*G1; %扰动法;
cov2 = 4*B2/cov_theta1*B2;
theta2 = (G2'/cov2*G2)\G2'/cov2*h2;%第二次WLS估计值Zp；

theta = sign(theta1(1:2,:)).*(abs(theta2).^(1/2));
theta = round(theta)';
end