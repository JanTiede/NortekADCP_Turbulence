 function [uw vw u2 v2 w2 q b1f b1m b2f b2m b3f b3m b4f b4m b5f b5m b10f uf vf uv um vm] = RS5(no,avt,th,phi1,phi2,phi3,b1,b2,b3,b4,b5,u,v)
%     Calculation of Reynolds stresses uw, vw from a 5 beam ADCP.
%     Radial beam velocities are decomposed into a mean and a flucutating
%     part. The mean part is calculated for periods of user specified amount
%     on minutes avt. 

%% Noise estimation and removal
noise_slanted1 = 0.0734
noise_slanted2 = 0.0734
noise_slanted3 = 0.0734
noise_slanted4 = 0.0734
noise_vert = 0.0734

%% Calculate mean and fluctuations for Beam 1.
phi1 = nanmean(phi1)
phi2 = nanmean(phi2)
phi3 = nanmean(phi3)
npoints = 8*60*avt
for i = 1: size(b1,2)
br1 = b1(:,i);
M  = size(br1, 1) - mod(size(br1,1), npoints);
y  = reshape(br1(1:M), npoints, []);
r = transpose(sum(y, 1) / npoints);
b1m(:,i) = repelem(r, npoints);
end
for i = 1: size(b2,2)
br2 = b2(:,i);
M  = size(br2, 1) - mod(size(br2,1), npoints);
y  = reshape(br2(1:M), npoints, []);
r = transpose(sum(y, 1) / npoints);
b2m(:,i) = repelem(r, npoints);
end
for i = 1: size(b3,2)
br3 = b3(:,i);
M  = size(br3, 1) - mod(size(br3,1), npoints);
y  = reshape(br3(1:M), npoints, []);
r = transpose(sum(y, 1) / npoints);
b3m(:,i) = repelem(r, npoints);
end
for i = 1: size(b4,2)
br4 = b4(:,i);
M  = size(br4, 1) - mod(size(br4,1), npoints);
y  = reshape(br4(1:M), npoints, []);
r = transpose(sum(y, 1) / npoints);
b4m(:,i) = repelem(r, npoints);
end
for i = 1: size(b5,2)
b5r = b5(:,i);
M  = size(b5r, 1) - mod(size(b5r,1), npoints);
y  = reshape(b5r(1:M), npoints, []);
r = transpose(sum(y, 1) / npoints);
b5m(:,i) = repelem(r, npoints);
end
for i = 1: size(u,2)
ur = u(:,i);
M  = size(ur, 1) - mod(size(ur,1), npoints);
y  = reshape(ur(1:M), npoints, []);
r = transpose(sum(y, 1) / npoints);
um(:,i) = repelem(r, npoints);
end
for i = 1: size(v,2)
vr = v(:,i);
M  = size(vr, 1) - mod(size(vr,1), npoints);
y  = reshape(vr(1:M), npoints, []);
r = transpose(sum(y, 1) / npoints);
vm(:,i) = repelem(r, npoints);
end
b1f = b1(1:M,:) - b1m;
b10f = b1f;
b2f = b2(1:M,:) - b2m;
b3f = b3(1:M,:) - b3m;
b4f = b4(1:M,:) - b4m;
b5f = b5(1:M,:) - b5m;
uf = u(1:M,:) - um;
vf = v(1:M,:) - vm;
uv = uf.^2 .* vf.^2;

%% Filter for waves and other
% f = 8;
% T = 1/f;
% % frequencies = [0.1 0.9];
% [b,a] = butter(5, 0.09./(f/2),'high');
% filterButter = tf(b,a,T);
% figure;
% bode(filterButter);
% 
% t = 1:length(b1f);
% for i = 1:26
% y = double(b1f(:,i));
% y(isnan(y)) = 0;
% x(:,i) = y;
% b1f(:,i) = filtfilt(b,a,x(:,i));
% end
% for i = 1:26
% y = double(b2f(:,i));
% y(isnan(y)) = 0;
% x(:,i) = y;
% b2f(:,i) = filtfilt(b,a,x(:,i));
% end
% for i = 1:26
% y = double(b3f(:,i));
% y(isnan(y)) = 0;
% x(:,i) = y;
% b3f(:,i) = filtfilt(b,a,x(:,i));
% end
% for i = 1:26
% y = double(b4f(:,i));
% y(isnan(y)) = 0;
% x(:,i) = y;
% b1f(:,i) = filtfilt(b,a,x(:,i));
% end
% for i = 1:26
% y = double(b5f(:,i));
% y(isnan(y)) = 0;
% x(:,i) = y;
% b5f(:,i) = filtfilt(b,a,x(:,i));
% end
%%
for i = 1: size(b1f,2)
br1 = b1f(:,i).^2;
M  = size(br1, 1) - mod(size(br1,1), npoints);
y  = reshape(br1(1:M), npoints, []);
r = transpose(sum(y, 1) / npoints);
b1fm(:,i) = repelem(r, npoints);
end
for i = 1: size(b2f,2)
br2 = b2f(:,i).^2;
M  = size(br2, 1) - mod(size(br2,1), npoints);
y  = reshape(br2(1:M), npoints, []);
r = transpose(sum(y, 1) / npoints);
b2fm(:,i) = repelem(r, npoints);
end
for i = 1: size(b3f,2)
br3 = b3f(:,i).^2;
M  = size(br3, 1) - mod(size(br3,1), npoints);
y  = reshape(br3(1:M), npoints, []);
r = transpose(sum(y, 1) / npoints);
b3fm(:,i) = repelem(r, npoints);
end
for i = 1: size(b4f,2)
br4 = b4f(:,i).^2;
M  = size(br4, 1) - mod(size(br4,1), npoints);
y  = reshape(br4(1:M), npoints, []);
r = transpose(sum(y, 1) / npoints);
b4fm(:,i) = repelem(r, npoints);
end
for i = 1: size(b5f,2)
b5r = b5f(:,i).^2;
M  = size(b5r, 1) - mod(size(b5r,1), npoints);
y  = reshape(b5r(1:M), npoints, []);
r = transpose(sum(y, 1) / npoints);
b5fm(:,i) = repelem(r, npoints);
end
for i = 1: size(uf,2)
ur = uf(:,i);
M  = size(ur, 1) - mod(size(ur,1), npoints);
y  = reshape(ur(1:M), npoints, []);
r = transpose(sum(y, 1) / npoints);
ufm(:,i) = repelem(r, npoints);
end
for i = 1: size(vf,2)
vr = vf(:,i);
M  = size(vr, 1) - mod(size(vr,1), npoints);
y  = reshape(vr(1:M), npoints, []);
r = transpose(sum(y, 1) / npoints);
vfm(:,i) = repelem(r, npoints);
end

%% Calculate uw and vw for ADCP without pitch and roll.
% uw = (1/(2*sin(2*th))).*(b2fm-b1fm);
% vw = (1/(2*sin(2*th))).*(b4fm-b3fm);

%% Calculate tangential and normal stress (Dewey & Stringer).
uw = (-1/(4*sin(th)^6*cos(th)^2))*sin(th)^5*cos(th)*(b3f.^2-b1f.^2) ...
      + 2*sin(th)^4*cos(th)^2*phi3*(b3f.^2+b1f.^2) - 4*sin(th)^4*cos(th)^2 ...
     *phi3*b5f.^2  - 4*sin(th)^6*cos(th)^2*phi2*uf.*vf;

vw = (-1/(4*sin(th)^6*cos(th)^2))*sin(th)^5*cos(th)*(b4f.^2-b2f.^2) ...
     - 2*sin(th)^4*cos(th)^2*phi2*(b4f.^2+b2f.^2) + 4*sin(th)^4*cos(th)^2 ...
     *phi2*b5f.^2 + 4*sin(th)^6*cos(th)^2*phi3*uf.*vf;

u2 = (-1/(4*sin(th)^6*cos(th)^2))*(-2*sin(th)^4*cos(th)^2*(b3f.^2+b1f.^2 ...
    -2*cos(th)^2*b5f.^2))+2*sin(th)^5*cos(th)*phi3*(b3f.^2-b1f.^2);

v2 = (-1/(4*sin(th)^6*cos(th)^2))*-2*sin(th)^4*cos(th)^2*(b4f.^2+b2f.^2 ...
    -2*cos(th)^2*b5f.^2)-2*sin(th)^4*cos(th)^2*phi3*(b3f.^2-b1f.^2) ...
    +2*sin(th)^3*cos(th)^3*phi3*(b3f.^2-b1f.^2)-2*sin(th)^5*cos(th)*phi2*(b4f.^2-b2f.^2);

w2 = (-1/(4*sin(th)^6*cos(th)^2))*-2*sin(th)^5*cos(th)*(b3f.^2-b1f.^2) ...
    +2*sin(th)^5*cos(th)*phi2*(b4f.^2-b2f.^2)-4*sin(th)^6*cos(th)^2*b5f.^2;

q = (((1/sin(4*sin(th)^2))*((b1f.^2+b2f.^2+b3f.^2+b4f.^2)-2*(2*cos(th)^2)*b5f.^2) - (cot(th)-1)*phi3*(b2f.^2-b1f.^2))*2);
q = (q);
end

