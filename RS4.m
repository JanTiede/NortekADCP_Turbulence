function [uw vw b1f b2f b3f b4f wf] = RS4(avt,th,phi1,phi2,phi3,b1,b2,b3,b4,u,v,w)
%4BEAM Calculation of Reynolds stresses uw, vw from a 4 beam ADCP.
%   Radial beam velocities are decomposed into a mean and a flucutating
%   part. The mean part is calculated for periods of user specified amount
%   on minutes avt. 
th = th
phi1 = nanmean(phi1);
phi2 = nanmean(phi2)
phi3 = nanmean(phi3)

%% Noise estimation and removal
noise_slanted1 = 0.0734
noise_slanted2 = 0.0734
noise_slanted3 = 0.0734
noise_slanted4 = 0.0734
b1 = b1 - noise_slanted1;
b2 = b2 - noise_slanted2;
b3 = b3 - noise_slanted3;
b4 = b4 - noise_slanted4;

%% Calculate mean and fluctuations for Beam 1.
npoints = 8*60*avt % 8 Hz * 60s * avt min
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
for i = 1: size(w,2)
brw = w(:,i);
M  = size(brw, 1) - mod(size(brw,1), npoints);
y  = reshape(brw(1:M), npoints, []);
r = transpose(sum(y, 1) / npoints);
wfm(:,i) = repelem(r, npoints);
end
for i = 1: size(u,2)
bru = u(:,i);
M  = size(bru, 1) - mod(size(bru,1), npoints);
y  = reshape(bru(1:M), npoints, []);
r = transpose(sum(y, 1) / npoints);
ufm(:,i) = repelem(r, npoints);
end
for i = 1: size(v,2)
brv = w(:,i);
M  = size(brv, 1) - mod(size(brv,1), npoints);
y  = reshape(brv(1:M), npoints, []);
r = transpose(sum(y, 1) / npoints);
vfm(:,i) = repelem(r, npoints);
end
b1f = b1(1:M,:) - b1m;
b2f = b2(1:M,:) - b2m;
b3f = b3(1:M,:) - b3m;
b4f = b4(1:M,:) - b4m;
uf = u(1:M,:) - ufm;
vf = v(1:M,:) - vfm;
wf = w(1:M,:) - wfm;
uv = uf.^2 .* vf.^2;

%% Calculate u'w' disregarding pitch and roll.
uw = -1/(2*sin(2*th)).*(b2f.^2-b1f.^2);
%% Calculate v'w' disregarding pitch and roll.
vw = -1/(2*sin(2*th)).*(b4f.^2-b3f.^2);
%% Calculate u'w' with small angle approx.
% uw =  -1/(2*sin(2*th)).*(b2f.^2-b1f.^2)+(phi3/sin(th)^2).*(0.5*(b2f.^2+b1f.^2)-wf.^2)-phi2*uv;

%% Calculate v'w' with small angle approx.
% vw =  -1/(2*sin(2*th)).*(b4f.^2-b3f.^2)+(phi2/sin(th)^2).*(0.5*(b4f.^2+b3f.^2)-wf.^2)+phi3*uv;

%% Calculate u'w'.
% uw = (-1/(4*sin(th)^6*cos(th)^2))*sin(th)^5*cos(th)*(b3f.^2-b1f.^2) ...
%       + 2*sin(th)^4*cos(th)^2*phi3*(b3f.^2+b1f.^2) - 4*sin(th)^4*cos(th)^2 ...
%      *phi3*wf.^2  - 4*sin(th)^6*cos(th)^2*phi2*uf.*vf;

%% Calculate v'w'.
% vw = (-1/(4*sin(th)^6*cos(th)^2))*sin(th)^5*cos(th)*(b4f.^2-b2f.^2) ...
%      - 2*sin(th)^4*cos(th)^2*phi2*(b4f.^2+b2f.^2) + 4*sin(th)^4*cos(th)^2 ...
%      *phi2*wf.^2 + 4*sin(th)^6*cos(th)^2*phi3*uf.*vf;

end

