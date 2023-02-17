clear
close all
clc

Delta = 0.023;
delta = 0.012; %[seconds]
Gamma = 2.675987E8; % rad/s/T
D0 = 3000; % [micrometer^2/seconds]
l_delta = sqrt(D0*delta);

bval = load('20211104_114911LRTE55bmax15000multipledirectionss004a001.bval');
bvec = load('20211104_114911LRTE55bmax15000multipledirectionss004a001.bvec');
S_img = load_untouch_nii('real_4_11_2021_multi_gibbsCorrSubVoxShift_driftCo_TED.nii.gz');
S_img = double(S_img.img);
S_img = S_img(end:-1:1,:,:,:);
S_img(52,36:38,38,:) = 0;

mask = load_untouch_nii('real_4_11_2021_multi_gibbsCorrSubVoxShift_driftCo_TED_brain_mask.nii.gz');
mask = double(mask.img);
mask = mask(end:-1:1,:,:,:);
se = strel('disk',2);        
mask = imerode(mask,se);

[b_sort, order] = sort(bval);
[C, ia, ic] = unique(b_sort);

q = sqrt(C/10^6/(Delta - delta/3));
l_g = (D0*delta./q).^(1/3);

p = length(C);
S0 = S_img(:,:,:,bval == 0);

S_shell = zeros(size(S_img,1),size(S_img,2),size(S_img,3),p-1);

S_img_sort = S_img(:,:,:,order);
for kk = 1:p-1
    S_shell(:,:,:,kk) = mean(S_img_sort(:,:,:,ia(kk):ia(kk+1)-1),4).*mask;
end
S_shell(:,:,:,p) = mean(S_img_sort(:,:,:,ia(p):end),4).*mask;
S_shell(S_shell<0) = 10^-3;

S_norm = zeros(size(S_shell,1),size(S_shell,2),size(S_shell,3),size(S_shell,4));

for k = 1:size(S_shell,4)
    S_norm(:,:,:,k) = S_shell(:,:,:,k).*mask./(S_shell(:,:,:,1));
end

S_norm(isnan(S_norm)) = 10^-3;
S_norm(S_norm == inf) = 10^-3;
S_norm(S_norm >1) = 1;
S_norm(S_norm <0) = 10^-3;

[s1, s2, s3, s4] = size(S_norm);

opts = optimset('Display','off');
lb = [0.01, 0.01,  0.01, 0.01, 0.01];
ub = [20, 20,  5, 5, 5];
x1 = [6, 3, 2,  2, 2];

xmap = zeros(s1,s2,s3,5);
for i = 1:s1
    for j = 1:s2
        for k = 1:s3
            if mask(i,j,k)>0
            signal = squeeze(S_norm(i,j,k,:));
xmap(i,j,k,:) = lsqcurvefit(@localization_code, x1, q', signal ,lb,ub,opts);
            end
        end
    end
    i
end

q1 = 0:0.01:1.36;

xdata(:,1) = q;
xdata(1,2) = Delta;
xdata(2,2) = delta;

opts = optimset('Display','off');
lb = [0.01, 0.01];
ub = [3, 2];
x1 = [2, 1];

map_d_k = zeros(s1,s2,s3,2);
for i = 1:s1
    for j = 1:s2
        for k = 1:s3
            if mask(i,j,k)>0
            signal = squeeze(S_norm(i,j,k,:));
map_d_k(i,j,k,:) = lsqcurvefit(@kurtosis_d, x1, xdata(1:5,:), signal(1:5,:) ,lb,ub,opts);

            end
        end
    end
    i
end

kk = 37;
figure,
ha = tight_subplot(3,5,[0.01 0.02],[0.01 0.05],[.01 .03]) ; % [between], [down, up], [left, right]
axes(ha(1))
imagesc(rot90(flipud(xmap(:,:,kk,1))),[0,13]), colormap('gray'), axis off, axis equal, colorbar
set(gca, 'FontSize', 36)
title( 'r [\mum]' , 'FontWeight','normal')

axes(ha(2))
imagesc(rot90(flipud(xmap(:,:,kk,2))),[0,4]), colormap('gray'), axis off, axis equal, colorbar
set(gca, 'FontSize', 36)
title( 'u [\mum]' , 'FontWeight','normal')

axes(ha(3))
imagesc(rot90(flipud(xmap(:,:,kk,3))),[0,3]), colormap('gray'), axis off, axis equal, colorbar
set(gca, 'FontSize', 36)
title( 'v [\mum]' , 'FontWeight','normal')

axes(ha(4))
imagesc(rot90(flipud(xmap(:,:,kk,4))),[0,3]), colormap('gray'), axis off, axis equal, colorbar
set(gca, 'FontSize', 36)
title( 'w [\mum]' , 'FontWeight','normal')

axes(ha(5))
imagesc(rot90(flipud(xmap(:,:,kk,5))),[0,3]), colormap('gray'), axis off, axis equal, colorbar
set(gca, 'FontSize', 36)
title( 'p [\mum]' , 'FontWeight','normal')

axes(ha(6))
axis off

axes(ha(7))
axis off

axes(ha(8))
axis off

axes(ha(9))
axis off

axes(ha(10))
axis off

r1 = xmap(:,:,kk,1);
u1 = xmap(:,:,kk,2);
v1 = xmap(:,:,kk,3);
w1 = xmap(:,:,kk,4);
p1 = xmap(:,:,kk,5);

r = xmap(:,:,:,1);
u = xmap(:,:,:,2);
v = xmap(:,:,:,3);
w = xmap(:,:,:,4);
p = xmap(:,:,:,5);

td = Delta - delta/3;
td = td*1000;
D_new = r.^(2/3) .* u.^(4/3) / td;

coeff_K = 1/6 * D_new.^2 .* td^2; 

K_new = r.^(2/3) .* v.^(10/3) ./ coeff_K;
K_new(isnan(K_new)) = 0;

figure,
kk = 37;
ha = tight_subplot(3,4,[0.07 0.01],[0.01 0.07],[.04 .01]) ; % [between], [down, up], [left, right]
axes(ha(1))
imagesc(rot90(flipud(map_d_k(:,:,kk,1))),[0,3]), colormap('gray'), axis off, axis equal, colorbar
set(gca, 'FontSize', 36)
title( 'D [\mum^2/ms]' , 'FontWeight','normal')

axes(ha(5))
imagesc(rot90(flipud(D_new(:,:,kk))),[0,3]), colormap('gray'), axis off, axis equal, colorbar
set(gca, 'FontSize', 36)

axes(ha(2))
imagesc(rot90(flipud(map_d_k(:,:,kk,2))),[0,2]), colormap('gray'), axis off, axis equal, colorbar
set(gca, 'FontSize', 36)
title( 'K' , 'FontWeight','normal')

axes(ha(6))
imagesc(rot90(flipud(K_new(:,:,kk))),[0,2]), colormap('gray'), axis off, axis equal, colorbar
set(gca, 'FontSize', 36)

axes(ha(7))
set(gca, 'FontSize', 36)
title( 'DKI' , 'FontWeight','normal')

axes(ha(8))
set(gca, 'FontSize', 36)
title( 'Proposal' , 'FontWeight','normal')

figure,
kk = 37;
ha = tight_subplot(3,4,[0.07 0.01],[0.01 0.07],[.04 .01]) ; % [between], [down, up], [left, right]
axes(ha(1))
imagesc(rot90(flipud(map_d_k(:,:,kk,1))-D_new(:,:,kk)),[-3,3]), axis off, axis equal, colorbar
set(gca, 'FontSize', 36)
addpath("./cbrewer");
colormap(brewermap([], '-RdBu'));

axes(ha(2))
imagesc(rot90(flipud(map_d_k(:,:,kk,2)-K_new(:,:,kk))),[-1,1]), axis off, axis equal, colorbar
set(gca, 'FontSize', 36)
addpath("./cbrewer");
colormap(brewermap([], '-RdBu'));

axes(ha(7))
set(gca, 'FontSize', 36)
title( 'DKI - Proposal' , 'FontWeight','normal')

D_dki = map_d_k(:,:,kk,1);
D_pro = D_new(:,:,kk);

K_dki = map_d_k(:,:,kk,2);
K_pro = K_new(:,:,kk);

figure,
kk = 37;
ha = tight_subplot(3,4,[0.07 0.01],[0.01 0.07],[.04 .01]) ; % [between], [down, up], [left, right]
axes(ha(1))
plot(D_dki(:),D_pro(:),'.'), axis equal, ylim([0,3]), xlim([0,3]), ylabel('Proposal') , xlabel('DKI')
set(gca, 'FontSize', 36)

axes(ha(2))
plot(K_dki(:),K_pro(:),'.'), axis equal, ylim([0,1]), xlim([0,1]),  xlabel('DKI')
set(gca, 'FontSize', 36)

axes(ha(5))
axis off

axes(ha(6))
axis off

figure,
ha = tight_subplot(2,5,[0.1 0.07],[0.01 0.01],[.05 .01]) ; % [between], [down, up], [left, right]
axes(ha(1))
plot(u1(:),r1(:),'.'), axis equal, ylim([0,13]), xlim([0,13]),  xlabel('u [\mum]'), ylabel('r [\mum]')
set(gca, 'FontSize', 36)

axes(ha(2))
plot(v1(:),r1(:),'.'), axis equal, ylim([0,13]), xlim([0,13]),  xlabel('v [\mum]'), ylabel('r [\mum]')
set(gca, 'FontSize', 36)

axes(ha(3))
plot(w1(:),r1(:),'.'), axis equal, ylim([0,13]), xlim([0,13]),  xlabel('w [\mum]'), ylabel('r [\mum]')
set(gca, 'FontSize', 36)

axes(ha(4))
plot(p1(:),r1(:),'.'), axis equal, ylim([0,13]), xlim([0,13]),  xlabel('p [\mum]'), ylabel('r [\mum]')
set(gca, 'FontSize', 36)

axes(ha(5))
plot(v1(:),u1(:),'.'), axis equal, ylim([0,13]), xlim([0,13]),  xlabel('v [\mum]'), ylabel('u [\mum]')
set(gca, 'FontSize', 36)

axes(ha(6))
plot(w1(:),u1(:),'.'), axis equal, ylim([0,13]), xlim([0,13]),  xlabel('w [\mum]'), ylabel('u [\mum]')
set(gca, 'FontSize', 36)

axes(ha(7))
plot(p1(:),u1(:),'.'), axis equal, ylim([0,13]), xlim([0,13]),  xlabel('p [\mum]'), ylabel('u [\mum]')
set(gca, 'FontSize', 36)

axes(ha(8))
plot(w1(:),v1(:),'.'), axis equal, ylim([0,5]), xlim([0,5]),  xlabel('w [\mum]'), ylabel('v [\mum]')
set(gca, 'FontSize', 36)

axes(ha(9))
plot(p1(:),v1(:),'.'), axis equal, ylim([0,5]), xlim([0,5]),  xlabel('p [\mum]'), ylabel('v [\mum]')
set(gca, 'FontSize', 36)

axes(ha(10))
plot(p1(:),w1(:),'.'), axis equal, ylim([0,5]), xlim([0,5]),  xlabel('p [\mum]'), ylabel('w [\mum]')
set(gca, 'FontSize', 36)
