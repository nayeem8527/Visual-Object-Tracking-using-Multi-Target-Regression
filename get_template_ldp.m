function [T, TW] = get_template_ldp(frame, bb, T, TW, err, trackpars)
%
img = frame;
[h, w] = size(img);
tmp = double(imresize(img(max(0,bb(1)):min(bb(3),h), max(0,bb(2)):min(bb(4),w)), trackpars.nsize));

options.scale         = [1,2,3];
options.deltax        = 12;
options.deltay        = 12;
options.n             = 2;
options.ldporient     = [0 , 1 , 2 , 3];
options.color         = 0;
options.patchsize     = 20;
options.rmextremebins = 1;
options.norm          = 4;
options.clamp         = 0.2;

% tmp = double(vl_hog(im2single(tmp), 8));
% tmp = extractLBPFeatures(uint8(tmp));
[tmp , ~] =  denseMBLDP(uint8(tmp) , options );
tmp = tmp';
[coeff,~] = pca(tmp);
coeff = coeff(:,1:300);
tmp = tmp*coeff;
tmp = sum(tmp,2);
tmp = (tmp-min(tmp))/(max(tmp)-min(tmp));
%     ind = find(hog==0);
%     hog(ind(1:size(ind,1)))=0.001;
% tmp = sum(tmp,2);
% tmp = (tmp-min(tmp))/(max(tmp)-min(tmp));
% ind = find(tmp==0);
% tmp(ind(1:size(ind,1)))=0.001;
T = [T  tmp(:) / sqrt(sum(tmp(:) .* tmp(:)))];
TW = [TW exp(-err)];
if size(T,2) > trackpars.lengthT
    T = T(:,2:end);
    TW = TW(2:end);
end
