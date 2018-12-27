function [T, TW] = get_template(frame, bb, T, TW, err, trackpars,ppp)
%
img = frame;
[h, w] = size(img);
% tmp = double(imresize(img(max(0,bb(1)):min(bb(3),h), max(0,bb(2)):min(bb(4),w)), trackpars.nsize));
% tmp = double(vl_hog(im2single(tmp), 8));

% tmp = double(imresize(imcrop(img,bb), trackpars.nsize));
% tmp = fhog(single(tmp));
% tt = tmp(:);
% tt = tt(1:496,:);
tt = feat_extractor_hog(img,bb,ppp);


% T = [T  tmp(:) / sqrt(sum(tmp(:) .* tmp(:)))];
T = [T  tt / sqrt(sum(tt .* tt))];
TW = [TW exp(-err)];
if size(T,2) > trackpars.lengthT
    T = T(:,2:end);
    TW = TW(2:end);
end
