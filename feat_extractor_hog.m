function [feats] = feat_extractor_hog(I,boxes)
feats = [];
for j=1:size(I,3)
    img = I(:,:,j);
for i = j%:size(boxes,1)
    sample = imcrop(img,boxes(i,:));
    sample = imresize(sample,[32,32]);
    %% hog features
    hog = fhog(single(sample));
%     hog = double(vl_hog(im2single(sample), 8));
%     feat.feaArr = [feat.feaArr double(hog)];
    temp = hog(:);
    temp = temp(1:496,:);
    feats = [feats temp];
end
end
end
