function [feats] = feat_extractor_hog(img,boxes,bbox)
feats = [];
for i = 1:size(boxes,1)
    sample = imcrop(img,boxes(i,:));
    sample = imresize(sample,[90,80]);
    %% hog features
    hog = fhog(single(sample));
%     hog = double(vl_hog(im2single(sample), 8));
%     feat.feaArr = [feat.feaArr double(hog)];
    temp = hog(:);
    temp = temp(1:3410,:);
    feats = [feats temp];
end
end
