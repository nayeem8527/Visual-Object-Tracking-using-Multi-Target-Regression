pos = boxes;
figure
imshow(img)
hold on
for i=1:size(pos,1)
    rectangle('Position',pos(i,:),'EdgeColor','g');
    hold on
    waitforbuttonpress    
end

% for i=1:size(neg,1)
%     rectangle('Position',neg(i,:),'EdgeColor','r');
%     hold on
% end