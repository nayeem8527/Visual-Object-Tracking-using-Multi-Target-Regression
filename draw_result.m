function drawopt = draw_result(drawopt, fno, frame, bbox)
if (isempty(drawopt))     
  figure('position',[50 100 size(frame,2) size(frame,1)]); clf;                               
  set(gcf,'DoubleBuffer','on','MenuBar','none');%,'Visible','off');

  drawopt.curaxis = [];
  drawopt.curaxis.frm  = axes('position', [0.00 0 1.00 1.0]);
end

% curaxis = drawopt.curaxis;
% axes(curaxis.frm);      
imagesc(frame, [0,255]); 
colormap(gray);  
hold on;     

% x = [bbox(2) bbox(4) bbox(4) bbox(2) bbox(2)];
% y = [bbox(1) bbox(1) bbox(3) bbox(3) bbox(1)];
% line(x, y, 'linewidth', 3, 'color', 'r');
rectangle('Position', bbox, 'linewidth', 3, 'EdgeColor', 'r');
text(10, 15, '#', 'Color','y', 'FontWeight','bold', 'FontSize',24);
text(30, 15, num2str(fno), 'Color','y', 'FontWeight','bold', 'FontSize',24);

% axis equal off;
hold off;
drawnow;    