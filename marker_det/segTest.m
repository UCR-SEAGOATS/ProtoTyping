files = dir('../dataSets/');
vid = VideoReader('../DataSets/marker_indoor1.mp4');
results = VideoWriter('results.avi');
close all
figure
j = 1;
i = 1;
open(results);

while hasFrame(vid)
    if j > 30
        frame = readFrame(vid);
        J = shapeSegmentation(frame);

        RGB = double(cat(3, J, J, J));
        writeVideo(results, [RGB]);
        j = 0;
        i = i+1
    end
    
    j = j + 1;
%     subplot(1,2,1)
%     imshow(frame)
%     subplot(1,2,2)
%     imshow(J)
%     drawnow
end

close(results)