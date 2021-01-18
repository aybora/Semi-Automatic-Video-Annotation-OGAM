function [evaluationResult,clickIndex] = SimulateEvaluateTrack(frameIndex,bbox,foldnumb,RedTxt)

orgIndex = frameIndex;
frameIndex = frameIndex-frameIndex(1)+1;
evaluationResult = zeros(1, frameIndex(end) - frameIndex(1) + 1);
clickIndex = 0;
foldernumber = foldnumb*100000;
deneme2 = RedTxt(:,1);

for p = 1:size(frameIndex,2)
    
    targetvalue1 = orgIndex(p)+foldernumber;
    firstloc2 = find((deneme2==targetvalue1),1,'first');
    
    temparr1 = bbox(orgIndex(p),:);
    if(foldnumb>=10)
        Iss = strcat(num2str(targetvalue1),'.png');
    else
        Iss = strcat('0',num2str(targetvalue1),'.png');
    end
    
    if(isfile(Iss))
        if(rectint(temparr1,temparr1)~=0)
            d31 = temparr1(1);
            d41 = temparr1(3);
            d11 = temparr1(2);
            d21 = temparr1(4);
            
            Is = imread(Iss);
            Is = Is(d11:d21,d31:d41,:);
            
            montagematrix(p) = {imresize(Is,[100 100])};
        end
    end
    
end

M = montage(montagematrix,'Size',[5 9]);
imgOut = M.CData;

boxSizex = round(size(imgOut,1)/5);
boxSizey = round(size(imgOut,2)/9);

figure = imshow(imgOut);

%% loop

while(true)
   
    [x,y,button] = ginput(1);
    x = ceil(x / boxSizex);
    y = ceil(y / boxSizey);
    id = (y - 1) * 9 + x;
    
    clickIndex = clickIndex+1;
    
    evaluationResult(frameIndex(id)) = button;
    if(id~=1)
        for i=(frameIndex(oldid)+1):(frameIndex(id)-1)
            evaluationResult(i) = oldbutton;
        end
    end
    
    oldid = id;
    oldbutton = button;
    
    if(id==size(frameIndex,2))
        break;
    end
    
end

end
