function [evaluationResult,clickIndex] = SimulateEvaluateTrack(frameIndex,bbox,orgdata,foldnumb,RedTxt)

orgIndex = frameIndex;
frameIndex = frameIndex-frameIndex(1)+1;
evaluationResult = zeros(1, frameIndex(end) - frameIndex(1) + 1);
clickIndex = 0;
foldernumber = foldnumb*100000;
oldbutton = 3;
deneme1 = orgdata(:,1);
deneme2 = RedTxt(:,1);

for p = 1:size(frameIndex,2)
    
    targetvalue1 = orgIndex(p)+foldernumber;
    firstloc1 = find((deneme1==targetvalue1),1,'first');
    firstloc2 = find((deneme2==targetvalue1),1,'first');
    
    if((orgdata(firstloc1,2)==0) & ((size(firstloc1,1)~=0) && (size(firstloc2,1)~=0)))
        evaluationResult(frameIndex(p)) = 3;
        if(p~=1)
            for i=(frameIndex(p-1)+1):(frameIndex(p)-1)
                evaluationResult(i) = oldbutton;
            end
            
            if((oldbutton ~= 3) && (p ~= size(frameIndex,2)))
                clickIndex = clickIndex+1;
            end
        end
        
        
    elseif(size(firstloc2,1)==0)
        evaluationResult(frameIndex(p)) = oldbutton;
        if(p~=1)
            for i=(frameIndex(p-1)+1):(frameIndex(p)-1)
                evaluationResult(i) = oldbutton;
            end
        end
        
    else
        matrixORG = orgdata(firstloc1,[3:6]);
        areaORG = rectint(matrixORG,matrixORG);
        temparr1 = bbox(orgIndex(p),:);
        d1 = temparr1(1);
        d3 = temparr1(3)-temparr1(1);
        d2 = temparr1(2);
        d4 = temparr1(4)-temparr1(2);
        
        matrixYOLO = [d1 d2 d3 d4];
        areaYOLO = rectint(matrixYOLO,matrixYOLO);
        
        if(size(firstloc1,1)==0)
            IOU = 0;
        else
            IOU = rectint(matrixORG,matrixYOLO);
            IOU = IOU/(areaORG+areaYOLO-IOU);
        end
        
        if(IOU>0.5)
            evaluationResult(frameIndex(p))=1;
            if(p~=1)
                for i=(frameIndex(p-1)+1):(frameIndex(p)-1)
                    evaluationResult(i) = oldbutton;
                end
                
                if((oldbutton ~= 1) && (p ~= size(frameIndex,2)))
                    clickIndex = clickIndex+1;
                end
            end
            
        elseif(IOU == 0)
            evaluationResult(frameIndex(p))= 3;
            if(p~=1)
                for i=(frameIndex(p-1)+1):(frameIndex(p)-1)
                    evaluationResult(i) = oldbutton;
                end
                
                if((oldbutton ~= 3) && (p ~= size(frameIndex,2)))
                    clickIndex = clickIndex+1;
                end
            end 
            
        elseif(IOU <= 0.5)
            evaluationResult(frameIndex(p))= 2;
            if(p~=1)
                for i=(frameIndex(p-1)+1):(frameIndex(p)-1)
                    evaluationResult(i) = oldbutton;
                end
                
                if((oldbutton ~= 2) && (p ~= size(frameIndex,2)))
                    clickIndex = clickIndex+1;
                end
            end
        end
        
        %
        
    end
    
    if((p == 1)||(p == size(frameIndex,2)))
        clickIndex = clickIndex+1;
    end
    
    oldbutton = evaluationResult(frameIndex(p));
    
end

end
