function [outputIndex] = PrepareTrackEvaluation(startFrame, endFrame, correlationScore, objectnessScore, avgObjectnessScore, displacement)

if((endFrame-startFrame)>42)

N = 8;

keyFrames = round(linspace(startFrame, endFrame, N));

outputIndex = zeros(1, (N - 1) * 5 + N);

for n = 2 : N
    
    partStart = keyFrames(n - 1);
    partEnd = keyFrames(n);
    
    % CC
    [~,minCorrIND] = min(correlationScore(partStart + 1 : partEnd - 1));
    
    % OBJ
    objectnessScore(partStart + minCorrIND) = inf;
    avgObjectnessScore(partStart + minCorrIND) = inf;
    displacement(partStart + minCorrIND) = 0;
    [~,minObjIND] = min(objectnessScore(partStart + 1 : partEnd - 1));
    
    % AVG OBJ
    avgObjectnessScore(partStart + minObjIND) = inf;
    displacement(partStart + minObjIND) = 0;
    [~,minAvgObjIND] = min(avgObjectnessScore(partStart + 1 : partEnd - 1));
    
    % DISP
    displacement(partStart + minAvgObjIND) = 0;
    [~,maxDispIND] = max(displacement(partStart + 1 : partEnd - 1));

    % INDEXES
    indexes = sort([partStart, partEnd, partStart + minCorrIND, partStart + minObjIND, partStart + minAvgObjIND, partStart + maxDispIND], 'ascend');
    
    indexDiff = diff(indexes);
    
    [maxVal, maxIND] = max(indexDiff);
    
    % correct indexes
    outputIndex((n - 2) * 6 + 1 : (n - 1)*6+1) = sort([indexes, round(indexes(maxIND) + maxVal / 2)], 'ascend');
     
end

elseif((endFrame-startFrame)<43)

    outputIndex = startFrame:endFrame;
    
end

end