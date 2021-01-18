AUTO CLICK

1. copy the 'reals.mat' file that contains ground truth object locations into the dataset folder.
2. copy the MATLAB files in either IoU0.5 or IoU0.2 depending on which you want to work with, into the same dataset folder.
3. Open 'FinalBoss' MATLAB file and then in the 5th line write the number of the video that you want to work with.
4. In the 6th line write the name of the batch of detections. This batch should also be located inside the folder.
5. Run the file

Outputs:
(for the example video number 2)
- Video2FinalRedText.txt contains the detections labeled as hit after the auto clicking has been done.
- Video2FinalMissText.txt contains the detections labeled as miss after the auto clicking has been done.
- Video2PerformanceText.txt contains the performance test done after the auto clicking operation has been done. This file also contains the number of clicks done during auto clicking.
- Video2BeforeClickPerformanceText.txt contains the performance test done before the auto clicking operation.


MANUAL CLICKING

1. copy the MATLAB files in Manual Click into the dataset folder.
2. Open 'FinalBoss' MATLAB file and then in the 5th line write the number of the video that you want to work with.
3. In the 6th line write the name of the batch of detections. This batch should also be located inside the folder.
4. Run the 'FinalBoss' MATLAB file

Operations:
- Each track will present you 43 of its detections (less if track size is below 43).
- Left click on the frame means it is a hit, right click on the frame means it is a miss, middle click on the frame means it is a weak hit.
- First and last frames must be clicked at all times.
- Example: 1st frame is left clicked, and then 20th frame is middle clicked, all frames from 1st frame to 20th frame will be counted as left clicked.

Outputs:
(for the example video number 2)
- Video2FinalRedText.txt contains the detections labeled as hit after the auto clicking has been done.
- Video2FinalMissText.txt contains the detections labeled as miss after the auto clicking has been done.