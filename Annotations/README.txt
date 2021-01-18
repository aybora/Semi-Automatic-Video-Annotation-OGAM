Complete uav_detection_2 dataset and not annotated videos of uav_detection dataset is annotated with the algorithm.
Annotation flags have a small difference than the original one: We have used flags instead of number_of_drones.

frame_number, flag, x_min, y_min, width, height

flag: 

0: There is one drone present in the frame.
-1: There is no drone present in the frame, you can feed it as a negative example to your network.
-2: Frame is blurry/drone is partially visible in the frame, you should not feed it during training, or count it while evaluating performance.

