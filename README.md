# Boulder Calculator
This matlab app provides position, size, and volume of boulders based on their contours detected by Detectron2.

1. Detect contours of boulders with Detectron2 and save the results as a JSON file in COCO format using the attached subapp.
2. Load the JSON file, the ortho image used for detection and the corresponding DSM. 
3. Calculate lengths of a axis and b axis, area, orientation, location and volumes.

A simple usage example can be found in 'simple_example.m'.
