# creates a video from the pictures in the Images_I folder

import cv2
import numpy as np
import glob

img_array = []
for filename in glob.glob('C:/Users/franz/fenicsproject1/Videomaker/Images_E/*.jpg'):
    img = cv2.imread(filename)
    height, width, layers = img.shape
    size = (width, height)
    img_array.append(img)

out = cv2.VideoWriter('switzerland.avi', cv2.VideoWriter_fourcc(*'DIVX'), 1, size)

for i in range(len(img_array)):
    out.write(img_array[i])
out.release()

print("Finished")
