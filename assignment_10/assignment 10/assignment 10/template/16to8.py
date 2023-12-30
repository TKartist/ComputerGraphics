import cv2
import numpy as np

image_16bits = cv2.imread("./terrains/Lugano.png", cv2.IMREAD_UNCHANGED)
normalized = cv2.normalize(image_16bits, None, 0, 255, cv2.NORM_MINMAX)
image_8bits = cv2.cvtColor(normalized, cv2.COLOR_GRAY2RGB)
cv2.imwrite("output_image.jpg", image_8bits)
