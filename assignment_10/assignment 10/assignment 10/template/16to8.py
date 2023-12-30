import cv2
import numpy as np

image_16bits = cv2.imread("./terrains/Lugano.png", cv2.IMREAD_UNCHANGED)
normalized = cv2.normalize(image_16bits, None, 0, 255, cv2.NORM_MINMAX)

image_8bits_3channels = np.zeros(
    (image_16bits.shape[0], image_16bits.shape[1], 3), dtype=np.uint8
)
image_8bits_3channels[:, :, 0] = (normalized & 0xFF).astype(np.uint8)
image_8bits_3channels[:, :, 1] = ((image_16bits >> 8) & 0xFF).astype(np.uint8)

cv2.imwrite("output_image.jpg", image_8bits_3channels)

print("Pixel Values of image_8bits_3channels:")
print(image_8bits_3channels)
print(normalized)
