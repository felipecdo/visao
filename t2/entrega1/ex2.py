import math
import cv2
import numpy as np
import matplotlib.pyplot as plt

import siamxt   

def image_area_open_filter(img_neg, neighborhood, area):
    mxt_neg = siamxt.MaxTreeAlpha(img_neg, neighborhood)
    mxt_neg.areaOpen(area)
    return mxt_neg.getImage()
    
original_img = cv2.imread('../EP2/fruit.png', cv2.IMREAD_GRAYSCALE)
original_img = original_img.astype(np.uint16)

fig, imgs = plt.subplots(1,1, figsize=(15,15))
imgs.imshow(original_img, cmap='gray')
plt.axis('off')
plt.title('Imagem Original')
plt.show()

# Structuring element with connectivity-8
neighborhood = np.ones((3,3),dtype = bool)

# Negating the image
img_max = original_img.max()
img_neg = img_max - original_img

areas = [10,100,1000,10000]

filtered_images = []

for area in areas:
    img_filtered = image_area_open_filter(img_neg, neighborhood, area)
    img_filtered = img_max - img_filtered
    filtered_images.append(img_filtered)
    
# Displaying all filtered image
fig, imgs = plt.subplots(2,2, figsize=(15,12))
plt.axis('off')
for i in range(len(filtered_images)):
    imgs[i/2,i%2].imshow(filtered_images[i][0:200,100:350], cmap='gray')
    imgs[i/2,i%2].axis('off')
    imgs[i/2,i%2].set_title('Poda com area = '+ str(areas[i]))

plt.show()

# Show the best parameter result
img_filtered =  img_filtered = image_area_open_filter(img_neg, neighborhood, 100)
img_filtered = img_max - img_filtered

#Displaying the filtered image
fig, imgs = plt.subplots(1,1, figsize=(15,15))
plt.axis('off')
plt.title('Imagem Resultado')
imgs.imshow(img_filtered, cmap='gray')

plt.show()