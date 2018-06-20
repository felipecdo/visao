
import math
import cv2
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage 

import siamxt   

original_img = cv2.imread('../EP2/knee.pgm')
original_img = original_img.astype(np.uint8)

fig, imgs = plt.subplots(1,1, figsize=(15,15))
imgs.imshow(original_img, cmap='gray')
plt.title('Imagem Original')
plt.show()

img_gray = cv2.cvtColor(original_img, cv2.COLOR_BGR2GRAY)

# Rounded Structuring Element
struct = np.array([[False, True, False],
                   [True, True, True],
                   [False, True, False]])

dilation = ndimage.grey_dilation(img_gray, structure=struct)
erosion = ndimage.grey_erosion(img_gray, structure=struct)

gradient_img = dilation - erosion

# Neighborhood connectivity-8
neighborhood = np.ones((3,3),dtype = bool)

# Number of leaves to be preserved
n = 8

max_from_gradient = gradient_img.max()
mxt_gradient = siamxt.MaxTreeAlpha(max_from_gradient - gradient_img, neighborhood)

# Computes volume attribute of the max-tree nodes
area = mxt_gradient.node_array[3,:]
height = mxt_gradient.computeHeight()
volume = area * height

# Apply extinct filter
volume_ext = mxt_gradient.computeExtinctionValues(volume,"volume")
mxt_gradient.extinctionFilter(volume_ext, n)

#Recovering the image 
img_gradient_filtered =  max_from_gradient - mxt_gradient.getImage()

#Displaying the filtered image
fig, imgs = plt.subplots(1,2, figsize=(15,15))
imgs[0].imshow(gradient_img, cmap='gray')
imgs[0].set_title('Imagem Resultado do Gradiente Morfologico')
imgs[1].imshow(img_gradient_filtered, cmap='gray')
imgs[1].set_title('Imagem filtrada por extincao pelo atributo de volume')
plt.show()


# Negating again the image (build max tree again)
max_gradient = img_gradient_filtered.max()
img_gradient_filtered_neg = max_gradient - img_gradient_filtered

mxt_neg = siamxt.MaxTreeAlpha(img_gradient_filtered_neg, neighborhood)

#Selecting nodes that fit the criteria (only leaves)
nodes = (mxt_neg.node_array[1,:] == 0)

#Filtering
mxt_neg.contractDR(nodes)

img_filtered = mxt_neg.getImage()

img_map = np.where(img_filtered > 0,1,0)
img_map = img_map.astype(np.uint8)

# Marker labelling
ret, markers = cv2.connectedComponents(img_map)

colors = [
    [0,0,0],
    [0,0,1],
    [0,1,0],
    [0,1,1],
    [1,0,0],
    [1,0,1],
    [1,1,0],
    [1,1,1],
    [0.5,0,0.5],
    [0,0.5,0.5],
]

markers_img = np.zeros(original_img.shape)

for i in range(0,len(markers_img)):
    for j in range(0,len(markers_img[0])):
        markers_img[i,j] = colors[markers[i,j]]



# Displaying the Mask image
fig, imgs = plt.subplots(1,2, figsize=(15,15))
imgs[0].imshow(img_map, cmap='gray')
imgs[0].set_title('Imagem mascara com as 8 folhas - sem marcacao')
imgs[1].imshow(markers_img)
imgs[1].set_title('Imagem mascara com as 8 folhas - com marcacao dos componentes conectados')
plt.axis('off')
plt.show()

markers = cv2.watershed(original_img, markers)

frontier_img = np.zeros(original_img.shape)

for i in range(0,frontier_img.shape[0]):
    for j in range(0,frontier_img.shape[1]):
        frontier_img[i,j] = colors[markers[i,j]]

fig, imgs = plt.subplots(1,2, figsize=(15,15))
imgs[0].imshow(markers_img)
imgs[0].set_title('Imagem mascara com as 8 folhas - com marcacao dos componentes conectados')
imgs[1].imshow(frontier_img, cmap='gray')
imgs[1].set_title('Imagem resultante do watershed')
plt.axis('off')
plt.show()

