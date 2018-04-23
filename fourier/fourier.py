# Based on MAC0417/5768 Visão e Processamento de Imagens
# from Prof. Dr. Paulo A. V. de Miranda
import math
import cv2
import numpy as np

def creatematrix(nrows, ncols):
    M = []
    for i in range(0,nrows):
        L = []
        for j in range(0,ncols):
            L.append(0)
        M.append(L)
    return M

def createcosineimage(nrows, ncols, theta, freq):
    M = creatematrix(nrows, ncols)
    dx = math.cos(theta)
    dy = math.sin(theta)
    for i in range(len(M)):
        for j in range(len(M[0])):
            #projecao vetorial:
            D = dx*j + dy*i
            M[i][j] = math.cos(D*2.0*math.pi*freq)/2.0+0.5
    A = np.array(M)
    A = np.floor(255*A)
    A = A.astype(np.uint8)
    return A


img = cv2.imread('quadrado3.png',0) #quadrado.png
h, w = img.shape
img = img.astype(np.float64)

cosine_img = createcosineimage(h, w, 0.0*math.pi/2.0, 0.05)
img = img + cosine_img

f = np.fft.fft2(img)
fshift = np.fft.fftshift(f)
#fshift = np.copy(f)
spectrum = 20*np.log(1+np.abs(fshift))


Imax = img.max()
img = np.floor(255*(img/float(Imax)))
img = img.astype(np.uint8)
cv2.imshow('Image', img)

cv2.imshow('Cosine', cosine_img)

Smax = spectrum.max()
spectrum = np.floor(255*(spectrum/Smax))
spectrum = spectrum.astype(np.uint8)
cv2.imshow('Spectrum', spectrum)
cv2.waitKey(0)
cv2.imwrite('spectrum.png', spectrum)

