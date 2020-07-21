import cv2
import numpy as np
from math import sqrt, atan, pi, floor, atan

from PIL import Image

def hsv2rgb(h, s, v):
    h = float(h)
    s = float(s)
    v = float(v)
    h60 = h / 60.0
    h60f = floor(h60)
    hi = int(h60f) % 6
    f = h60 - h60f
    p = v * (1 - s)
    q = v * (1 - f * s)
    t = v * (1 - (1 - f) * s)
    r, g, b = 0, 0, 0
    if hi == 0: r, g, b = v, t, p
    elif hi == 1: r, g, b = q, v, p
    elif hi == 2: r, g, b = p, v, t
    elif hi == 3: r, g, b = p, q, v
    elif hi == 4: r, g, b = t, p, v
    elif hi == 5: r, g, b = v, p, q
    r, g, b = int(r * 255), int(g * 255), int(b * 255)
    return r, g, b


def main(imgname):
	img_color = cv2.imread(imgname)
	img = cv2.imread(imgname, 0)
	# gray_img = (img-np.min(img))/(np.max(img)-np.min(img))*255
	# gray_img[gray_img>30] = 255
	gray_img = img_color[:, :, 2]

	img = cv2.GaussianBlur(img,(5,5),0)
	PI_4 = pi/4
	PI_2 = pi/2

	# calculate gx and gy
	sobelx = cv2.Sobel(img, cv2.CV_64F, 1, 0, ksize=3)
	sobely = cv2.Sobel(img, cv2.CV_64F, 0, 1, ksize=3)
	cv2.imwrite("gx.png", sobelx)
	cv2.imwrite("gy.png", sobely)

	sobelx = sobelx.astype(np.float32)
	sobely = sobely.astype(np.float32)

	WINDOWRADIUS = 5
	WINDOWSIZE = 2*WINDOWRADIUS+1
	gx = cv2.copyMakeBorder(sobelx, WINDOWRADIUS, WINDOWRADIUS, WINDOWRADIUS, WINDOWRADIUS, cv2.BORDER_REFLECT)
	gy = cv2.copyMakeBorder(sobely, WINDOWRADIUS, WINDOWRADIUS, WINDOWRADIUS, WINDOWRADIUS, cv2.BORDER_REFLECT)

	gx_sqr = gx*gx
	gy_sqr = gy*gy
	gx_gy = gx*gy
	gx_sqr_gausse = cv2.GaussianBlur(gx_sqr, (WINDOWSIZE, WINDOWSIZE), 0)
	gy_sqr_gausse = cv2.GaussianBlur(gy_sqr, (WINDOWSIZE, WINDOWSIZE), 0)
	gx_gy_gausse = cv2.GaussianBlur(gx_gy, (WINDOWSIZE, WINDOWSIZE), 0)

	theta = np.zeros((img.shape[0], img.shape[1]))
	confidence = np.zeros((img.shape[0], img.shape[1]))
	
	for i in range(img.shape[0]):
		for j in range(img.shape[1]):
			matrix = np.array([[gx_sqr_gausse[i+WINDOWRADIUS][j+WINDOWRADIUS], gx_gy_gausse[i+WINDOWRADIUS][j+WINDOWRADIUS] ], [gx_gy_gausse[i+WINDOWRADIUS][j+WINDOWRADIUS], gy_sqr_gausse[i+WINDOWRADIUS][j+WINDOWRADIUS] ] ])
			a,b=np.linalg.eig(matrix) 
			max_a = max(a[0], a[1])
			min_a = min(a[0], a[1])
			if(max_a+min_a != 0):
				confidence[i][j] = (max_a-min_a)/(max_a+min_a)

			if(gy_sqr_gausse[i+WINDOWRADIUS][j+WINDOWRADIUS] - gx_sqr_gausse[i+WINDOWRADIUS][j+WINDOWRADIUS] == 0):
				if(gx_gy_gausse[i+WINDOWRADIUS][j+WINDOWRADIUS]>0):
					theta[i][j] = pi/4
				else:
					theta[i][j] = -pi/4
			else:
				theta[i][j] = 1.0 / 2 * atan(2 * gx_gy_gausse[i+WINDOWRADIUS][j+WINDOWRADIUS] / (gy_sqr_gausse[i+WINDOWRADIUS][j+WINDOWRADIUS] - gx_sqr_gausse[i+WINDOWRADIUS][j+WINDOWRADIUS]))
			if(gy_sqr_gausse[i+WINDOWRADIUS][j+WINDOWRADIUS] < gx_sqr_gausse[i+WINDOWRADIUS][j+WINDOWRADIUS]):
				theta[i][j] += pi/2

	theta = theta / pi * 180 + 90
	confidence = (confidence-np.min(confidence))/(np.max(confidence)-np.min(confidence))*255

	# write hsv img: h:gradient_d, s:confidence, v:gray_img
	img_hsv = np.ones((img.shape[0], img.shape[1], 3))
	result_img = np.ones((img.shape[0], img.shape[1], 3))

	img_hsv[:, :, 0] = theta #gradient_d
	img_hsv[:, :, 1] = confidence
	img_hsv[:, :, 2] = gray_img

	img_hsv = img_hsv.astype(np.uint8)
	result_img = cv2.cvtColor(img_hsv, cv2.COLOR_HSV2BGR)

	#print(np.max(result_img), np.min(result_img))
	#result_img = (result_img-np.min(result_img))/(np.max(result_img)-np.min(result_img))*255
	
	#result_img = cv2.applyColorMap(gradient_direction.astype(np.uint8), cv2.COLORMAP_HSV)
	#print(np.max(gradient_direction), np.min(gradient_direction))
	#result_img = (gradient_direction-np.min(gradient_direction))/(np.max(gradient_direction)-np.min(gradient_direction))*255

	cv2.imwrite("result_hsv.png", result_img)

if __name__ == '__main__':
	imgname = "./F-actin/00-1-2-f.tif"
	#imgname = "./chirp.tif"
	main(imgname)