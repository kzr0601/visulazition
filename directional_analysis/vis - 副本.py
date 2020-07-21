import cv2
import numpy as np
from math import sqrt, atan, pi

from PIL import Image


def main(imgname):
	img = cv2.imread(imgname, 0)
	PI_4 = pi/4
	PI_2 = pi/2

	# calculate gx and gy
	sobelx = cv2.Sobel(img, cv2.CV_64F, 1, 0, ksize=3)
	sobely = cv2.Sobel(img, cv2.CV_64F, 0, 1, ksize=3)
	cv2.imwrite("gx.png", sobelx)
	cv2.imwrite("gy.png", sobely)

	sobelx = sobelx.astype(np.float32)
	sobely = sobely.astype(np.float32)

	gradient_value = np.zeros((img.shape[0], img.shape[1]))
	gradient_direction = np.zeros((img.shape[0], img.shape[1]))
	
	# calculate gradient direction and value
	THRE = 20
	mask = np.ones((img.shape[0], img.shape[1], 3))
	for i in range(img.shape[0]):
		for j in range(img.shape[1]):
			gradient_value[i][j] = sqrt(pow(sobelx[i][j], 2)+ pow(sobely[i][j], 2))
			if(gradient_value[i][j] < THRE):
				mask[i, j, :] = 0
			# if(sobelx[i][j] == 0):
			# 	gradient_direction[i][j] = 0
			# else:
			# 	gradient_direction[i][j] = atan(sobely[i][j] / sobelx[i][j])
	gradient_direction = cv2.phase(sobelx, sobely, angleInDegrees=True)	# True in degree, False in radius
	# gradient_direction = 255*(gradient_direction-np.min(gradient_direction))/(np.max(gradient_direction)-np.min(gradient_direction))
	# gradient_value = cv2.addWeighted(sobelx, 0.5, sobely, 0.5, 0)
	
	# do color map
	gradient_value_color = cv2.applyColorMap(gradient_value.astype(np.uint8), cv2.COLORMAP_JET)
	gradient_direction = cv2.applyColorMap(gradient_direction.astype(np.uint8), cv2.COLORMAP_JET)


	# gd = np.zeros((img.shape[0], img.shape[1], 3))
	# for i in range(img.shape[0]):
	# 	for j in range(img.shape[1]):
	# 		if(abs(gradient_direction[i][j]) < PI_4 ):
	# 			gd[i][j][1] = 255# * abs(gradient_direction[i][j]) / PI_4
	# 		else:
	# 			gd[i][j][2] = 255# * (abs(gradient_direction[i][j])-PI_4)/PI_4

	gd = gradient_direction
	gd = gd * mask

	# write data
	cv2.imwrite("gradient_value.png", gradient_value_color)
	# cv2.imwrite("gradient_value_gray.png", gradient_value)
	cv2.imwrite("gradient_direction.png", gd)

	gradient_value = cv2.imread("gradient_value.png", 0)
	img = cv2.imread("gradient_direction.png")
	# print(np.max(gradient_value), np.min(gradient_value))
	gradient_value = (gradient_value-np.min(gradient_value))/(np.max(gradient_value)-np.min(gradient_value))
	# print(np.max(gradient_value), np.min(gradient_value))
	for i in range(gradient_value.shape[0]):
		for j in range(gradient_value.shape[1]):
			gradient_value[i][j] = min(1, gradient_value[i][j]+0.3)
	for i in range(3):
		img[:, :, i] = img[:, :, i]*gradient_value
	cv2.imwrite("gradient_direction_2.png", img)

	# b, g, r = cv2.split(img)
	# # alpha_channel = gradient_value.astype(np.uint8)
	# alpha_channel = np.ones(b.shape, dtype=b.dtype) * 255
	# alpha_channel[:, :int(b.shape[0] / 2)] = 100
	# img_RGBA = cv2.merge((b, g, r, alpha_channel))
	# cv2.imwrite("gradient_direction_rgba.png", img_RGBA)

	# img = Image.open("gradient_direction.png")
	# img = img.convert('RGBA')
	# # r, g, b, alpha = img.split()
	# # a = Image.open("gradient_value_gray.png")
	# # if(a.mode!="L" or a.mode!="1"):
	# # 	a = a.convert("L")
	# # img.putalpha(a)
	# img_blender = Image.new('RGBA', img.size, (0,0,0,0))
	# img = Image.blend(img_blender, img, 1)
	# img.save("gradient_direction_2.png")
	# # print(np.max(alpha), np.min(alpha))

if __name__ == '__main__':
	imgname = "./F-actin/00-1-2-f.tif"
	main(imgname)