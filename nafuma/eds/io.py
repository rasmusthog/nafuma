from PIL import Image
import numpy as np
import cv2

def read_image(path, weight=None, colour=None, crop=None, resize=None, brightness=None):

    img = np.array(Image.open(path))

    if colour is not None:
        img = change_colour(img, colour)

    if brightness is not None:
        img = increase_brightness(img, increase=brightness)

    if crop is not None:
        img = crop_image(img, crop)

    if resize is not None:
        img = resize_image(img, resize)

    if weight is not None:
        img = scale_image(img, weight)

    return img


def scale_image(image, factor):

    for i in range(0,image.shape[0]):
        for j in range(0, image.shape[1]):
            image[i][j][0] = image[i][j][0]*factor
            image[i][j][1] = image[i][j][1]*factor
            image[i][j][2] = image[i][j][2]*factor

    return image


def crop_image(image, factor):

    y, x = image.shape[0:2]

    new_y, new_x = int(y*factor), int(x*factor)

    image = image[:new_y, :new_x]

    res = cv2.resize(image, dsize=(x, y), interpolation=cv2.INTER_CUBIC)

    return res


def resize_image(image, factor):

    y, x = image.shape[0:2]

    new_y, new_x = int(y*factor), int(x*factor)

    res = cv2.resize(image, dsize=(new_x, new_y), interpolation=cv2.INTER_CUBIC)

    return res


def increase_brightness(image, brightness):

    for i in range(0,image.shape[0]):
        for j in range(0, image.shape[1]):
            image[i][j][0] = image[i][j][0]+brightness
            image[i][j][1] = image[i][j][1]+brightness
            image[i][j][2] = image[i][j][2]+brightness


    return image


def add_images(image1, image2):

    assert image1.shape == image2.shape

    compound_image = np.zeros((image1.shape[0], image1.shape[1], image1.shape[2]))
    for i in range(image1.shape[0]):
        for j in range(image1.shape[1]):
            compound_image[i][j] = [0, 0, 0]

            compound_image[i][j][0] = int(int(image1[i][j][0]) + int(image2[i][j][0]))
            compound_image[i][j][1] = int(int(image1[i][j][1]) + int(image2[i][j][1]))
            compound_image[i][j][2] = int(int(image1[i][j][2]) + int(image2[i][j][2]))



    return compound_image




def get_colour(image):


    colour = [0, 0, 0]
    for i in range(image.shape[0]):
        for j in range(image.shape[1]):
            if image[i][j][0] > colour[0]:
                colour[0] = image[i][j][0]
            
            if image[i][j][1] > colour[1]:
                colour[1] = image[i][j][1]
            
            if image[i][j][2] > colour[2]:
                colour[2] = image[i][j][2]

    colour = np.array(colour)

    return colour
            

def change_colour(image, new_colour):

    new_colour = np.array(new_colour)

    old_colour = get_colour(image)


    for i in range(image.shape[0]):
        for j in range(image.shape[1]):
            factor = max(image[i][j]) / max(old_colour)
            image[i][j] = new_colour.astype(float) * factor


    return image
