import numpy as np
import cv2


def conditional_median_filter(img: np.ndarray, kernel_size: int = 3, th_factor: float = 3.0):
    pad_img = np.pad(img, kernel_size//2, mode='reflect')
    N, M = pad_img.shape
    new_img = np.zeros_like(img)
    
    for i in np.arange(kernel_size//2, N - kernel_size//2):
        for j in np.arange(kernel_size//2, M - kernel_size//2):
            sub_img = pad_img[i-kernel_size//2:i+kernel_size//2+1, j-kernel_size//2:j+kernel_size//2+1]
            sub_img_vec = np.delete(np.hstack(sub_img), kernel_size**2//2)
            sub_img_std = np.std(sub_img_vec)
            sub_img_mean = np.mean(sub_img_vec)
            if np.abs(sub_img[kernel_size//2,kernel_size//2] - sub_img_mean) > th_factor * sub_img_std:
                new_img[i-kernel_size//2, j-kernel_size//2] = np.median(sub_img_vec)
            else:
                new_img[i - kernel_size // 2, j - kernel_size // 2] = pad_img[i,j]

    return new_img


def non_bone_voxel_filter(img: np.ndarray, th_val: float, std_th: float, kernel_size: int = 3):
    pad_img = np.pad(img, kernel_size // 2, mode='edge')
    N, M, H = pad_img.shape
    new_img = np.zeros_like(img)
    # pad_img = np.sum(cv2.GaussianBlur(pad_img, (kernel_size, kernel_size), 7*np.std(img)), axis=2) / H


    for k in np.arange(kernel_size // 2, H - kernel_size // 2):
        pad_img = cv2.bilateralFilter(img[:, :, k - kernel_size // 2], 9, 0.75, 0.75)
        for i in np.arange(kernel_size // 2, N - kernel_size // 2):
            for j in np.arange(kernel_size // 2, M - kernel_size // 2):
                img_val = img[i - kernel_size // 2, j - kernel_size // 2, k - kernel_size // 2]
                if img_val > th_val:
                    sub_img = pad_img[i - kernel_size // 2:i + kernel_size // 2 + 1,
                              j - kernel_size // 2:j + kernel_size // 2 + 1]
                    sub_img_vec = np.hstack(sub_img)
                    sub_img_std = np.std(sub_img_vec)
                    if sub_img_std < std_th*((img_val - th_val) / th_val):
                        new_img[i - kernel_size // 2, j - kernel_size // 2, k - kernel_size // 2] = 1

    return new_img

def cond_median_filter(img: np.ndarray, th_val: float, kernel_size: int = 3):
    pad_img = np.pad(img, kernel_size // 2, mode='edge')
    N, M, H = pad_img.shape
    new_img = np.zeros_like(img)
    # pad_img = np.sum(cv2.GaussianBlur(pad_img, (kernel_size, kernel_size), 7*np.std(img)), axis=2) / H
    for k in np.arange(kernel_size // 2, H - kernel_size // 2):
        for i in np.arange(kernel_size // 2, N - kernel_size // 2):
            for j in np.arange(kernel_size // 2, M - kernel_size // 2):
                img_val = img[i - kernel_size // 2, j - kernel_size // 2, k - kernel_size // 2]
                if img_val == 1:
                    sub_img = pad_img[i - kernel_size // 2:i + kernel_size // 2 + 1,
                              j - kernel_size // 2:j + kernel_size // 2 + 1,
                              k - kernel_size // 2:k + kernel_size // 2 + 1]
                    sub_img_vec = np.hstack(sub_img)
                    sub_img_sum = np.sum(sub_img_vec)
                    if sub_img_sum > th_val:
                        new_img[i - kernel_size // 2, j - kernel_size // 2, k - kernel_size // 2] = 1

    return new_img