import os

import numpy as np


def regresion_lineal(X, Y) -> tuple([float, float]):
    """
    Simple linear regression
    Parameters
    ----------
    X : array like
    Y : array like
    returns: m y b de la linea mx + b
    """
    X, Y = np.array(X), np.array(Y)
    N = len(X)
    Delta = N * sum(X ** 2) - (sum(X)) ** 2
    m = (N * sum(X * Y) - sum(X) * sum(Y)) / Delta
    b = (sum(X ** 2) * sum(Y) - sum(X) * sum(X * Y)) / Delta

    return m, b


def img2bw(image, low_th=None, upp_th=None) -> np.array:
    """
    Binarize the given image according to the given threshold. The lower and
    the upper threshold are included as True's in the returned matrix:
    e.g: 
        low_th = 4, upp_th = 6 -> 3 = False, 4, 5, and 6 = True, 7 = False and
        so on.
    Parameters
    ----------
    image : numpy array
        Numpy array with the image data.
    lowel_trhesh : int - float, optional
        Minimal threshold for binarizing the image.
        If none given, default values is None.
    upper_trhesh : int - float, optional
        Maximal threshold for binarizing the image.
        If none given, default value is None
    ----------
    Returns: numpy array of bools
        Numpy array with the same dimensionality as the input, but with all
        boolean values.
    """
    if not low_th:
        low_th = 0
    if not upp_th:
        upp_th = np.max(image) + 1
    return (image >= low_th) & (image <= upp_th)


def fits_finder(path=".", extension=".fits") -> list:
    """
    Finds all .fits images in the current working tree
    Parameters
    ----------
    path : string, optional
        The default is "." The main path where to search files with given
        extension.
    extension : string, optional
        The default is ".png". The files' extension to be search in
        subdirectories.
        e.g: 
            ".png", ".jpg"
    Returns: file paths
    """
    file_paths = []
    for root, dirs, files in os.walk(path):
        for name in files:
            path = os.path.join(root, name)
            size = os.stat(path).st_size / 1000
            if (
                path.endswith(extension) and size < 500
            ):  # agregar que busque las que arrancan con proc
                file_paths.append(path)
    return file_paths
