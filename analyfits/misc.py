import os

import numpy as np


def regresion_lineal(X, Y):
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


def img2bw(image, low_th=None, upp_th=None):
    """
    Binarize the given image according to the given threshold. The lower and
    the upper threshold are included as True's in the returned matrix:
    e.g: 
        low_th = 4, upp_th = 6 -> 3 = False, 4, 5, and 6 = True, 7 = False and
        so on.
    Parameters
    ----------
    image : numpy array
        Array de numpy con los datos de la imagen
    lowel_trhesh : int - float, optional
        threshold minimo para el binarizado de la imagen.
        Si nada es ingresado, entonces se toma como nulo
    upper_trhesh : int - float, optional
        threshold máximo para el binarizado de la imagen.
        Si nada es ingresado, entonces se toma como nulo
    ----------
    Returns: numpy array of bools
        array de numpy con la misma dimensionalidad que el
        input, pero con todos valores booleanos.
    """
    if not low_th:
        low_th = 0
    if not upp_th:
        upp_th = np.max(image) + 1
    return (image >= low_th) & (image <= upp_th)


def fits_finder(directorio=".", extension=".fits"):
    """
    Finds all .fits images in the current working tree
    Parameters
    ----------
    directorio : string, optional
        The default is "." El directorio "madre" del cual se quieren buscar
        las imágenes .png de todos sus subdirectorios
    extension : string, optional
        The default is ".png". La extensión del tipo de archivo que se quiere
        buscar en los subdirectorios. Es necesario usar ".extension": EJ:
        ".png", ".jpg"
    Returns: string de directorios y archivos .png
    """
    archivos = []
    for root, dirs, files in os.walk(directorio):
        for name in files:
            path = os.path.join(root, name)
            size = os.stat(path).st_size / 1000
            if (
                path.endswith(extension) and size < 500
            ):  # agregar que busque las que arrancan con proc
                archivos.append(path)
    return archivos
