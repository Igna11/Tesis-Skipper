import numpy as np
import skimage.measure
from scipy import ndimage as ndi
from scipy.stats import poisson
import pandas as pd

from analyfits.manipufits import ADU2e
from analyfits.manipufits import single_fits2double_fits
from analyfits.misc import img2bw


def _regionprops_data(image, labels):
    """
    This function must only be called by another functions:
    Parameters
    ----------
    image : numpy array
        Input image array
    labels : numpy array
        Labels of the regions of interes for the input image

    Returns
    -------
        rps : regionprops object
    """

    def image_intensity(region, intensities):
        """
        Es una función interna de get_cluster_info
        Suma las intensidades de una region en unidades de electrones
        """
        return np.sum(intensities[region])

    rps = skimage.measure.regionprops(
        labels,
        intensity_image=image,
        cache=False,
        extra_properties=[image_intensity],
    )
    return rps


def get_cluster_info(image, labels) -> dict:
    """
    Parameters
    ----------
    image : numpy array
        Array de numpy con los datos de la imagen
    labels : numpy array
        Array de numpy de la misma estructura dimensional que el array
        imagen. Contiene 0 donde no hay features y contiene enteros
        crecientes para etiquetar (labelear) los features.
    ----------
    Returns: dict
        Diccionario con las siguientes keys:
            areas,
            centros,
            coordenadas,
            perimetro,
            energias,
            box
    """

    rps = _regionprops_data(image, labels)
    areas = [r.area for r in rps]
    energy = [r.image_intensity for r in rps]
    centros = [r.centroid for r in rps]
    coords = [r.coords for r in rps]
    perimetro = [r.perimeter for r in rps]

    dic_props = {
        "areas": areas,
        "centros": centros,
        "coordenadas": coords,
        "perimetro": perimetro,
        "energias": energy,
    }
    return dic_props


def cluster_info_extract(
    image: str, low_th=2, upp_th=300, ohdu=0, img="o", mu_bkg=0.1858
):
    """
    Extracts the cluster information from image given the lower threshold,
    upper threshold ohdu, the kind of data desired to generate (img='o',
    'l' or 's') and mu_bkg

    Parameters
    ----------
    image: str
        The string path of the image to process.
    low_th: int, optional
        Lower threshol to remove events with less or equal charge than low_th.
    upp_th: int, optional
        Upper threshol to remove events with more or equal than upp_th.
    ohdu: int, optional
        Select the cuadrant of the sensor to analyze the image. 0, 1, 2, 3.
    img: str, optional
        Type of output data:
            "o" for original clusters.
            "l" for simulated noise with original labels <- deprecated.
            "s" for simulated noise and new labels <- deprecated.
    mu_bkg: float, optional
        Expected value from the poissonian distribution of the background
        noise. mu_bkg = 0.1858 default value.
    Returns
    -------
    dic: dict
        dictionary with cluster features of the image
    n_features: int
        number of features found in a given image
    """
    error_message: str = "Invalid string, try 'o', 'l' or 's'"
    image_o = ADU2e(image, ohdu=ohdu)[0]

    # Imagen con clústers unicamente
    image_s = single_fits2double_fits(image, ohdu=ohdu)[1]
    # Le agrego ruido a la imagen con clusters unicamente
    image_s += poisson.rvs(mu_bkg, size=(50, 493))

    # Binarizo ambas imagenes:
    image_o_bw = img2bw(image_o, lower_thresh=low_th, upper_thresh=upp_th)
    image_s_bw = img2bw(image_s, lower_thresh=low_th, upper_thresh=upp_th)

    # Genero las labels y también guardo el número de features
    label_im_o, n_features_o = ndi.label(
        image_o_bw, structure=[[0, 1, 0], [1, 1, 1], [0, 1, 0]]
    )
    label_im_s, n_features_s = ndi.label(
        image_s_bw, structure=[[0, 1, 0], [1, 1, 1], [0, 1, 0]]
    )

    # Armo el diccionario con información que devuelve get_cluster_size_etc
    if img == "o":
        dic = get_cluster_info(image_o, label_im_o)
        n_features = n_features_o
    elif img == "l":
        dic = get_cluster_info(image_s, label_im_o)
        n_features = n_features_o
    elif img == "s":
        dic = get_cluster_info(image_s, label_im_s)
        n_features = n_features_s
    else:
        raise ValueError(error_message)

    return dic, n_features


def image2cluster_info(
    fits_imgs_list, low_th=2, upp_th=300, ohdu=0, img="o", mu_bkg=0.1858
):
    """
    Extrae la información de los clusters de todas las imágenes y lo guarda en
    listas
    Parameters
    ----------
    image: str
        The string path of the image to process.
    low_th: int, optional
        Lower threshol to remove events with less or equal charge than low_th.
    upp_th: int, optional
        Upper threshol to remove events with more or equal than upp_th.
    ohdu: int, optional
        Select the cuadrant of the sensor to analyze the image. 0, 1, 2, 3.
    img: str, optional
        Type of output data:
            "o" for original clusters.
            "l" for simulated noise with original labels <- deprecated.
            "s" for simulated noise and new labels <- deprecated.
    mu_bkg: float, optional
        Expected value from the poissonian distribution of the background
        noise. mu_bkg = 0.1858 default value.
    
    Returns 
    -------
    lista_dics: list
        list with dictionaryies of all images in list_imgs_fits.
    lista_img_idxs: list
        List of numpy arrays with of lenght of the total features per image,
        with the image index en each element, e.g:
            lista_img_idxs = [array(1,1,...,1),
                              array(2,2,..,2),
                              ...,
                              array(N,N,...,N)]
    """

    # Armo listas donde guardo la info
    lista_dics = []
    lista_n_features = []
    lista_img_idxs = []

    # recorro las imagenes que quiero usar
    for i, image in enumerate(fits_imgs_list):

        dic, n_features = cluster_info_extract(
            image,
            low_th=low_th,
            upp_th=upp_th,
            ohdu=ohdu,
            img=img,
            mu_bkg=mu_bkg,
        )
        lista_dics.append(dic)
        lista_n_features.append(n_features)
        lista_img_idxs.append(np.ones(n_features) * i)
        # Printeo el progreso del ciclo for
        print("\r%.2f%%" % (100 * (i + 1) / 925), end="")

    return lista_dics, lista_img_idxs


def df_gen(lista_dics, lista_img_idxs):
    """
    Esta función arma a partir de la listas de información un DataFrame de
    pandas para manipular la información más comodamente
    Parameters
    ----------
    lista_dics : list
        lista de diccionarios
    lista_img_idxs : list
        lista de índices para identificar una misma imagen
    ----------
    Returns: pd.DataFrame
        DataFrame con toda la información de las imágenes
    """

    # armo listas para meter los datos
    energias = []
    areas = []
    centros = []
    perimetros = []

    # meto los datos en las listas
    for dics in lista_dics:
        energias += dics["energias"]
        areas += dics["areas"]
        centros += dics["centros"]
        perimetros += dics["perimetro"]

    # concateno los índices para que sea un solo array
    lista_img_idxs = np.concatenate(lista_img_idxs)

    # Armo el dicionario
    diccionario = {}
    diccionario["img_idx"] = lista_img_idxs
    diccionario["energia"] = energias
    diccionario["area"] = areas
    diccionario["centro"] = centros
    diccionario["perimetro"] = perimetros

    return pd.DataFrame.from_dict(diccionario)


def border_events(img_path, low_th=2, ohdu=0, verbose=False):
    """
    Counts the total charge in the border of all clusters in a given
    image, as well as the total area of the sum of the borders (pixels)

    Parameters
    ----------
    img_path : str path
        pathname or name of the image
    low_th : int, optional
        Lower threshold, equivalent to EPIX, e.g. epix=1.5 == low_th=2.
    verbose: bool, optional
        Set true to see errors exceptions
    Returns
    -------
    charg : int
        Total charge of the clusters of the image
    pixels : int
        Total area (number of pixels) of the border of the clusters
    n_clusters : int
        Total number of clusters
    """
    # Cargo la imagen
    img = ADU2e(img_path, ohdu=ohdu)[0]
    # binarizo con thresh=2 equivalente a epix 1.5
    img_bw = img2bw(img, lower_thresh=2)
    # Structura
    struc = [[1, 1, 1], [1, 1, 1], [1, 1, 1]]

    # Le busco los labels
    img_label, count = ndi.label(img_bw, structure=struc)

    # Usando los labels, voy a buscar los clusters que cumplan que
    # tienen entre 160 y 200 electrones

    rps = _regionprops_data(img, img_label)

    # Aprovecho y cuento la cantidad de clusters de 180 por imagen
    cluster_index = []
    n_clusters = 0
    for r in rps:
        intens = r.image_intensity
        if intens >= 160 and intens <= 200:
            cluster_index.append(1)
            n_clusters += 1
        else:
            cluster_index.append(0)

    # Ahora, buscando en qué índice de la lista cluster_index hay 1's,
    # puedo matchearlo con el label del mismo valor que el número del
    # índice. Esos labels van a sobrevivir, los demás van a ser 0.

    # para eso, transformo la lista en array de numpy
    cluster_index_np = np.array(cluster_index)

    # Lista de vals/indices donde los clusters tienen ese rango de carga
    values = list(np.where(cluster_index_np == 1)[0])

    # Me armo una lista con condiciones (arrays True/False)
    conditions = [img_label == value for value in values]

    # Ahora me armo la máscara solo en los clusters deseados, para eso
    # redefino la img_label. Equivale a haber encontrado solo clusters con
    # carga entre 160 y 200 electrones
    try:
        img_label = np.select(conditions, values)
        # Ahora dilato esos clusters
        img_dilat = ndi.morphology.binary_dilation(img_label, structure=struc)

        # Transformo a 0's y 1's la img_dilatada
        img_dilat = img_dilat * 1  # np.where(img_dilat!=0, 1, 0)

        # Transformo a 0's y 1's la img_label
        img_label = np.where(img_label != 0, 1, 0)

        # Hago la máscara usando la imagen dilatada y img_label
        mask = img_dilat - img_label

        # Se la aplico a la imagen original y sumo la carga
        img_bordes = img * mask
        charg = np.sum(img_bordes)

        # Se la aplico a la binarizada
        pixels = np.sum(mask)

        return charg, pixels, n_clusters
    except ValueError as e:
        if verbose:
            print("imagen sin clusters entre 160 y 200 electrones")
        return 0, 0, 0
