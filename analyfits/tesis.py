import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


def regresion_lineal(X, Y):
    """
    Simple linear regression
    Parameters
    ----------
    X : array like
        DESCRIPTION
    Y : array like
        DESCRIPTION
    returns: m y b de la linea mx + b
    """
    N = len(X)
    Delta = N * sum(X ** 2) - (sum(X)) ** 2
    m = (N * sum(X * Y) - sum(X) * sum(Y)) / Delta
    b = (sum(X ** 2) * sum(Y) - sum(X) * sum(X * Y)) / Delta

    return m, b


def fits_finder(directorio=".", extension=".fits"):
    """
    Finds all .fits images in the current working tree
    Parameters
    ----------
    directorio : string, optional
        DESCRIPTION. The default is "."
        El directorio "madre" del cual se quieren buscar las imágenes .png
        de todos sus subdirectorios
    extension : string, optional
        DESCRIPTION. The default is ".png"
        La extensión del tipo de archivo que se quiere buscar en los
        subdirectorios. Es necesario usar ".extension": EJ: ".png", ".jpg"
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


def ADU2e(
    src_path,
    alpha=1.9848e-3,
    beta=-9.77376e-11,
    gamma=1.87747e-15,
    delta=-7.08404e-21,
    ohdu=0,
    save=False,
):
    """
    From a .fits image, and using the polynomial fit for ADU -> electron
    returns both image data:
        image in ADU's
        image in e- units.

    Parameters
    ----------
    src_path : string
        Directory path of the file
    alpha : float, optional
        By default, alpha = 1.9848e-3
        ref: ske_Calibration_T123k.txt
        ref: Tesis Kevin Anderson, tabla 3.1: alpha = 1.99e-3
    beta : float, optional
        By default, beta = -9.77376e-11
        ref: ske_Calibration_T123k.txt
        ref: Tesis Kevin Anderson, tabla 3.1: beta = 8.61e-11
    gamma : float, optional
        By default, beta = 1.87747e-15
        ref: ske_Calibration_T123k.txt
        ref: Tesis Kevin Anderson, tabla 3.1: gamma = -6.31e-17
    delta : float, optional
        By default, beta = -7.08404e-21
        ref: ske_Calibration_T123k.txt
        ref: Tesis Kevin Anderson, tabla 3.1: delta = 2.73e-23
        the values used here differ from the ones in the Kevin Anderson's tesis
    ohdu : int, optional
        0, 1, 2, 3 allowed (each sensor of the total sensor)
    save : bool, optional
        When set True, a new fits file in electron units is saved
    """
    with fits.open(src_path) as fits_img:
        ADU_img_data = fits_img[ohdu].data
        e_img_data = np.round(
            ADU_img_data * alpha
            + ADU_img_data ** 2 * beta
            + ADU_img_data ** 3 * gamma
            + ADU_img_data ** 4 * delta
        )
        # =====================================
        # Los datos vienen con valores negativos:
        # filtro los valores negativos y los transformo en 0
        ADU_img_data[ADU_img_data < 0] = 0
        e_img_data[e_img_data < 0] = 0

        if save is True:
            fits_img[ohdu].data = e_img_data
            tgt_path = src_path.strip(".fits") + "electron_units.fits"
            fits_img.writeto(tgt_path)
    return e_img_data, ADU_img_data


def single_fits2double_fits(
    src_path, tgt_path="", threshold=1, ohdu=0, save=False
):
    """
    From one .fits file generates 2 new fits files.

    Parameters
    ----------
    src_path : string
        Path to the original .fits file or filename
    tgt_path : string
        path to the 2 new .fits files or filename WITHOUT EXTENSION.
        e.g /media/usr/.../image_name <-- correct
        e.g /media/usr/.../image_name.fits <-- incorrect
    threshold : int, optional
        Threshold to start counting electrons. The data will be divided into
        two diferents arrays, one with all pixels having a number of electrons
        greater than threshold and the other one having the rest.
        The default is 1.
    ohdu : int, optional.
        0,1,2,3 allowed (each sensor of the total sensor)
    save : bool, optional.
        If set True, then generates the 2 .fits files with the threshold
        applied
    Returns
    -------
    None.

    """
    electron_fits = ADU2e(src_path, ohdu=ohdu)[0]

    # make some copies of the original data to modify it avoiding aliasing
    elec_minor = electron_fits.copy()
    elec_major = electron_fits.copy()

    # All pixels with a number of electrons > threshold are set to 0
    elec_minor[elec_minor > threshold] = 0
    # All pixels with a number of electrons < threshold are set to 0
    elec_major[elec_major <= threshold] = 0

    # differenciate the minor and major paths
    tgt_path_minor = tgt_path + "minor.fits"
    tgt_path_major = tgt_path + "major.fits"
    if save is True:
        with fits.open(src_path) as fits_img:
            # save elec_minor
            fits_img[ohdu].data = elec_minor
            try:
                fits_img.writeto(tgt_path_minor)
            except OSError as e:
                print(e)
                print("A file with the same name already exists.")
            # save elec_major
            fits_img[ohdu].data = elec_major
            try:
                fits_img.writeto(tgt_path_major)
            except OSError as e:
                print(e)
                print("A file with the same name already exists.")
    return elec_minor, elec_major


def get_cluster_info(image, labels):
    """
    Parameters
    ----------
    image : numpy array
        DESCRIPTION:
        Array de numpy con los datos de la imagen
    labels : numpy array
        DESCRIPTION.
        Array de numpy de la misma estructura dimensional que el array
        imagen. Contiene 0 donde no hay features y contiene enteros
        crecientes para etiquetar (labelear) los features.
    ----------
    Returns: dict
        DESCRIPTION:
        Diccionario con las siguientes keys:
            areas,
            centros,
            coordenadas,
            perimetro,
            energias,
            box
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


def img2bw(image, lower_thresh=None, upper_thresh=None):
    """
    Binariza la imagen ingresada según el threshold dado:
    Parameters
    ----------
    image : numpy array
        DESCRIPTION:
        Array de numpy con los datos de la imagen
    lowel_trhesh : int - float, optional
        DESCRIPTION.
        threshold minimo para el binarizado de la imagen.
        Si nada es ingresado, entonces se toma como nulo
    upper_trhesh : int - float, optional
        DESCRIPTION.
        threshold máximo para el binarizado de la imagen.
        Si nada es ingresado, entonces se toma como nulo
    ----------
    Returns: numpy array of bools
        DESCRIPTION:
        array de numpy con la misma dimensionalidad que el
        input, pero con todos valores booleanos.
    """
    if lower_thresh is None:
        lower_thresh = 0
    if upper_thresh is None:
        upper_thresh = np.max(image) + 1
    return (image >= lower_thresh) & (image <= upper_thresh)


def image2cluster_info(fits_imgs_list, low_th=2, upp_th=300, ohdu=0, img="o"):
    """
    Extrae la información de los clusters de todas las imágenes y lo guarda en
    listas
    Parameters
    ----------0.81916
    fits_imgs_list : list
        DESCRIPTION.
        lista de directorios donde están las imágenes .fits
    img : str - optional
        DESCRIPTION.
        para indicarle si va a buscar clusters:
            "o": Busca los features en la imagen original.
            "l": Busca los features en la imagen simulada, pero usa los labels
            para de la imagen original.
            "s": Busca los features de la iamgen simulada.
    low_th : int - optional
        DESCRIPTION.
        lower detection threshold
    upp_th : int - optional
        DESCRIPTION.
        Upper detection threshold
    ohdu : int - optional
        DESCRIPTION.
        The HDU data: 0, 1, 2, 3 available
    ----------
    Returns: lista_dics, lista_b_features, lista_img_idxs
        DESCRIPTION:
        lista con la información
    """

    # Armo listas donde guardo la info
    lista_dics = []
    lista_n_features = []
    lista_img_idxs = []

    # Hardcodeo el mu de la poissoniana
    mu = 0.2243

    # recorro las imagenes que quiero usar
    for i, image in enumerate(fits_imgs_list):
        # Imagen completa
        image_o = ADU2e(image, ohdu=ohdu)[0]

        # Imagen con clústers unicamente
        image_s = single_fits2double_fits(image, ohdu=ohdu)[1]
        # Le agrego ruido a la imagen con clusters unicamente
        image_s += poisson.rvs(mu, size=(50, 493))

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
            print("opción incorrecta: Solo 'o', 'l', 's' son válidos")
            return 0
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
        DESCRIPTION:
        lista de diccionarios
    lista_img_idxs : list
        DESCRIPTION.
        lista de índices para identificar una misma imagen
    ----------
    Returns: pd.DataFrame
        DESCRIPTION:
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


def border_events(img_path, low_th=2, ohdu=0):
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

    def image_intensity(region, intensities):
        """
        Es una función interna de get_cluster_info
        Suma las intensidades de una region en unidades de electrones
        """
        return np.sum(intensities[region])

    rps = skimage.measure.regionprops(
        img_label,
        intensity_image=img,
        cache=False,
        extra_properties=[image_intensity],
    )

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
        carga = np.sum(img_bordes)

        # Se la aplico a la binarizada
        pixeles = np.sum(mask)

        return carga, pixeles, n_clusters
    except ValueError as e:
        print("imagen %i sin clusters entre 160 y 200 electrones" % i)
        return 0, 0, 0
