import numpy as np
from astropy.io import fits

from analyfits import calibration


def _fits_saver(fits_image, tgt_path="", ohdu=0):
    """
    Saves the fits_image data into the given tgt_path, with the given ohdu
    Parameters
    ----------
    fits_image :
        The data to save into .fits image
    tgt_path : str, optional
        pathname for the new file
    ohdu : int, optional
        0, 1, 2, 3 available.
    """
    with fits.open(src_path) as fits_img:
        fits_img[ohdu].data = fits_image
        try:
            fits_img.writeto(tgt_path)
        except OSError as e:
            print(e)
            print("A file with the same name already exists.")


def ADU2e(
    src_path,
    alpha=calibration.ALPHA,
    beta=calibration.BETA,
    gamma=calibration.GAMMA,
    delta=calibration.DELTA,
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
        Directory path of the file.
    alpha : float, optional
        By default, alpha = 1.9848e-3. Coefficient of the polynomial fit.
    beta : float, optional
        By default, beta = -9.77376e-11 Coefficient of the polynomial fit.
    gamma : float, optional
        By default, beta = 1.87747e-15 Coefficient of the polynomial fit.
    delta : float, optional
        By default, beta = -7.08404e-21 Coefficient of the polynomial fit.
    ohdu : int, optional
        0, 1, 2, 3 allowed (each sensor of the total sensor)
    save : bool, optional
        When set True, a new fits file in electron units is saved
    """
    with fits.open(src_path) as fits_img:
        ADU_img_data = fits_img[ohdu].data
        e_original = (
            ADU_img_data * alpha
            + ADU_img_data ** 2 * beta
            + ADU_img_data ** 3 * gamma
            + ADU_img_data ** 4 * delta
        )
        #e_img_data = e_original
        e_img_data = np.round(e_original)
        # =====================================
        # Los datos vienen con valores negativos:
        # filtro los valores negativos y los transformo en 0
        ADU_img_data[ADU_img_data < 0] = 0
        e_img_data[e_img_data < 0] = 0

        if save:
            tgt_path = src_path.strip(".fits") + "electron_units.fits"
            _fits_saver(e_img_data, tgt_path=tgt_path, ohdu=ohdu)

    return e_img_data, ADU_img_data


def single_fits2double_fits(
    src_path, tgt_path="", threshold=1, ohdu=0, save=False
):
    """
    From one .fits image generates 2 different images and optionally generates
    2 new fits files to be saved.

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
    if save:
        # differenciate the minor and major paths
        tgt_path_minor = tgt_path + "minor.fits"
        tgt_path_major = tgt_path + "major.fits"
        _fits_saver(elec_minor, tgt_path=tgt_path_minor, ohdu=ohdu)
        _fits_saver(elec_major, tgt_path=tgt_path_major, ohdu=ohdu)

    return elec_minor, elec_major
