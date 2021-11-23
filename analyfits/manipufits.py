import numpy as np
from astropy.io import fits

from analyfits import calibration
from analyfits.calibration import Calibration


class ManipulateFits:
    def __init__(self):
        """
        Initialices the calibration for the OHDU 1 (numered as 0 in code)
        which is the most frequently used calibration
        """
        self._ohdu = 0
        self._alpha = 1.9848e-3
        self._beta = -9.77376e-11
        self._gamma = 1.87747e-15
        self._delta = -7.08404e-21

    def set_calibration(self, alpha, beta, gamma, delta):
        """Sets a custom calibration if needed"""
        self._alpha = alpha
        self._beta = beta
        self._gamma = gamma
        self._delta = delta

    def set_ohdu(self, ohdu):
        """
        Sets the desired OHDU to use in the analysis and picks the correct
        calibration.
        """
        self._ohdu = ohdu
        if ohdu == 0:
            OHDUObject = calibration.OHDU_1
        elif ohdu == 1:
            OHDUObject = calibration.OHDU_2
        elif ohdu == 2:
            OHDUObject = calibration.OHDU_3
        elif ohdu == 3:
            OHDUObject = calibration.OHDU_4
        cal = Calibration(OHDUObject)
        self._alpha = cal.ALPHA
        self._beta = cal.BETA
        self._gamma = cal.GAMMA
        self._delta = cal.DELTA

    def ADU2e(self, src_path, ADU=False):
        """
        From a .fits image, and using the polynomial fit for ADU -> electron
        returns either both image data:
            image in ADU's
            image in e- units.
        Parameters
        ----------
        src_path : string
            Path to the original .fits file or filename
        ADU : Bool, optional
            Default value False, if true is given then returns the numpy array
            in units of ADU's, otherwise in units of electrons.
        Returns
        -------
        e_img_data
        """
        with fits.open(src_path) as fits_img:
            self.ADU_img_data = fits_img[self._ohdu].data
            e_original = (
                self.ADU_img_data * self._alpha
                + self.ADU_img_data ** 2 * self._beta
                + self.ADU_img_data ** 3 * self._gamma
                + self.ADU_img_data ** 4 * self._delta
            )
            # e_img_data = e_original
            self.e_img_data = np.round(e_original)

            # cleaning of negative values
            if ADU:
                self.ADU_img_data[self.ADU_img_data < 0] = 0
                return self.ADU_img_data

            self.e_img_data[self.e_img_data < 0] = 0

        return self.e_img_data

    def single_fits2double_fits(self, src_path, threshold=1):
        """
        From one .fits image generates 2 different images and optionally
        generates 2 new fits files to be saved.

        Parameters
        ----------
        src_path : string
            Path to the original .fits file or filename
        threshold : int, optional
            Threshold to start counting electrons. The data is divided into
            two diferents arrays, one with all pixels having a number of
            electrons greater than threshold and the other one having the rest.
            The default is 1.
        Returns
        -------
        elec_minor : np.array
            The matrix array of the image with minor number of events per pixel
        elec_major : np.array
            The matrix array of the image with major number of events per pixel
        """
        electron_fits = self.ADU2e(src_path)

        # make some copies of the original data to modify it avoiding aliasing
        self.elec_minor = electron_fits.copy()
        self.elec_major = electron_fits.copy()

        # All pixels with a number of electrons > threshold are set to 0
        self.elec_minor[self.elec_minor > threshold] = 0
        # All pixels with a number of electrons < threshold are set to 0
        self.elec_major[self.elec_major <= threshold] = 0

        return self.elec_minor, self.elec_major

    def save_fits(self, src_path, tgt_path):
        """
        Saves the fits_image data in units of electrons only, into the
        given tgt_path
        Parameters
        ----------
        tgt_path : str
            pathname for the new file
        """
        with fits.open(src_path) as fits_img:
            fits_img[self._ohdu].data = self.e_img_data
            try:
                fits_img.writeto(tgt_path)
            except OSError as e:
                print(e)
                print("A file with the same name already exists.")

    def easy_plot(self, cmap="hot", vmin=0, vmax=10):
        """Shortcut to watch the original image"""
        import matplotlib.pyplot as plt

        plt.figure(figsize=(17, 2))
        plt.imshow(self.e_img_data, cmap=cmap, vmin=vmin, vmax=vmax)
