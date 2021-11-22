"""
Calibration of the polynomial coefficients for fitting ADU's vs Electrons
ref: ske_Calibration_T123k.txt
ref: Tesis Kevin Anderson, tabla 3.1: alpha = 1.99e-3
version: 03/10/2021
"""


class OHDU_1:
    ALPHA = 1.9848e-3
    BETA = -9.77376e-11
    GAMMA = 1.87747e-15
    DELTA = -7.08404e-21


class OHDU_2:
    ALPHA = 2.0577e-3
    BETA = -1.88216e-09
    GAMMA = 1.77606e-14
    DELTA = -7.29647e-20


class OHDU_3:
    ALPHA = 1.65799e-3
    BETA = -1.26944e-10
    GAMMA = 1.19425e-15
    DELTA = -3.83054e-21


class OHDU_4:
    ALPHA = 2.14969e-3
    BETA = -1.09097e-08
    GAMMA = 1.13386e-13
    DELTA = -4.02504e-19


class Calibration:
    def __init__(self, OHDU_object):
        calibration = OHDU_object
        self.ALPHA = calibration.ALPHA
        self.BETA = calibration.BETA
        self.GAMMA = calibration.GAMMA
        self.DELTA = calibration.DELTA
