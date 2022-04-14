"""
Calibration of the polynomial coefficients for fitting ADU's vs Electrons
ref: ske_Calibration_T123k.txt
ref: Tesis Kevin Anderson, tabla 3.1: alpha = 1.99e-3
version: 03/10/2021
"""


class OHDU_1:
    ALPHA = 2.02239e-3
    BETA = 5.09616e-11
    GAMMA = 2.37466e-16
    DELTA = -8.59986e-22


class OHDU_2:
    ALPHA = 1.95224e-3
    BETA = 3.39492e-11
    GAMMA = -3.02458e-16
    DELTA = 4.15541e-22


class OHDU_3:
    ALPHA = 1.68488e-3
    BETA = 5.82793e-11
    GAMMA = -2.94793e-16
    DELTA = 5.3502e-22


class OHDU_4:
    ALPHA = 1.78221e-3
    BETA = -8.09364e-11
    GAMMA = 1.02599e-16
    DELTA = 3.55681e-22


class Calibration:
    def __init__(self, OHDU_object):
        calibration = OHDU_object
        self.ALPHA = calibration.ALPHA
        self.BETA = calibration.BETA
        self.GAMMA = calibration.GAMMA
        self.DELTA = calibration.DELTA
