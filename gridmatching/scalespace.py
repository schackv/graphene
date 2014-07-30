import numpy as np
from numpy import pi, exp

"""Module kindly provided by abll @ github.com/andersbll"""


class ScaleSpace:
    def __init__(self, img_shape, sigmas, dys, dxs):
        ''' Compute the scale-space of an image.
        Upon initialization, this class precomputes the Gaussian windows used
        to smooth images of a fixed shape to save the computations at later
        points.
        '''
        assert(len(sigmas) == len(dys) == len(dxs))
        h, w = img_shape
        g_y, g_x = np.mgrid[-.5+.5/h:.5:1./h, -.5+.5/w :.5: 1./w]
        self.filters = []
        for sigma, dy, dx in zip(sigmas, dys, dxs):
            g = exp(- (g_x**2 + g_y**2) * (pi*2*sigma)**2 / 2.)
            g = np.fft.fftshift(g)
            if dy > 0 or dx > 0:
                #TODO change list(range to np.linspace or similar
                dg_y = np.array((list(range(0, h//2))+list(range(-h//2, 0))),
                                dtype=float, ndmin=2) / h
                dg_x = np.array((list(range(0, w//2))+list(range(-w//2, 0))),
                                dtype=float, ndmin=2) / w
                dg = (dg_y.T ** dy) * (dg_x ** dx) * (1j * 2 * pi) ** (dy + dx)
                g = np.multiply(g, dg)
            self.filters.append(g)

    def compute_f(self, img_f):
        ''' Compute the scale space of an image in the fourier domain.'''
        return [np.multiply(img_f, f) for f in self.filters]

    def compute(self, img):
        ''' Compute the scale space of an image.'''
        img_f = np.fft.fft2(img)
        return [np.fft.ifft2(np.multiply(img_f, f)).real for f in self.filters]


def scale(img, sigma, dy=0, dx=0):
    '''Compute the scale-space of an image. sigma is the scale parameter. dx
    and dy specify the differentiation order along the x and y axis
    respectively.'''
    ss = ScaleSpace(img.shape, [sigma], [dy], [dx])
    return ss.compute(img)[0]


def gradient_orientation(img, sigma):
    '''Calculate gradient orientations at scale sigma.'''
    Ly = scale(img, sigma, dy=1, dx=0)
    Lx = scale(img, sigma, dy=0, dx=1)
    go = np.arctan2(Ly, Lx)
    go_m = np.sqrt(Lx**2+Ly**2)
    return go, go_m


def shape_index(img, sigma):
    '''Calculate the shape index at scale sigma.'''
    Lyy = scale(img, sigma, dy=2, dx=0)
    Lxy = scale(img, sigma, dy=1, dx=1)
    Lxx = scale(img, sigma, dy=0, dx=2)

    si = 2/np.pi * np.arctan((-Lxx-Lyy) / np.sqrt( 4*Lxy**2 + (Lxx - Lyy)**2))
    si_c = .5*np.sqrt(Lxx**2 + 2*Lxy**2 + Lyy**2)
    return si, si_c
