#handles instrument functions and noise simulation
#baseline, offset, wavenumber shifts?

# TODO: Add support for polarized radiation

import hapi
import numpy as np
from arsca_datatf import wl2wn, wn2wl, perp_norm_vec, arb_rotation
from arsca_utilities import linear_to_multi_idx

def convolve(wl,spec,reso):
    wn = wl2wn(wl)
    spec_wn = spec[::-1]
    reso_wn = wl2wn((np.mean(wl) - reso / 2.0)) - wl2wn((np.mean(wl) + reso / 2.0))
    wing_width_wn = reso_wn
    wn_,spec_wn_,i,j,slit = hapi.convolveSpectrum(wn,spec_wn,Resolution=reso_wn,AF_wing=wing_width_wn,
                                                SlitFunction=hapi.SLIT_GAUSSIAN)
    wl_ = wn2wl(wn_)
    spec_ = spec_wn_[::-1]

    return (wl_,spec_)

def camera_fov(image_size,fov_angle,view_center,up_vec):
    if fov_angle[0] > 2 * np.pi or fov_angle[1] > 2 * np.pi:
        print("FoV angles larger than full circle! Check if they are in radians!")
    n_fov = image_size[0] * image_size[1]
    #image coordinates are the traditional ones with
    vert_ang_step = fov_angle[0] / image_size[0]
    horz_ang_step = fov_angle[1] / image_size[1]
    view_vecs = np.zeros((n_fov,3))
    side_vec = perp_norm_vec(view_center,up_vec)

    # after this, start_vec will point to the pixel in the upper left corner of
    # the image
    start_vec_ = arb_rotation(view_center, -fov_angle[0] / 2 + vert_ang_step / 2, side_vec)
    start_vec = arb_rotation(start_vec_, -fov_angle[1] / 2 + horz_ang_step / 2, up_vec)

    for i_vec in range(n_fov):
        # get the image coordinates
        # at this point we decide how the image gets indexed
        # The image indices are row-major ordered:
        # (0,0) = 0, (0,1) = 1, ..., (1,0) = image_size[1] and so on.
        (i,j) = linear_to_multi_idx(i_vec, image_size)
        vert_ang = i * vert_ang_step
        horz_ang = j * horz_ang_step
        curr_vec_ = arb_rotation(start_vec, horz_ang, up_vec)
        curr_vec = arb_rotation(curr_vec_, vert_ang, side_vec)
        view_vecs[i_vec,:] = curr_vec
    return view_vecs
