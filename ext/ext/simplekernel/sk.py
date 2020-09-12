#!/usr/bin/env python


import yaml
import sys


kernal_body = R'''
template <TPARAMS>
__global__
void KERNEL_NAME(count_buffer restrict nebuf
, energy_buffer restrict ebuf, virial_buffer restrict vbuf
, grad_prec* restrict gx, grad_prec* restrict gy
, grad_prec* restrict gz, TINKER_IMAGE_PARAMS, real cut, real off
, const real* restrict x, const real* restrict y, const real* restrict z
, int n, const Spatial::SortedAtom* restrict sorted
, int nakpl, const int* restrict iakpl
, int niak, const int* restrict iak, const int* restrict lst
, int nexclude, const int (*restrict exclude)[2]
, const real* restrict exclude_scale, Spatial2::ScaleInfo info
GLOBAL_VARIABLES
GLOBAL_ARRAYS
PER_ATOM_ARRAYS
)
{
    constexpr bool do_e = Ver::e;
    constexpr bool do_a = Ver::a;
    constexpr bool do_g = Ver::g;
    constexpr bool do_v = Ver::v;


    const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
    const int iwarp = ithread / WARP_SIZE;
    const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
    const int ilane = threadIdx.x & (WARP_SIZE - 1);


    using ebuf_prec = energy_buffer_traits::type;
    using vbuf_prec = virial_buffer_traits::type;
    ebuf_prec etl;
    vbuf_prec vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz;
    if CONSTEXPR (do_e) {
        etl = 0;
    }
    if CONSTEXPR (do_v) {
        vtlxx = 0;
        vtlyx = 0;
        vtlzx = 0;
        vtlyy = 0;
        vtlzy = 0;
        vtlzz = 0;
    }
    int ctl;
    if CONSTEXPR (do_a) {
        ctl = 0;
    }


    //* /
    // exclude
    for (int ii = ithread; ii < nexclude; ii += blockDim.x * gridDim.x) {
        int i = exclude[ii][0];
        int k = exclude[ii][1];
        real scale = exclude_scale[ii];


        LOAD_ATOM_I_PARAMS
        real xi = x[i];
        real yi = y[i];
        real zi = z[i];


        LOAD_ATOM_K_PARAMS
        real xr = xi - x[k];
        real yr = yi - y[k];
        real zr = zi - z[k];


        PAIRWISE_SCALED


        if CONSTEXPR (do_e) {
            etl += cvt_to<ebuf_prec>(e);
            if CONSTEXPR (do_a) {
                if (e != 0) {
                    ctl += 1;
                }
            }
        }
        if CONSTEXPR (do_g) {
            real dedx, dedy, dedz;
            de *= invr;
            dedx = de * xr;
            dedy = de * yr;
            dedz = de * zr;
            atomic_add(dedx, gx, i);
            atomic_add(dedy, gy, i);
            atomic_add(dedz, gz, i);
            atomic_add(-dedx, gx, k);
            atomic_add(-dedy, gy, k);
            atomic_add(-dedz, gz, k);
            if CONSTEXPR (do_v) {
                vtlxx += cvt_to<vbuf_prec>(xr * dedx);
                vtlyx += cvt_to<vbuf_prec>(yr * dedx);
                vtlzx += cvt_to<vbuf_prec>(zr * dedx);
                vtlyy += cvt_to<vbuf_prec>(yr * dedy);
                vtlzy += cvt_to<vbuf_prec>(zr * dedy);
                vtlzz += cvt_to<vbuf_prec>(zr * dedz);
            }
        }
    }
    // */


    //* /
    // block pairs that have scale factors
    for (int iw = iwarp; iw < nakpl; iw += nwarp) {
        real fix, fiy, fiz;
        real fkx, fky, fkz;
        if CONSTEXPR (do_g) {
            fix = 0;
            fiy = 0;
            fiz = 0;
            fkx = 0;
            fky = 0;
            fkz = 0;
        }


        int tri, tx, ty;
        tri = iakpl[iw];
        tri_to_xy(tri, tx, ty);


        int shiid = ty * WARP_SIZE + ilane;
        int shatomi = min(shiid, n - 1);
        real shxi = sorted[shatomi].x;
        real shyi = sorted[shatomi].y;
        real shzi = sorted[shatomi].z;
        int shi = sorted[shatomi].unsorted;
        LOAD_ATOM_SHI_PARAMS


        int kid = tx * WARP_SIZE + ilane;
        int atomk = min(kid, n - 1);
        real xk = sorted[atomk].x;
        real yk = sorted[atomk].y;
        real zk = sorted[atomk].z;
        int k = sorted[atomk].unsorted;
        LOAD_ATOM_K_PARAMS


        int bit0 = info.bit0[iw * WARP_SIZE + ilane];
        for (int j = 0; j < WARP_SIZE; ++j) {
            int srclane = (ilane + j) & (WARP_SIZE - 1);
            int srcmask = 1 << srclane;
            int bit = bit0 & srcmask;
            LOAD_SHFL_PARAMS
            int iid = shiid;
            real xr = shxi - xk;
            real yr = shyi - yk;
            real zr = shzi - zk;


            bool incl = iid < kid and kid < n and bit == 0;
            PAIRWISE_CODE


            if CONSTEXPR (do_e) {
                etl += incl ? cvt_to<ebuf_prec>(e) : 0;
                if CONSTEXPR (do_a) {
                    if (incl and e != 0) {
                        ctl += 1;
                    }
                }
            }
            if CONSTEXPR (do_g) {
                real dedx, dedy, dedz;
                de = incl ? de * invr : 0;
                dedx = de * xr;
                dedy = de * yr;
                dedz = de * zr;
                fix += dedx;
                fiy += dedy;
                fiz += dedz;
                fkx -= dedx;
                fky -= dedy;
                fkz -= dedz;
                if CONSTEXPR (do_v) {
                    vtlxx += cvt_to<vbuf_prec>(xr * dedx);
                    vtlyx += cvt_to<vbuf_prec>(yr * dedx);
                    vtlzx += cvt_to<vbuf_prec>(zr * dedx);
                    vtlyy += cvt_to<vbuf_prec>(yr * dedy);
                    vtlzy += cvt_to<vbuf_prec>(zr * dedy);
                    vtlzz += cvt_to<vbuf_prec>(zr * dedz);
                }
            }


            SHFL_ATOM_I_PARAMS
            shiid = __shfl_sync(ALL_LANES, shiid, ilane + 1);
            shxi = __shfl_sync(ALL_LANES, shxi, ilane + 1);
            shyi = __shfl_sync(ALL_LANES, shyi, ilane + 1);
            shzi = __shfl_sync(ALL_LANES, shzi, ilane + 1);
            fix = __shfl_sync(ALL_LANES, fix, ilane + 1);
            fiy = __shfl_sync(ALL_LANES, fiy, ilane + 1);
            fiz = __shfl_sync(ALL_LANES, fiz, ilane + 1);
        }


        if CONSTEXPR (do_g) {
            atomic_add(fix, gx, shi);
            atomic_add(fiy, gy, shi);
            atomic_add(fiz, gz, shi);
            atomic_add(fkx, gx, k);
            atomic_add(fky, gy, k);
            atomic_add(fkz, gz, k);
        }
    } // end loop block pairs
    // */


    //* /
    // block-atoms
    for (int iw = iwarp; iw < niak; iw += nwarp) {
        real fix, fiy, fiz;
        real fkx, fky, fkz;
        if CONSTEXPR (do_g) {
            fix = 0;
            fiy = 0;
            fiz = 0;
            fkx = 0;
            fky = 0;
            fkz = 0;
        }


        int ty = iak[iw];
        int shatomi = ty * WARP_SIZE + ilane;
        real shxi = sorted[shatomi].x;
        real shyi = sorted[shatomi].y;
        real shzi = sorted[shatomi].z;
        int shi = sorted[shatomi].unsorted;
        LOAD_ATOM_SHI_PARAMS


        int atomk = lst[iw * WARP_SIZE + ilane];
        real xk = sorted[atomk].x;
        real yk = sorted[atomk].y;
        real zk = sorted[atomk].z;
        int k = sorted[atomk].unsorted;
        LOAD_ATOM_K_PARAMS


        for (int j = 0; j < WARP_SIZE; ++j) {
            LOAD_SHFL_PARAMS
            real xr = shxi - xk;
            real yr = shyi - yk;
            real zr = shzi - zk;


            bool incl = atomk > 0;
            PAIRWISE_CODE


            if CONSTEXPR (do_e) {
                etl += incl ? cvt_to<ebuf_prec>(e) : 0;
                if CONSTEXPR (do_a) {
                    if (incl and e != 0) {
                        ctl += 1;
                    }
                }
            }
            if CONSTEXPR (do_g) {
                real dedx, dedy, dedz;
                de = incl ? de * invr : 0;
                dedx = de * xr;
                dedy = de * yr;
                dedz = de * zr;
                fix += dedx;
                fiy += dedy;
                fiz += dedz;
                fkx -= dedx;
                fky -= dedy;
                fkz -= dedz;
                if CONSTEXPR (do_v) {
                    vtlxx += cvt_to<vbuf_prec>(xr * dedx);
                    vtlyx += cvt_to<vbuf_prec>(yr * dedx);
                    vtlzx += cvt_to<vbuf_prec>(zr * dedx);
                    vtlyy += cvt_to<vbuf_prec>(yr * dedy);
                    vtlzy += cvt_to<vbuf_prec>(zr * dedy);
                    vtlzz += cvt_to<vbuf_prec>(zr * dedz);
                }
            }


            SHFL_ATOM_I_PARAMS
            shxi = __shfl_sync(ALL_LANES, shxi, ilane + 1);
            shyi = __shfl_sync(ALL_LANES, shyi, ilane + 1);
            shzi = __shfl_sync(ALL_LANES, shzi, ilane + 1);
            fix = __shfl_sync(ALL_LANES, fix, ilane + 1);
            fiy = __shfl_sync(ALL_LANES, fiy, ilane + 1);
            fiz = __shfl_sync(ALL_LANES, fiz, ilane + 1);
        }


        if CONSTEXPR (do_g) {
            atomic_add(fix, gx, shi);
            atomic_add(fiy, gy, shi);
            atomic_add(fiz, gz, shi);
            atomic_add(fkx, gx, k);
            atomic_add(fky, gy, k);
            atomic_add(fkz, gz, k);
        }
    } // end loop block-atoms
    // */


    if CONSTEXPR (do_a) {
        atomic_add(ctl, nebuf, ithread);
    }
    if CONSTEXPR (do_e) {
        atomic_add(etl, ebuf, ithread);
    }
    if CONSTEXPR (do_v) {
        atomic_add(vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz, vbuf, ithread);
    }
}
'''


def shft_per_atom_params(prefix, d):
    out = ''
    for k in sorted(d):
        e = d[k]
        out = out + '%s%s = __shfl_sync(ALL_LANES, %s%s, ilane + 1);' % (prefix, e['name'],prefix, e['name'])
    return out

def load_shfl_params(newprefix, oldprefix, d):
    out = ''
    for k in sorted(d):
        e = d[k]
        out = out + '%s %s%s = %s%s;' % (e['type'], newprefix, e['name'], oldprefix, e['name'])
    return out

def load_per_atom_params(prefix, idx, d):
    ''' type prefix_stem = array[idx];'''
    out = ''
    for k in sorted(d):
        e = d[k]
        out = out + '%s %s%s = %s[%s];' % (e['type'], prefix, e['name'], e['name'], idx)
    return out


def generate(yaml_file):
    with open(yaml_file) as file:
        config = yaml.full_load(file)

    kernal_name = config['KERNEL_NAME']

    global_variables = ''
    d_gvar = config['GLOBAL_VARIABLES']
    for k in sorted(d_gvar):
        e = d_gvar[k]
        global_variables = global_variables + ', %s %s' % (e['param_type'], e['name'])

    global_arrays = ''
    d_gar = config['GLOBAL_ARRAYS']
    for k in sorted(d_gar):
        e = d_gar[k]
        global_arrays = global_arrays + ', %s %s' % (e['param_type'], e['name'])

    per_atom_arrays = ''
    d_paa = config['PER_ATOM_ARRAYS']
    for k in sorted(d_paa):
        e = d_paa[k]
        per_atom_arrays = per_atom_arrays + ', %s %s' % (e['param_type'], e['name'])

    output = kernal_body\
        .replace('TPARAMS', config['TPARAMS'])\
        .replace('KERNEL_NAME', kernal_name)\
        .replace('GLOBAL_VARIABLES', global_variables)\
        .replace('GLOBAL_ARRAYS', global_arrays)\
        .replace('PER_ATOM_ARRAYS', per_atom_arrays)\
        .replace('LOAD_ATOM_I_PARAMS', load_per_atom_params('i', 'i', d_paa))\
        .replace('LOAD_ATOM_K_PARAMS', load_per_atom_params('k', 'k', d_paa))\
        .replace('LOAD_ATOM_SHI_PARAMS', load_per_atom_params('shi', 'shi', d_paa))\
        .replace('LOAD_SHFL_PARAMS', load_shfl_params('i', 'shi', d_paa))\
        .replace('SHFL_ATOM_I_PARAMS', shft_per_atom_params('shi', d_paa))\
        .replace('PAIRWISE_CODE', config['PAIRWISE_CODE'])\
        .replace('PAIRWISE_SCALED', config['PAIRWISE_SCALED'])

    print(output)


if __name__ == '__main__':
    generate(sys.argv[1])
