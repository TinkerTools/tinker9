import yaml
import sys
import re


class yamlkey:
    kernel_name = 'KERNEL_NAME'
    template_params = 'TEMPLATE_PARAMS'
    constexpr_flags = 'CONSTEXPR_FLAGS'
    i_force = 'I_FORCE'
    k_force = 'K_FORCE'
    i_position = 'I_POSITION'
    k_position = 'K_POSITION'
    i_variables = 'I_VARIABLES'
    k_variables = 'K_VARIABLES'

    count = 'COUNT'
    energy = 'ENERGY'
    virial = 'VIRIAL'
    gradient = 'GRADIENT'

    cut_distance = 'CUT_DISTANCE'
    off_distance = 'OFF_DISTANCE'

    exclude_info = 'EXCLUDE_INFO'
    scale_1x_type = 'SCALE_1X_TYPE'

    extra_params = 'EXTRA_PARAMS'

    scaled_pairwise = 'SCALED_PAIRWISE_INTERACTION'
    full_pairwise = 'FULL_PAIRWISE_INTERACTION'


alphabets = {
    '0': 'a', '1': 'b', '2': 'c', '3': 'd', '4': 'e',
    '5': 'f', '6': 'g', '7': 'h', '8': 'i', '9': 'j',
    '10': 'k', '11': 'l', '12': 'm', '13': 'n', '14': 'o', '15': 'p'}


def function_parameter(ptype, pname):
    if ptype in ['int', 'real']:
        return ', {} {}'.format(ptype, pname)
    elif ptype == 'int_const_array':
        return ', const int* restrict {}'.format(pname)
    elif ptype == 'real_const_array':
        return ', const real* restrict {}'.format(pname)
    elif 'int2_const_array' in ptype:
        return ', const int (*restrict {})[2]'.format(pname)
    elif 'real2_const_array' in ptype:
        return ', const real (*restrict {})[2]'.format(pname)
    elif 'real3_const_array' in ptype:
        return ', const real (*restrict {})[3]'.format(pname)
    elif 'real4_const_array' in ptype:
        return ', const real (*restrict {})[4]'.format(pname)
    else:
        assert(False), 'Do not know how to parse type: {}'.format(ptype)


def load_scale_parameter(ptype, stem, input, separate_scaled_pairwise):
    if ptype == 'real_const_array':
        v = ''
        if input is None:
            if not separate_scaled_pairwise:
                v = 'real {}a = 1;'.format(stem)
        else:
            v = 'real {}a = {}[ii];'.format(stem, input)
        return v
    elif ptype == 'real2_const_array':
        v = ''
        if input is None:
            if not separate_scaled_pairwise:
                v = v + 'real {}a = 1;'.format(stem)
                v = v + 'real {}b = 1;'.format(stem)
        else:
            v = v + 'real {}a = {}[ii][0];'.format(stem, input)
            v = v + 'real {}b = {}[ii][1];'.format(stem, input)
        return v
    elif ptype == 'real3_const_array':
        v = ''
        if input is None:
            if not separate_scaled_pairwise:
                v = v + 'real {}a = 1;'.format(stem)
                v = v + 'real {}b = 1;'.format(stem)
                v = v + 'real {}c = 1;'.format(stem)
        else:
            v = v + 'real {}a = {}[ii][0];'.format(stem, input)
            v = v + 'real {}b = {}[ii][1];'.format(stem, input)
            v = v + 'real {}c = {}[ii][2];'.format(stem, input)
        return v
    elif ptype == 'real4_const_array':
        v = ''
        if input is None:
            if not separate_scaled_pairwise:
                v = v + 'real {}a = 1;'.format(stem)
                v = v + 'real {}b = 1;'.format(stem)
                v = v + 'real {}c = 1;'.format(stem)
                v = v + 'real {}d = 1;'.format(stem)
        else:
            v = v + 'real {}a = {}[ii][0];'.format(stem, input)
            v = v + 'real {}b = {}[ii][1];'.format(stem, input)
            v = v + 'real {}c = {}[ii][2];'.format(stem, input)
            v = v + 'real {}d = {}[ii][3];'.format(stem, input)
        return v
    else:
        pattern = re.compile('(\\D*)(\\d*)_const_array,')
        match = re.match(pattern, ptype)
        if match and len(match.groups()) == 2:
            # real4_const_array ==> (real) (4)
            t = match.group(1)
            # dim = match.group(2)
            ss = ptype.split(',')
            v = ''
            for i in range(1,len(ss)):
                idx = ss[i]
                al = alphabets[idx]
                if input is None:
                    if not separate_scaled_pairwise:
                        v = v + '{} {}{} = 1;'.format(t, stem, al)
                else:
                    v = v + '{} {}{} = {}[ii][{}];'.format(t, stem, al, input, idx)
            return v
        assert(False), 'Do not know how to load type: {}'.format(ptype)


def get_src(src, index):
    if ',' in src:
        vs = src.split(',')
        return '{}[{}][{}]'.format(vs[0], index, vs[1])
    else:
        return '{}[{}]'.format(src, index)


def ik_force(iork, op, lst):
    v = ''
    for it in lst:
        prefix = ''
        if iork == 'i':
            # prefix = 'sh'
            prefix = ''
        suffix = ''
        name = it['name']
        shared = it['location'] == 'shared'
        if op == 'decl':
            t = it['type']
            if shared:
                v = v + '__shared__ {} {}{}[BLOCK_DIM];'.format(t, prefix, name)
            else:
                v = v + '{} {}{};'.format(t, prefix, name)
        elif op == 'zero':
            if shared:
                suffix = '[threadIdx.x]'
            v = v + '{}{}{} = 0;'.format(prefix, name, suffix)
        elif op == 'init':
            if shared:
                suffix = '[threadIdx.x]'
            src = it['from']
            if ',' not in src:
                v = v + '{}{}{} = {}[{}{}];'.format(prefix, name, suffix, src, prefix, iork)
            else:
                vs = src.split(',')
                v = v + '{}{}{} = {}[{}{}][{}];'.format(prefix, name, suffix, vs[0], prefix, iork, vs[1])
        elif op == 'init_exclude_xyz':
            src = it['from']
            if iork == 'i':
                if shared:
                    v = v + '{}[klane] = {}[{}{}];'.format(name, src, prefix, iork)
                else:
                    v = v + '{} = {}[{}{}];'.format(name, src, prefix, iork)
            else:
                v = v + '{} = {}[{}{}];'.format(name, src, prefix, iork)
        elif op == 'init_exclude':
            src = it['from']
            if iork == 'i':
                if shared:
                    v = v + '{}[klane] = {};'.format(name, get_src(src, '{}{}'.format(prefix, iork)))
                else:
                    v = v + '{} = {};'.format(name, get_src(src, '{}{}'.format(prefix, iork)))
            else:
                v = v + '{} = {};'.format(name, get_src(src, '{}{}'.format(prefix, iork)))
        elif op == 'init_xyz':
            if shared:
                suffix = '[threadIdx.x]'
            v = v + '{}{}{} = sorted[{}atom{}].{};'.format(prefix, name, suffix, prefix, iork, it['from'])
        elif op == 'load':
            t = it['type']
            if shared:
                suffix = '[klane]'
            v = v + '{} {} = {}{}{};'.format(t, name, prefix, name, suffix)
        elif op == 'shfl':
            if not shared:
                v = v + '{0:}{1:} = __shfl_sync(ALL_LANES, {0:}{1:}, ilane + 1);'.format(prefix, name)
        elif op == 'save':
            if shared:
                suffix = '[threadIdx.x]'
            dst = it['addto']
            if ',' not in dst:
                v = v + 'atomic_add({}{}{}, {}, {}{});'.format(prefix, name, suffix, dst, prefix, iork)
            else:
                vs = dst.split(',')
                v = v + 'atomic_add({}{}{}, &{}[{}{}][{}]);'.format(prefix, name, suffix, vs[0], prefix, iork, vs[1])
    return v


def replace(s, d):
    output = s
    for k in d.keys():
        v = d[k]
        if v == None:
            v = ''
        output = output.replace(k, v)
    return output


def replace_i_vars(s, cfg):
    lst = []
    d = {}
    for kcfg in [yamlkey.i_position, yamlkey.i_force, yamlkey.i_variables]:
        if kcfg in cfg.keys():
            vcfg = cfg[kcfg]
            for l in vcfg:
                lst.append(l)
    for l in lst:
        n = l['name']
        shared = l['location'] == 'shared'
        old_name = '@' + n + '@'
        new_name = n
        if shared:
            new_name = n + '[klane]'
        d[old_name] = new_name
    output = replace(s, d)
    return output


kernel_template_a = '''
TEMPLATE_PARAMS
__global__
void KERNEL_NAMEa(TINKER_IMAGE_PARAMS \
    COUNT_KERNEL_PARAMS ENERGY_KERNEL_PARAMS \
    VIRIAL_KERNEL_PARAMS GRADIENT_KERNEL_PARAMS \
    CUT_KERNEL_PARAMS OFF_KERNEL_PARAMS \
    EXTRA_KERNEL_PARAMS \
    , const Spatial::SortedAtom* restrict sorted \
    , int niak, const int* restrict iak, const int* restrict lst)
{
    KERNEL_CONSTEXPR_FLAGS \
    const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
    const int iwarp = ithread / WARP_SIZE;
    const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
    const int ilane = threadIdx.x & (WARP_SIZE - 1);


    DECLARE_ZERO_LOCAL_COUNT \
    DECLARE_ZERO_LOCAL_ENERGY \
    DECLARE_ZERO_LOCAL_VIRIAL \
    DECLARE_POSITION_I_AND_K \
    DECLARE_FORCE_I_AND_K \
    DECLARE_PARAMS_I_AND_K


    for (int iw = iwarp; iw < niak; iw += nwarp) {
        KERNEL_ZERO_LOCAL_FORCE


        int ty = iak[iw];
        int atomi = ty * WARP_SIZE + ilane;
        int i = sorted[atomi].unsorted;
        int atomk = lst[iw * WARP_SIZE + ilane];
        int k = sorted[atomk].unsorted; \
        KERNEL_INIT_XYZ


        KERNEL_INIT_PARAMS_I_AND_K


        for (int j = 0; j < WARP_SIZE; ++j) {
            KERNEL_KLANE2 \
            bool incl = atomk > 0; \
            KERNEL_SCALE_1 \
            KERNEL_FULL_PAIRWISE_INTERACTION


            KERNEL_SHUFFLE_XYZ_I \
            KERNEL_SHUFFLE_PARAMS_I \
            KERNEL_SHUFFLE_LOCAL_FORCE
        }


        KERNEL_SAVE_LOCAL_FORCE
    }


    KERNEL_SUM_COUNT KERNEL_SUM_ENERGY KERNEL_SUM_VIRIAL
}
'''


kernel_template_b = '''
TEMPLATE_PARAMS
__global__
void KERNEL_NAMEb(TINKER_IMAGE_PARAMS \
    COUNT_KERNEL_PARAMS ENERGY_KERNEL_PARAMS \
    VIRIAL_KERNEL_PARAMS GRADIENT_KERNEL_PARAMS \
    CUT_KERNEL_PARAMS OFF_KERNEL_PARAMS \
    EXTRA_KERNEL_PARAMS \
    EXCLUDE_INFO_KERNEL_PARAMS \
    , const Spatial::SortedAtom* restrict sorted \
    , int n, int nakpl, const int* restrict iakpl)
{
    KERNEL_CONSTEXPR_FLAGS \
    const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
    const int iwarp = ithread / WARP_SIZE;
    const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
    const int ilane = threadIdx.x & (WARP_SIZE - 1);


    DECLARE_ZERO_LOCAL_COUNT \
    DECLARE_ZERO_LOCAL_ENERGY \
    DECLARE_ZERO_LOCAL_VIRIAL \
    DECLARE_POSITION_I_AND_K \
    DECLARE_FORCE_I_AND_K \
    DECLARE_PARAMS_I_AND_K


    for (int iw = iwarp; iw < nakpl; iw += nwarp) {
        KERNEL_ZERO_LOCAL_FORCE


        int tri, tx, ty;
        tri = iakpl[iw];
        tri_to_xy(tri, tx, ty);


        int iid = ty * WARP_SIZE + ilane;
        int atomi = min(iid, n - 1);
        int i = sorted[atomi].unsorted;
        int kid = tx * WARP_SIZE + ilane;
        int atomk = min(kid, n - 1);
        int k = sorted[atomk].unsorted;
        KERNEL_INIT_XYZ


        KERNEL_INIT_PARAMS_I_AND_K


        KERNEL_LOAD_INFO_VARIABLES \
        for (int j = 0; j < WARP_SIZE; ++j) {
            int srclane = (ilane + j) & (WARP_SIZE - 1);
            int srcmask = 1 << srclane;
            KERNEL_KLANE1 \
            bool incl = iid < kid and kid < n; \
            KERNEL_EXCLUDE_BIT \
            KERNEL_SCALE_1 \
            KERNEL_FULL_PAIRWISE_INTERACTION


            iid = __shfl_sync(ALL_LANES, iid, ilane + 1); \
            KERNEL_SHUFFLE_XYZ_I \
            KERNEL_SHUFFLE_PARAMS_I \
            KERNEL_SHUFFLE_LOCAL_FORCE
        }


        KERNEL_SAVE_LOCAL_FORCE
    }


    KERNEL_SUM_COUNT KERNEL_SUM_ENERGY KERNEL_SUM_VIRIAL
}
'''


kernel_template_c = '''
TEMPLATE_PARAMS
__global__
void KERNEL_NAMEc(TINKER_IMAGE_PARAMS \
    COUNT_KERNEL_PARAMS ENERGY_KERNEL_PARAMS \
    VIRIAL_KERNEL_PARAMS GRADIENT_KERNEL_PARAMS \
    CUT_KERNEL_PARAMS OFF_KERNEL_PARAMS \
    EXTRA_KERNEL_PARAMS \
    EXCLUDE_SCALE_KERNEL_PARAMS)
{
    KERNEL_CONSTEXPR_FLAGS \
    const int ithread = threadIdx.x + blockIdx.x * blockDim.x;


    DECLARE_ZERO_LOCAL_COUNT \
    DECLARE_ZERO_LOCAL_ENERGY \
    DECLARE_ZERO_LOCAL_VIRIAL \
    DECLARE_POSITION_I_AND_K \
    DECLARE_FORCE_I_AND_K \
    DECLARE_PARAMS_I_AND_K


    for (int ii = ithread; ii < nexclude; ii += blockDim.x * gridDim.x) {
        KERNEL_SCALED_KLANE \
        KERNEL_ZERO_LOCAL_FORCE


        int i = exclude[ii][0];
        int k = exclude[ii][1];
        KERNEL_LOAD_1X_SCALES


        KERNEL_INIT_EXCLUDE_XYZ
        KERNEL_INIT_EXCLUDE_PARAMS_I_AND_K


        constexpr bool incl = true;
        KERNEL_SCALED_PAIRWISE_INTERACTION


        KERNEL_SAVE_LOCAL_FORCE
    }


    KERNEL_SUM_COUNT KERNEL_SUM_ENERGY KERNEL_SUM_VIRIAL
}
'''


if __name__ == '__main__':
    yaml_file = sys.argv[1]
    use_klane = False
    with open(yaml_file) as input_file:
        config = yaml.full_load(input_file)
        for k in config.keys():
            lst = config[k]
            if isinstance(lst, list):
                for d in lst:
                    if isinstance(d, dict):
                        if 'def' in d.keys():
                            s = d['def']
                            vs = s.split()
                            d['location'] = vs[0]
                            d['type'] = vs[1]
                            d['name'] = vs[2]
                            for i in range(3, len(vs)):
                                s2 = vs[i]
                                if s2.startswith('from:'):
                                    s3 = s2.replace('from:', '')
                                    d['from'] = s3
                                elif s2.startswith('addto:'):
                                    s3 = s2.replace('addto:', '')
                                    d['addto'] = s3
        k = 'I_FORCE'
        if k in config.keys():
            lst = config[k]
            for d in lst:
                k2 = 'location'
                if k2 in d.keys():
                    v2 = d[k2]
                    if v2 == 'shared':
                        use_klane = True


    d = {}


    k, v = 'KERNEL_NAME', config[yamlkey.kernel_name]
    d[k] = v


    k, v = 'TEMPLATE_PARAMS', ''
    kcfg = yamlkey.template_params
    if kcfg in config.keys():
        v = config[kcfg]
    d[k] = v


    k, v = 'KERNEL_CONSTEXPR_FLAGS', ''
    kcfg = yamlkey.constexpr_flags
    if kcfg in config.keys():
        v = config[kcfg]
    d[k] = v


    k, v = 'EXTRA_KERNEL_PARAMS', ''
    kcfg = yamlkey.extra_params
    if kcfg in config.keys():
        v = config[kcfg]
    d[k] = v


    # klane
    k, v = 'KERNEL_KLANE1', ''
    k2, v2 = 'KERNEL_KLANE2', ''
    k3, v3 = 'KERNEL_SCALED_KLANE', ''
    if use_klane:
        v = v + 'int klane = srclane + threadIdx.x - ilane;'
        v2 = v2 + 'int srclane = (ilane + j) & (WARP_SIZE - 1);'
        v2 = v2 + 'int klane = srclane + threadIdx.x - ilane;'
        v3 = 'const int klane = threadIdx.x;'
    d[k], d[k2], d[k3] = v, v2, v3


    # count
    k, v = 'COUNT_KERNEL_PARAMS', ''
    k2, v2 = 'DECLARE_ZERO_LOCAL_COUNT', ''
    k3, v3 = 'KERNEL_SUM_COUNT', ''
    kcfg = yamlkey.count
    if kcfg in config.keys():
        vcfg, decl, zero, total = config[kcfg], '', '', ''
        for t in vcfg:
            v = v + ', count_buffer restrict {}'.format(t)
            decl = decl + 'int {}tl;'.format(t)
            zero = zero + '{}tl = 0;'.format(t)
            total = total + 'atomic_add({}tl, {}, ithread);'.format(t, t)
        v2 = '%s if CONSTEXPR (do_a) {%s}' % (decl, zero)
        v3 = 'if CONSTEXPR (do_a) {%s}' % total
    d[k] = v
    d[k2] = v2
    d[k3] = v3


    # energy
    k, v = 'ENERGY_KERNEL_PARAMS', ''
    k2, v2 = 'DECLARE_ZERO_LOCAL_ENERGY', ''
    k3, v3 = 'KERNEL_SUM_ENERGY', ''
    kcfg = yamlkey.energy
    if kcfg in config.keys():
        vcfg, decl, zero, total = config[kcfg], '', '', ''
        for t in vcfg:
            v = v + ', energy_buffer restrict {}'.format(t)
            decl = decl + 'ebuf_prec {}tl;'.format(t)
            zero = zero + '{}tl = 0;'.format(t)
            total = total + 'atomic_add({}tl, {}, ithread);'.format(t, t)
        v2 = 'using ebuf_prec = energy_buffer_traits::type;\n'
        v2 = v2 + '%s if CONSTEXPR (do_e) {%s}' % (decl, zero)
        v3 = 'if CONSTEXPR (do_e) {%s}' % total
    d[k] = v
    d[k2] = v2
    d[k3] = v3


    # virial
    k, v = 'VIRIAL_KERNEL_PARAMS', ''
    k2, v2 = 'DECLARE_ZERO_LOCAL_VIRIAL', ''
    k3, v3 = 'KERNEL_SUM_VIRIAL', ''
    kcfg = yamlkey.virial
    if kcfg in config.keys():
        vcfg, decl, zero, total = config[kcfg], '', '', ''
        for t in vcfg:
            v = v + ', virial_buffer restrict {}'.format(t)
            decl = decl + 'vbuf_prec {}tlxx'.format(t)
            zero = zero + '{}tlxx = 0;'.format(t)
            total = total + 'atomic_add({}tlxx'.format(t)
            for sufx in ['yx', 'zx', 'yy', 'zy', 'zz']:
                decl = decl + ', {}tl{}'.format(t, sufx)
                zero = zero + '{}tl{} = 0;'.format(t, sufx)
                total = total + ', {}tl{}'.format(t, sufx)
            decl = decl + ';'
            zero = zero + '\n'
            total = total + ', {}, ithread);'.format(t)
        v2 = 'using vbuf_prec = virial_buffer_traits::type;\n'
        v2 = v2 + '%s if CONSTEXPR (do_v) {%s}' % (decl, zero)
        v3 = 'if CONSTEXPR (do_v) {%s}' % total
    d[k] = v
    d[k2] = v2
    d[k3] = v3


    # gradient
    k, v = 'GRADIENT_KERNEL_PARAMS', ''
    kcfg = yamlkey.gradient
    if kcfg in config.keys():
        vcfg = config[kcfg]
        for t in vcfg:
            v = v + ', grad_prec* restrict {}'.format(t)
    k0, v0 = 'DECLARE_FORCE_I_AND_K', ''
    k1, v1 = 'KERNEL_ZERO_LOCAL_FORCE', ''
    k2, v2 = 'KERNEL_SAVE_LOCAL_FORCE', ''
    k3, v3 = 'KERNEL_SHUFFLE_LOCAL_FORCE', ''
    kcfg = yamlkey.i_force
    if kcfg in config.keys():
        vcfg = config[kcfg]
        v0 = v0 + ik_force('i', 'decl', vcfg)
        v1 = v1 + ik_force('i', 'zero', vcfg)
        v2 = v2 + ik_force('i', 'save', vcfg)
        v3 = v3 + ik_force('i', 'shfl', vcfg)
    kcfg = yamlkey.k_force
    if kcfg in config.keys():
        vcfg = config[kcfg]
        v0 = v0 + ik_force('k', 'decl', vcfg)
        v1 = v1 + ik_force('k', 'zero', vcfg)
        v2 = v2 + ik_force('k', 'save', vcfg)
    kcfg = yamlkey.constexpr_flags
    if kcfg in config.keys():
        vcfg = config[kcfg]
        if 'constexpr bool do_g =' in vcfg:
            v1 = 'if CONSTEXPR (do_g) {%s}' % v1
            v2 = 'if CONSTEXPR (do_g) {%s}' % v2
    if v3 != '':
        v3 = 'if CONSTEXPR (do_g) {%s}' % v3
    d[k], d[k0], d[k1], d[k2], d[k3] = v, v0, v1, v2, v3


    # cut, off
    k, v = 'CUT_KERNEL_PARAMS', ''
    kcfg = yamlkey.cut_distance
    if kcfg in config.keys():
        vcfg = config[kcfg]
        for t in vcfg:
            v = v + ', real {}'.format(t)
    d[k] = v


    k, v = 'OFF_KERNEL_PARAMS', ''
    kcfg = yamlkey.off_distance
    if kcfg in config.keys():
        vcfg = config[kcfg]
        for t in vcfg:
            v = v + ', real {}'.format(t)
    d[k] = v


    # i and k pos
    k, v = 'DECLARE_POSITION_I_AND_K', ''
    k1, v1 = 'KERNEL_INIT_XYZ', ''
    k2, v2 = 'KERNEL_INIT_EXCLUDE_XYZ', ''
    # k3, v3 = 'KERNEL_LOAD_XYZ_I', ''
    k4, v4 = 'KERNEL_SHUFFLE_XYZ_I', ''
    kcfg = yamlkey.i_position
    if kcfg in config.keys():
        vcfg = config[kcfg]
        v = v + ik_force('i', 'decl', vcfg)
        v1 = v1 + ik_force('i', 'init_xyz', vcfg)
        v2 = v2 + ik_force('i', 'init_exclude_xyz', vcfg)
        # v3 = v3 + ik_force('i', 'load', vcfg)
        v4 = v4 + ik_force('i', 'shfl', vcfg)
    kcfg = yamlkey.k_position
    if kcfg in config.keys():
        vcfg = config[kcfg]
        v = v + ik_force('k', 'decl', vcfg)
        v1 = v1 + ik_force('k', 'init_xyz', vcfg)
        v2 = v2 + ik_force('k', 'init_exclude_xyz', vcfg)
    # d[k], d[k1], d[k2], d[k3], d[k4] = v, v1, v2, v3, v4
    d[k], d[k1], d[k2], d[k4] = v, v1, v2, v4


    # i and k variable
    k, v = 'DECLARE_PARAMS_I_AND_K', ''
    k2, v2 = 'KERNEL_INIT_PARAMS_I_AND_K', ''
    # k3, v3 = 'KERNEL_LOAD_PARAMS_I', ''
    k4, v4 = 'KERNEL_SHUFFLE_PARAMS_I', ''
    k5, v5 = 'KERNEL_INIT_EXCLUDE_PARAMS_I_AND_K', ''
    kcfg = yamlkey.i_variables
    if kcfg in config.keys():
        vcfg = config[kcfg]
        v = v + ik_force('i', 'decl', vcfg)
        v2 = v2 + ik_force('i', 'init', vcfg)
        # v3 = v3 + ik_force('i', 'load', vcfg)
        v4 = v4 + ik_force('i', 'shfl', vcfg)
        v5 = v5 + ik_force('i', 'init_exclude', vcfg)
    kcfg = yamlkey.k_variables
    if kcfg in config.keys():
        vcfg = config[kcfg]
        v = v + ik_force('k', 'decl', vcfg)
        v2 = v2 + ik_force('k', 'init', vcfg)
        v5 = v5 + ik_force('k', 'init_exclude', vcfg)
    # d[k], d[k2], d[k3], d[k4], d[k5] = v, v2, v3, v4, v5
    d[k], d[k2], d[k4], d[k5] = v, v2, v4, v5


    ############################################################################
    # exclude


    k, v = 'EXCLUDE_INFO_KERNEL_PARAMS', ''
    k2, v2 = 'KERNEL_LOAD_INFO_VARIABLES', ''
    k3, v3 = 'KERNEL_EXCLUDE_BIT', ''
    kcfg = yamlkey.exclude_info
    if kcfg in config.keys():
        vcfg = config[kcfg]
        for t in vcfg:
            v = v + ', const unsigned* restrict {}'.format(t)
            v2 = v2 + 'unsigned int {0:}0 = {0:}[iw * WARP_SIZE + ilane];'.format(t)
            v3 = v3 + ' and ({}0 & srcmask) == 0'.format(t)
    if v3 != 0:
        v3 = 'incl = incl {};'.format(v3)
    d[k], d[k2], d[k3] = v, v2, v3


    k, v = 'EXCLUDE_SCALE_KERNEL_PARAMS', ''
    k1, v1 = 'KERNEL_LOAD_1X_SCALES', ''
    k3, v3 = 'KERNEL_SCALE_1', ''
    kcfg = yamlkey.scale_1x_type
    separate_scaled_pairwise = yamlkey.scaled_pairwise in config.keys()
    if kcfg in config.keys():
        vcfg = config[kcfg]
        v = ', int nexclude, const int (*restrict exclude)[2]'
        v = v + function_parameter(vcfg, 'exclude_scale')
        v1 = load_scale_parameter(vcfg, 'scale', 'exclude_scale', separate_scaled_pairwise)
        v3 = load_scale_parameter(vcfg, 'scale', None, separate_scaled_pairwise)
    v = v + ', const real* restrict x'
    v = v + ', const real* restrict y'
    v = v + ', const real* restrict z'
    d[k], d[k1], d[k3] = v, v1, v3


    k, v = 'KERNEL_FULL_PAIRWISE_INTERACTION', ''
    k1, v1 = 'KERNEL_SCALED_PAIRWISE_INTERACTION', ''
    kcfg = yamlkey.full_pairwise
    if kcfg in config.keys():
        vcfg = config[kcfg]
        v = v + vcfg
        v1 = v
    kcfg = yamlkey.scaled_pairwise
    if kcfg in config.keys():
        vcfg = config[kcfg]
        v1 = vcfg
    d[k], d[k1] = v, v1


    ############################################################################
    # output


    print('// ck.py Version 2.0.0')
    print('\n\n')


    if yamlkey.scale_1x_type in config.keys():
        output = replace(kernel_template_c, d)
        output = replace_i_vars(output, config)
        print(output)
        print('\n\n')


    output = replace(kernel_template_b, d)
    output = replace_i_vars(output, config)
    print(output)
    print('\n\n')


    output = replace(kernel_template_a, d)
    output = replace_i_vars(output, config)
    print(output)
    print('\n\n')
