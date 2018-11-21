//
// Created by peter on 11/20/18.
//


typedef struct {
    unsigned short int __control_word;
    unsigned short int __glibc_reserved1;
    unsigned short int __status_word;
    unsigned short int __glibc_reserved2;
    unsigned short int __tags;
    unsigned short int __glibc_reserved3;
    unsigned int __eip;
    unsigned short int __cs_selector;
    unsigned int __opcode:11;
    unsigned int __glibc_reserved4:5;
    unsigned int __data_offset;
    unsigned short int __data_selector;
    unsigned short int __glibc_reserved5;

    unsigned int __mxcsr;

} fenv_t;

typedef int int4;
typedef union {
    int4 i[2];
    double x;
} mynumber;


static __inline __attribute__ ((__always_inline__)) void
libc_feholdexcept_setround_sse(fenv_t *e, int r) {
    unsigned int mxcsr;
    asm ("stmxcsr" " %0" : "=m" (*&mxcsr));
    e->__mxcsr = mxcsr;
    mxcsr = ((mxcsr | 0x1f80) & ~0x603f) | (r << 3);
    asm volatile ("ldmxcsr" " %0" : : "m" (*&mxcsr));
}

static __inline __attribute__ ((__always_inline__)) void
libc_fesetenv_sse(fenv_t *e) {
    asm volatile ("ldmxcsr" " %0" : : "m" (e->__mxcsr));
}


static __inline __attribute__ ((__always_inline__)) int
fegetround(void) {
#ifdef FE_TONEAREST
    return FE_TONEAREST;
#else
    return 0;
#endif
}

static const double inroot[128] = {
        1.40872145012100, 1.39792649065766, 1.38737595123859, 1.37706074531819,
        1.36697225234682, 1.35710228748795, 1.34744307370643, 1.33798721601135,
        1.32872767765984, 1.31965775814772, 1.31077107283046, 1.30206153403386,
        1.29352333352711, 1.28515092624400, 1.27693901514820, 1.26888253714903,
        1.26097664998256, 1.25321671998073, 1.24559831065844, 1.23811717205462,
        1.23076923076923, 1.22355058064300, 1.21645747403153, 1.20948631362953,
        1.20263364480453, 1.19589614840310, 1.18927063399547, 1.18275403352732,
        1.17634339535009, 1.17003587860341, 1.16382874792529, 1.15771936846787,
        1.15170520119791, 1.14578379846309, 1.13995279980655, 1.13420992801334,
        1.12855298537376, 1.12297985014975, 1.11748847323133, 1.11207687497107,
        1.10674314218572, 1.10148542531442, 1.09630193572405, 1.09119094315276,
        1.08615077328341, 1.08117980543918, 1.07627647039410, 1.07143924829188,
        1.06666666666667, 1.06195729855996, 1.05730976072814, 1.05272271193563,
        1.04819485132867, 1.04372491688551, 1.03931168393861, 1.03495396376504,
        1.03065060224133, 1.02640047855933, 1.02220250399990, 1.01805562076124,
        1.01395880083916, 1.00991104495649, 1.00591138153909, 1.00195886573624,
        0.99611649018350, 0.98848330114434, 0.98102294317595, 0.97372899112030,
        0.96659534932828, 0.95961623024651, 0.95278613468066, 0.94609983358253,
        0.93955235122353, 0.93313894963169, 0.92685511418159, 0.92069654023750,
        0.91465912076005, 0.90873893479530, 0.90293223677296, 0.89723544654727,
        0.89164514012056, 0.88615804099474, 0.88077101210109, 0.87548104826333,
        0.87028526915267, 0.86518091269740, 0.86016532891275, 0.85523597411976,
        0.85039040552437, 0.84562627613070, 0.84094132996422, 0.83633339758291,
        0.83180039185606, 0.82734030399203, 0.82295119979782, 0.81863121615464,
        0.81437855769486, 0.81019149366693, 0.80606835497581, 0.80200753138734,
        0.79800746888611, 0.79406666717674, 0.79018367731967, 0.78635709949278,
        0.78258558087123, 0.77886781361798, 0.77520253297841, 0.77158851547266,
        0.76802457717971, 0.76450957210799, 0.76104239064719, 0.75762195809661,
        0.75424723326565, 0.75091720714229, 0.74763090162560, 0.74438736831878,
        0.74118568737933, 0.73802496642311, 0.73490433947940, 0.73182296599416,
        0.72878002987884, 0.72577473860242, 0.72280632232420, 0.71987403306536,
        0.71697714391715, 0.71411494828392, 0.71128675915902, 0.70849190843208};


double ieee754_sqrt(double x) {
    static const double
            rt0 = 9.99999999859990725855365213134618E-01,
            rt1 = 4.99999999495955425917856814202739E-01,
            rt2 = 3.75017500867345182581453026130850E-01,
            rt3 = 3.12523626554518656309172508769531E-01;
    static const double big = 134217728.0;
    double y, t, del, res, res1, hy, z, zz, p, hx, tx, ty, s;
    mynumber a, c = {{0, 0}};
    int4 k;

    a.x = x;
    k = a.i[1];
    a.i[1] = (k & 0x001fffff) | 0x3fe00000;
    t = inroot[(k & 0x001fffff) >> 14];
    s = a.x;

    if (k > 0x000fffff && k < 0x7ff00000) {
        int rm = fegetround();
        fenv_t env;
        libc_feholdexcept_setround_sse(&env, 0);
        double ret;
        y = 1.0 - t * (t * s);
        t = t * (rt0 + y * (rt1 + y * (rt2 + y * rt3)));
        c.i[1] = 0x20000000 + ((k & 0x7fe00000) >> 1);
        y = t * s;
        hy = (y + big) - big;
        del = 0.5 * t * ((s - hy * hy) - (y - hy) * (y + hy));
        res = y + del;
        if (res == (res + 1.002 * ((y - res) + del)))
            ret = res * c.x;
        else {
            res1 = res + 1.5 * ((y - res) + del);
            p = 134217729.0 * (res);
            hx = ((res) - p) + p;
            tx = (res) - hx;
            p = 134217729.0 * (res1);
            hy = ((res1) - p) + p;
            ty = (res1) - hy;
            z = (res) * (res1);
            zz = (((hx * hy - z) + hx * ty) + tx * hy) + tx * ty;;
            res = ((((z - s) + zz) < 0) ? (((res1) > (res)) ? (res1) : (res)) :
                   (((res1) < (res)) ? (res1) : (res)));
            ret = res * c.x;
        }
        do {
            if (sizeof(ret) <= sizeof(double) || __builtin_types_compatible_p (__typeof(ret), _Float128))
                    __asm __volatile ("" : : "x" (ret));
            else __asm __volatile ("" : : "f" (ret));
        }
        while (0);
        libc_fesetenv_sse(&env);
        double dret = x / ret;
        if (dret != ret) {
            double force_inexact = 1.0 / 3.0;
            do {
                if (sizeof(force_inexact) <= sizeof(double) ||
                    __builtin_types_compatible_p (__typeof(force_inexact), _Float128))
                        __asm __volatile ("" : : "x" (force_inexact));
                else __asm __volatile ("" : : "f" (force_inexact));
            }
            while (0);


            switch (rm) {

                case 0x800:
                    if (dret > ret)
                        ret = (res + 0x1p-1022) * c.x;
                    break;


                case 0x400:


                case 0xc00:


                    if (dret < ret)
                        ret = (res - 0x1p-1022) * c.x;
                    break;


                default:
                    break;
            }
        }


        return ret;
    } else {
        if ((k & 0x7ff00000) == 0x7ff00000)
            return x * x + x;
        if (x == 0)
            return x;
        if (k < 0)
            return (x - x) / (x - x);
        return 0x1p-256 * ieee754_sqrt(x * 0x1p512);
    }
}
