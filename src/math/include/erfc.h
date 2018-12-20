//
// Created by peter on 12/19/18.
//

#ifndef PRICER_ERFC_H
#define PRICER_ERFC_H


typedef struct {
    double x, y;
} Sleef_double2;

typedef long long int64_t;

static inline int64_t __attribute__((always_inline)) doubleToRawLongBits(double d) {
    union {
        double f;
        int64_t i;
    } tmp;
    tmp.f = d;
    return tmp.i;
}

static inline double __attribute__((always_inline)) longBitsToDouble(int64_t i) {
    union {
        double f;
        int64_t i;
    } tmp;
    tmp.i = i;
    return tmp.f;
}

static inline Sleef_double2 __attribute__((always_inline)) dd(double h, double l) {
    Sleef_double2 ret;
    ret.x = h;
    ret.y = l;
    return ret;
}

static inline double __attribute__((always_inline)) fabsk(double x) {
    return longBitsToDouble(0x7fffffffffffffffLL & doubleToRawLongBits(x));
}

static inline double __attribute__((always_inline)) mla(double x, double y, double z) { return x * y + z; }

static inline double __attribute__((always_inline)) rintk(double x) {
    return x < 0 ? (int) (x - 0.5) : (int) (x + 0.5);
}

static inline int __attribute__((always_inline)) xisnan(double x) { return x != x; }

static inline double __attribute__((always_inline)) upper(double d) {
    return longBitsToDouble(doubleToRawLongBits(d) & 0xfffffffff8000000LL);
}

static inline Sleef_double2 __attribute__((always_inline)) ddmul_d2_d_d(double x, double y) {
    double xh = upper(x), xl = x - xh;
    double yh = upper(y), yl = y - yh;
    Sleef_double2 r;

    r.x = x * y;
    r.y = xh * yh - r.x + xl * yh + xh * yl + xl * yl;

    return r;
}

static inline Sleef_double2 __attribute__((always_inline)) dddiv_d2_d2_d2(Sleef_double2 n, Sleef_double2 d) {
    double t = 1.0 / d.x;
    double dh = upper(d.x), dl = d.x - dh;
    double th = upper(t), tl = t - th;
    double nhh = upper(n.x), nhl = n.x - nhh;

    Sleef_double2 q;

    q.x = n.x * t;

    double u = -q.x + nhh * th + nhh * tl + nhl * th + nhl * tl +
               q.x * (1 - dh * th - dh * tl - dl * th - dl * tl);

    q.y = t * (n.y - q.x * d.y) + u;

    return q;
}

static inline Sleef_double2 __attribute__((always_inline)) ddmul_d2_d2_d(Sleef_double2 x, double y) {
    double xh = upper(x.x), xl = x.x - xh;
    double yh = upper(y), yl = y - yh;
    Sleef_double2 r;

    r.x = x.x * y;
    r.y = xh * yh - r.x + xl * yh + xh * yl + xl * yl + x.y * y;

    return r;
}

static inline Sleef_double2 __attribute__((always_inline)) ddmul_d2_d2_d2(Sleef_double2 x, Sleef_double2 y) {
    double xh = upper(x.x), xl = x.x - xh;
    double yh = upper(y.x), yl = y.x - yh;
    Sleef_double2 r;

    r.x = x.x * y.x;
    r.y = xh * yh - r.x + xl * yh + xh * yl + xl * yl + x.x * y.y + x.y * y.x;

    return r;
}


static inline Sleef_double2 __attribute__((always_inline)) ddadd2_d2_d2_d2(Sleef_double2 x, Sleef_double2 y) {
    Sleef_double2 r;

    r.x = x.x + y.x;
    double v = r.x - x.x;
    r.y = (x.x - (r.x - v)) + (y.x - v);
    r.y += x.y + y.y;

    return r;
}

static inline Sleef_double2 __attribute__((always_inline)) ddadd2_d2_d2_d(Sleef_double2 x, double y) {
    Sleef_double2 r;

    r.x = x.x + y;
    double v = r.x - x.x;
    r.y = (x.x - (r.x - v)) + (y - v);
    r.y += x.y;

    return r;
}

static inline Sleef_double2 __attribute__((always_inline)) ddadd2_d2_d_d2(double x, Sleef_double2 y) {
    Sleef_double2 r;

    r.x = x + y.x;
    double v = r.x - x;
    r.y = (x - (r.x - v)) + (y.x - v) + y.y;

    return r;
}


static inline Sleef_double2 __attribute__((always_inline)) ddsub_d2_d2_d2(Sleef_double2 x, Sleef_double2 y) {
    // |x| >= |y|

    Sleef_double2 r;

    r.x = x.x - y.x;
    r.y = x.x - r.x - y.x + x.y - y.y;

    return r;
}

static inline Sleef_double2 __attribute__((always_inline)) ddsqu_d2_d2(Sleef_double2 x) {
    double xh = upper(x.x), xl = x.x - xh;
    Sleef_double2 r;

    r.x = x.x * x.x;
    r.y = xh * xh - r.x + (xh + xh) * xl + xl * xl + x.x * (x.y + x.y);

    return r;
}

static inline double __attribute__((always_inline)) pow2i(int q) {
    return longBitsToDouble(((int64_t) (q + 0x3ff)) << 52);
}

static inline double __attribute__((always_inline)) ldexp2k(double d, int e) { // faster than ldexpk, short reach
    return d * pow2i(e >> 1) * pow2i(e - (e >> 1));
}

#define L2U .69314718055966295651160180568695068359375
#define L2L .28235290563031577122588448175013436025525412068e-12
#define R_LN2 1.442695040888963407359924681001892137426645954152985934135449406931
#define SLEEF_NAN __builtin_nan("")

static inline Sleef_double2 __attribute__((always_inline)) expk2(Sleef_double2 d) {
    int q = (int) rintk((d.x + d.y) * R_LN2);
    Sleef_double2 s, t;
    double u;

    s = ddadd2_d2_d2_d(d, q * -L2U);
    s = ddadd2_d2_d2_d(s, q * -L2L);

    u = +0.1602472219709932072e-9;
    u = mla(u, s.x, +0.2092255183563157007e-8);
    u = mla(u, s.x, +0.2505230023782644465e-7);
    u = mla(u, s.x, +0.2755724800902135303e-6);
    u = mla(u, s.x, +0.2755731892386044373e-5);
    u = mla(u, s.x, +0.2480158735605815065e-4);
    u = mla(u, s.x, +0.1984126984148071858e-3);
    u = mla(u, s.x, +0.1388888888886763255e-2);
    u = mla(u, s.x, +0.8333333333333347095e-2);
    u = mla(u, s.x, +0.4166666666666669905e-1);

    t = ddadd2_d2_d2_d(ddmul_d2_d2_d(s, u), +0.1666666666666666574e+0);
    t = ddadd2_d2_d2_d(ddmul_d2_d2_d2(s, t), 0.5);
    t = ddadd2_d2_d2_d2(s, ddmul_d2_d2_d2(ddsqu_d2_d2(s), t));

    t = ddadd2_d2_d_d2(1, t);

    t.x = ldexp2k(t.x, q);
    t.y = ldexp2k(t.y, q);

    return d.x < -1000 ? dd(0, 0) : t;
}

static inline double __attribute__((always_inline)) xerfc_u15(double a) {
    double s = a, r = 0, t;
    Sleef_double2 u, d, x;
    a = fabsk(a);
    int o0 = a < 1.0, o1 = a < 2.2, o2 = a < 4.2, o3 = a < 27.3;
    u = o0 ? ddmul_d2_d_d(a, a) : o1 ? dd(a, 0) : dddiv_d2_d2_d2(dd(1, 0), dd(a, 0));

    t = o0 ? +0.6801072401395386139e-20 : o1 ? +0.3438010341362585303e-12 : o2 ? -0.5757819536420710449e+2
                                                                               : +0.2334249729638701319e+5;
    t = mla(t, u.x, o0 ? -0.2161766247570055669e-18 : o1 ? -0.1237021188160598264e-10 : o2 ? +0.4669289654498104483e+3
                                                                                           : -0.4695661044933107769e+5);
    t = mla(t, u.x, o0 ? +0.4695919173301595670e-17 : o1 ? +0.2117985839877627852e-09 : o2 ? -0.1796329879461355858e+4
                                                                                           : +0.3173403108748643353e+5);
    t = mla(t, u.x, o0 ? -0.9049140419888007122e-16 : o1 ? -0.2290560929177369506e-08 : o2 ? +0.4355892193699575728e+4
                                                                                           : +0.3242982786959573787e+4);
    t = mla(t, u.x, o0 ? +0.1634018903557410728e-14 : o1 ? +0.1748931621698149538e-07 : o2 ? -0.7456258884965764992e+4
                                                                                           : -0.2014717999760347811e+5);
    t = mla(t, u.x, o0 ? -0.2783485786333451745e-13 : o1 ? -0.9956602606623249195e-07 : o2 ? +0.9553977358167021521e+4
                                                                                           : +0.1554006970967118286e+5);
    t = mla(t, u.x, o0 ? +0.4463221276786415752e-12 : o1 ? +0.4330010240640327080e-06 : o2 ? -0.9470019905444229153e+4
                                                                                           : -0.6150874190563554293e+4);
    t = mla(t, u.x, o0 ? -0.6711366622850136563e-11 : o1 ? -0.1435050600991763331e-05 : o2 ? +0.7387344321849855078e+4
                                                                                           : +0.1240047765634815732e+4);
    t = mla(t, u.x, o0 ? +0.9422759050232662223e-10 : o1 ? +0.3460139479650695662e-05 : o2 ? -0.4557713054166382790e+4
                                                                                           : -0.8210325475752699731e+2);
    t = mla(t, u.x, o0 ? -0.1229055530100229098e-08 : o1 ? -0.4988908180632898173e-05 : o2 ? +0.2207866967354055305e+4
                                                                                           : +0.3242443880839930870e+2);
    t = mla(t, u.x, o0 ? +0.1480719281585086512e-07 : o1 ? -0.1308775976326352012e-05 : o2 ? -0.8217975658621754746e+3
                                                                                           : -0.2923418863833160586e+2);
    t = mla(t, u.x, o0 ? -0.1636584469123399803e-06 : o1 ? +0.2825086540850310103e-04 : o2 ? +0.2268659483507917400e+3
                                                                                           : +0.3457461732814383071e+0);
    t = mla(t, u.x, o0 ? +0.1646211436588923575e-05 : o1 ? -0.6393913713069986071e-04 : o2 ? -0.4633361260318560682e+2
                                                                                           : +0.5489730155952392998e+1);
    t = mla(t, u.x, o0 ? -0.1492565035840623511e-04 : o1 ? -0.2566436514695078926e-04 : o2 ? +0.9557380123733945965e+1
                                                                                           : +0.1559934132251294134e-2);
    t = mla(t, u.x, o0 ? +0.1205533298178967851e-03 : o1 ? +0.5895792375659440364e-03 : o2 ? -0.2958429331939661289e+1
                                                                                           : -0.1541741566831520638e+1);
    t = mla(t, u.x, o0 ? -0.8548327023450850081e-03 : o1 ? -0.1695715579163588598e-02 : o2 ? +0.1670329508092765480e+0
                                                                                           : +0.2823152230558364186e-5);
    t = mla(t, u.x, o0 ? +0.5223977625442187932e-02 : o1 ? +0.2089116434918055149e-03 : o2 ? +0.6096615680115419211e+0
                                                                                           : +0.6249999184195342838e+0);
    t = mla(t, u.x, o0 ? -0.2686617064513125222e-01 : o1 ? +0.1912855949584917753e-01 : o2 ? +0.1059212443193543585e-2
                                                                                           : +0.1741749416408701288e-8);

    d = ddmul_d2_d2_d(u, t);
    d = ddadd2_d2_d2_d2(d, o0 ? dd(0.11283791670955126141, -4.0175691625932118483e-18) :
                           o1 ? dd(-0.10277263343147646779, -6.2338714083404900225e-18) :
                           o2 ? dd(-0.50005180473999022439, 2.6362140569041995803e-17) :
                           dd(-0.5000000000258444377, -4.0074044712386992281e-17));
    d = ddmul_d2_d2_d2(d, u);
    d = ddadd2_d2_d2_d2(d, o0 ? dd(-0.37612638903183753802, 1.3391897206042552387e-17) :
                           o1 ? dd(-0.63661976742916359662, 7.6321019159085724662e-18) :
                           o2 ? dd(1.601106273924963368e-06, 1.1974001857764476775e-23) :
                           dd(2.3761973137523364792e-13, -1.1670076950531026582e-29));
    d = ddmul_d2_d2_d2(d, u);
    d = ddadd2_d2_d2_d2(d, o0 ? dd(1.1283791670955125586, 1.5335459613165822674e-17) :
                           o1 ? dd(-1.1283791674717296161, 8.0896847755965377194e-17) :
                           o2 ? dd(-0.57236496645145429341, 3.0704553245872027258e-17) :
                           dd(-0.57236494292470108114, -2.3984352208056898003e-17));

    x = ddmul_d2_d2_d(o1 ? d : dd(-a, 0), a);
    x = o1 ? x : ddadd2_d2_d2_d2(x, d);
    x = o0 ? ddsub_d2_d2_d2(dd(1, 0), x) : expk2(x);
    x = o1 ? x : ddmul_d2_d2_d2(x, u);

    r = o3 ? (x.x + x.y) : 0;
    if (s < 0) r = 2 - r;
    r = xisnan(s) ? SLEEF_NAN : r;
    return r;
}


#endif //PRICER_ERFC_H
