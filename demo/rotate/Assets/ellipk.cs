using UnityEngine;


//taken from ellpk.c of the ellf.tgz on http://www.netlib.org/cephes/
//translated 09.01.2022 by Sebastian Jandl (to c# and float types)

public class ellipk
{
    public static float PIO2 = 1.57079632679489661923f;       /* pi/2 */
    public static float PI = 3.14159265358979323846f;       /* pi */
    public static float MAXLOG = 7.08396418532264106224E2f;     /* log 2**1022 */
    public static float LOG2E = 1.4426950408889634073599f;     /* 1/log(2) */
    public static float LOGE2 = 6.93147180559945309417E-1f;    /* log(2) */
    public static float MINLOG = -7.08396418532264106224E2f;     /* log 2**-1022 */
    public static float C1 = 1.3862943611198906188E0f; /* log(4) */

    public static float[] P =
    {
     1.37982864606273237150E-4f,
     2.28025724005875567385E-3f,
     7.97404013220415179367E-3f,
     9.85821379021226008714E-3f,
     6.87489687449949877925E-3f,
     6.18901033637687613229E-3f,
     8.79078273952743772254E-3f,
     1.49380448916805252718E-2f,
     3.08851465246711995998E-2f,
     9.65735902811690126535E-2f,
     1.38629436111989062502E0f
    };

    public static float[] Q =
    {
     2.94078955048598507511E-5f,
     9.14184723865917226571E-4f,
     5.94058303753167793257E-3f,
     1.54850516649762399335E-2f,
     2.39089602715924892727E-2f,
     3.01204715227604046988E-2f,
     3.73774314173823228969E-2f,
     4.88280347570998239232E-2f,
     7.03124996963957469739E-2f,
     1.24999999999870820058E-1f,
     4.99999999999999999821E-1f
    };
    public static bool isNan(float x)
    {
        return x == Mathf.Infinity || x == Mathf.NegativeInfinity;
    }

    public static float ellik(float phi, float m)
    {
        float a, b, c, e, temp, t, K;
        int d, mod, sign, npio2;

        if (m == 0.0)
            return (phi);
        a = 1.0f - m;
        if (a == 0.0)
        {
            if (Mathf.Abs(phi) >= PIO2)
            {
                Debug.Log("math error");
                return (Mathf.Infinity);
            }
            return Mathf.Log((Mathf.Tan((PIO2 + phi) / 2.0f)));
        }
        npio2 = (int)Mathf.Floor(phi / PIO2);
        if (npio2 % 2 > 0)
            npio2 += 1;
        if (npio2 > 0)
        {
            K = ellpk(a);
            phi = phi - npio2 * PIO2;
        }
        else
            K = 0.0f;
        if (phi < 0.0)
        {
            phi = -phi;
            sign = -1;
        }
        else
            sign = 0;
        b = Mathf.Sqrt(a);
        t = Mathf.Tan(phi);
        if (Mathf.Abs(t) > 10.0)
        {
            /* Transform the amplitude */
            e = 1.0f / (b * t);
            /* ... but avoid multiple recursions.  */
            if (Mathf.Abs(e) < 10.0)
            {
                e = Mathf.Atan(e);
                if (npio2 == 0)
                    K = ellpk(a);
                temp = K - ellik(e, m);
                goto done;
            }
        }
        a = 1.0f;
        c = Mathf.Sqrt(m);
        d = 1;
        mod = 0;

        while (Mathf.Abs(c / a) > 1e-6f /*Mathf.Epsilon*/)
        {
            temp = b / a;
            phi = phi + Mathf.Atan(t * temp) + mod * PI;
            mod = (int)((phi + PIO2) / Mathf.PI);
            t = t * (1.0f + temp) / (1.0f - temp * t * t);
            c = (a - b) / 2.0f;
            temp = Mathf.Sqrt(a * b);
            a = (a + b) / 2.0f;
            b = temp;
            d += d;
        }

        temp = (Mathf.Atan(t) + mod * PI) / (d * a);

    done:
        if (sign < 0)
            temp = -temp;
        temp += npio2 * K;
        return (temp);
    }

    public static float ellpk(float x)
    {

        if ((x < 0.0) || (x > 1.0))
        {
            Debug.Log("math error");
            return (0.0f);
        }

        if (x > Mathf.Epsilon)
        {
            return (polevl(x, P, 10) - Mathf.Log(x) * polevl(x, Q, 10));
        }
        else
        {
            if (x == 0.0)
            {
                Debug.Log("math error");
                return (Mathf.Infinity);
            }
            else
            {
                return (C1 - 0.5f * Mathf.Log(x));
            }
        }
    }

    public static float polevl(float x, float[] coef, int N)
    {
        float ans;
        int i;
        ans = coef[0];
        i = N;

        do
            ans = ans * x + coef[N - i + 1];
        while (--i > 0);

        return (ans);
    }

    public static float p1evl(float x, float[] coef, int N)
    {
        float ans;
        int i;

        ans = x + coef[0];
        i = N - 1;

        do
            ans = ans * x + coef[N - i + 1];
        while (--i > 0);

        return (ans);
    }


    public static float cosh(float x)
    {
        float y;

        if (isNan(x))
            return (x);
        if (x < 0)

            x = -x;
        if (x > (MAXLOG + LOGE2))
        {
            Debug.Log("math error");
            return (Mathf.Infinity);
        }
        if (x >= (MAXLOG - LOGE2))
        {
            y = Mathf.Exp(0.5f * x);
            y = (0.5f * y) * y;
            return (y);
        }
        y = Mathf.Exp(x);
        y = 0.5f * (y + 1.0f / y);
        return (y);
    }


    public static float sinh(float x)
    {
        float a;

        if (x == 0.0)
            return (x);

        a = Mathf.Abs(x);
        if ((x > (MAXLOG + LOGE2)) || (x > -(MINLOG - LOGE2)))
        {
            Debug.Log("math error");
            if (x > 0)
                return (Mathf.Infinity);
            else
                return (Mathf.NegativeInfinity);
        }
        if (a > 1.0)
        {
            if (a >= (MAXLOG - LOGE2))
            {
                a = Mathf.Exp(0.5f * a);
                a = (0.5f * a) * a;
                if (x < 0)
                    a = -a;
                return (a);
            }
            a = Mathf.Exp(a);
            a = 0.5f * a - (0.5f / a);
            if (x < 0)
                a = -a;
            return (a);
        }

        a *= a;
        return (x + x * a * (polevl(a, P, 3) / p1evl(a, Q, 3)));
    }


    public static float tanh(float x)
    {
        float s, z;

        if (x == 0.0)
            return (x);

        z = Mathf.Abs(x);
        if (z > 0.5 * MAXLOG)
        {
            if (x > 0)
                return (1.0f);
            else
                return (-1.0f);
        }
        if (z >= 0.625)
        {
            s = Mathf.Exp(2.0f * z);
            z = 1.0f - 2.0f / (s + 1.0f);
            if (x < 0)
                z = -z;
        }
        else
        {
            if (x == 0.0)
                return (x);
            s = x * x;
            z = polevl(s, P, 2) / p1evl(s, Q, 3);
            z = x * s * z;
            z = x + z;
        }
        return (z);
    }

    public static int ellpj(float u, float m, ref float sn, ref float cn, ref float dn, ref float ph)
    {
        float ai, b, phi, t, twon;
        float[] a = new float[9];
        float[] c = new float[9];
        int i;


        /* Check for special cases */

        if (m < 0.0 || m > 1.0)
        {
            Debug.Log("math error");
            sn = 0.0f;
            cn = 0.0f;
            ph = 0.0f;
            dn = 0.0f;
            return -1;
        }
        if (m < 1.0e-9)
        {
            t = Mathf.Sin(u);
            b = Mathf.Cos(u);
            ai = 0.25f * m * (u - t * b);
            sn = t - ai * b;
            cn = b + ai * t;
            ph = u - ai;
            dn = 1.0f - 0.5f * m * t * t;
            return 0;
        }

        if (m >= 0.9999999999f)
        {
            ai = 0.25f * (1.0f - m);
            b = cosh(u);
            t = tanh(u);
            phi = 1.0f / b;
            twon = b * sinh(u);
            sn = t + ai * (twon - u) / (b * b);
            ph = 2.0f * Mathf.Atan(Mathf.Exp(u)) - PIO2 + ai * (twon - u) / b;
            ai *= t * phi;
            cn = phi - ai * (twon - u);
            dn = phi + ai * (twon + u);
            return (0);
        }


        /*	A. G. M. scale		*/
        a[0] = 1.0f;
        b = Mathf.Sqrt(1.0f - m);
        c[0] = Mathf.Sqrt(m);
        twon = 1.0f;
        i = 0;

        while (Mathf.Abs(c[i] / a[i]) > 1e-20f /*Mathf.Epsilon*/)
        {
            if (i > 7)
            {
                Debug.Log("math error");
                break;
            }
            ai = a[i];
            ++i;
            c[i] = (ai - b) / 2.0f;
            t = Mathf.Sqrt(ai * b);
            a[i] = (ai + b) / 2.0f;
            b = t;
            twon *= 2.0f;
        }

        /* backward recurrence */
        phi = twon * a[i] * u;
        do
        {
            t = c[i] * Mathf.Sin(phi) / a[i];
            b = phi;
            phi = (Mathf.Asin(t) + phi) / 2.0f;
        }
        while (--i > 0);

        t = Mathf.Sin(phi);
        sn = t;
        cn = Mathf.Cos(phi);
        /* Thanks to Hartmut Henkel for reporting a bug here:  */
        dn = Mathf.Sqrt(1.0f - m * t * t);
        ph = phi;
        return (0);
    }
}
