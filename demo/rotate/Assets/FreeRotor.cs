using System.Collections;
using System.Collections.Generic;
using UnityEngine;


/// <summary>
/// Implementation of the exact rigid body rotation dynamic from the paper:
///   "Numerical implementation of the exact dynamics of free rigid bodies"
///   from Ramses van Zon, Jeremy Schofield
/// https://doi.org/10.1016/j.jcp.2006.11.019
/// 
/// by Sebastian Jandl, 09.01.2022
/// </summary>

public class FreeRotor : MonoBehaviour
{
    //settings
    public Vector3 I;           //body frame axis Inertias
    public Vector3 w_initial;   //body frame initial angular velocities

    [Range(0.0f, 5.0f)]
    public float t_mult;

    //debug oberservers
    public Vector3 omegas;      //developed ang-velo. during precession cycle
    public float E_kin;         //kinetic rotation energy from "omegas"
    public float E_kin_star;    //kinetic rotation energy from prev-now/rotation
    public float E_kin_comb;    //kinetic rotation energy mean
    public Vector3 L_before;    //angular momentum from intial conditions
    public Vector3 L_past;      //angular momentum from rotation development
    Quaternion u_t;      
    Quaternion a;
    public float omega_p;   
    public float cn, sn, dn, ph;
    [Range(0.0f, Mathf.PI*2)]
    public float psi;

    //variables set in initialization
    Quaternion initial_rotation;
    Quaternion q_b;
    Vector3 I_s;
    public bool orderflag;
    public Vector3 wm_s;
    [Range(0.0f, 20.0f)]
    public float epsilon;
    float m , A_1, A_2;
    int nt;
    float[] c_r, c_i;
    public float K;
    public float u;
    float L_2, L;
    Vector3 L_vec_0;

    Quaternion q_T10_y;
    Quaternion q_T10_x;
    Quaternion q_T10_;
    Quaternion q_A_0;
    Matrix4x4 U, V;

    float old_t, old_psi;
    public float delta_psi;
    public float angle;
    public Vector3 axis;

    //function to convert a rotation Matrix to a Quaternion, !probably not correct with unitys left hand axis system(?)!
    public static Quaternion QuaternionFromMatrix(Matrix4x4 m)
    {
        // Adapted from: http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm
        Quaternion q = new Quaternion();
        q.w  = Mathf.Sqrt(Mathf.Max(0, 1 + m[0, 0] + m[1, 1] + m[2, 2])) / 2;
        q.x  = Mathf.Sqrt(Mathf.Max(0, 1 + m[0, 0] - m[1, 1] - m[2, 2])) / 2;
        q.y  = Mathf.Sqrt(Mathf.Max(0, 1 - m[0, 0] + m[1, 1] - m[2, 2])) / 2;
        q.z  = Mathf.Sqrt(Mathf.Max(0, 1 - m[0, 0] - m[1, 1] + m[2, 2])) / 2;
        q.x *= Mathf.Sign(q.x * (m[2, 1] - m[1, 2]));
        q.y *= Mathf.Sign(q.y * (m[0, 2] - m[2, 0]));
        q.z *= Mathf.Sign(q.z * (m[1, 0] - m[0, 1]));
        return q;
    }

    public static void Swap<T>(ref T lhs, ref T rhs)
    {
        T temp;
        temp = lhs;
        lhs = rhs;
        rhs = temp;
    }

    public static void swap_mat(ref Matrix4x4 mat, int a, int b, int c)
    {
        Vector4 temp = mat.GetRow(a);
        mat.SetRow(a, mat.GetRow(b));
        mat.SetRow(b, temp);
        mat.SetRow(c, -mat.GetRow(c));
    }

    //function to convert the provided x=sin(a) argument to the a=asin(x) argument supported by the lib
    public static float ellipk_sin_arg(float x, float m)
    {
        //valid range [-1,1]
        if (Mathf.Abs(x) > 1.0f)
        {
            Debug.Log("math error");
            return x;
        }

        bool is_neg = x < 0.0f;

        //somehow the ellik is wrong when using negitiv values, use the sin/asin mirror symmetry here
        if(is_neg)
            return -ellipk.ellik(Mathf.Asin(-x), m);
        else 
            return ellipk.ellik(Mathf.Asin(x), m);
    }

    // Start is called before the first frame update
    void Start()
    {
        //create types
        Vector3 w0_s    = new Vector3();    //vector of initial angular velocities (reordered if necessary)
        omegas          = new Vector3();    
        I_s             = new Vector3();
        wm_s            = new Vector3();

        float K_d;
        float q;
        float etha;
        float Xi; 
        float dA2;
        int n,k;
        float r_0;
        float i_0; 

        //TODO: sort and arrange Intertias/omegas so it is usable in all configurations
        //I.y must be middle component for now
        if (!(I.y > Mathf.Min(I.x,I.z) && I.y < Mathf.Max(I.x,I.z)))
        { Debug.Log("I.y must be middle component!"); }

        //get initial rotation 
        q_A_0 = Quaternion.Inverse(transform.rotation);
        old_t = 0.0f;

        //==============================BEGIN===============================================
        // initialize routine according to paper
        //==================================================================================
        L_2 = Vector3.Scale(I, w_initial).sqrMagnitude; // I.x * I.x * w_initial.x * w_initial.x + I.y * I.y * w_initial.y * w_initial.y + I.z * I.z * w_initial.z * w_initial.z;
        float two_E = I.x * w_initial.x * w_initial.x + I.y * w_initial.y * w_initial.y + I.z * w_initial.z * w_initial.z;
        L = Mathf.Sqrt(L_2);
        float L_n;

        V = Matrix4x4.identity;
        U = Matrix4x4.identity;

        //initial ordering to bring L1 < L2 < L3

        //bring middle to middel
        int[] idx = new int[3];
        float[] I_temp = new float[3];

        I_temp[0] = I.x;
        I_temp[1] = I.y;
        I_temp[2] = I.z;
        idx[0] = 0;
        idx[1] = 1;
        idx[2] = 2;

        if (I_temp[0] > I_temp[1])
        { 
            Swap(ref I_temp[0], ref I_temp[1]); Swap(ref idx[0], ref idx[1]);
            swap_mat(ref V, 0, 1, 2);
        }

        if (I_temp[1] > I_temp[2])
        { 
            Swap(ref I_temp[1], ref I_temp[2]); Swap(ref idx[1], ref idx[2]);
            swap_mat(ref V, 1, 2,0);
        }

        if (I_temp[0] > I_temp[1])
        { 
            Swap(ref I_temp[0], ref I_temp[1]); Swap(ref idx[0], ref idx[1]);
            swap_mat(ref V, 0, 1, 2);
        }

        //second ording to care about  energy
        //orderflag = (two_E > L_2 / I.y && I.x < I.z) || (two_E < L_2 / I.y && I.x > I.z);
        orderflag = (two_E > L_2 / I_temp[1] && I_temp[0] < I_temp[2]) || (two_E < L_2 / I_temp[1] && I_temp[0] > I_temp[2]);
        if (orderflag)
        {
            U.m00 = 0;
            U.m02 = 1;
            U.m11 = -1;
            U.m20 = 1;
            U.m22 = 0;
        }

        I_s = U * V * I;
        I_s.x = Mathf.Abs(I_s.x); //correct Inertia always positiv
        I_s.y = Mathf.Abs(I_s.y); //correct Inertia always positiv
        I_s.z = Mathf.Abs(I_s.z); //correct Inertia always positiv


        w0_s = U * V * w_initial;

        Vector3 L_vec = Vector3.Scale(I_s, w0_s);
        L_vec_0 = L_vec;
        L_n = Mathf.Sqrt(L_vec.x * L_vec.x + L_vec.y * L_vec.y);

        //rotate to bring y comp to 0
        float angle_z = Mathf.Atan2((-L_vec.y / L_n) , (L_vec.x / L_n));
        q_T10_y = Quaternion.Euler(0,0, angle_z * Mathf.Rad2Deg);

        //rotate to bring x comp to 0
        float angle_y = Mathf.Atan2((-L_n / L) , (L_vec.z / L));
        q_T10_x = Quaternion.Euler(0, angle_y * Mathf.Rad2Deg, 0);


        L_before = Vector3.Scale(I_s, w0_s);

        //B is the reusable rotation Matrix to combine the: initial rotation (A_0),
        //  switch of components for jacobi ordering (U_star), and z-axis alignment (T)

        wm_s.x  =  Mathf.Sign(w0_s.x) * Mathf.Sqrt((L_2 - two_E * I_s.z) / (I_s.x * (I_s.x - I_s.z)));
        wm_s.y  = -Mathf.Sign(w0_s.x) * Mathf.Sqrt((L_2 - two_E * I_s.z) / (I_s.y * (I_s.y - I_s.z)));
        wm_s.z  =  Mathf.Sign(w0_s.z) * Mathf.Sqrt((L_2 - two_E * I_s.x) / (I_s.z * (I_s.z - I_s.x)));
        omega_p =  Mathf.Sign(I_s.y - I_s.z) * Mathf.Sign(w0_s.z) * Mathf.Sqrt((L_2 - two_E * I_s.x) * (I_s.z - I_s.y) / (I_s.x * I_s.y * I_s.z));

        //inputs for the elliptic functions are calculated
        m       = (L_2 - two_E * I_s.z) * (I_s.x - I_s.y) / ((L_2 - two_E * I_s.x) * (I_s.z - I_s.y));
        epsilon = -4.65792f;//
        epsilon = ellipk_sin_arg(w0_s.y / wm_s.y, m);
        K       = ellipk_sin_arg(1.0f , m);
        K_d     = ellipk_sin_arg(1.0f, (1.0f - m));
        q       = Mathf.Exp(-Mathf.PI * K_d / K);
        etha    = Mathf.Sign(w0_s.z) * K_d - ellipk_sin_arg(I_s.z * wm_s.z / L, (1.0f - m));
        Xi      = Mathf.Exp(Mathf.PI * etha / K);
        A_2     = L / I_s.x + Mathf.PI * omega_p * (Xi + 1.0f) / (2 * K * (Xi - 1.0f));

        n = 1;
        do
        {
            float q2n = Mathf.Pow(q, 2 * n);
            dA2 = -(Mathf.PI * omega_p / K) * (q2n / (1.0f - q2n)) * (Mathf.Pow(Xi, n) - Mathf.Pow(Xi, -n));
            A_2 += dA2;
            n++;
        } while (Mathf.Abs(dA2) > Mathf.Epsilon);
        nt = (int)(Mathf.Log(Mathf.Epsilon) / Mathf.Log(q));

        r_0 = 0;
        i_0 = 0;
        c_r = new float[nt + 1];
        c_i = new float[nt + 1];
        for (n = 0; n <= nt; n++)
        {
            float temp = 2 * Mathf.Pow(q, n * (n + 1) + 1.0f / 4.0f);
            c_r[n] = Mathf.Pow(-1.0f, n)    * temp * ellipk.cosh((2 * n + 1.0f) * Mathf.PI * etha / (2.0f * K));
            c_i[n] = Mathf.Pow(-1.0f, n+1)  * temp * ellipk.sinh((2 * n + 1.0f) * Mathf.PI * etha / (2.0f * K));
            r_0 += c_r[n] * Mathf.Sin((2 * n + 1.0f) * Mathf.PI * epsilon / (2.0f * K));
            i_0 += c_i[n] * Mathf.Cos((2 * n + 1.0f) * Mathf.PI * epsilon / (2.0f * K));

            //added this to exit before "nan" values emerge
            if (c_r[n] == 0.0f && c_i[n] == 0.0f)
                break;
        }

        k = r_0 > 0 ? 0 : (int)Mathf.Sign(i_0);
        A_1 = Mathf.Atan(i_0 / r_0) + k * Mathf.PI;
        //===============================END================================================
    }


    // Update is called once per frame
    void Update()
    {
        //==============================BEGIN===============================================
        // evolve routine according to paper
        //==================================================================================

        //note that this uses an absolute time - no step integrations required
        float t = Time.time * t_mult; 

        cn = 0;
        sn = 0;
        dn = 0;
        ph = 0;

        u = omega_p * t + epsilon;
        ellipk.ellpj(u, m, ref sn, ref cn, ref dn, ref ph);

        omegas.x = wm_s.x * cn;
        omegas.y = wm_s.y * sn;
        omegas.z = wm_s.z * dn;

        float real_V1 = 0;
        float img_V1 = 0;

        for(int n = 0; n <= nt; n++)
        {
            real_V1 += c_r[n] * Mathf.Sin((2 * n + 1) * Mathf.PI * (omega_p * t + epsilon) / (2 * K));
            img_V1  += c_i[n] * Mathf.Cos((2 * n + 1) * Mathf.PI * (omega_p * t + epsilon) / (2 * K));
        }
        
        float C = Mathf.Cos(A_1 + A_2 * t);
        float S = Mathf.Sin(A_1 + A_2 * t);

        float cos_psi = (C * real_V1 + S * img_V1) / Mathf.Sqrt(real_V1 * real_V1 + img_V1 * img_V1);
        float sin_psi = (S * real_V1 - C * img_V1) / Mathf.Sqrt(real_V1 * real_V1 + img_V1 * img_V1);

        psi = Mathf.Atan2(sin_psi,cos_psi);
        delta_psi = psi < old_psi ? psi + Mathf.PI*2 - old_psi : psi - old_psi;
        delta_psi /= (t - old_t);
        old_psi = psi;

        Vector3 L_vec = Vector3.Scale(I_s, omegas);
        float L_n = Mathf.Sqrt(L_vec.x * L_vec.x + L_vec.y * L_vec.y);

        //rotate to bring y comp to 0
        float angle_z = Mathf.Atan2((-L_vec.y / L_n) , (L_vec.x / L_n));
        Quaternion q_T11_y = Quaternion.Euler(0, 0, angle_z * Mathf.Rad2Deg);

        //rotate to bring x comp to 0
        float angle_y = Mathf.Atan2((-L_n / L) , (L_vec.z / L));
        Quaternion q_T11_x = Quaternion.Euler(0, angle_y * Mathf.Rad2Deg, 0);

        omegas = V.transpose * U.transpose * omegas;
        Quaternion rot_q_ = Quaternion.Euler(0,0, psi * Mathf.Rad2Deg); //usually the rotation angle should be set negative, but I think due to unity -z axis it needs to use it as is.

        //total
        Quaternion q_all = Quaternion.Inverse(QuaternionFromMatrix(V)) * Quaternion.Inverse(QuaternionFromMatrix(U)) * Quaternion.Inverse(q_T11_y) * Quaternion.Inverse(q_T11_x) * rot_q_ * q_T10_x * q_T10_y * QuaternionFromMatrix(U) * QuaternionFromMatrix(V) * q_A_0;

        //I think q_all is how the space around should turn, to rotate the object itself we have to invert the rotation
        Quaternion new_rotation = Quaternion.Inverse(q_all);

        //skalar moment of inertia
        Vector3 some_rot = transform.rotation * Vector3.forward;
        Vector3 some_rot2 = new_rotation * Vector3.forward;

        (new_rotation * Quaternion.Inverse(transform.rotation)).ToAngleAxis(out angle, out axis);
        angle = Quaternion.Angle(transform.rotation, new_rotation);
        float rad_vel = angle * Mathf.Deg2Rad / (t - old_t);
        axis = axis.normalized;
            
        float I = Vector3.Dot(axis, Vector3.Scale(I_s, axis));
        E_kin_star = 0.5f * I * rad_vel * rad_vel;

        I = Vector3.Dot(omegas.normalized, Vector3.Scale(I_s, omegas.normalized)); // skalar I = n * I * n^T, we have I as the diagonal vector so n.dot(scale(I,n)) is fine
        E_kin = 0.5f * I * omegas.sqrMagnitude;

        transform.rotation = new_rotation;
        old_t = t;
    }
}
