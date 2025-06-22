using System;
using System.IO.Compression;
using System.Numerics;
using System.Runtime.InteropServices;

namespace TransformApp
{
    class Program
    {
        static void Main(string[] args)
        {
            // Parameter initialization
            Vector3 u1 = new Vector3(1, 0, 0);
            double theta1 = 30.0;
            Vector3 t1 = new Vector3(0.2, 0, 0);
            Vector3 u2 = new Vector3(0, 0, 1);
            double theate2 = 45.0;
            Vector3 t2 = new Vector3(0, 0, 10);

            // H1 and H2

            double[,] H1 = CreateHomegeneousMatrix(u1, theta1, t1);
            double[,] H2 = CreateHomegeneousMatrix(u2, theate2, t2);

            // Hfinal

            double[,] Hfinal = MultiplyMatrices(H2, H1);

            //Transform the point
            double[] P_scs = new double[] { 1, 2, 3, 1 };
            double[] P_wcs = MultiplyMatrixByVector(Hfinal, P_scs);

            PrintMatrix("H1 (SCS → PCS)", H1);
            PrintMatrix("H2 (PCS → WCS)", H2);
            PrintMatrix("H_total", Hfinal);
            PrintVector("P_wcs", P_wcs);

        }
        static double[,] CreateHomegeneousMatrix(Vector3 axis, double angleDeg, Vector3 translation)
        {
            Quaternion q = Quaternion.FromAxisAngle(axis, angleDeg);
            double[,] R = q.ToRotationMatrix();
            double[,] H = new double[4, 4];

            //Rotation part
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    H[i, j] = R[i, j];
                }
            }
            //Translation part
            H[0, 3] = translation.X;
            H[1, 3] = translation.Y;
            H[2, 3] = translation.Z;

            H[3, 0] = H[3, 1] = H[3, 2] = 0;
            H[3, 3] = 1;

            return H;
        }

        // Matrix multiplication
        static double[,] MultiplyMatrices(double[,] A, double[,] B)
        {
            var C = new double[4, 4];
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    for (int k = 0; k < 4; k++)
                    {
                        C[i, j] += A[i, k] * B[k, j];

                    }
                }
            }
            return C;
        }
        static double[] MultiplyMatrixByVector(double[,] matrix, double[] vector)
        {
            var r = new double[4];
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    r[i] += matrix[i, j] * vector[j];
                }
            }
            return r;

        }
                static void PrintMatrix(string name, double[,] M)
        {
            Console.WriteLine(name + ":");
            for (int i = 0; i < 4; i++)
            {
                Console.Write("  | ");
                for (int j = 0; j < 4; j++)
                    Console.Write($"{M[i,j],8:F4} ");
                Console.WriteLine("|");
            }
            Console.WriteLine();
        }
        static void PrintVector(string name, double[] v)
        {
            Console.Write($"{name}: [");
            for (int i = 0; i < v.Length; i++)
                Console.Write($"{v[i]:F6}" + (i<v.Length-1?", ":""));
            Console.WriteLine("]");
        }
    }

        struct Vector3
        {
            public double X, Y, Z;

            public Vector3(double x, double y, double z) { X = x; Y = y; Z = z; }
        }

        struct Quaternion
        {
            public double X, Y, Z, W;
            public Quaternion(double x, double y, double z, double w) { X = x; Y = y; Z = z; W = w; }

            public static Quaternion FromAxisAngle(Vector3 axis, double angleDeg)
            {
                double rad = Math.PI * angleDeg / 180;
                double s = Math.Sin(rad / 2);
                double c = Math.Cos(rad / 2);
                double norm = Math.Sqrt(axis.X*axis.X+axis.Y*axis.Y+axis.Z*axis.Z);
                return new Quaternion(axis.X/norm * s, axis.Y/norm * s, axis.Z/norm * s, c);
            }
            public double[,] ToRotationMatrix()
            {
                double[,] R = new double[3, 3];
                R[0, 0] = 1 - 2 * (Y * Y + Z * Z);
                R[0, 1] = 2 * (X * Y - W * Z);
                R[0, 2] = 2 * (X * Z + W * Y);
                R[1, 0] = 2 * (X * Y + W * Z);
                R[1, 1] = 1 - 2 * (X * X + Z * Z);
                R[1, 2] = 2 * (Y * Z - W * X);
                R[2, 0] = 2 * (X * Z - W * Y);
                R[2, 1] = 2 * (Y * Z + W * X);
                R[2, 2] = 1 - 2 * (X * X + Y * Y);

                return R;
            } 
            
        }
    }
