using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Threading;

using System.Globalization;
using System.IO.Ports;
using System.Runtime.InteropServices;
using System.Text.RegularExpressions;


namespace Denoising
{
    class Program
    {
        public static double t0 = 0;
        public static List<string> ResultSt = new List<string>();

        static void Main(string[] args)
        {
            //double [] date1 = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
            List<double> data1 = new List<double>();
            data1 = Load("data.txt");
            int M = data1.Count;
            double J = Math.Log(2, M);
            //double sigm = 8.5E-7; //сигма для 3D
            //double sigm = 9.05E-08; //сигма для 3D
            double sigm = 5.05E-06; //сигма для 3D
            t0 = Math.Sqrt(2 * Math.Log(M) / M) * 2 * sigm;  // часть отсечки. ее ещу нужно добножить на Math.Pow(2, (J-j)/2)) на каждом шаге 
            // j - шаг итерации. С большим шагом t будет все меньше и меньше. 

            List<double> data2 = new List<double>();
            List<double> datapconvAp1 = new List<double>();
            List<double> datapconvDec1 = new List<double>();
            List<double> datapconvAp2 = new List<double>();
            List<double> datapconvDec2 = new List<double>();
            List<double> datapconvAp3 = new List<double>();
            List<double> datapconvDec3 = new List<double>();
            List<double> datapconvAp4 = new List<double>();
            List<double> datapconvDec4 = new List<double>();
            List<double> datapconvAp5 = new List<double>();
            List<double> datapconvDec5 = new List<double>();
            List<double> datapconvAp6 = new List<double>();
            List<double> datapconvDec6 = new List<double>();
            List<double> datapconvAp7 = new List<double>();
            List<double> datapconvDec7 = new List<double>();
            List<double> datapconv3 = new List<double>();

            double[] CL = new double[] {    (1 + Math.Sqrt(3)) / (4 * Math.Sqrt(2)),
                                            (3 + Math.Sqrt(3)) / (4 * Math.Sqrt(2)),
                                            (3 - Math.Sqrt(3)) / (4 * Math.Sqrt(2)),
                                            (1 - Math.Sqrt(3)) / (4 * Math.Sqrt(2)) }; // D4

            /*double[] CL = new double[] {  (1 + Math.Sqrt(3)) / (4 * Math.Sqrt(2)),
                                            (3 + Math.Sqrt(3)) / (4 * Math.Sqrt(2)),
                                            (3 - Math.Sqrt(3)) / (4 * Math.Sqrt(2)),
                                            (1 - Math.Sqrt(3)) / (4 * Math.Sqrt(2)) }; // D4
            */
            //double[] CL = new double[] { 1, 1 }; // Хаара
            //double[] CL = new double[] { 0.0, 1.0, 2.0, 3.0 }; 


            int N = CL.Length;
            double[] CH = new double[N];
            for (int i = 0; i < N; i++)
            {
                CH[i] = Math.Pow(-1, i) * CL[N - i - 1];
            }

            printAr(CL, "CL");
            printAr(CH, "CH");
            var result = icoeffs(CL, CH);
            printAr(result.Item1, "icoeffsv 1");
            printAr(result.Item2, "icoeffsv 2");

            var resultWave1 = pconv(data1, CL, CH, 0, J, 0);
            datapconvAp1 = resultWave1.Item1;
            datapconvDec1 = resultWave1.Item2;

            var resultWave2 = pconv(datapconvAp1, CL, CH, 0, J, 1);
            datapconvAp2 = resultWave2.Item1;
            datapconvDec2 = resultWave2.Item2;
            var resultWave3 = pconv(datapconvAp2, CL, CH, 0, J, 2);
            datapconvAp3 = resultWave3.Item1;
            datapconvDec3 = resultWave3.Item2;
            var resultWave4 = pconv(datapconvAp3, CL, CH, 0, J, 3);
            datapconvAp4 = resultWave4.Item1;
            datapconvDec4 = resultWave4.Item2;
            var resultWave5 = pconv(datapconvAp4, CL, CH, 0, J, 4);
            datapconvAp5 = resultWave5.Item1;
            datapconvDec5 = resultWave5.Item2;
            var resultWave6 = pconv(datapconvAp5, CL, CH, 0, J, 5);
            datapconvAp6 = resultWave6.Item1;
            datapconvDec6 = resultWave6.Item2;
            var resultWave7 = pconv(datapconvAp6, CL, CH, 0, J, 6);
            datapconvAp7 = resultWave7.Item1;
            datapconvDec7 = resultWave7.Item2;


            AddColSt(datapconvAp1, "datapconvAp1");
            AddColSt(datapconvDec1, "datapconvDec1");
            AddColSt(datapconvAp2, "datapconvAp2 ");
            AddColSt(datapconvDec2, "datapconvDec2 ");
            AddColSt(datapconvAp3, "datapconvAp3 ");
            AddColSt(datapconvDec3, "datapconvDec3 ");
            AddColSt(datapconvAp4, "datapconvAp4 ");
            AddColSt(datapconvDec4, "datapconvDec4 ");
            AddColSt(datapconvAp5, "datapconvAp5 ");
            AddColSt(datapconvDec5, "datapconvDec5 ");
            AddColSt(datapconvAp6, "datapconvAp6 ");
            AddColSt(datapconvDec7, "datapconvDec7 ");

            var resultWaveback7 = pconvBack(ListAdd(datapconvAp7, datapconvDec7), result.Item1, result.Item2, 2);
            AddColSt(ListAdd(resultWaveback7.Item1, resultWaveback7.Item2), "back7");
            var resultWaveback6 = pconvBack(ListAdd(datapconvAp6, datapconvDec6), result.Item1, result.Item2, 2);
            AddColSt(ListAdd(resultWaveback6.Item1, resultWaveback6.Item2), "back6");
            var resultWaveback5 = pconvBack(ListAdd(datapconvAp5, datapconvDec5), result.Item1, result.Item2, 2);
            AddColSt(ListAdd(resultWaveback5.Item1, resultWaveback5.Item2), "back5");

            var resultWaveback4 = pconvBack(ListAdd(datapconvAp4, datapconvDec4), result.Item1, result.Item2, 2);
            AddColSt(ListAdd(resultWaveback4.Item1, resultWaveback4.Item2), "back4");
            var resultWaveback3 = pconvBack(ListAdd(datapconvAp3, datapconvDec3), result.Item1, result.Item2, 2);
            AddColSt(ListAdd(resultWaveback3.Item1, resultWaveback3.Item2), "back3");

            var resultWaveback2 = pconvBack(ListAdd(datapconvAp2, datapconvDec2), result.Item1, result.Item2, 2);
            AddColSt(ListAdd(resultWaveback2.Item1, resultWaveback2.Item2), "back2");
            var resultWaveback1 = pconvBack(ListAdd(ListAdd(resultWaveback2.Item1, resultWaveback2.Item2), datapconvDec1), result.Item1, result.Item2, 2);
            AddColSt(ListAdd(resultWaveback1.Item1, resultWaveback1.Item2), "back1");


            SaveListSt(ResultSt);
            //SaveWavelet(datapconv);
            //Save(datapconv2);
            //printAr(date1, "date1");
            //printAr(datepconv, "datepconv");



            //data2 = pconv(datapconv, result.Item1, result.Item2, 2);
            //Save(data2);
            //printAr(date2, "date2");
            Console.ReadKey();

        }


        static void printAr(double[] arr, string st)
        {
            Console.WriteLine("\n" + st);
            for (int i = 0; i < arr.Length; i++)
            {
                Console.Write(arr[i] + "\t");
            }
            Console.WriteLine();

        }


        static List<double> ListAdd(List<double> ap, List<double> dec)
        {
            List<double> outList = new List<double>();

            for (int k = 0; k < ap.Count; k++) // 
            {
                outList.Add(ap[k]); // Низкочастотный коэффициент
                outList.Add(dec[k]);          // Высокочастотный коэффициент

            }

            return outList;
        }


        static Tuple<List<double>, List<double>> pconv(List<double> data, double[] CL, double[] CH, int delta, double J, int step)
        {
            int N = CL.Length;
            int M = data.Count;
            double sL = 0, sH = 0;
            //double t=t0*Math.Pow(2, (J-step)/2);
            double t = t0 * Math.Pow(2, (J) / 2);
            List<double> outsL = new List<double>();
            List<double> outhL = new List<double>();

            for (int k = 0; k < M; k += 2) // Перебираем числа 0, 2, 4…
            {
                sL = 0;                        // Низкочастотный коэффициент
                sH = 0;                       // Высокочастотный коэффициент
                for (int i = 0; i < N; i++) // Находим сами взвешенные суммы
                {
                    sL += data[(M + k + i - delta) % M] * CL[i];
                    sH += data[(M + k + i - delta) % M] * CH[i];
                }


                if (Math.Abs(sH) < t)
                    sH = 0;
                else
                    sH = Math.Sign(sH) * (Math.Abs(sH) - t);

                outsL.Add(sL);                 // Добавляем коэффициенты в список
                outhL.Add(sH);
            }



            var tuple = new Tuple<List<double>, List<double>>(outsL, outhL);
            return tuple;
        }

        static Tuple<List<double>, List<double>> pconvBack(List<double> data, double[] CL, double[] CH, int delta)
        {
            int N = CL.Length;
            int M = data.Count;
            double sL = 0, sH = 0;

            List<double> outsL = new List<double>();
            List<double> outhL = new List<double>();

            for (int k = 0; k < M; k += 2) // Перебираем числа 0, 2, 4…
            {
                sL = 0;                        // Низкочастотный коэффициент
                sH = 0;                       // Высокочастотный коэффициент
                for (int i = 0; i < N; i++) // Находим сами взвешенные суммы
                {
                    sL += data[(M + k + i - delta) % M] * CL[i];
                    sH += data[(M + k + i - delta) % M] * CH[i];
                }

                outsL.Add(sL);                 // Добавляем коэффициенты в список
                outhL.Add(sH);
            }



            var tuple = new Tuple<List<double>, List<double>>(outsL, outhL);
            return tuple;
        }


        static Tuple<double[], double[]> icoeffs(double[] CL, double[] CH)
        {
            int N = CL.Length;
            double[] iCL = new double[N];  // Коэффициенты первой строки
            double[] iCH = new double[N];  // Коэффициенты второй строки

            for (int k = 0; k < N; k += 2) // Перебираем числа 0, 2, 4…
            {
                //iCL[k] = CL[k-2];
                //iCL[k+1] = CH[k-2];
                //iCH[k] = CL[k-1];
                //iCH[k+1] = CH[k-1];   
                iCL[k] = CL[(k + N - 2) % N];
                iCL[k + 1] = CH[(k + N - 2) % N];
                iCH[k] = CL[(k + N - 1) % N];
                iCH[k + 1] = CH[(k + N - 1) % N];
            }
            var tuple = new Tuple<double[], double[]>(iCL, iCH);
            return tuple;
        }

        static public void Save(List<double> list, string st)
        {
            //string File_name = DateTime1.ToString() + ".txt";
            using (StreamWriter sw = new StreamWriter(st + ".txt", true))
            {


                foreach (double s in list)
                {
                    sw.WriteLine(s.ToString());

                }
            }
        }
        static public void SaveWavelet(List<double> list)
        {
            //string File_name = DateTime1.ToString() + ".txt";
            using (StreamWriter sw = new StreamWriter("temp.txt", true))
            {


                for (int i = 0; i < (list.Count - 1); i += 2)

                    sw.WriteLine(list[i].ToString() + "\t" + list[i + 1].ToString());


            }
        }


        static public void SaveListSt(List<string> list)
        {
            //string File_name = DateTime1.ToString() + ".txt";
            using (StreamWriter sw = new StreamWriter("ResultSt.txt", true))
            {


                foreach (string line in list)
                {
                    sw.WriteLine(line);

                }

            }
        }

        static public List<double> Load(string file)
        {
            List<double> dataFile = new List<double>();
            ResultSt.Add("data" + "\t");
            //string File_name = DateTime1.ToString() + ".txt";
            foreach (string line in File.ReadLines(file))
            {
                dataFile.Add(double.Parse(line));
                ResultSt.Add(line + "\t");
            }
            ResultSt.Add("\t");
            return dataFile;

        }

        static public void AddColSt(List<double> list, string col)
        {
            //string File_name = DateTime1.ToString() + ".txt";
            ResultSt[0] += col + "\t";
            int i = 1;
            foreach (double line in list)
            {
                ResultSt[i] += line.ToString() + "\t";
                i++;

            }
            for (int j = i; j < ResultSt.Count; j++)
                ResultSt[j] += "\t";


        }


    }
}
