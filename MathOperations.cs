using System;
using System.Numerics;
using System.Collections.Generic;
using System.Linq;
using System.Diagnostics;

namespace SpeechEndpointDetection
{
    public class MathOperations
    {
        /// <summary>
        /// Вычисление поворачивающего модуля e^(-i*2*PI*k/N)
        /// </summary>
        /// <param name="k"></param>
        /// <param name="N"></param>
        /// <returns></returns>
        private static Complex w(int k, int N)
        {
            double arg = -2 * Math.PI * k / N;
            return new Complex(Math.Cos(arg), Math.Sin(arg));
        }

        /// <summary>
        /// Возвращает спектр сигнала
        /// </summary>
        /// <param name="x">Массив значений сигнала. Количество значений должно быть степенью 2</param>
        /// <returns>Массив со значениями спектра сигнала</returns>
        public static Complex[] FastFourierTransform(double[] x)
        {
            int length = GetNewLength(x.Length);
            Complex[] newX = new Complex[length];
            for (int i = 0; i < x.Length; i++)
            {
                newX[i] = new Complex(x[i], 0);
            }
            for (int i = x.Length; i < length; i++)
            {
                newX[i] = new Complex(0, 0);
            }
            return FastFourierTransform(newX);
        }

        /// <summary>
        /// Возвращает спектр сигнала
        /// </summary>
        /// <param name="x">Массив значений сигнала. Количество значений должно быть степенью 2</param>
        /// <returns>Массив со значениями магнитуд спектра сигнала</returns>
        public static double[] FastFourierTransformMagn(double[] x)
        {
            Complex[] fft = FastFourierTransform(x);
            double[] result = new double[fft.Length];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = fft[i].Magnitude;
            }
            return result;
        }

        /// <summary>
        /// Расчёт длины массива со значениями спектра
        /// </summary>
        /// <param name="length">Длина массива с отсчётами сигнала</param>
        /// <returns></returns>
        private static int GetNewLength(int length)
        {
            if (!IsPowerOfTwo((uint)length))
            {
                int newLength = 2;
                while (newLength < length)
                    newLength *= 2;
                return newLength;
            }
            else
            {
                return length;
            }
        }

        /// <summary>
        /// Проверка, является ли число степенью двойки
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        private static bool IsPowerOfTwo (uint x)
        {
            return ((x != 0) && ((x & (~x + 1)) == x));
        }
        /// <summary>
        /// Возвращает спектр сигнала
        /// </summary>
        /// <param name="x">Массив значений сигнала. Количество значений должно быть степенью 2</param>
        /// <returns>Массив со значениями спектра сигнала</returns>
        public static Complex[] FastFourierTransform(Complex[] x)
        {
            Complex[] X;
            int N = x.Length;
            if (N == 2)
            {
                X = new Complex[2];
                X[0] = x[0] + x[1];
                X[1] = x[0] - x[1];
            }
            else
            {
                Complex[] x_even = new Complex[N / 2];
                Complex[] x_odd = new Complex[N / 2];
                for (int i = 0; i < N / 2; i++)
                {
                    x_even[i] = x[2 * i];
                    x_odd[i] = x[2 * i + 1];
                }
                Complex[] X_even = FastFourierTransform(x_even);
                Complex[] X_odd = FastFourierTransform(x_odd);
                X = new Complex[N];
                for (int i = 0; i < N / 2; i++)
                {
                    X[i] = X_even[i] + w(i, N) * X_odd[i];
                    X[i + N / 2] = X_even[i] - w(i, N) * X_odd[i];
                }
            }
            return X;
        }

        /// <summary>
        /// Считает среднее геометрическое амплитуд последовательности комплексных чисел
        /// </summary>
        /// <param name="s">Последовательность комплексных чисел</param>
        /// <returns>Среднее геометрическое последовательности комплексных чисел</returns>
        public static double GeometricMean(double[] s)
        {
            double result = 1.0;
            int state = 0;

            double av = s.Average();
            for (int i = 0; i < s.Length; i++)
            {
                result = result * s[i] / av;
                if (double.IsInfinity(result))
                {
                    state = 1;
                    break;
                }
                if (result == 0)
                {
                    state = 2;
                    break;
                }
            }
            if (state == 1)
            {
                result = 1.0;
                for (int i = 0; i < s.Length; i++)
                {
                    result = result * s[i] / 10.0;
                }
                result = Math.Pow(result, 1.0 / (s.Length)) * 10.0;
            }
            if (state == 2)
            {
                int numplus = 0;
                int numminus = 0;
                result = 1.0;
                for (int i = 0; i < s.Length; i++)
                {
                    if (s[i] != 0)
                        result = result * s[i] / av;
                    if (result < 0.001)
                    {
                        result *= 1000;
                        numplus++;
                    }
                    if (result > 1000)
                    {
                        result /= 1000;
                        numminus++;
                    }
                }
                result = Math.Pow(result, 1.0 / (s.Length)) * av * Math.Pow(1000.0, -(double)(numplus - numminus)/(double)(s.Length));
            }
            if (state == 0)
                result = Math.Pow(result, 1.0 / (s.Length)) * av;
            return result;
        }

        /// <summary>
        /// Считает среднее арифметическое амплитуд последовательности комплексных чисел
        /// </summary>
        /// <param name="s">Последовательность комплексных чисел</param>
        /// <returns>Среднее арифметическое последовательности комплексных чисел</returns>
        public static double ArithmeticMean(double[] s)
        {
            double result = 0;
            for (int i = 0; i < s.Length; i++)
            {
                result = result + s[i];
            }
            result = result / (s.Length);
            return result;
        }

        /// <summary>
        /// Получает значение SFM для спектра сигнала
        /// </summary>
        /// <param name="s">спектр сигнала</param>
        /// <returns>SFM</returns>
        public static double GetSFM(double[] s)
        {
            double g = MathOperations.GeometricMean(s);
            double a = MathOperations.ArithmeticMean(s);
            double log = Math.Log10(g / a);
            return 10 * log;
        }

        /// <summary>
        /// Вычисление поворачивающего модуля e^(i*2*PI*k/N)
        /// </summary>
        /// <param name="k"></param>
        /// <param name="N"></param>
        /// <returns></returns>
        private static Complex iw(int k, int N)
        {
            if (k % N == 0) return 1;
            double arg = 2 * Math.PI * k / N;
            return new Complex(Math.Cos(arg), Math.Sin(arg));
        }

        /// <summary>
        /// Нахождение среднего арифметического значений спектра
        /// </summary>
        /// <param name="s">Спектр</param>
        /// <param name="wk">Частота</param>
        /// <param name="m">Конечный отсчёт сигнала</param>
        /// <param name="R">Количество отсчётов, для которых считается среднее арифметическое</param>
        /// <returns></returns>
        public static double ArithmeticMean(double[][] s, int wk, int m, int R)
        {
            double result = 0;
            for (int i = m - R; i < m; i++)
            {
                result = result + s[i][wk] * s[i][wk];
            }
            result = result / R;
            return result;
        }

        /// <summary>
        /// Нахождение среднего геометрического значений спектра
        /// </summary>
        /// <param name="s">Спектр</param>
        /// <param name="wk">Частота</param>
        /// <param name="m">Конечный отсчёт сигнала</param>
        /// <param name="R">Количество отсчётов, для которых считается среднее геометрическое</param>
        /// <returns></returns>
        public static double GeometricMean(double[][] s, int wk, int m, int R)
        {
            double result = 1.0;
            int state = 0;
            for (int i = m - R; i < m; i++)
            {
                result = result * s[i][wk] * s[i][wk];
                if (double.IsInfinity(result))
                {
                    state = 1;
                    break;
                }
                if (result == 0)
                {
                    state = 2;
                    break;
                }
            }
            if (state == 1)
            {
                result = 1.0;
                for (int i = m - R; i < m; i++)
                {
                    result = result * s[i][wk] * s[i][wk] / 10.0;
                }
                result = Math.Pow(result, 1.0 / R) * 10.0;
            }
            if (state == 2)
            {
                result = 1.0;
                for (int i = m - R; i < m; i++)
                {
                    result = result * s[i][wk] * s[i][wk] * 10.0;
                }
                result = Math.Pow(result, 1.0 / R) / 10.0;
            }
            if (state == 0)
                result = Math.Pow(result, 1.0 / R);
            return result;
        }
    }
}
