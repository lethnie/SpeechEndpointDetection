using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.Diagnostics;

namespace SpeechEndpointDetection.Algorithms
{
    public class Parameters
    {
        /// <summary>
        /// длительность аудио в секундах
        /// </summary>
        private double duration = 0;
        /// <summary>
        /// частота дискретизации
        /// </summary>
        private int sampleRate = 0;
        /// <summary>
        /// полный сигнал, не разбитый на кадры
        /// </summary>
        private double[] fullSignal;
        /// <summary>
        /// разбитый на кадры сигнал с наложенным окном
        /// </summary>
        private double[][] wSignal;
        /// <summary>
        /// разбитый на кадры сигнал
        /// </summary>
        private double[][] signal;

        /// <summary>
        /// количество отсчётов в миллисекунде
        /// </summary>
        private int samplesInMs;
        /// <summary>
        /// количество отсчётов в кадре
        /// </summary>
        private int sizeOfFrame;
        /// <summary>
        /// количество отсчётов в сдвиге
        /// </summary>
        private int sizeOfShift;
        /// <summary>
        /// длина кадра в миллисекундах
        /// </summary>
        private int frSize = 20;
        /// <summary>
        /// длина сдвига в миллисекундах
        /// </summary>
        private int shSize = 10;
        /// <summary>
        /// количество кадров в сигнале
        /// </summary>
        private int numOfFrames;
        private int M = 3;
        private int R = 30;

        private double[] sfm;
        private double[] f;
        private double[][] spectrum;

        /// <summary>
        /// По умолчанию: размер кадра 20 мс, размер смещения - 10 мс, окно - прямоугольное, M = 3, R = 30
        /// </summary>
        /// <param name="sign">Исходный сигнал</param>
        /// <param name="dur">Длительность аудио (в секундах)</param>
        /// <param name="sRate">Частота дискретизации</param>
        public Parameters(double[] sign, double dur, int sRate)
        {
            fullSignal = sign;
            duration = dur;
            sampleRate = sRate;
            samplesInMs = (int)(fullSignal.Length / (duration * 1000));
            sizeOfFrame = frSize * samplesInMs;
            sizeOfShift = shSize * samplesInMs;
            numOfFrames = (int)Math.Ceiling((double)sign.Length / sizeOfShift) - 1;
            wSignal = GetSignal(-1);
            signal = GetSignal(-1);
            FFT();
        }

        /// <summary>
        /// По умолчанию: размер кадра 20 мс, размер смещения - 10 мс, M = 3, R = 30
        /// </summary>
        /// <param name="sign">Исходный сигнал</param>
        /// <param name="dur">Длительность аудио (в секундах)</param>
        /// <param name="sRate">Частота дискретизации</param>
        /// <param name="wType">Вид используемого окна: 0 - окно Ханна, 1 - Хэмминга, -1 - без окна</param>
        public Parameters(double[] sign, double dur, int sRate, int wType)
        {
            fullSignal = sign;
            duration = dur;
            sampleRate = sRate;
            samplesInMs = (int)(fullSignal.Length / (duration * 1000));
            sizeOfFrame = frSize * samplesInMs;
            sizeOfShift = shSize * samplesInMs;
            numOfFrames = (int)Math.Ceiling((double)sign.Length / sizeOfShift) - 1;
            wSignal = GetSignal(wType);
            signal = GetSignal(-1);
            FFT();
        }

        /// <summary>
        /// По умолчанию: M = 3, R = 30
        /// </summary>
        /// <param name="sign">Исходный сигнал</param>
        /// <param name="dur">Длительность аудио (в секундах)</param>
        /// <param name="sRate">Частота дискретизации</param>
        /// <param name="wType">Вид используемого окна: 0 - окно Ханна, 1 - Хэмминга, -1 - без окна</param>
        /// <param name="fSize">Размер кадра (в миллисекундах)</param>
        /// <param name="sSize">Размер смещения (в миллисекундах)</param>
        public Parameters(double[] sign, double dur, int sRate, int wType, int fSize, int sSize)
        {
            fullSignal = sign;
            duration = dur;
            sampleRate = sRate;
            frSize = fSize;
            shSize = sSize;
            samplesInMs = (int)(fullSignal.Length / (duration * 1000));
            sizeOfFrame = frSize * samplesInMs;
            sizeOfShift = shSize * samplesInMs;
            numOfFrames = (int)Math.Ceiling((double)sign.Length / sizeOfShift) - 1;
            wSignal = GetSignal(wType);
            signal = GetSignal(-1);
            FFT();
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="sign">Исходный сигнал</param>
        /// <param name="dur">Длительность аудио (в секундах)</param>
        /// <param name="sRate">Частота дискретизации</param>
        /// <param name="wType">Вид используемого окна: 0 - окно Ханна, 1 - Хэмминга, -1 - без окна</param>
        /// <param name="fSize">Размер кадра (в миллисекундах)</param>
        /// <param name="sSize">Размер смещения (в миллисекундах)</param>
        /// <param name="m">Размер окна для преобразования Уэлча-Бартлетта</param>
        /// <param name="r">Размер окна для долгосрочных характеристик</param>
        public Parameters(double[] sign, double dur, int sRate, int wType, int fSize, int sSize, int m, int r)
        {
            fullSignal = sign;
            duration = dur;
            sampleRate = sRate;
            frSize = fSize;
            shSize = sSize;
            M = m;
            R = r;
            samplesInMs = (int)(fullSignal.Length / (duration * 1000));
            sizeOfFrame = frSize * samplesInMs;
            sizeOfShift = shSize * samplesInMs;
            numOfFrames = (int)Math.Ceiling((double)sign.Length / sizeOfShift) - 1;
            wSignal = GetSignal(wType);
            signal = GetSignal(-1);
            FFT();
        }

        /// <summary>
        /// Получение количества отсчётов в одном кадре
        /// </summary>
        /// <returns></returns>
        public int GetFrameSizeInSamples()
        {
            return sizeOfFrame;
        }

        /// <summary>
        /// Получение количества отсчётов в одном сдвиге
        /// </summary>
        /// <returns></returns>
        public int GetShiftSizeInSamples()
        {
            return sizeOfShift;
        }

        /// <summary>
        /// Получение количества миллисекунд в одном кадре
        /// </summary>
        /// <returns></returns>
        public int GetFrameSizeInMs()
        {
            return frSize;
        }

        /// <summary>
        /// Получение количества миллисекунд в одном сдвиге
        /// </summary>
        /// <returns></returns>
        public int GetShiftSizeInMs()
        {
            return shSize;
        }

        /// <summary>
        /// Получение разделённого по кадрам сигнала
        /// </summary>
        /// <param name="wtype">Вид используемого окна: 0 - окно Ханна, 1 - Хэмминга, -1 - без окна</param>
        /// <returns></returns>
        private double[][] GetSignal(int wtype)
        {         
            double[][] signal = new double[numOfFrames][];
            double[] w = GetWindow(sizeOfFrame, wtype);
            for (int i = 0; i < numOfFrames; i++)
            {
                signal[i] = new double[sizeOfFrame];
                for (int j = 0; (j < sizeOfFrame); j++)
                {
                    if (i * sizeOfShift + j >= fullSignal.Length)
                        signal[i][j] = 0;
                    else
                        signal[i][j] = fullSignal[i * sizeOfShift + j] * w[j];
                }
            }
            return signal;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="N"></param>
        /// <param name="type">0 - Hanning, 1 - Hamming</param>
        /// <returns></returns>
        public static double[] GetWindow(int N, int type)
        {
            double[] result = new double[N];
            switch (type)
            {
                case 0:
                    {
                        for (int n = 0; n < N; n++)
                        {
                            result[n] = 0.5 * (1.0 - Math.Cos(2 * Math.PI * n / (N - 1)));
                        }
                        break;
                    }
                case 1:
                    {
                        for (int n = 0; n < N; n++)
                        {
                            result[n] = 0.53836 - 0.46164 * Math.Cos(2 * Math.PI * n / (N - 1));
                        }
                        break;
                    }
                default:
                    {
                        for (int n = 0; n < N; n++)
                        {
                            result[n] = 1;
                        }
                        break;
                    }
            }
            return result;
        }

        /// <summary>
        /// Получение исходного сигнала, не разделённого на кадры
        /// </summary>
        /// <returns></returns>
        public double[] GetSignal()
        {
            double max = 0;
            for (int i = 0; i < fullSignal.Length; i++)
            {
                if (Math.Abs(fullSignal[i]) > max)
                    max = Math.Abs(fullSignal[i]);
            }
            for (int i = 0; i < fullSignal.Length; i++)
            {
                fullSignal[i] = fullSignal[i] / max;
            }
            return fullSignal;
        }
      
        /// <summary>
        /// Получение энергии сигнала для каждого кадра
        /// </summary>
        /// <returns></returns>
        public double[] GetEnergy()
        {
            double[] energy = new double[numOfFrames];
            for (int i = 0; i < numOfFrames; i++)
            {
                double sum = 0;
                for (int j = 0; j < signal[i].Length; j++)
                {
                    sum += signal[i][j] * signal[i][j];
                }
                energy[i] = sum;
            }
            return energy;
        }

        /// <summary>
        /// Получение среднего числа пересечений нуля сигналом для каждого кадра
        /// </summary>
        /// <returns></returns>
        public double[] GetZeroCrossingRate()
        {
            double[] zcr = new double[numOfFrames];
            for (int i = 0; i < numOfFrames; i++)
            {
                double sum = 0;
                for (int j = 1; j < signal[i].Length; j++)
                {
                    sum += Math.Abs(Math.Sign(signal[i][j]) - Math.Sign(signal[i][j - 1]));
                }
                zcr[i] = sum / (2 * sizeOfFrame);
            }
            return zcr;
        }

        /// <summary>
        /// Выполнение быстрого преобразования Фурье для сигнала, получение меры спектральной плоскостности и доминирующих частот
        /// </summary>
        private void FFT()
        {
            sfm = new double[numOfFrames];
            f = new double[numOfFrames];
            spectrum = new double[numOfFrames][];        

            for (long i = 0; i < numOfFrames; i++)
            {
                spectrum[i] = MathOperations.FastFourierTransformMagn(wSignal[i]);
                sfm[i] = MathOperations.GetSFM(spectrum[i]);
                long maxL = 0;
                for (int j = 1; j < spectrum[i].Length / 2; j++)
                {
                    if ((spectrum[i][j] > spectrum[i][maxL]))
                        maxL = j;
                }
                f[i] = (maxL / ((double)frSize / 1000.0)) * GetFrequencyCorrection(spectrum[i].Length);
            }
        }

        /// <summary>
        /// Получение коэффициента для перевода порядкового номера частоты в спектре в реальную частоту
        /// </summary>
        /// <param name="spectrumLength">Длина спектра</param>
        /// <returns></returns>
        private double GetFrequencyCorrection(int spectrumLength)
        {
            return (double)sizeOfFrame / (double)spectrumLength;
        }

        /// <summary>
        /// Получение спектра сигнала
        /// </summary>
        /// <returns></returns>
        public double[][] GetSpectrum()
        {
            return spectrum;
        }

        /// <summary>
        /// Получение меры спектральной плоскостности для каждого кадра
        /// </summary>
        /// <returns></returns>
        public double[] GetSpectralFlatnessMeasure()
        {
            return sfm;
        }

        /// <summary>
        /// Преобразование Бартлетта-Уэлча
        /// </summary>
        /// <param name="signal"></param>
        /// <returns></returns>
        private double[][] WelchBartlett(double[][] signal)
        {
            double[][] s = new double[numOfFrames][];
            for (int n = 0; n < numOfFrames; n++)
            {
                s[n] = new double[spectrum[n].Length];
                for (int k = 0; k < s[n].Length; k++)
                {
                    double sum = 0;
                    int p = n;
                    if (n < M)
                        p = 0;
                    else p = n - M;
                    for (; p < n; p++)
                    {
                        sum += (spectrum[p][k] * spectrum[p][k]);
                    }
                    sum = sum / M;
                    s[n][k] = sum;
                }
            }
            return s;
        }

        /// <summary>
        /// Получение значений для характеристики LSFM
        /// </summary>
        /// <param name="m"></param>
        /// <param name="K"></param>
        /// <param name="s"></param>
        /// <returns></returns>
        private double LFunction(int m, int K, double[][] s)
        {
            double L = 0.0;
            int mm = 0;

            if (m <= R)
                mm = 0;
            else mm = m - R;

            for (int k0 = 0; k0 < K; k0++)
            {
                double means = MathOperations.GeometricMean(s, k0, m, m - mm) / MathOperations.ArithmeticMean(s, k0, m, m - mm);
                if (means == 0)
                    L = 0;
                else
                    L = L + Math.Log10(means);
            }
            return L;
        }

        /// <summary>
        /// Получение долгосрочной меры спектральной плоскостности для каждого кадра
        /// </summary>
        /// <returns></returns>
        public double[] GetLongTermSpectralFlatnessMeasure()
        {
            double[][] s = WelchBartlett(wSignal);
            double[] l = new double[numOfFrames];

            for (int m = 0; m < numOfFrames; m++)
            {
                l[m] = LFunction(m, s[m].Length, s);
            }
            return l;
        }

        /// <summary>
        /// Получение значений для характеристики LTSV
        /// </summary>
        /// <param name="m"></param>
        /// <param name="k"></param>
        /// <param name="s"></param>
        /// <returns></returns>
        private double KsiFunction(int m, int k, double[][] s)
        {
            double ksi = 0;
            int mm = 0;
            if (m < R)
                mm = 0;
            else mm = m - R;
            for (int m2 = mm; m2 < m; m2++)
            {
                double denom = 0;
                for (int m3 = mm; m3 < m; m3++)
                {
                    denom += s[m3][k];
                }
                if (m < R)
                    denom *= R / (m + 1);
                if (s[m2][k] == 0)
                    ksi += 0;
                else
                    ksi += (s[m2][k] / denom) * Math.Log10(s[m2][k] / denom);
            }
            if (m < R)
                ksi *= R / (m + 1);
            return -ksi;
        }

        /// <summary>
        /// Получение долгосрочной меры изменчивости сигнала
        /// </summary>
        /// <returns></returns>
        public double[] GetLongTermSignalVariability()
        {
            double[][] s = WelchBartlett(wSignal);
            double[] l = new double[numOfFrames];

            for (int m = 0; m < numOfFrames; m++)
            {
                double sum = 0;
                double[] ksi = new double[s[m].Length];
                for (int k = 0; k < s[m].Length; k++)
                {
                    ksi[k] = KsiFunction(m, k, s);
                }
                double mKsi = ksi.Average();
                for (int k = 0; k < s[m].Length; k++)
                {
                    sum += (ksi[k] - mKsi) * (ksi[k] - mKsi);
                }
                sum = sum / s[m].Length;
                l[m] = sum;
            }
            return l;
        }
        
        /// <summary>
        /// Получение характеристики Energy Entropy
        /// </summary>
        /// <returns></returns>
        public double[] GetEnergyEntropy()
        {
            double[] result = new double[numOfFrames];

            for (int i = 0; i < numOfFrames; i++)
            {
                double energy = 0;
                for (int j = 0; j < sizeOfFrame; j++)
                {
                    energy += wSignal[i][j] * wSignal[i][j];
                }

                double entropy = 0;
                int K = spectrum[i].Length / 2 + 1;
                double[] p = new double[K];
                double denom = 0;
                for (int k = 0; k < K; k++)
                {
                    denom += spectrum[i][k] * spectrum[i][k];
                }
                for (int k = 0; k < K; k++)
                {
                    p[k] = spectrum[i][k] * spectrum[i][k] / denom;
                    entropy -= p[k] * Math.Log10(p[k]);
                }

                result[i] = Math.Sqrt(1 + Math.Abs(energy * entropy));
            }
            return result;
        }

        /// <summary>
        /// Энергия Тиджера для каждого кадра
        /// </summary>
        /// <returns></returns>
        public double[] GetTeagerEnergy()
        {
            double[] result = new double[numOfFrames];

            for (int i = 0; i < numOfFrames; i++)
            {
                result[i] = 0;
                double df = sampleRate / spectrum[i].Length;
                int K = spectrum[i].Length / 2 + 1;
                int k = (int)(300 / df);
                for (; k * df < 4000; k++)
                {
                    result[i] += ((k) * df) * ((k) * df) * Math.Abs(spectrum[i][k]) * Math.Abs(spectrum[i][k]);
                }
                result[i] = Math.Sqrt(result[i]);
            }
            return result;
        }

        /// <summary>
        /// Получение характеристики Mean Delta
        /// </summary>
        /// <returns></returns>
        public double[] GetMeanDelta()
        {
            double[] result = new double[numOfFrames];

            for (int i = 0; i < numOfFrames; i++)
            {
                int K = spectrum[i].Length / 2;
                int L = K / 2;
                double[] R = new double[L];
                for (int l = 0; l < L; l++)
                {
                    R[l] = 0;
                    for (int k = 0; k < K - l; k++)
                    {
                        R[l] += spectrum[i][k] * spectrum[i][k] * spectrum[i][k + l] * spectrum[i][k + l];
                    }
                }
                double[] deltaR = new double[L];
                int Q = 3;
                int sumQ = 0;
                for (int q = -Q; q <= Q; q++)
                {
                    sumQ += q * q;
                }
                for (int l = 0; l < Q; l++)
                {
                    deltaR[l] = 0;
                    for (int q = 0; q <= Q; q++)
                    {
                        deltaR[l] += q * R[l + q] / sumQ;
                    }
                    deltaR[l] *= 2;
                }
                for (int l = L - Q; l < L; l++)
                {
                    deltaR[l] = 0;
                    for (int q = -Q; q <= 0; q++)
                    {
                        deltaR[l] += q * R[l + q] / sumQ;
                    }
                    deltaR[l] *= 2;
                }
                for (int l = Q; l < L - Q; l++)
                {
                    deltaR[l] = 0;
                    for (int q = -Q; q <= Q; q++)
                    {
                        deltaR[l] += q * R[l + q] / sumQ;
                    }
                }

                double[] deltaRS = new double[L];
                int J = 3;
                for (int l = 0; l < L; l++)
                {
                    double max = deltaR[l];
                    int left = -J;
                    int right = J;
                    if (l - J < 0)
                        left = -l;
                    if (l + J >= L)
                        right = (L - l - 1);
                    for (int z = left; z <= right; z++)
                    {
                        if (deltaR[l + z] > max)
                            max = deltaR[l + z];
                    }
                    deltaRS[l] = max;
                }
                double sum = 0;
                for (int l = 0; l < L; l++)
                {
                    sum += Math.Abs(deltaRS[l]);
                }
                result[i] = Math.Sqrt(sum);
            }
            return result;
        }

        /// <summary>
        /// Доминирующая частота для каждого кадра
        /// </summary>
        /// <returns></returns>
        public double[] GetFrequency()
        {
            return f;
        }
    }
}
