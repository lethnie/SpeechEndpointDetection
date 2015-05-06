using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SpeechEndpointDetection.Algorithms
{
    public class Detection
    {
        private const double frequencyThreshFloor = 300;
        private const double frequencyThreshCeil = 4000;

        /// <summary>
        /// Минимальное число кадров, помеченных неречевыми
        /// </summary>
        private int minFalseFrames = 25;
        /// <summary>
        /// Минимальное число кадров, помеченных речевыми
        /// </summary>
        private int minTrueFrames = 20;
        /// <summary>
        /// Количество начальных кадров, считающихся неречевыми, используемых для определения пороговых значений характеристик
        /// </summary>
        private int silentFrames = 30;
        /// <summary>
        /// Количество миллисекунд в одном кадре
        /// </summary>
        private int frameSize = 0;
        /// <summary>
        /// Количество миллисекунд в одном сдвиге
        /// </summary>
        private int shiftSize = 0;
        /// <summary>
        /// Параметры и характеристики выбранного аудио
        /// </summary>
        private Parameters parameters;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="par">Параметры и характеристики выбранного аудио</param>
        public Detection(Parameters par)
        {
            parameters = par;
            frameSize = par.GetFrameSizeInMs();
            shiftSize = par.GetShiftSizeInMs();
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="par">Параметры и характеристики выбранного аудио</param>
        /// <param name="silentFrames">Количество начальных кадров, считающихся неречевыми, используемых для определения пороговых значений характеристик</param>
        public Detection(Parameters par, int silentFrames)
        {
            parameters = par;
            frameSize = par.GetFrameSizeInMs();
            shiftSize = par.GetShiftSizeInMs();
            this.silentFrames = silentFrames;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="par">Параметры и характеристики выбранного аудио</param>
        /// <param name="minFalseFrames">Минимальное число кадров, помеченных неречевыми</param>
        /// <param name="minTrueFrames">Минимальное число кадров, помеченных речевыми</param>
        public Detection(Parameters par, int minFalseFrames, int minTrueFrames)
        {
            parameters = par;
            frameSize = par.GetFrameSizeInMs();
            shiftSize = par.GetShiftSizeInMs();
            this.minFalseFrames = minFalseFrames;
            this.minTrueFrames = minTrueFrames;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="par">Параметры и характеристики выбранного аудио</param>
        /// <param name="silentFrames"> Количество начальных кадров, считающихся неречевыми, используемых для определения пороговых значений характеристик</param>
        /// <param name="minFalseFrames">Минимальное число кадров, помеченных неречевыми</param>
        /// <param name="minTrueFrames">Минимальное число кадров, помеченных речевыми</param>
        public Detection(Parameters par, int silentFrames, int minFalseFrames, int minTrueFrames)
        {
            parameters = par;
            frameSize = par.GetFrameSizeInMs();
            shiftSize = par.GetShiftSizeInMs();
            this.minFalseFrames = minFalseFrames;
            this.minTrueFrames = minTrueFrames;
            this.silentFrames = silentFrames;
        }

        /// <summary>
        /// Корректировка границ выделенных речевых фрагментов
        /// </summary>
        /// <param name="frames"></param>
        private void BorderCorrection(bool[] frames)
        {
            bool detected = false;
            int count = 0;
            for (int i = 0; i < frames.Length - 2; i++)
            {
                if (detected)
                {
                    if (frames[i])
                        count++;
                    else
                    {
                        if (count <= 2)
                        {
                            for (int j = i - 1; j >= i - count; j--)
                            {
                                frames[j] = false;
                            }
                        }
                        detected = false;
                    }
                }
                else
                {
                    if (frames[i])
                    {
                        detected = true;
                        count = 1;
                    }
                }
            }
            detected = false;
            for (int i = 0; i < frames.Length - minFalseFrames; i++)
            {
                if (detected)
                {
                    if (!frames[i])
                    {
                        int j = i + 1;
                        for (; j < i + minFalseFrames; j++)
                        {
                            if (frames[j])
                            {
                                for (int k = i; k < j; k++)
                                    frames[k] = true;
                                i = j;
                                break;
                            }
                        }
                        if (j == i + minFalseFrames)
                            detected = false;
                    }
                }
                else
                {
                    if (frames[i])
                    {
                        detected = true;
                    }
                }
            }
            detected = false;
            count = 0;
            for (int i = 0; i < frames.Length - minTrueFrames; i++)
            {
                if (detected)
                {
                    if (frames[i])
                        count++;
                    else
                    {
                        if (count < minTrueFrames)
                        {
                            for (int j = i - 1; j >= i - count; j--)
                            {
                                frames[j] = false;
                            }
                        }
                        detected = false;
                    }
                }
                else
                {
                    if (frames[i])
                    {
                        detected = true;
                        count = 1;
                    }
                }
            }
        }

        /// <summary>
        /// Преобразование массива, определяющего принадлежность кадров к речи, к набору границ речи
        /// </summary>
        /// <param name="speech">Соответствие кадров речи</param>
        /// <returns>Границы речи (в миллисекундах)</returns>
        private List<int[]> GetBounds(bool[] speech)
        {
            List<int[]> result = new List<int[]>();

            for (int i = 0; i < speech.Length; i++)
            {
                if (speech[i])
                {
                    int[] res = new int[2];
                    res[0] = i * shiftSize;
                    while ((i < speech.Length) && (speech[i]))
                        i++;
                    res[1] = i * shiftSize + frameSize - shiftSize;
                    result.Add(res);
                }
            }
            return result;
        }

        /// <summary>
        /// Алгоритм определения фрагментов речи на основе энергии сигнала и среднего числа пересечений нуля
        /// </summary>
        /// <returns>Границы речи (в миллисекундах)</returns>
        public List<int[]> AlgorithmEnergyZCR()
        {
            return GetBounds(AlgorithmEnergyZCR_());
        }

        private bool[] AlgorithmEnergyZCR_()
        {
            double[] energy = parameters.GetEnergy();
            double[] zcr = parameters.GetZeroCrossingRate();
            double[] noiseBuff = new double[silentFrames];
            int zcrCount = 26;
            double m = 0;
            double D = 0;
            for (int i = 0; i < silentFrames; i++)
            {
                m += zcr[i];
            }
            m = m / silentFrames;
            for (int i = 0; i < silentFrames; i++)
            {
                D += (zcr[i] - m) * (zcr[i] - m);
            }
            D = D / (silentFrames - 1);
            double zcrThresh = Math.Min(25.0 / frameSize, m + 2 * Math.Sqrt(D));
            m = 0;
            D = 0;
            for (int i = 0; i < silentFrames; i++)
            {
                noiseBuff[i] = energy[i];
                m += energy[i];
            }
            m = m / silentFrames;
            for (int i = 0; i < silentFrames; i++)
            {
                D += (energy[i] - m) * (energy[i] - m);
            }
            D = D / (silentFrames - 1);
            double max = 0;
            for (int i = 0; i < energy.Length; i++)
            {
                if ((energy[i] > max) && (zcr[i] <= zcrThresh))
                    max = energy[i];
            }
            double l1 = 0.03 * (max - m) + m;
            double l2 = 4 * m;
            double ITL = Math.Min(l1, l2);
            double ITU = 5 * ITL;

            bool[] result = new bool[energy.Length];
            for (int i = 0; i < result.Length; i++)
                result[i] = false;

            int start = -1;
            int finish = -1;
            double p = 0.1;
            double sigmaOld = Math.Sqrt(D);
            int buffInd = 0;

            for (int i = silentFrames; i < energy.Length; i++)
            {
                if (energy[i] >= ITL)
                {
                    start = i;
                    while ((start != -1) && (finish == -1))
                    {
                        if (i >= energy.Length)
                            break;
                        if ((energy[i] >= ITU) && (zcr[i] <= zcrThresh))
                            finish = i;
                        if (energy[i] < ITL)
                            start = -1;
                        i++;
                    }
                    if (finish != -1)
                    {
                        while ((i < energy.Length) && (energy[i] >= ITL))
                        {
                            i++;
                        }
                        finish = i - 1;
                        for (int j = start; j <= finish; j++)
                        {
                            result[j] = true;
                        }
                    }
                    start = -1;
                    finish = -1;
                }
                else
                {
                    noiseBuff[buffInd] = energy[i];
                    buffInd++;
                    if (buffInd >= 30)
                        buffInd = 0;
                    double avg = noiseBuff.Average();
                    double sumsq = noiseBuff.Select(val => (val - avg) * (val - avg)).Sum();
                    double sigmaNew = Math.Sqrt(sumsq / (noiseBuff.Length - 1));
                    if (sigmaNew / sigmaOld >= 1.25)
                        p = 0.25;
                    else
                    {
                        if (sigmaNew / sigmaOld >= 1.1)
                            p = 0.2;
                        else
                        {
                            if (sigmaNew / sigmaOld >= 1.0)
                                p = 0.15;
                            else
                                p = 0.1;
                        }
                    }
                    l1 = 0.03 * (max - avg) + avg;
                    l2 = 4 * avg;
                    ITL = (1 - p) * ITL + p * Math.Min(l1, l2);
                    ITU = 5 * ITL;
                    sigmaOld = sigmaNew;
                }
            }
            BorderCorrection(result);

            for (int i = 1; i < result.Length - 1; i++)
            {
                int count = 0;
                int last = -1;
                if ((result[i]) && (!result[i - 1]))
                {
                    for (int j = i - 1; (j > i - zcrCount) && (j >= 0); j--)
                    {
                        if (zcr[j] >= zcrThresh)
                        {
                            count++;
                            last = j;
                        }
                    }
                    if (count >= 3)
                    {
                        for (int j = i - 1; j >= last; j--)
                        {
                            result[j] = true;
                        }
                    }
                }
                count = 0;
                last = -1;
                if ((result[i]) && (!result[i + 1]))
                {
                    for (int j = i + 1; (j < i + zcrCount) && (j < result.Length); j++)
                    {
                        if (zcr[j] >= zcrThresh)
                        {
                            count++;
                            last = j;
                        }
                    }
                    if (count >= 3)
                    {
                        for (int j = i + 1; j <= last; j++)
                        {
                            result[j] = true;
                        }
                    }
                }
            }
            return result;
        }

        /// <summary>
        /// Алгоритм определения фрагментов речи на основе энергии сигнала, доминирующей частоты и меры спектральной плоскостности
        /// </summary>
        /// <returns>Границы речи (в миллисекундах)</returns>
        public List<int[]> AlgorithmEnergySFMFrequency()
        {
            return GetBounds(AlgorithmEnergySFMFrequency_());
        }

        private bool[] AlgorithmEnergySFMFrequency_()
        {
            double[] energy = parameters.GetEnergy();
            double[] sfm = parameters.GetSpectralFlatnessMeasure();
            double[] frequency = parameters.GetFrequency();

            bool[] sfmSpeech = new bool[sfm.Length];
            bool[] energySpeech = new bool[energy.Length];
            double m = 0;
            double D = 0;
            double[] noiseBuff = new double[silentFrames];
            for (int i = 0; i < silentFrames; i++)
            {
                m += energy[i];
                noiseBuff[i] = energy[i];
            }
            m = m / silentFrames;
            for (int i = 0; i < silentFrames; i++)
            {
                D += (energy[i] - m) * (energy[i] - m);
            }
            D = D / (silentFrames - 1);
            double max = 0;
            for (int i = 0; i < energy.Length; i++)
            {
                if (energy[i] > max)
                    max = energy[i];
            }
            double l1 = 0.03 * (max - m) + m;
            double l2 = 4 * m;
            double ITL = Math.Min(l1, l2);
            double ITU = 5 * ITL;

            int start = -1;
            int finish = -1;
            double p = 0.1;
            double sigmaOld = Math.Sqrt(D);
            int buffInd = 0;

            for (int i = silentFrames; i < energy.Length; i++)
            {
                if (energy[i] >= ITL)
                {
                    start = i;
                    while ((start != -1) && (finish == -1))
                    {
                        if (i >= energy.Length)
                            break;
                        if (energy[i] >= ITU)
                            finish = i;
                        if (energy[i] < ITL)
                            start = -1;
                        i++;
                    }
                    if (finish != -1)
                    {
                        while ((i < energy.Length) && (energy[i] >= ITL))
                        {
                            i++;
                        }
                        finish = i - 1;
                        for (int j = start; j <= finish; j++)
                        {
                            energySpeech[j] = true;
                        }
                    }
                    start = -1;
                    finish = -1;
                }
                else
                {
                    noiseBuff[buffInd] = energy[i];
                    buffInd++;
                    if (buffInd >= 30)
                        buffInd = 0;
                    double avg = noiseBuff.Average();
                    double sumsq = noiseBuff.Select(val => (val - avg) * (val - avg)).Sum();
                    double sigmaNew = Math.Sqrt(sumsq / (noiseBuff.Length - 1));
                    if (sigmaNew / sigmaOld >= 1.25)
                        p = 0.25;
                    else
                    {
                        if (sigmaNew / sigmaOld >= 1.1)
                            p = 0.2;
                        else
                        {
                            if (sigmaNew / sigmaOld >= 1.0)
                                p = 0.15;
                            else
                                p = 0.1;
                        }
                    }
                    l1 = 0.03 * (max - avg) + avg;
                    l2 = 4 * avg;
                    ITL = (1 - p) * ITL + p * Math.Min(l1, l2);
                    ITU = 5 * ITL;
                    sigmaOld = sigmaNew;
                }
            }
            //BorderCorrection(energySpeech);

            m = 0;
            D = 0;
            noiseBuff = new double[silentFrames];
            buffInd = 0;
            for (int i = 0; i < silentFrames; i++)
            {
                m += sfm[i];
                noiseBuff[i] = sfm[i];
            }
            m = m / silentFrames;
            for (int i = 0; i < silentFrames; i++)
            {
                D += (sfm[i] - m) * (sfm[i] - m);
            }
            D = D / (silentFrames - 1);
            max = 0;
            for (int i = 0; i < sfm.Length; i++)
            {
                if (sfm[i] < max)
                    max = sfm[i];
            }
            l1 = 0.03 * (max - m) + m;
            l2 = 4 * m;
            ITL = Math.Max(l1, l2);
            ITU = ITL - 3;

            start = -1;
            finish = -1;
            p = 0.1;
            sigmaOld = Math.Sqrt(D);

            for (int i = silentFrames; i < sfm.Length; i++)
            {
                if (sfm[i] <= ITL)
                {
                    start = i;
                    int count = 0;
                    while ((start != -1) && (finish == -1))
                    {
                        if (i >= sfm.Length)
                            break;
                        if (sfm[i] <= ITU)
                            finish = i;
                        if (sfm[i] > ITL)
                            count++;
                        i++;
                        if (count > 2)
                            start = -1;
                    }
                    if (finish != -1)
                    {
                        count = 0;
                        while ((i < sfm.Length - 1) && (count < 3))
                        {
                            i++;
                            if (sfm[i] > ITL)
                                count++;
                        }
                        finish = i - 1;
                        for (int j = start; j <= finish; j++)
                        {
                            sfmSpeech[j] = true;
                        }
                    }
                    start = -1;
                    finish = -1;
                }
                else
                {
                    noiseBuff[buffInd] = sfm[i];
                    buffInd++;
                    if (buffInd >= 30)
                        buffInd = 0;
                    double avg = noiseBuff.Average();
                    double sumsq = noiseBuff.Select(val => (val - avg) * (val - avg)).Sum();
                    double sigmaNew = Math.Sqrt(sumsq / (noiseBuff.Length - 1));
                    if (sigmaNew / sigmaOld >= 1.25)
                        p = 0.25;
                    else
                    {
                        if (sigmaNew / sigmaOld >= 1.1)
                            p = 0.2;
                        else
                        {
                            if (sigmaNew / sigmaOld >= 1.0)
                                p = 0.15;
                            else
                                p = 0.1;
                        }
                    }
                    l1 = 0.03 * (max - avg) + avg;
                    l2 = 4 * avg;
                    ITL = (1 - p) * ITL + p * Math.Max(l1, l2);
                    sigmaOld = sigmaNew;
                }
            }
            //BorderCorrection(sfmSpeech);

            bool[] result = new bool[energy.Length];
            for (int i = 0; i < result.Length; i++)
                result[i] = false;

            for (int i = silentFrames; i < energy.Length; i++)
            {
                int count = 0;
                if (energySpeech[i])
                    count++;
                if ((frequency[i] >= frequencyThreshFloor) && (frequency[i] <= frequencyThreshCeil))
                    count++;
                if (sfmSpeech[i])
                    count++;
                if (count > 1)
                    result[i] = true;
            }

            BorderCorrection(result);
            return result;
        }

        /// <summary>
        /// Алгоритм определения фрагментов речи на основе меры изменчивости сигнала
        /// </summary>
        /// <returns>Границы речи (в миллисекундах)</returns>
        public List<int[]> AlgorithmLTSV()
        {
            return GetBounds(AlgorithmLTSV_());
        }
        
        private bool[] AlgorithmLTSV_()
        {
            double[] ltsv = parameters.GetLongTermSignalVariability();
            bool[] result = new bool[ltsv.Length];

            double[] gamma = new double[ltsv.Length];
            int numOfBuff = 60;
            int nn = numOfBuff / silentFrames - 1;

            double[] nPsi = new double[numOfBuff];
            double[] sPsi = new double[numOfBuff];
            for (int i = 0; i < nn; i++)
            {
                for (int j = 0; j < silentFrames; j++)
                {
                    nPsi[i * silentFrames + j] = ltsv[j];
                }
            }
            for (int j = nn * silentFrames; j < nPsi.Length; j++)
            {
                nPsi[j] = ltsv[j - nn * silentFrames];
            }
            for (int i = 0; i < numOfBuff; i++)
            {
                sPsi[i] = 0;
            }
            double alpha = 0.3;
            for (int i = 0; i < silentFrames; i++)
            {
                gamma[i] = nPsi.Max();
            }

            int R = 30;
            for (int i = silentFrames; i < ltsv.Length - R; i++)
            {
                double count = 0;
                for (int j = i; j < i + R; j++)
                {
                    if (ltsv[j] > gamma[i - 1])
                        count = count + 1;
                }
                count = (double)count / (double)R;
                if (count >= 0.8)
                {
                    result[i] = true;
                }
                else
                {
                    result[i] = false;
                }

                if (ltsv[i] > gamma[i - 1])
                {
                    if (sPsi[1] == 0)
                    {
                        for (int j = 0; j < numOfBuff; j++)
                        {
                            sPsi[j] = ltsv[i];
                        }
                    }
                    else
                    {
                        for (int j = numOfBuff - 1; j > 0; j--)
                        {
                            sPsi[j] = sPsi[j - 1];
                        }
                        sPsi[0] = ltsv[i];
                    }
                }
                else
                {
                    for (int j = numOfBuff - 1; j > 0; j--)
                    {
                        nPsi[j] = nPsi[j - 1];
                    }
                    nPsi[0] = ltsv[i];
                }
                if (sPsi[0] == 0)
                {
                    gamma[i] = nPsi.Max();
                }
                else
                {
                    gamma[i] = alpha * sPsi.Min() + (1 - alpha) * nPsi.Max();
                }
            }
            BorderCorrection(result);
            return result;
        }

        /// <summary>
        /// Алгоритм определения фрагментов речи на основе долгосрочной меры спектральной плоскостности
        /// </summary>
        /// <returns>Границы речи (в миллисекундах)</returns>
        public List<int[]> AlgorithmLSFM()
        {
            return GetBounds(AlgorithmLSFM_());
        }
        
        private bool[] AlgorithmLSFM_()
        {
            double[] lsfm = parameters.GetLongTermSpectralFlatnessMeasure();
            bool[] result = new bool[lsfm.Length];

            double[] gamma = new double[lsfm.Length];
            int numOfBuff = 60;
            int nn = numOfBuff / silentFrames - 1;

            double[] nPsi = new double[numOfBuff];
            double[] sPsi = new double[numOfBuff];
            for (int i = 0; i < nn; i++)
            {
                for (int j = 0; j < silentFrames; j++)
                {
                    nPsi[i * silentFrames + j] = lsfm[j];
                }
            }
            for (int j = nn * silentFrames; j < nPsi.Length; j++)
            {
                nPsi[j] = lsfm[j - nn * silentFrames];
            }
            for (int i = 0; i < numOfBuff; i++)
            {
                sPsi[i] = 0;
            }
            double alpha = 0.3;
            for (int i = 0; i < silentFrames; i++)
            {
                gamma[i] = nPsi.Min();
            }

            int R = 30;
            for (int i = silentFrames; i < lsfm.Length - R; i++)
            {
                double count = 0;
                for (int j = i; j < i + R; j++)
                {
                    if (lsfm[j] < gamma[i - 1])
                        count = count + 1;
                }
                count = (double)count / (double)R;
                if (count >= 0.8)
                {
                    result[i] = true;
                }
                else
                {
                    result[i] = false;
                }

                if (lsfm[i] < gamma[i - 1])
                {
                    if (sPsi[1] == 0)
                    {
                        for (int j = 0; j < numOfBuff; j++)
                        {
                            sPsi[j] = lsfm[i];
                        }
                    }
                    else
                    {
                        for (int j = numOfBuff - 1; j > 0; j--)
                        {
                            sPsi[j] = sPsi[j - 1];
                        }
                        sPsi[0] = lsfm[i];
                    }
                }
                else
                {
                    for (int j = numOfBuff - 1; j > 0; j--)
                    {
                        nPsi[j] = nPsi[j - 1];
                    }
                    nPsi[0] = lsfm[i];
                }
                if (sPsi[0] == 0)
                {
                    gamma[i] = nPsi.Min();
                }
                else
                {
                    gamma[i] = alpha * sPsi.Max() + (1 - alpha) * nPsi.Min();
                }
            }
            BorderCorrection(result);
            return result;
        }

        /// <summary>
        /// Алгоритм определения фрагментов речи на основе усреднённой спектральной дельта-функции автокорреляции (Mean-Delta)
        /// </summary>
        /// <returns>Границы речи (в миллисекундах)</returns>
        public List<int[]> AlgorithmMeanDelta()
        {
            return GetBounds(AlgorithmMeanDelta_());
        }
        
        private bool[] AlgorithmMeanDelta_()
        {
            double[] meanD = parameters.GetMeanDelta();
            bool[] result = new bool[meanD.Length];

            double[][] gamma = new double[meanD.Length][];
            int numOfBuff = 60;
            int nn = numOfBuff / silentFrames - 1;

            double[] noiseBuff = new double[numOfBuff];
            double[] speechBuff = new double[numOfBuff];
            for (int i = 0; i < nn; i++)
            {
                for (int j = 0; j < silentFrames; j++)
                {
                    noiseBuff[i * silentFrames + j] = meanD[j];
                }
            }
            for (int j = nn * silentFrames; j < noiseBuff.Length; j++)
            {
                noiseBuff[j] = meanD[j - nn * silentFrames];
            }
            for (int i = 0; i < numOfBuff; i++)
            {
                speechBuff[i] = 0;
            }
            double[] g = GetThresholds(noiseBuff, speechBuff);
            for (int i = 0; i < silentFrames; i++)
            {
                gamma[i] = g;
            }

            for (int i = silentFrames; i < meanD.Length; i++)
            {
                if (meanD[i] >= gamma[i - 1][0])
                {
                    int first = i;
                    int count = 0;
                    while ((i < meanD.Length) && (meanD[i] < gamma[i - 1][1]) && (count < 4))
                    {
                        if (meanD[i] < gamma[i - 1][0])
                            count++;
                        gamma[i] = gamma[i - 1];
                        i++;
                    }

                    if ((count >= 4) || (i >= meanD.Length))
                    {
                        if (i - first < numOfBuff)
                        {
                            for (int jj = numOfBuff - 1; jj >= i - first; jj--)
                            {
                                noiseBuff[jj] = noiseBuff[jj - (i - first)];
                            }
                            for (int jj = i - first - 1; jj >= 0; jj--)
                            {
                                noiseBuff[jj] = meanD[jj + first];
                            }
                        }
                        else
                        {
                            for (int jj = 0; jj < numOfBuff; jj++)
                            {
                                noiseBuff[jj] = meanD[i - jj - 1];
                            }
                        }
                    }
                    else
                    {
                        for (int j = first; j <= i; j++)
                        {
                            result[j] = true;
                        }
                        count = 0;
                        while ((i < meanD.Length) && (meanD[i] >= gamma[i - 1][0]) && (count < 4))
                        {
                            if (meanD[i] < gamma[i - 1][0])
                                count++;
                            result[i] = true;
                            gamma[i] = gamma[i - 1];
                            i++;
                        }
                        if (i - first < numOfBuff)
                        {
                            for (int jj = numOfBuff - 1; jj >= i - first; jj--)
                            {
                                speechBuff[jj] = speechBuff[jj - (i - first)];
                            }
                            for (int jj = i - first - 1; jj >= 0; jj--)
                            {
                                speechBuff[jj] = meanD[jj + first];
                            }
                        }
                        else
                        {
                            for (int jj = 0; jj < numOfBuff; jj++)
                            {
                                speechBuff[jj] = meanD[i - jj - 1];
                            }
                        }
                    }
                }
                else
                {
                    for (int j = numOfBuff - 1; j > 0; j--)
                    {
                        noiseBuff[j] = noiseBuff[j - 1];
                    }
                    noiseBuff[0] = meanD[i];
                }
                if (i < meanD.Length)
                    gamma[i] = GetThresholds(noiseBuff, speechBuff);
            }
            BorderCorrection(result);
            return result;
        }

        /// <summary>
        /// Алгоритм определения фрагментов речи на основе энергии Тиджера
        /// </summary>
        /// <returns>Границы речи (в миллисекундах)</returns>
        public List<int[]> AlgorithmTeagerEnergy()
        {
            return GetBounds(AlgorithmTeagerEnergy_());
        }
        
        private bool[] AlgorithmTeagerEnergy_()
        {
            double[] tEnergy = parameters.GetTeagerEnergy();
            bool[] result = new bool[tEnergy.Length];

            double[][] gamma = new double[tEnergy.Length][];
            int numOfBuff = 60;
            int nn = numOfBuff / silentFrames - 1;

            double[] noiseBuff = new double[numOfBuff];
            double[] speechBuff = new double[numOfBuff];
            for (int i = 0; i < nn; i++)
            {
                for (int j = 0; j < silentFrames; j++)
                {
                    noiseBuff[i * silentFrames + j] = tEnergy[j];
                }
            }
            for (int j = nn * silentFrames; j < noiseBuff.Length; j++)
            {
                noiseBuff[j] = tEnergy[j - nn * silentFrames];
            }
            for (int i = 0; i < numOfBuff; i++)
            {
                speechBuff[i] = 0;
            }
            double[] g = GetThresholds(noiseBuff, speechBuff);
            for (int i = 0; i < silentFrames; i++)
            {
                gamma[i] = g;
            }

            for (int i = silentFrames; i < tEnergy.Length; i++)
            {
                if (tEnergy[i] >= gamma[i - 1][0])
                {
                    int first = i;
                    int count = 0;
                    while ((i < tEnergy.Length) && (tEnergy[i] < gamma[i - 1][1]) && (count < 4))
                    {
                        if (tEnergy[i] < gamma[i - 1][0])
                            count++;
                        gamma[i] = gamma[i - 1];
                        i++;
                    }

                    if ((count >= 4) || (i >= tEnergy.Length))
                    {
                        if (i - first < numOfBuff)
                        {
                            for (int jj = numOfBuff - 1; jj >= i - first; jj--)
                            {
                                noiseBuff[jj] = noiseBuff[jj - (i - first)];
                            }
                            for (int jj = i - first - 1; jj >= 0; jj--)
                            {
                                noiseBuff[jj] = tEnergy[jj + first];
                            }
                        }
                        else
                        {
                            for (int jj = 0; jj < numOfBuff; jj++)
                            {
                                noiseBuff[jj] = tEnergy[i - jj - 1];
                            }
                        }
                    }
                    else
                    {
                        for (int j = first; j <= i; j++)
                        {
                            result[j] = true;
                        }
                        count = 0;
                        while ((i < tEnergy.Length) && (tEnergy[i] >= gamma[i - 1][0]) && (count < 4))
                        {
                            if (tEnergy[i] < gamma[i - 1][0])
                                count++;
                            result[i] = true;
                            gamma[i] = gamma[i - 1];
                            i++;
                        }
                        if (i - first < numOfBuff)
                        {
                            for (int jj = numOfBuff - 1; jj >= i - first; jj--)
                            {
                                speechBuff[jj] = speechBuff[jj - (i - first)];
                            }
                            for (int jj = i - first - 1; jj >= 0; jj--)
                            {
                                speechBuff[jj] = tEnergy[jj + first];
                            }
                        }
                        else
                        {
                            for (int jj = 0; jj < numOfBuff; jj++)
                            {
                                speechBuff[jj] = tEnergy[i - jj - 1];
                            }
                        }
                    }
                }
                else
                {
                    for (int j = numOfBuff - 1; j > 0; j--)
                    {
                        noiseBuff[j] = noiseBuff[j - 1];
                    }
                    noiseBuff[0] = tEnergy[i];
                }
                if (i < tEnergy.Length)
                    gamma[i] = GetThresholds(noiseBuff, speechBuff);
            }
            BorderCorrection(result);
            return result;
        }

        /// <summary>
        /// Алгоритм определения фрагментов речи на основе энергии Тиджера и среднего числа пересечений нуля
        /// </summary>
        /// <returns>Границы речи (в миллисекундах)</returns>
        public List<int[]> AlgorithmTeagerEnergyZCR()
        {
            return GetBounds(AlgorithmTeagerEnergyZCR_());
        }
        
        private bool[] AlgorithmTeagerEnergyZCR_()
        {
            double[] tEnergy = parameters.GetTeagerEnergy();
            double[] zcr = parameters.GetZeroCrossingRate();

            bool[] result = new bool[tEnergy.Length];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = false;
            }
            double[][] gamma = new double[tEnergy.Length][];
            int numOfBuff = 60;
            int nn = numOfBuff / silentFrames - 1;

            double[] noiseBuff = new double[numOfBuff];
            double[] speechBuff = new double[numOfBuff];
            for (int i = 0; i < nn; i++)
            {
                for (int j = 0; j < silentFrames; j++)
                {
                    noiseBuff[i * silentFrames + j] = tEnergy[j];
                }
            }
            for (int j = nn * silentFrames; j < noiseBuff.Length; j++)
            {
                noiseBuff[j] = tEnergy[j - nn * silentFrames];
            }
            for (int i = 0; i < numOfBuff; i++)
            {
                speechBuff[i] = 0;
            }
            double[] g = GetThresholds(noiseBuff, speechBuff);
            for (int i = 0; i < silentFrames; i++)
            {
                gamma[i] = g;
            }

            for (int i = silentFrames; i < tEnergy.Length; i++)
            {
                if (tEnergy[i] >= gamma[i - 1][0])
                {
                    int first = i;
                    int count = 0;
                    while ((i < tEnergy.Length) && (tEnergy[i] < gamma[i - 1][1]) && (count < 4))
                    {
                        if (tEnergy[i] < gamma[i - 1][0])
                            count++;
                        gamma[i] = gamma[i - 1];
                        i++;
                    }

                    if ((count >= 4) || (i >= tEnergy.Length))
                    {
                        if (i - first < numOfBuff)
                        {
                            for (int jj = numOfBuff - 1; jj >= i - first; jj--)
                            {
                                noiseBuff[jj] = noiseBuff[jj - (i - first)];
                            }
                            for (int jj = i - first - 1; jj >= 0; jj--)
                            {
                                noiseBuff[jj] = tEnergy[jj + first];
                            }
                        }
                        else
                        {
                            for (int jj = 0; jj < numOfBuff; jj++)
                            {
                                noiseBuff[jj] = tEnergy[i - jj - 1];
                            }
                        }
                    }
                    else
                    {
                        for (int j = first; j <= i; j++)
                        {
                            result[j] = true;
                        }
                        count = 0;
                        while ((i < tEnergy.Length) && (tEnergy[i] >= gamma[i - 1][0]) && (count < 4))
                        {
                            if (tEnergy[i] < gamma[i - 1][0])
                                count++;
                            result[i] = true;
                            gamma[i] = gamma[i - 1];
                            i++;
                        }
                        if (i - first < numOfBuff)
                        {
                            for (int jj = numOfBuff - 1; jj >= i - first; jj--)
                            {
                                speechBuff[jj] = speechBuff[jj - (i - first)];
                            }
                            for (int jj = i - first - 1; jj >= 0; jj--)
                            {
                                speechBuff[jj] = tEnergy[jj + first];
                            }
                        }
                        else
                        {
                            for (int jj = 0; jj < numOfBuff; jj++)
                            {
                                speechBuff[jj] = tEnergy[i - jj - 1];
                            }
                        }
                    }
                }
                else
                {
                    for (int j = numOfBuff - 1; j > 0; j--)
                    {
                        noiseBuff[j] = noiseBuff[j - 1];
                    }
                    noiseBuff[0] = tEnergy[i];
                }
                if (i < tEnergy.Length)
                    gamma[i] = GetThresholds(noiseBuff, speechBuff);
            }
            BorderCorrection(result);

            double m = 0;
            double D = 0;
            for (int i = 0; i < silentFrames; i++)
            {
                m += zcr[i];
            }
            m = m / silentFrames;
            for (int i = 0; i < silentFrames; i++)
            {
                D += (zcr[i] - m) * (zcr[i] - m);
            }
            D = D / (silentFrames - 1);
            double zcrThresh = Math.Min(25.0 / frameSize, m + 2 * Math.Sqrt(D));

            for (int i = 1; i < zcr.Length - 1; i++)
            {
                int count = 0;
                int last = -1;
                if ((result[i]) && (!result[i - 1]))
                {
                    for (int j = i - 1; (j > i - 10) && (j >= 0); j--)
                    {
                        if (zcr[j] >= zcrThresh)
                        {
                            count++;
                            last = j;
                        }
                    }
                    if (count >= 3)
                    {
                        for (int j = i - 1; j >= last; j--)
                        {
                            result[j] = true;
                        }
                    }
                }
                count = 0;
                last = -1;
                if ((result[i]) && (!result[i + 1]))
                {
                    for (int j = i + 1; (j < i + 10) && (j < zcr.Length); j++)
                    {
                        if (zcr[j] >= zcrThresh)
                        {
                            count++;
                            last = j;
                        }
                    }
                    if (count >= 3)
                    {
                        for (int j = i + 1; j <= last; j++)
                        {
                            result[j] = true;
                        }
                    }
                }
            }
            return result;
        }

        /// <summary>
        /// Алгоритм определения фрагментов речи на основе энергии Тиджера, меры спектральной плоскостности, усреднённой спектральной дельта-функции автокорреляции (Mean-Delta) и среднего числа пересечений нуля
        /// </summary>
        /// <returns>Границы речи (в миллисекундах)</returns>
        public List<int[]> AlgorithmTeagerEnergySFMMeanDeltaZCR()
        {
            return GetBounds(AlgorithmTeagerEnergySFMMeanDeltaZCR_());
        }
        
        private bool[] AlgorithmTeagerEnergySFMMeanDeltaZCR_()
        {
            double[] tEnergy = parameters.GetTeagerEnergy();
            double[] zcr = parameters.GetZeroCrossingRate();
            double[] sfm = parameters.GetSpectralFlatnessMeasure();

            double varSfm = 0;
            double meanSfm = sfm.Average();
            for (int i = 0; i < sfm.Length; i++)
            {
                varSfm += (sfm[i] - meanSfm) * (sfm[i] - meanSfm);
            }
            varSfm /= (sfm.Length - 1);

            bool[] result = new bool[tEnergy.Length];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = false;
            }
            double[][] gamma = new double[tEnergy.Length][];
            int numOfBuff = 60;
            int nn = numOfBuff / silentFrames - 1;

            double[] noiseBuff = new double[numOfBuff];
            double[] speechBuff = new double[numOfBuff];
            for (int i = 0; i < nn; i++)
            {
                for (int j = 0; j < silentFrames; j++)
                {
                    noiseBuff[i * silentFrames + j] = tEnergy[j];
                }
            }
            for (int j = nn * silentFrames; j < noiseBuff.Length; j++)
            {
                noiseBuff[j] = tEnergy[j - nn * silentFrames];
            }
            for (int i = 0; i < numOfBuff; i++)
            {
                speechBuff[i] = 0;
            }
            double[] g = GetThresholds(noiseBuff, speechBuff);
            for (int i = 0; i < silentFrames; i++)
            {
                gamma[i] = g;
            }

            {
                for (int i = silentFrames; i < tEnergy.Length; i++)
                {
                    if (tEnergy[i] >= gamma[i - 1][0])
                    {
                        int first = i;
                        int count = 0;
                        while ((i < tEnergy.Length) && (tEnergy[i] < gamma[i - 1][1]) && (count < 4))
                        {
                            if (tEnergy[i] < gamma[i - 1][0])
                                count++;
                            gamma[i] = gamma[i - 1];
                            i++;
                        }

                        if ((count >= 4) || (i >= tEnergy.Length))
                        {
                            if (i - first < numOfBuff)
                            {
                                for (int jj = numOfBuff - 1; jj >= i - first; jj--)
                                {
                                    noiseBuff[jj] = noiseBuff[jj - (i - first)];
                                }
                                for (int jj = i - first - 1; jj >= 0; jj--)
                                {
                                    noiseBuff[jj] = tEnergy[jj + first];
                                }
                            }
                            else
                            {
                                for (int jj = 0; jj < numOfBuff; jj++)
                                {
                                    noiseBuff[jj] = tEnergy[i - jj - 1];
                                }
                            }
                        }
                        else
                        {
                            for (int j = first; j <= i; j++)
                            {
                                result[j] = true;
                            }
                            count = 0;
                            while ((i < tEnergy.Length) && (tEnergy[i] >= gamma[i - 1][0]) && (count < 4))
                            {
                                if (tEnergy[i] < gamma[i - 1][0])
                                    count++;
                                result[i] = true;
                                gamma[i] = gamma[i - 1];
                                i++;
                            }
                            if (i - first < numOfBuff)
                            {
                                for (int jj = numOfBuff - 1; jj >= i - first; jj--)
                                {
                                    speechBuff[jj] = speechBuff[jj - (i - first)];
                                }
                                for (int jj = i - first - 1; jj >= 0; jj--)
                                {
                                    speechBuff[jj] = tEnergy[jj + first];
                                }
                            }
                            else
                            {
                                for (int jj = 0; jj < numOfBuff; jj++)
                                {
                                    speechBuff[jj] = tEnergy[i - jj - 1];
                                }
                            }
                        }
                    }
                    else
                    {
                        for (int j = numOfBuff - 1; j > 0; j--)
                        {
                            noiseBuff[j] = noiseBuff[j - 1];
                        }
                        noiseBuff[0] = tEnergy[i];
                    }
                    if (i < tEnergy.Length)
                        gamma[i] = GetThresholds(noiseBuff, speechBuff);
                }
            }

            if ((varSfm < 0.13) || (varSfm > 1))
            {
                double[] meanD = parameters.GetMeanDelta();
                bool[] result2 = new bool[meanD.Length];

                double[][] gamma2 = new double[meanD.Length][];

                double[] noiseBuff2 = new double[numOfBuff];
                double[] speechBuff2 = new double[numOfBuff];
                for (int i = 0; i < nn; i++)
                {
                    for (int j = 0; j < silentFrames; j++)
                    {
                        noiseBuff2[i * silentFrames + j] = meanD[j];
                    }
                }
                for (int j = nn * silentFrames; j < noiseBuff2.Length; j++)
                {
                    noiseBuff2[j] = meanD[j - nn * silentFrames];
                }
                for (int i = 0; i < numOfBuff; i++)
                {
                    speechBuff2[i] = 0;
                }
                double[] g2 = GetThresholds(noiseBuff2, speechBuff2);
                for (int i = 0; i < silentFrames; i++)
                {
                    gamma2[i] = g2;
                }

                for (int i = silentFrames; i < meanD.Length; i++)
                {
                    if (meanD[i] >= gamma2[i - 1][0])
                    {
                        int first = i;
                        int count = 0;
                        while ((i < meanD.Length) && (meanD[i] < gamma2[i - 1][1]) && (count < 4))
                        {
                            if (meanD[i] < gamma2[i - 1][0])
                                count++;
                            gamma2[i] = gamma2[i - 1];
                            i++;
                        }

                        if ((count >= 4) || (i >= meanD.Length))
                        {
                            if (i - first < numOfBuff)
                            {
                                for (int jj = numOfBuff - 1; jj >= i - first; jj--)
                                {
                                    noiseBuff2[jj] = noiseBuff2[jj - (i - first)];
                                }
                                for (int jj = i - first - 1; jj >= 0; jj--)
                                {
                                    noiseBuff2[jj] = meanD[jj + first];
                                }
                            }
                            else
                            {
                                for (int jj = 0; jj < numOfBuff; jj++)
                                {
                                    noiseBuff2[jj] = meanD[i - jj - 1];
                                }
                            }
                        }
                        else
                        {
                            for (int j = first; j <= i; j++)
                            {
                                result2[j] = true;
                            }
                            count = 0;
                            while ((i < meanD.Length) && (meanD[i] >= gamma2[i - 1][0]) && (count < 4))
                            {
                                if (meanD[i] < gamma2[i - 1][0])
                                    count++;
                                result2[i] = true;
                                gamma2[i] = gamma2[i - 1];
                                i++;
                            }
                            if (i - first < numOfBuff)
                            {
                                for (int jj = numOfBuff - 1; jj >= i - first; jj--)
                                {
                                    speechBuff2[jj] = speechBuff2[jj - (i - first)];
                                }
                                for (int jj = i - first - 1; jj >= 0; jj--)
                                {
                                    speechBuff2[jj] = meanD[jj + first];
                                }
                            }
                            else
                            {
                                for (int jj = 0; jj < numOfBuff; jj++)
                                {
                                    speechBuff2[jj] = meanD[i - jj - 1];
                                }
                            }
                        }
                    }
                    else
                    {
                        for (int j = numOfBuff - 1; j > 0; j--)
                        {
                            noiseBuff2[j] = noiseBuff2[j - 1];
                        }
                        noiseBuff2[0] = meanD[i];
                    }
                    if (i < meanD.Length)
                        gamma2[i] = GetThresholds(noiseBuff2, speechBuff2);
                }
                BorderCorrection(result2);

                for (int i = 0; i < result2.Length; i++)
                {
                    if (result2[i])
                        result[i] = true;
                }
            }
            BorderCorrection(result);

            double m = 0;
            double D = 0;
            for (int i = 0; i < silentFrames; i++)
            {
                m += zcr[i];
            }
            m = m / silentFrames;
            for (int i = 0; i < silentFrames; i++)
            {
                D += (zcr[i] - m) * (zcr[i] - m);
            }
            D = D / (silentFrames - 1);
            double zcrThresh = Math.Min(25.0 / frameSize, m + 2 * Math.Sqrt(D));

            for (int i = 1; i < zcr.Length - 1; i++)
            {
                int count = 0;
                int last = -1;
                if ((result[i]) && (!result[i - 1]))
                {
                    for (int j = i - 1; (j > i - 10) && (j >= 0); j--)
                    {
                        if (zcr[j] >= zcrThresh)
                        {
                            count++;
                            last = j;
                        }
                    }
                    if (count >= 3)
                    {
                        for (int j = i - 1; j >= last; j--)
                        {
                            result[j] = true;
                        }
                    }
                }
                count = 0;
                last = -1;
                if ((result[i]) && (!result[i + 1]))
                {
                    for (int j = i + 1; (j < i + 10) && (j < zcr.Length); j++)
                    {
                        if (zcr[j] >= zcrThresh)
                        {
                            count++;
                            last = j;
                        }
                    }
                    if (count >= 3)
                    {
                        for (int j = i + 1; j <= last; j++)
                        {
                            result[j] = true;
                        }
                    }
                }
            }
            return result;
        }

        /// <summary>
        /// Алгоритм определения фрагментов речи на основе энергии и энтропии
        /// </summary>
        /// <returns>Границы речи (в миллисекундах)</returns>
        public List<int[]> AlgorithmEnergyEntropy()
        {
            return GetBounds(AlgorithmEnergyEntropy_());
        }
        
        private bool[] AlgorithmEnergyEntropy_()
        {
            double[] enen = parameters.GetEnergyEntropy();
            bool[] result = new bool[enen.Length];

            double[][] gamma = new double[enen.Length][];
            int numOfBuff = 60;
            int nn = numOfBuff / silentFrames - 1;

            double[] noiseBuff = new double[numOfBuff];
            double[] speechBuff = new double[numOfBuff];
            for (int i = 0; i < nn; i++)
            {
                for (int j = 0; j < silentFrames; j++)
                {
                    noiseBuff[i * silentFrames + j] = enen[j];
                }
            }
            for (int j = nn * silentFrames; j < noiseBuff.Length; j++)
            {
                noiseBuff[j] = enen[j - nn * silentFrames];
            }
            for (int i = 0; i < numOfBuff; i++)
            {
                speechBuff[i] = 0;
            }
            double[] g = GetThresholds(noiseBuff, speechBuff);
            for (int i = 0; i < silentFrames; i++)
            {
                gamma[i] = g;
            }

            for (int i = silentFrames; i < enen.Length; i++)
            {
                if (enen[i] >= gamma[i - 1][0])
                {
                    int first = i;
                    int count = 0;
                    while ((i < enen.Length) && (enen[i] < gamma[i - 1][1]) && (count < 4))
                    {
                        if (enen[i] < gamma[i - 1][0])
                            count++;
                        gamma[i] = gamma[i - 1];
                        i++;
                    }

                    if ((count >= 4) || (i >= enen.Length))
                    {
                        if (i - first < numOfBuff)
                        {
                            for (int jj = numOfBuff - 1; jj >= i - first; jj--)
                            {
                                noiseBuff[jj] = noiseBuff[jj - (i - first)];
                            }
                            for (int jj = i - first - 1; jj >= 0; jj--)
                            {
                                noiseBuff[jj] = enen[jj + first];
                            }
                        }
                        else
                        {
                            for (int jj = 0; jj < numOfBuff; jj++)
                            {
                                noiseBuff[jj] = enen[i - jj - 1];
                            }
                        }
                    }
                    else
                    {
                        for (int j = first; j <= i; j++)
                        {
                            result[j] = true;
                        }
                        count = 0;
                        while ((i < enen.Length) && (enen[i] >= gamma[i - 1][0]) && (count < 4))
                        {
                            if (enen[i] < gamma[i - 1][0])
                                count++;
                            result[i] = true;
                            gamma[i] = gamma[i - 1];
                            i++;
                        }
                        if (i - first < numOfBuff)
                        {
                            for (int jj = numOfBuff - 1; jj >= i - first; jj--)
                            {
                                speechBuff[jj] = speechBuff[jj - (i - first)];
                            }
                            for (int jj = i - first - 1; jj >= 0; jj--)
                            {
                                speechBuff[jj] = enen[jj + first];
                            }
                        }
                        else
                        {
                            for (int jj = 0; jj < numOfBuff; jj++)
                            {
                                speechBuff[jj] = enen[i - jj - 1];
                            }
                        }
                    }
                }
                else
                {
                    for (int j = numOfBuff - 1; j > 0; j--)
                    {
                        noiseBuff[j] = noiseBuff[j - 1];
                    }
                    noiseBuff[0] = enen[i];
                }
                if (i < enen.Length)
                    gamma[i] = GetThresholds(noiseBuff, speechBuff);
            }
            BorderCorrection(result);
            return result;
        }

        /// <summary>
        /// Алгоритм определения фрагментов речи на основе энергии Тиджера, долгосрочной меры спектральной плоскостности и доминирующей частоты
        /// </summary>
        /// <returns>Границы речи (в миллисекундах)</returns>
        public List<int[]> AlgorithmTeagerEnergyLSFMFrequency()
        {
            return GetBounds(AlgorithmTeagerEnergyLSFMFrequency_());
        }

        private bool[] AlgorithmTeagerEnergyLSFMFrequency_()
        {
            bool[] energySpeech = AlgorithmTeagerEnergyZCR_();
            bool[] sfmSpeech = AlgorithmLSFM_();
            double[] frequency = parameters.GetFrequency();
            bool[] result = new bool[energySpeech.Length];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = false;
            }
            for (int i = silentFrames; i < energySpeech.Length; i++)
            {
                int count = 0;
                if (energySpeech[i])
                    count++;
                if ((frequency[i] >= frequencyThreshFloor) && (frequency[i] <= frequencyThreshCeil))
                    count++;
                if (sfmSpeech[i])
                    count++;
                if (count > 1)
                    result[i] = true;
            }

            BorderCorrection(result);
            return result;
        }

        /// <summary>
        /// Алгоритм определения фрагментов речи на основе энергии Тиджера, долгосрочной меры спектральной плоскостности, доминирующей частоты, меры изменчивости сигнала и характеристики Mean Delta
        /// </summary>
        /// <returns>Границы речи (в миллисекундах)</returns>
        public List<int[]> AlgorithmTeagerEnergyLSFMFrequencyLTSVMeanDelta()
        {
            return GetBounds(AlgorithmTeagerEnergyLSFMFrequencyLTSVMeanDelta_());
        }

        private bool[] AlgorithmTeagerEnergyLSFMFrequencyLTSVMeanDelta_()
        {
            bool[] energySpeech = AlgorithmTeagerEnergyZCR_();
            bool[] sfmSpeech = AlgorithmLSFM_();
            bool[] ltsvSpeech = AlgorithmLTSV_();
            bool[] mdSpeech = AlgorithmMeanDelta_();

            double[] frequency = parameters.GetFrequency();
            bool[] result = new bool[energySpeech.Length];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = false;
            }
            for (int i = silentFrames; i < energySpeech.Length; i++)
            {
                int count = 0;
                if (energySpeech[i])
                    count++;
                if ((frequency[i] >= frequencyThreshFloor) && (frequency[i] <= frequencyThreshCeil))
                    count++;
                if (sfmSpeech[i])
                    count++;
                if (ltsvSpeech[i])
                    count++;
                if (mdSpeech[i])
                    count++;
                if (count > 2)
                    result[i] = true;
            }

            BorderCorrection(result);
            return result;
        }

        /// <summary>
        /// Алгоритм определения фрагментов речи на основе энергии Тиджера, долгосрочной меры спектральной плоскостности, доминирующей частоты и меры изменчивости сигнала
        /// </summary>
        /// <returns>Границы речи (в миллисекундах)</returns>
        public List<int[]> AlgorithmTeagerEnergyLSFMFrequencyLTSV()
        {
            return GetBounds(AlgorithmTeagerEnergyLSFMFrequencyLTSV_());
        }

        private bool[] AlgorithmTeagerEnergyLSFMFrequencyLTSV_()
        {
            bool[] energySpeech = AlgorithmTeagerEnergyZCR_();
            bool[] sfmSpeech = AlgorithmLSFM_();
            bool[] ltsvSpeech = AlgorithmLTSV_();

            double[] frequency = parameters.GetFrequency();
            bool[] result = new bool[energySpeech.Length];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = false;
            }
            for (int i = silentFrames; i < energySpeech.Length; i++)
            {
                int count = 0;
                if (energySpeech[i])
                    count++;
                if ((frequency[i] >= frequencyThreshFloor) && (frequency[i] <= frequencyThreshCeil))
                    count++;
                if (sfmSpeech[i])
                    count++;
                if (ltsvSpeech[i])
                    count++;
                if (count > 2)
                    result[i] = true;
            }

            BorderCorrection(result);
            return result;
        }

        /// <summary>
        /// Обновление пороговых значений
        /// </summary>
        /// <param name="noise">Буфер шума</param>
        /// <param name="speech">Буфер речи</param>
        /// <returns>Новые значения верхнего и нижнего порогов</returns>
        private double[] GetThresholds(double[] noise, double[] speech)
        {
            double[] thresholds = new double[2];
            double alpha = 0.03;
            double beta = 1.5;
            double gamma = 0.05;

            double muDown = noise.Average();
            double sigma = 0;
            for (int i = 0; i < noise.Length; i++)
            {
                sigma += (noise[i] - muDown) * (noise[i] - muDown);
            }
            sigma = Math.Sqrt(sigma / (noise.Length - 1));
            muDown += sigma;

            if (speech.Length == 0)
            {
                thresholds[0] = muDown;
                thresholds[1] = thresholds[0] + (beta) * (thresholds[0] - noise.Average());//beta * thresholds[0];
                return thresholds;
            }
            double muUp = speech.Average();
            if (muDown / muUp < gamma)
                muDown = gamma * muUp;

            thresholds[0] = muDown + alpha * (muUp - muDown);
            thresholds[1] = thresholds[0] + (beta) * (thresholds[0] - noise.Average());//beta * thresholds[0];//
            return thresholds;
        }
    }
}
