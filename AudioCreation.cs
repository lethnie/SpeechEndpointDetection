using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using SpeechEndpointDetection.Noises;
using System.Numerics;
using System.Diagnostics;

namespace SpeechEndpointDetection
{
    public static class AudioCreation
    {
        /// <summary>
        /// Получение нового сигнала на основе существующего с добавлением шума
        /// </summary>
        /// <param name="signal">Текущий сигнал</param>
        /// <param name="snr">Требуемое отношение мощности текущего сигнала к мощности шума</param>
        /// <param name="noiseType">Тип шума</param>
        /// <returns></returns>
        public static double[] GetNewSignal(double[] signal, int snr, Noises.Noises noiseType)
        {
            double coeff = 0;
            double[] result = new double[signal.Length];

            double sigRMS = 0;
            double maxSig = signal[0];
            for (int i = 0; i < signal.Length; i++)
            {
                if (signal[i] >= maxSig)
                {
                    sigRMS = signal[i] * signal[i];
                    maxSig = signal[i];
                }
            }
            double noiseRMS = 0;
            double[] noise = GetSignal(NoisesInfo.GetNoiseStream(noiseType), out noiseRMS);
            if (sigRMS != 0)
            {
                coeff = Math.Sqrt(sigRMS / (Math.Pow(10.0, snr / 10.0)) / noiseRMS);
            }
            else
            {
                coeff = 1;
            }
            double max = 0;
            int noiseNum = 0;
            for (int k = 0; k < signal.Length; k++)
            {
                if (noiseNum >= noise.Length)
                    noiseNum = 0;
                double sum = signal[k] + coeff * noise[noiseNum];
                result[k] = sum / 2;
                noiseNum++;
                if (Math.Abs(result[k]) > max)
                {
                    max = result[k];
                }
            }
            for (int k = 0; k < signal.Length; k++)
            {
                result[k] = result[k] / max;
            }
            return result;
        }

        /// <summary>
        /// Получение нового сигнала на основе существующего с добавлением другого сигнала
        /// </summary>
        /// <param name="signal">Текущий сигнал</param>
        /// <param name="snr">Требуемое отношение мощности текущего сигнала к мощности второго сигнала</param>
        /// <param name="newSignal">Добавляемый сигнал</param>
        /// <returns></returns>
        public static double[] GetNewSignal(double[] signal, int snr, Stream newSignal)
        {
            double coeff = 0;
            double[] result = new double[signal.Length];

            double sigRMS = 0;
            double maxSig = signal[0];
            for (int i = 0; i < signal.Length; i++)
            {
                if (signal[i] >= maxSig)
                {
                    sigRMS = signal[i] * signal[i];
                    maxSig = signal[i];
                }
            }
            double noiseRMS = 0;
            double[] noise = GetSignal(newSignal, out noiseRMS);
            if (sigRMS != 0)
            {
                coeff = Math.Sqrt(sigRMS / (Math.Pow(10, snr / 10)) / noiseRMS);
            }
            else
            {
                coeff = 1;
            }
            double max = 0;
            int noiseNum = 0;
            for (int k = 0; k < signal.Length; k++)
            {
                if (noiseNum >= noise.Length)
                    noiseNum = 0;
                double sum = signal[k] + coeff * noise[noiseNum];
                result[k] = sum / 2;
                noiseNum++;
                if (Math.Abs(result[k]) > max)
                {
                    max = result[k];
                }
            }
            for (int k = 0; k < signal.Length; k++)
            {
                result[k] = result[k] / max;
            }
            return result;
        }

        /// <summary>
        /// Чтение wav файла
        /// </summary>
        /// <param name="file"></param>
        /// <param name="rms">Мощность возвращаемого сигнала</param>
        /// <returns>Сигнал</returns>
        private static double[] GetSignal(Stream file, out double rms)
        {
            BinaryReader reader = new BinaryReader(file);
            int chunkID = reader.ReadInt32();
            int fileSize = reader.ReadInt32();
            int riffType = reader.ReadInt32();
            int fmtID = reader.ReadInt32();
            int fmtSize = reader.ReadInt32();
            int fmtCode = reader.ReadInt16();
            int channels = reader.ReadInt16();
            int sampleRate = reader.ReadInt32();
            int fmtAvgBPS = reader.ReadInt32();
            int fmtBlockAlign = reader.ReadInt16();
            int bitDepth = reader.ReadInt16();
            double r = Math.Pow(2.0, bitDepth - 1);
            if (fmtSize >= 18)
            {
                int fmtExtraSize = reader.ReadInt16();
                byte[] extraData = reader.ReadBytes(fmtExtraSize);
            }
            reader.BaseStream.Seek(fmtSize + 20, SeekOrigin.Begin);
            int dataID = reader.ReadInt32();
            int dataSize = reader.ReadInt32();
            byte[] data = reader.ReadBytes(dataSize);

            var channelLength = data.Length / 2;
            if (channels >= 2)
                channelLength = channelLength / channels;
            double[][] sChannel = new double[channels][];
            for (int j = 0; j < sChannel.Length; j++)
            {
                sChannel[j] = new double[channelLength];
            }
            {
                var x = 0;
                var length = data.Length;
                var maxS = length - 2 * channels + 1;
                for (int s = 0; s < maxS; )
                {
                    for (int j = 0; j < sChannel.Length; j++)
                    {
                        sChannel[j][x] = bytesToDouble(data[s], data[s + 1], r);
                        s += 2;
                    }
                    x++;
                }
            }
            reader.Close();

            double[] signal = new double[sChannel[0].Length];
            double max = 0;
            rms = 0;
            for (int i = 0; i < sChannel[0].Length; i++)
            {
                double sum = 0;
                for (int j = 0; j < sChannel.Length; j++)
                {
                    sum += sChannel[j][i];
                }
                signal[i] = sum / sChannel.Length;
                if (Math.Abs(signal[i]) > max)
                {
                    max = signal[i];
                }
            }
            for (int i = 0; i < signal.Length; i++)
            {
                signal[i] = signal[i] / max;
                rms = rms + signal[i] * signal[i];
            }
            rms = rms / signal.Length;
            return signal;
        }

        private static double bytesToDouble(byte firstByte, byte secondByte, double r)
        {
            short s = (short)((secondByte << 8) | firstByte);
            return s / r;
        }
    }
}
