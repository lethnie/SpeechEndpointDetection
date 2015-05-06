using System;
using System.IO;
using System.Numerics;
using System.Collections.Generic;
using SpeechEndpointDetection;

namespace SpeechEndpointDetection
{
    /// <summary>
    /// Чтение аудиодорожки, определение участков с речью
    /// </summary>
    public class WaveReader
    {
        //Параметры файла
        /// <summary>
        /// Содержит символы "RIFF". Является началом RIFF-цепочки
        /// </summary>
        private int chunkID;
        /// <summary>
        /// Oставшийся размер цепочки, начиная с этой позиции
        /// </summary>
        private int fileSize;
        /// <summary>
        /// Содержит символы "WAVE"
        /// </summary>
        private int riffType;
        /// <summary>
        /// Содержит символы "fmt "
        /// </summary>
        private int fmtID;
        /// <summary>
        /// Oставшийся размер подцепочки, начиная с этой позиции
        /// </summary>
        private int fmtSize;
        /// <summary>
        /// Аудиоформат
        /// </summary>
        private int fmtCode;
        /// <summary>
        /// Количество каналов. Моно = 1, Стерео = 2
        /// </summary>
        private int channels;
        /// <summary>
        /// Частота дискретизации. 8000 Гц, 44100 Гц и т.д.
        /// </summary>
        private int sampleRate;
        /// <summary>
        /// Количество байт, переданных за секунду воспроизведения
        /// </summary>
        private int fmtAvgBPS;
        /// <summary>
        /// Количество байт для одного сэмпла, включая все каналы
        /// </summary>
        private int fmtBlockAlign;
        /// <summary>
        /// Количество бит в сэмпле
        /// </summary>
        private int bitDepth;
        /// <summary>
        /// Дополнительные поля
        /// </summary>
        private int fmtExtraSize;
        /// <summary>
        /// Содержит символы "data" 
        /// </summary>
        private int dataID;
        /// <summary>
        /// Количество байт в области данных
        /// </summary>
        private int dataSize;
        /// <summary>
        /// Непосредственно WAV-данные
        /// </summary>
        private byte[] data;
        private byte[] extraData;
        /// <summary>
        /// Каналы
        /// </summary>
        private double[][] sChannel;
        /// <summary>
        /// Продолжительность аудио в секундах
        /// </summary>
        private double duration;
        /// <summary>
        /// Имя файла
        /// </summary>
        private string filename;
        /// <summary>
        /// Заголовок для речевых файлов
        /// </summary>
        private byte[] header = new byte[44];

        public WaveReader(string fileName)
        {
            filename = fileName;
            ReadWavfile();
            SetHeader();
            duration = (double)dataSize / (double)fmtAvgBPS;
        }

        /// <summary>
        /// Продолжительность аудио в секундах
        /// </summary>
        /// <returns></returns>
        public double GetDuration()
        {
            return duration;
        }

        /// <summary>
        /// Частота дискретизации
        /// </summary>
        /// <returns></returns>
        public int GetSampleRate()
        {
            return sampleRate;
        }

        /// <summary>
        /// Получение сигнала заданного аудио
        /// </summary>
        /// <returns>Сигнал</returns>
        public double[] GetSignal()
        {
            double[] signal = new double[sChannel[0].Length];
            double max = 0;
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
            }

            return signal;
        }

        private double r;

        private double bytesToDouble(byte firstByte, byte secondByte)
        {
            short s = (short) ((secondByte << 8) | firstByte);
            return s / r;
        }

        /// <summary>
        /// Создание нового аудиофайла-фрагмента текущего
        /// </summary>
        /// <param name="outFile">Имя создаваемого файла</param>
        /// <param name="start">Время начала фрагмента (мс)</param>
        /// <param name="finish">Время окончания фрагмента (мс)</param>
        public void CreateWaveFile(string outFile, int start, int finish)
        {
            int dSize = (int)Math.Ceiling((finish - start)*(fmtAvgBPS/1000.0)/8.0);
            byte[] result = GetHeaderByteArray(dSize);
            int st = 44;
            if (fmtSize >= 18)
            {
                st += 2 + fmtExtraSize;
            }
            
            int j = st;
            for (long bound = (long)(start * (fmtAvgBPS / 1000.0)); bound < (long)(finish * (fmtAvgBPS / 1000.0)); )
            {
                result[j] = data[bound];
                result[j + 1] = data[bound + 1];
                j += 2;
                bound += 2;
            }
            File.WriteAllBytes(outFile, result);
        }

        /// <summary>
        /// Чтение исходного файла аудио
        /// </summary>
        private void ReadWavfile()
        {
            BinaryReader reader = new BinaryReader(File.OpenRead(filename));
            chunkID = reader.ReadInt32();
            fileSize = reader.ReadInt32();
            riffType = reader.ReadInt32();
            fmtID = reader.ReadInt32();
            fmtSize = reader.ReadInt32();
            fmtCode = reader.ReadInt16();
            channels = reader.ReadInt16();
            sampleRate = reader.ReadInt32();
            fmtAvgBPS = reader.ReadInt32();
            fmtBlockAlign = reader.ReadInt16();
            bitDepth = reader.ReadInt16();
            r = Math.Pow(2.0, bitDepth - 1);
            if (fmtSize >= 18)
            {
                fmtExtraSize = reader.ReadInt16();
                extraData = reader.ReadBytes(fmtExtraSize);
            }
            reader.BaseStream.Seek(fmtSize + 20, SeekOrigin.Begin);
            dataID = reader.ReadInt32();
            dataSize = reader.ReadInt32();
            data = reader.ReadBytes(dataSize);

            var channelLength = data.Length / 2;
            if (channels >= 2)
                channelLength = channelLength / channels;
            sChannel = new double[channels][];
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
                        sChannel[j][x] = bytesToDouble(data[s], data[s + 1]);
                        s += 2;
                    }
                    x++;
                }
            }
            reader.Close();
        }

        /// <summary>
        /// Заготовка заголовка для речевых файлов
        /// </summary>
        private void SetHeader()
        {
            int num = 0;
            byte[] res = BitConverter.GetBytes(chunkID);
            for (int i = 0; i < res.Length; i++)
            {
                header[num] = res[i];
                num++;
            }
            num = 8;
            res = BitConverter.GetBytes(riffType);
            for (int i = 0; i < res.Length; i++)
            {
                header[num] = res[i];
                num++;
            }
            num = 12;
            res = BitConverter.GetBytes(fmtID);
            for (int i = 0; i < res.Length; i++)
            {
                header[num] = res[i];
                num++;
            }
            num = 16;
            res = BitConverter.GetBytes(fmtSize);
            for (int i = 0; i < res.Length; i++)
            {
                header[num] = res[i];
                num++;
            }
            num = 20;
            res = BitConverter.GetBytes(fmtCode);
            for (int i = 0; i < res.Length; i++)
            {
                header[num] = res[i];
                num++;
            }
            num = 22;
            res = BitConverter.GetBytes(channels);
            for (int i = 0; i < res.Length; i++)
            {
                header[num] = res[i];
                num++;
            }
            num = 24;
            res = BitConverter.GetBytes(sampleRate);
            for (int i = 0; i < res.Length; i++)
            {
                header[num] = res[i];
                num++;
            }
            num = 28;
            res = BitConverter.GetBytes(fmtAvgBPS);
            for (int i = 0; i < res.Length; i++)
            {
                header[num] = res[i];
                num++;
            }
            num = 32;
            res = BitConverter.GetBytes(fmtBlockAlign);
            for (int i = 0; i < res.Length; i++)
            {
                header[num] = res[i];
                num++;
            }
            num = 34;
            res = BitConverter.GetBytes(bitDepth);
            for (int i = 0; i < res.Length; i++)
            {
                header[num] = res[i];
                num++;
            }
            num = 36;
            if (fmtSize == 18)
            {
                res = BitConverter.GetBytes(fmtExtraSize);
                for (int i = 0; i < res.Length; i++)
                {
                    header[num] = res[i];
                    num++;
                }
                for (int i = 0; i < extraData.Length; i++)
                {
                    header[num] = extraData[i];
                    num++;
                }
                num = 38 + extraData.Length;
            }
            res = BitConverter.GetBytes(dataID);
            for (int i = 0; i < res.Length; i++)
            {
                header[num] = res[i];
                num++;
            }
        }

        /// <summary>
        /// Получение заголовка wav файла
        /// </summary>
        /// <param name="dSize">Размер области данных</param>
        /// <returns>Массив байтов, соответствующий заголовку</returns>
        private byte[] GetHeaderByteArray(int dSize)
        {
            int headerLength = 44;
            if (fmtSize == 18)
            {
                headerLength += 2 + fmtExtraSize;
            }
            int length = (headerLength + dSize) * 8; ;
            byte[] result = new byte[length];
            for (int i = 0; i < headerLength; i++)
            {
                result[i] = header[i];
                if (i == 4)
                {
                    byte[] res = BitConverter.GetBytes((dSize + 44 - 8) * 8);
                    for (int j = 0; j < res.Length; j++)
                    {
                        result[i] = res[j];
                        i++;
                    }
                    i--;
                }
                if (i == headerLength - 4)
                {
                    byte[] res = BitConverter.GetBytes(dSize * 8);
                    for (int j = 0; j < res.Length; j++)
                    {
                        result[i] = res[j];
                        i++;
                    }
                    i--;
                }
            }
            return result;
        }
    }
}
