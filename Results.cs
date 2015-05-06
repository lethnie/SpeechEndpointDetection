using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SpeechEndpointDetection
{
    /// <summary>
    /// Характеристики точности полученных результатов
    /// </summary>
    public class Results
    {
        private int duration;
        private bool[] real;
        private bool[] current;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="duration">Продолжительность аудио (в секундах)</param>
        /// <param name="realSpeech">Границы реальной речи</param>
        /// <param name="currentSpeech">Полученные в результате выполнения алгоритма границы речи</param>
        public Results(double duration, List<int[]> realSpeech, List<int[]> currentSpeech)
        {
            this.duration = (int)(1000*duration);
            real = new bool[(int)Math.Ceiling((double)this.duration / 10.0)];
            current = new bool[(int)Math.Ceiling((double)this.duration / 10.0)];
            for (int i = 0; i < real.Length; i++)
            {
                real[i] = false;
                current[i] = false;
            }
            for (int i = 0; i < realSpeech.Count; i++)
            {
                for (int j = realSpeech[i][0] / 10; (j < realSpeech[i][1] / 10) && (j < real.Length); j++)
                {
                    real[j] = true;
                }
            }
            
            for (int i = 0; i < currentSpeech.Count; i++)
            {
                for (int j = currentSpeech[i][0] / 10; (j < currentSpeech[i][1] / 10) && (j < current.Length); j++)
                    {
                        current[j] = true;
                    }
            }
        }

        /// <summary>
        /// Абсолютная точность: отношение количества правильных решений к количеству всех принимаемых решений
        /// </summary>
        /// <returns></returns>
        public double GetFullAccuracy()
        {
            double length = (double)real.Length;
            double count = length;
            for (int i = 0; i < real.Length; i++)
            {
                if (real[i] != current[i])
                    count--;
            }
            return count / length;
        }       

        /// <summary>
        /// Точность принятия решений в пользу отсутствия речи. Отношение правильных решений об отсутствии речи к необходимому числу таких решений
        /// </summary>
        /// <returns></returns>
        public double GetNonSpeechHitRate()
        {
            double length = 0.0;
            double count = 0.0;
            for (int i = 0; i < real.Length; i++)
            {
                if (!real[i])
                    length = length + 1;
                if (!real[i] && !current[i])
                    count++;
            }
            return count / length;
        }

        /// <summary>
        /// Точность принятия решений в пользу присутствия речи. Отношение правильных решений о наличии речи к необходимому числу таких решений
        /// </summary>
        /// <returns></returns>
        public double GetSpeechHitRate()
        {
            double length = 0.0;
            double count = 0.0;
            for (int i = 0; i < real.Length; i++)
            {
                if (real[i])
                    length = length + 1;
                if (real[i] && current[i])
                    count++;
            }
            return count / length;
        }

        private int boundaryFrames = 15;
        
        /// <summary>
        /// Ошибки в начале речи
        /// </summary>
        /// <returns></returns>
        public double GetBeginningErrors()
        {
            double length = 0;
            double count = 0.0;
            bool detected = false;
            for (int i = 0; i < real.Length; i++)
            {
                if (real[i] != current[i])
                    length++;
                if ((!detected) && (real[i]))
                {
                    int j = i;
                    while ((i < real.Length - 1) && (Math.Abs(i - j) < boundaryFrames))
                    {
                        i++;
                        if (real[i] != current[i])
                        {
                            count++;
                            length++;
                        }
                    }
                    
                    while ((j > 0) && (Math.Abs(i - j) < 2*boundaryFrames))
                    {
                        j--;
                        if (real[j] != current[j])
                            count++;
                    }
                    detected = true;
                    i--;
                }           
                if (!real[i])
                    detected = false;
            }
            return count / length;
        }

        /// <summary>
        /// Ошибки в конце речи
        /// </summary>
        /// <returns></returns>
        public double GetEndErrors()
        {
            double length = 0;
            double count = 0.0;
            bool detected = false;
            for (int i = real.Length - 1; i >= 0; i--)
            {
                if (real[i] != current[i])
                    length++;
                if ((!detected) && (real[i]))
                {
                    int j = i;
                    while ((i > 0) && (Math.Abs(i - j) < boundaryFrames))
                    {
                        i--;
                        if (real[i] != current[i])
                        {
                            count++;
                            length++;
                        }
                    }
                    while ((j < real.Length - 1) && (Math.Abs(i - j) < 2 * boundaryFrames))
                    {
                        j++;
                        if (real[j] != current[j])
                            count++;
                    }
                    detected = true;
                    i++;
                }
                if (!real[i])
                    detected = false;
            }
            return count / length;
        }

        /// <summary>
        /// Ошибки первого рода: неверное отнесение к речи
        /// </summary>
        /// <returns></returns>
        public double TypeIErrors()
        {
            double length = 0;
            double count = 0;
            for (int i = 0; i < real.Length; i++)
            {
                if (!real[i])
                    length++;
                if (!real[i] && current[i])
                    count++;
            }
            return count / length;
        }

        /// <summary>
        /// Ошибки второго рода: неверное отнесение к не-речи
        /// </summary>
        /// <returns></returns>
        public double TypeIIErrors()
        {
            double length = 0;
            double count = 0;
            for (int i = 0; i < real.Length; i++)
            {
                if (real[i])
                    length++;
                if (real[i] && !current[i])
                    count++;
            }
            return count / length;
        }
    }
}
