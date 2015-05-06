using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Resources;
using System.Reflection;
using System.IO;

namespace SpeechEndpointDetection.Noises
{
    public enum Noises
    {
        Blue,
        Brown,
        Pink,
        White,
        MachineGun,
        Grey,
        Violet,
        WhiteGaussian,
        Bus,
        Cafeteria,
        Kids
    }

    public static class NoisesInfo
    {
        /// <summary>
        /// Возвращает аудио с шумом
        /// </summary>
        /// <param name="noise">Тип шума</param>
        /// <returns></returns>
        public static Stream GetNoiseStream(Noises noise)
        {
            Assembly assembly;
            assembly = Assembly.GetExecutingAssembly();
            switch (noise)
            {
                case Noises.Blue:
                    {
                        return assembly.GetManifestResourceStream("SpeechEndpointDetection.Noises.blue.wav");
                    }
                case Noises.Brown:
                    {
                        return assembly.GetManifestResourceStream("SpeechEndpointDetection.Noises.brown.wav");
                    }
                case Noises.Bus:
                    {
                        return assembly.GetManifestResourceStream("SpeechEndpointDetection.Noises.bus.wav");
                    }
                case Noises.Cafeteria:
                    {
                        return assembly.GetManifestResourceStream("SpeechEndpointDetection.Noises.cafeteria.wav");
                    }
                case Noises.Grey:
                    {
                        return assembly.GetManifestResourceStream("SpeechEndpointDetection.Noises.grey.wav");
                    }
                case Noises.Kids:
                    {
                        return assembly.GetManifestResourceStream("SpeechEndpointDetection.Noises.kids.wav");
                    }
                case Noises.MachineGun:
                    {
                        return assembly.GetManifestResourceStream("SpeechEndpointDetection.Noises.machinegun.wav");
                    }
                case Noises.Pink:
                    {
                        return assembly.GetManifestResourceStream("SpeechEndpointDetection.Noises.pink.wav");
                    }
                case Noises.Violet:
                    {
                        return assembly.GetManifestResourceStream("SpeechEndpointDetection.Noises.violet.wav");
                    }
                case Noises.White:
                    {
                        return assembly.GetManifestResourceStream("SpeechEndpointDetection.Noises.white.wav");
                    }
                case Noises.WhiteGaussian:
                    {
                        return assembly.GetManifestResourceStream("SpeechEndpointDetection.Noises.whitegaussian.wav");
                    }
                default:
                    {
                        return assembly.GetManifestResourceStream("SpeechEndpointDetection.Noises.white.wav");
                    }
            }
        }
    }
}
