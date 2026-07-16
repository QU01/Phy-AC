//
// Casing.cs — helpers for the STATIC parts: per-stage casing rings
// (shell segment + its stator vanes are unioned by AxialCompressor),
// bolted end flanges and the compressor bleed port boss.
//

using System;
using System.Numerics;
using PicoGK;

namespace AxialCompressorDesigner
{
    public static class Casing
    {
        internal const float fLipMm = 4f;   // axial lip beyond first/last row

        public static Voxels voxShellSegment(AxialCompressorParameters p,
                                             float fZ0, float fZ1)
        {
            float[][] tip = p.TipLine;
            Func<float, Vector2> fnOuter = t =>
            {
                float fZ = fZ0 + t * (fZ1 - fZ0);
                return new Vector2(
                    fRTipAt(tip, fZ) + p.TipClearanceMm + p.CasingWallMm, fZ);
            };
            Voxels vox = LatticeUtils.voxRevolveZ(fnOuter, 0f, 64);
            Func<float, Vector2> fnInner = t =>
            {
                float fZ = (fZ0 - 2f) + t * ((fZ1 + 2f) - (fZ0 - 2f));
                return new Vector2(fRTipAt(tip, fZ) + p.TipClearanceMm, fZ);
            };
            return vox - LatticeUtils.voxRevolveZ(fnInner, 0f, 64);
        }

        /// <summary>Bolted flange ring at axial position fZ (thickness
        /// towards +z if bForward). Welded onto the shell by union.</summary>
        public static Voxels voxFlange(AxialCompressorParameters p, float fZ,
                                       bool bForward)
        {
            float[][] tip = p.TipLine;
            float fRO = fRTipAt(tip, fZ) + p.TipClearanceMm + p.CasingWallMm;
            float fRF = fRO + p.FlangeWMm;
            float fZBase = bForward ? fZ : fZ - p.FlangeTMm;
            Voxels vox = LatticeUtils.voxCylinderZ(
                new Vector3(0f, 0f, fZBase), fRF, p.FlangeTMm);
            // anillo: hueco al radio del annulus + holgura
            vox -= LatticeUtils.voxCylinderZ(
                new Vector3(0f, 0f, fZBase - 1f),
                fRTipAt(tip, fZ) + p.TipClearanceMm, p.FlangeTMm + 2f);
            // barrenos en el círculo de pernos
            float fRB = fRO + 0.5f * p.FlangeWMm;
            for (int k = 0; k < p.FlangeBoltCount; k++)
            {
                float fA = 2f * MathF.PI * k / p.FlangeBoltCount;
                vox -= LatticeUtils.voxCylinderZ(
                    new Vector3(fRB * MathF.Cos(fA), fRB * MathF.Sin(fA),
                                fZBase - 1f),
                    0.5f * p.FlangeBoltDMm, p.FlangeTMm + 2f);
            }
            return vox;
        }

        /// <summary>Bleed port: radial boss with a through hole into the
        /// annulus at axial position fZ, angular position 0 (+x).</summary>
        public static Voxels voxBleedBoss(AxialCompressorParameters p,
                                          float fZ)
        {
            float[][] tip = p.TipLine;
            float fRIn = fRTipAt(tip, fZ) + p.TipClearanceMm;
            float fROut = fRIn + p.CasingWallMm + p.BleedBossHMm;
            Lattice lat = new Lattice();
            lat.AddBeam(new Vector3(fRIn + 1f, 0f, fZ),
                        0.5f * p.BleedBossDMm,
                        new Vector3(fROut, 0f, fZ),
                        0.5f * p.BleedBossDMm, false);
            return new Voxels(lat);
        }

        public static Voxels voxBleedHole(AxialCompressorParameters p,
                                          float fZ)
        {
            float[][] tip = p.TipLine;
            float fRIn = fRTipAt(tip, fZ) + p.TipClearanceMm;
            float fROut = fRIn + p.CasingWallMm + p.BleedBossHMm;
            Lattice lat = new Lattice();
            lat.AddBeam(new Vector3(fRIn - 2f, 0f, fZ),
                        0.5f * p.BleedHoleDMm,
                        new Vector3(fROut + 2f, 0f, fZ),
                        0.5f * p.BleedHoleDMm, false);
            return new Voxels(lat);
        }

        internal static float fRTipAt(float[][] line, float fZ)
            => RotorDrum.fRHubAt(line, fZ);   // same clamped interpolation
    }
}
