//
// RotorDrum.cs — helpers for the ROTATING parts (turbojet construction):
//
//   * voxShaft: drive shaft with mounting stubs and center bore.
//   * voxHubShellSegment: a z-segment of the hub flow-path drum (hollow,
//     HubWallMm wall) — each per-stage disc part owns its segment plus
//     the spacer ring under the following stator, like real bolted discs.
//   * voxDisc: the disc under one rotor row — bore fit over the shaft,
//     web, rim following the hub line, the broached fir-tree slots and
//     the tie-bolt circle. The disc is NOT the shaft: it has its own
//     bore and slides onto it.
//
// The disc geometry is a port of geometry_generator._cq_disc; the two
// routes must build the same disc. Three things that phase 10's
// interference check found in the CadQuery route are already avoided
// here by construction:
//
//   * the rim follows the HUB LINE (with RimReliefMm below the blade
//     platform), not a cylinder at the centre radius — the hub line
//     climbs, and a cylinder pokes into the blade upstream of centre;
//   * the tie-bolt circle is COMMON to the whole stack — computed per
//     disc it drifts outwards stage by stage and the bolts cross the
//     webs of the forward discs;
//   * the drum shell stops at the disc band instead of running straight
//     through the rim.
//

using System;
using System.Collections.Generic;
using System.Numerics;
using PicoGK;

namespace AxialCompressorDesigner
{
    public static class RotorDrum
    {
        /// <summary>Rim thickness as a multiple of the web (same 1.6 the
        /// CadQuery route uses).</summary>
        const float fRimOverWeb = 1.6f;

        public static Voxels voxShaft(AxialCompressorParameters p)
        {
            float[][] hub = p.HubLine;
            float fZ0 = hub[0][0] - p.ShaftStubMm;
            float fL = (hub[hub.Length - 1][0] - hub[0][0])
                      + 2f * p.ShaftStubMm;
            Voxels vox = LatticeUtils.voxCylinderZ(
                new Vector3(0f, 0f, fZ0), p.DrumInnerRadiusMm, fL);
            float fBore = p.ShaftBoreFrac * p.DrumInnerRadiusMm;
            if (fBore > 0.5f)
                vox -= LatticeUtils.voxCylinderZ(
                    new Vector3(0f, 0f, fZ0 - 1f), fBore, fL + 2f);
            return vox;
        }

        public static Voxels voxHubShellSegment(AxialCompressorParameters p,
                                                float fZ0, float fZ1)
        {
            float[][] hub = p.HubLine;
            Func<float, Vector2> fnOuter = t =>
            {
                float fZ = fZ0 + t * (fZ1 - fZ0);
                return new Vector2(fRHubAt(hub, fZ), fZ);
            };
            Voxels vox = LatticeUtils.voxRevolveZ(fnOuter, 0f, 64);
            Func<float, Vector2> fnInner = t =>
            {
                float fZ = (fZ0 - 2f) + t * ((fZ1 + 2f) - (fZ0 - 2f));
                float fR = MathF.Max(fRHubAt(hub, fZ) - p.HubWallMm,
                                     p.DrumInnerRadiusMm + 1f);
                return new Vector2(fR, fZ);
            };
            return vox - LatticeUtils.voxRevolveZ(fnInner, 0f, 64);
        }

        /// <summary>Web thickness of the disc under one row.</summary>
        public static float fWebMm(AxialCompressorParameters p, RowParams row)
            => Math.Clamp(p.DiskWebFrac * (row.ZTeMm - row.ZLeMm),
                          p.DiskWebMinMm, p.DiskWebMaxMm);

        /// <summary>Bore radius and inner rim radius of a disc. Split out
        /// because the tie-bolt circle must be known for the WHOLE stack
        /// before any disc is built.</summary>
        public static void DiscRadii(AxialCompressorParameters p,
                                     float fRHubC, out float fRBore,
                                     out float fRRimI)
        {
            fRBore = p.DrumInnerRadiusMm + p.FirTree.ClearanceMm;
            fRRimI = MathF.Max(fRHubC - 1.9f * p.FirTree.DepthMm,
                               fRBore + 6f);
        }

        /// <summary>Tie-bolt circle radius, COMMON to the whole stack:
        /// each disc has its own inner rim radius (it grows with the hub
        /// line), so a circle computed disc by disc leaves the bolts
        /// crossing the webs of the forward discs.</summary>
        public static float fTieBoltRadius(AxialCompressorParameters p,
                                           List<RowParams> aRotors)
        {
            float fRRimMin = float.MaxValue, fRBore = 0f;
            foreach (RowParams row in aRotors)
            {
                DiscRadii(p, fRHubAt(p.HubLine, row.ZCenterMm),
                          out fRBore, out float fRRimI);
                fRRimMin = MathF.Min(fRRimMin, fRRimI);
            }
            return 0.5f * (fRBore + fRRimMin);
        }

        /// <summary>Axial band the disc and its blade roots occupy. The
        /// drum shell must stop here; inside it the disc rim closes the
        /// flow path.</summary>
        public static void DiscBand(AxialCompressorParameters p,
                                    RowParams row, out float fZLo,
                                    out float fZHi)
        {
            FirTree.RootAxialAndSlope(p, row, out float fAx, out _);
            float fHalf = 0.5f * MathF.Max(
                p.FirTree.Enabled ? FirTree.fSlotAxialFactor * fAx : 0f,
                fRimOverWeb * fWebMm(p, row));
            fZLo = row.ZCenterMm - fHalf;
            fZHi = row.ZCenterMm + fHalf;
        }

        /// <summary>The disc under one rotor row.</summary>
        public static Voxels voxDisc(AxialCompressorParameters p,
                                     RowParams row, float fRTie)
        {
            float fWeb = fWebMm(p, row);
            float fRimT = fRimOverWeb * fWeb;
            float fZc = row.ZCenterMm;
            float fRHubC = fRHubAt(p.HubLine, fZc);
            DiscRadii(p, fRHubC, out float fRBore, out float fRRimI);

            // web: annulus bore → inner rim, over the web thickness
            Voxels vox = LatticeUtils.voxCylinderZ(
                new Vector3(0f, 0f, fZc - 0.5f * fWeb), fRRimI, fWeb);

            // rim: annulus inner rim → hub line minus the relief, over the
            // (thicker) rim band. The relief is what keeps the rim out of
            // the blade: the platform forms the gas path, not the rim.
            float fZ0r = fZc - 0.5f * fRimT;
            float fZ1r = fZc + 0.5f * fRimT;
            Voxels voxRim = LatticeUtils.voxRevolveZ(t =>
            {
                float fZ = fZ0r + t * (fZ1r - fZ0r);
                return new Vector2(fRHubAt(p.HubLine, fZ) - p.RimReliefMm,
                                   fZ);
            }, 0f, 32);
            voxRim -= LatticeUtils.voxCylinderZ(
                new Vector3(0f, 0f, fZ0r - 1f), fRRimI, fRimT + 2f);
            vox += voxRim;

            // bore: the disc slides onto the shaft with the assembly fit
            vox -= LatticeUtils.voxCylinderZ(
                new Vector3(0f, 0f, fZ0r - 1f), fRBore, fRimT + 2f);

            // tie-bolt holes through the web
            for (int k = 0; k < p.TieBoltCount; k++)
            {
                float fA = 2f * MathF.PI * k / p.TieBoltCount;
                vox -= LatticeUtils.voxCylinderZ(
                    new Vector3(fRTie * MathF.Cos(fA), fRTie * MathF.Sin(fA),
                                fZc - fWeb),
                    0.5f * p.TieBoltDMm, 2f * fWeb);
            }

            // broached fir-tree slots, one per blade
            if (p.FirTree.Enabled && row.BladeCount > 0)
                vox -= FirTree.voxSlotRing(p, row, fRHubC);
            return vox;
        }

        /// <summary>Hub radius at axial position (clamped interpolation).</summary>
        internal static float fRHubAt(float[][] line, float fZ)
        {
            if (fZ <= line[0][0])
                return line[0][1];
            for (int i = 1; i < line.Length; i++)
            {
                if (fZ <= line[i][0])
                {
                    float t = (fZ - line[i - 1][0])
                            / MathF.Max(line[i][0] - line[i - 1][0], 1e-6f);
                    return line[i - 1][1] + t * (line[i][1] - line[i - 1][1]);
                }
            }
            return line[line.Length - 1][1];
        }
    }
}
