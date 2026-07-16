//
// RotorDrum.cs — helpers for the ROTATING parts (turbojet construction):
//
//   * voxShaft: drive shaft with mounting stubs and center bore.
//   * voxHubShellSegment: a z-segment of the hub flow-path drum (hollow,
//     HubWallMm wall) — each per-stage disc part owns its segment plus
//     the spacer ring under the following stator, like real bolted discs.
//   * voxDiscWeb: the disc web connecting shaft and hub shell under one
//     rotor row (the same disc model structures_core sizes at L0s).
//

using System;
using System.Numerics;
using PicoGK;

namespace AxialCompressorDesigner
{
    public static class RotorDrum
    {
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

        public static Voxels voxDiscWeb(AxialCompressorParameters p,
                                        RowParams row)
        {
            float fWeb = Math.Clamp(
                p.DiskWebFrac * (row.ZTeMm - row.ZLeMm),
                p.DiskWebMinMm, p.DiskWebMaxMm);
            float fRHub = fRHubAt(p.HubLine, row.ZCenterMm);
            Voxels vox = LatticeUtils.voxCylinderZ(
                new Vector3(0f, 0f, row.ZCenterMm - 0.5f * fWeb),
                fRHub, fWeb);
            // el disco calza sobre el eje: barreno al radio del eje
            vox -= LatticeUtils.voxCylinderZ(
                new Vector3(0f, 0f, row.ZCenterMm - 0.5f * fWeb - 1f),
                p.DrumInnerRadiusMm, fWeb + 2f);
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
