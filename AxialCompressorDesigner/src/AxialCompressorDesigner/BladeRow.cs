//
// BladeRow.cs — axial blade rows as distance-field-thickened camber sheets.
//
// Pattern inherited from Phy-CC Blades.cs v3 (see its HISTORY header):
// building closed pressure+suction meshes self-intersects wherever the
// local camber radius drops below half the thickness, and OpenVDB renders
// the fold-over as loose shards. The robust recipe is to mesh only the
// camber MID-SURFACE (one sheet of quads per blade) and thicken it with
// `Voxels.voxMeshShell`, which offsets in the signed-distance domain and
// cannot self-intersect. LE/TE and tip come out rounded at ½·thickness.
//
// Geometry: each section arrives as a camber polyline in its chord frame
// (x chordwise, y normal, centroid at the origin, mm) plus a SIGNED
// stagger (rotor +, stator −). For each span rib the polyline is rotated
// by the stagger, then wrapped onto the cylinder of its stacking radius:
//
//   phi = phi0 + y'/r      (arc length preserved)
//   P   = (r·cos phi, r·sin phi, zCenter + x')
//
// Rotor rows root at the hub (root rib sunk into the drum, tip rib pulled
// in by clearance + shell radius); stator rows hang from the casing
// (mirrored treatment).
//

using System;
using System.Numerics;
using PicoGK;

namespace AxialCompressorDesigner
{
    public static class BladeRow
    {
        /// <summary>Interpolated ribs inserted between adjacent sections.
        /// The EXACT free-vortex twist comes from Python (N_SPAN=13
        /// stations — with few stations the piecewise-linear span left
        /// visible creases at the original ribs); this subdivision only
        /// smooths the ruled surface between them.</summary>
        const int nSpanSubdiv = 2;

        /// <summary>
        /// Builds all blades of one row as a single voxel field.
        /// </summary>
        /// <param name="row">Row parameters (sections hub→tip).</param>
        /// <param name="fTipClearanceMm">Running clearance at the free end.</param>
        /// <param name="fRootSinkMm">How deep the root rib sinks into its body.</param>
        /// <param name="fMinThicknessMm">Lower bound on shell thickness
        /// (keep >= 2 voxels so thin rear-stage blades survive).</param>
        public static Voxels voxBuildRow(RowParams row, float fTipClearanceMm,
                                         float fRootSinkMm,
                                         float fMinThicknessMm)
        {
            float fThick = MathF.Max(row.MeanThicknessMm, fMinThicknessMm);
            float fHalf = 0.5f * fThick;

            int nSecs = row.Sections.Count;
            int nPts = row.Sections[0].CamberPoints.Length;

            // Stacking radii, adjusted at root (sink) and free end (pull in
            // by clearance + shell inflation so voxMeshShell does not eat
            // the running gap).
            var afR = new float[nSecs];
            for (int i = 0; i < nSecs; i++)
                afR[i] = row.Sections[i].RMm;
            if (row.Rotating)
            {
                afR[0] -= fRootSinkMm;                       // into the drum
                afR[nSecs - 1] -= fTipClearanceMm + fHalf;   // tip gap
            }
            else
            {
                afR[0] += fTipClearanceMm + fHalf;           // hub gap
                afR[nSecs - 1] += fRootSinkMm;               // into casing
            }

            // densify the span: interpolated ribs (r, stagger, camber)
            var aRibs = DensifyRibs(row, afR);

            Mesh msh = new Mesh();
            for (int nBlade = 0; nBlade < row.BladeCount; nBlade++)
            {
                float fPhi0 = 2f * MathF.PI * nBlade / row.BladeCount;
                AddBladeSheet(msh, aRibs, row.ZCenterMm, fPhi0);
            }
            return Voxels.voxMeshShell(msh, fHalf);
        }

        struct Rib
        {
            public float fR;
            public float fStaggerDeg;
            public float[][] afPts;
        }

        static Rib[] DensifyRibs(RowParams row, float[] afR)
        {
            int nSecs = row.Sections.Count;
            int nPts = row.Sections[0].CamberPoints.Length;
            int nRibs = (nSecs - 1) * nSpanSubdiv + 1;
            var aRibs = new Rib[nRibs];
            int k = 0;
            for (int i = 0; i < nSecs - 1; i++)
            {
                SectionParams a = row.Sections[i];
                SectionParams b = row.Sections[i + 1];
                for (int s = 0; s < nSpanSubdiv; s++)
                {
                    float t = s / (float)nSpanSubdiv;
                    var pts = new float[nPts][];
                    for (int j = 0; j < nPts; j++)
                        pts[j] = new float[]
                        {
                            a.CamberPoints[j][0] * (1 - t) + b.CamberPoints[j][0] * t,
                            a.CamberPoints[j][1] * (1 - t) + b.CamberPoints[j][1] * t,
                        };
                    aRibs[k++] = new Rib
                    {
                        fR = afR[i] * (1 - t) + afR[i + 1] * t,
                        fStaggerDeg = a.StaggerDeg * (1 - t) + b.StaggerDeg * t,
                        afPts = pts,
                    };
                }
            }
            SectionParams last = row.Sections[nSecs - 1];
            aRibs[k] = new Rib
            {
                fR = afR[nSecs - 1],
                fStaggerDeg = last.StaggerDeg,
                afPts = last.CamberPoints,
            };
            return aRibs;
        }

        static void AddBladeSheet(Mesh msh, Rib[] aRibs, float fZCenter,
                                  float fPhi0)
        {
            int nRibs = aRibs.Length;
            int nPts = aRibs[0].afPts.Length;
            var an = new int[nRibs, nPts];

            for (int i = 0; i < nRibs; i++)
            {
                Rib rib = aRibs[i];
                float fR = MathF.Max(rib.fR, 1f);
                float fG = rib.fStaggerDeg * MathF.PI / 180f;
                float fCosG = MathF.Cos(fG), fSinG = MathF.Sin(fG);
                for (int j = 0; j < nPts; j++)
                {
                    float x = rib.afPts[j][0];
                    float y = rib.afPts[j][1];
                    float xr = x * fCosG - y * fSinG;    // axial (staggered)
                    float yr = x * fSinG + y * fCosG;    // tangential
                    float fPhi = fPhi0 + yr / fR;        // wrap on cylinder
                    an[i, j] = msh.nAddVertex(new Vector3(
                        fR * MathF.Cos(fPhi),
                        fR * MathF.Sin(fPhi),
                        fZCenter + xr));
                }
            }

            // single sheet of quads (winding irrelevant for voxMeshShell)
            for (int i = 0; i < nRibs - 1; i++)
                for (int j = 0; j < nPts - 1; j++)
                    msh.AddQuad(an[i, j], an[i + 1, j],
                                an[i + 1, j + 1], an[i, j + 1]);
        }
    }
}
