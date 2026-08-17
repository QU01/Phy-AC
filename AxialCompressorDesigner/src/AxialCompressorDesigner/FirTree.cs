//
// FirTree.cs — fir-tree blade retention, layer 5c.
//
// This is a LINE-BY-LINE port of geometry_generator.firtree_profile and
// _cq_firtree_solid. The two must stay identical: the CadQuery path is
// the re-CAD route and this one is the manufacturing route, and until
// phase 10 they described different machines — the C# sank the blade
// root into the body (RootSinkMm) and knew nothing about the tree, so
// any structural margin computed on the fir-tree fillet (the Peterson
// K_t in structures_core already uses its radius) did not correspond to
// the part that gets printed.
//
// Convention, same as Python: u = tangential (mm, centred on 0),
// v = DEPTH inwards from the platform (0 at the platform, +depth at the
// bottom). The returned contour is closed and counter-clockwise.
//
// Whoever edits one side edits the other. test_phyac.py compares the two
// profiles point by point (block T20) and fails if they drift.
//

using System;
using System.Collections.Generic;
using System.Numerics;
using PicoGK;

namespace AxialCompressorDesigner
{
    public static class FirTree
    {
        /// <summary>Sections of the root/slot loft. With only the two end
        /// faces the ruled surface is a CHORD of the cylindrical wrap and
        /// its sagitta grows with axial length: the slot, 5% longer than
        /// the root, sank ~0.02 mm deeper and the root edge bit into the
        /// disc rim. Must match geometry_generator.N_FIRTREE_SECTIONS.
        /// </summary>
        public const int nSections = 5;

        /// <summary>Axial over-length of the broached slot vs the root.
        /// Must match the 1.05 factor in geometry_generator._cq_disc.
        /// </summary>
        public const float fSlotAxialFactor = 1.05f;

        /// <summary>
        /// Fir-tree profile in the (tangential, radial-depth) plane.
        /// Port of geometry_generator.firtree_profile — same formulas,
        /// same order, same closed CCW contour.
        /// </summary>
        /// <param name="fWTopMm">total width at the platform.</param>
        /// <param name="fDepthMm">radial depth of the root.</param>
        /// <param name="nTeeth">tooth pairs (2-3 in axial compressors).</param>
        /// <param name="fTaperDeg">half-angle of the tree.</param>
        /// <param name="fLobeFrac">neck depth as a fraction of the local
        /// half-width.</param>
        /// <param name="fClearanceMm">uniform oversize — use &gt; 0 for the
        /// DISC SLOT, which must be larger than the root to broach and to
        /// leave assembly fit.</param>
        public static float[][] afProfile(float fWTopMm, float fDepthMm,
                                          int nTeeth = 3,
                                          float fTaperDeg = 8f,
                                          float fLobeFrac = 0.22f,
                                          float fClearanceMm = 0f)
        {
            int n = Math.Max(nTeeth, 1);
            float fDepth = MathF.Max(fDepthMm, 1e-3f);
            float fP = fDepth / n;                    // radial pitch per tooth
            float fTanT = MathF.Tan(MathF.Max(fTaperDeg, 0f)
                                    * MathF.PI / 180f);
            float fW0 = 0.5f * MathF.Max(fWTopMm, 1e-3f) + fClearanceMm;

            float fHalfW(float v) => MathF.Max(fW0 - v * fTanT, 0.25f * fW0);

            // ABSOLUTE lobe: subtracting a fraction of the LOCAL width
            // leaves the lower teeth (where the tree is already thin)
            // unable to re-open, and the "fir tree" degenerates into a
            // stepped cone with no bearing flanks. The floor guarantees
            // the lobe always beats the taper drop between neck and the
            // next crest.
            float fLobeAbs = MathF.Max(fLobeFrac * fW0,
                                       1.4f * fTanT * 0.45f * fP);

            var aRight = new List<float[]> { new[] { fW0, 0f } };
            for (int k = 0; k < n; k++)
            {
                float fVCrest = k * fP;
                float fVNeck = fVCrest + 0.55f * fP;
                float fWc = fHalfW(fVCrest);
                // the neck is measured from the ENVELOPE at its own depth
                float fWn = MathF.Max(fHalfW(fVNeck) - fLobeAbs,
                                      0.18f * fW0);
                aRight.Add(new[] { fWc, fVCrest + 0.12f * fP });
                aRight.Add(new[] { fWn, fVNeck });
                if (k < n - 1)
                    aRight.Add(new[] { fHalfW(fVCrest + fP),
                                       fVCrest + 0.88f * fP });
            }
            float fWBot = aRight[aRight.Count - 1][0];
            aRight.Add(new[] { 0.85f * fWBot, fDepth });   // chamfered bottom

            var aOut = new float[2 * aRight.Count][];
            for (int i = 0; i < aRight.Count; i++)
                aOut[i] = aRight[i];
            for (int i = 0; i < aRight.Count; i++)
            {
                float[] pt = aRight[aRight.Count - 1 - i];
                aOut[aRight.Count + i] = new[] { -pt[0], pt[1] };
            }
            return aOut;
        }

        /// <summary>Axial length of the root and the dr/dz SLOPE of its
        /// platform, for one row. Port of
        /// geometry_generator._root_axial_and_slope: the platform is not a
        /// cylinder, it follows the hub line, which climbs stage by stage.
        /// Root, slot and disc rim share this same linear law, so the
        /// broaching fit is the only gap between them.</summary>
        public static void RootAxialAndSlope(AxialCompressorParameters p,
                                             RowParams row,
                                             out float fAxialMm,
                                             out float fSlope)
        {
            fAxialMm = p.FirTree.AxialFrac * (row.ZTeMm - row.ZLeMm);
            float fDz = 0.5f * MathF.Max(fAxialMm, 1e-6f);
            float fZc = row.ZCenterMm;
            fSlope = (RotorDrum.fRHubAt(p.HubLine, fZc + fDz)
                      - RotorDrum.fRHubAt(p.HubLine, fZc - fDz))
                     / MathF.Max(2f * fDz, 1e-6f);
        }

        /// <summary>
        /// One fir-tree solid (a root, or the tool that broaches its slot)
        /// in machine coordinates, at angular position fPhi0.
        ///
        /// Port of _cq_firtree_solid: the profile is WRAPPED on the
        /// cylinder of each radius — exactly like the blade sections — and
        /// lofted between nSections faces displaced tangentially by
        /// ±(L/2)·tan(gamma), which reproduces the real INCLINED BROACHING
        /// (the slot follows the blade chord, not the axis).
        /// </summary>
        public static Mesh mshSolid(float[][] afProfileUV, float fRPlatform,
                                    float fZCenter, float fAxialLen,
                                    float fStaggerDeg, float fRSlope,
                                    float fPhi0)
        {
            var msh = new Mesh();
            AddSolid(msh, afProfileUV, fRPlatform, fZCenter, fAxialLen,
                     fStaggerDeg, fRSlope, fPhi0);
            return msh;
        }

        /// <summary>Adds one fir-tree solid to an existing mesh (so a whole
        /// ring of slots becomes ONE voxel field instead of N).</summary>
        public static void AddSolid(Mesh msh, float[][] afProfileUV,
                                    float fRPlatform, float fZCenter,
                                    float fAxialLen, float fStaggerDeg,
                                    float fRSlope, float fPhi0)
        {
            int nPts = afProfileUV.Length;
            float fDz = 0.5f * fAxialLen;
            float fTanG = MathF.Tan(fStaggerDeg * MathF.PI / 180f);
            var an = new int[nSections, nPts];

            for (int i = 0; i < nSections; i++)
            {
                float s = -1f + 2f * i / (nSections - 1f);
                float fRp = fRPlatform + s * fDz * fRSlope;
                float fDu = s * fDz * fTanG;
                float fZ = fZCenter + s * fDz;
                for (int j = 0; j < nPts; j++)
                {
                    float u = afProfileUV[j][0];
                    float v = afProfileUV[j][1];
                    float fR = MathF.Max(fRp - v, 1e-3f);
                    float fPhi = fPhi0 + (u + fDu) / fR;
                    an[i, j] = msh.nAddVertex(new Vector3(
                        fR * MathF.Cos(fPhi), fR * MathF.Sin(fPhi), fZ));
                }
            }

            // ruled side walls around the closed contour
            for (int i = 0; i < nSections - 1; i++)
                for (int j = 0; j < nPts; j++)
                {
                    int j1 = (j + 1) % nPts;
                    msh.AddQuad(an[i, j], an[i + 1, j],
                                an[i + 1, j1], an[i, j1]);
                }

            // End caps. The contour is built as right + mirrored(right),
            // so index k on the right chain pairs with nPts-1-k on the
            // left one at the SAME depth v: zipping the two halves closes
            // the face without inventing a centroid vertex.
            AddCapZip(msh, an, 0, nPts);
            AddCapZip(msh, an, nSections - 1, nPts);
        }

        static void AddCapZip(Mesh msh, int[,] an, int iSec, int nPts)
        {
            int nHalf = nPts / 2;
            for (int k = 0; k + 1 < nHalf; k++)
                msh.AddQuad(an[iSec, k], an[iSec, k + 1],
                            an[iSec, nPts - 2 - k], an[iSec, nPts - 1 - k]);
        }

        /// <summary>Every root of one row as a SINGLE voxel field. One
        /// voxelization instead of BladeCount of them: building them one
        /// by one and unioning cost ~160 s per row at voxel 2 mm.
        /// </summary>
        public static Voxels voxRootRing(AxialCompressorParameters p,
                                         RowParams row, float fRPlatform)
        {
            RootAxialAndSlope(p, row, out float fAx, out float fSlope);
            float[][] prof = afProfile(p.FirTree.WTopMm, p.FirTree.DepthMm,
                                       p.FirTree.NTeeth, p.FirTree.TaperDeg,
                                       p.FirTree.LobeFrac, 0f);
            var msh = new Mesh();
            for (int k = 0; k < row.BladeCount; k++)
                AddSolid(msh, prof, fRPlatform, row.ZCenterMm, fAx,
                         row.Sections[0].StaggerDeg, fSlope,
                         2f * MathF.PI * k / row.BladeCount);
            return new Voxels(msh);
        }

        /// <summary>Every broached slot of one disc as a SINGLE voxel
        /// field: one boolean subtraction instead of BladeCount of them.
        /// </summary>
        public static Voxels voxSlotRing(AxialCompressorParameters p,
                                         RowParams row, float fRPlatform)
        {
            RootAxialAndSlope(p, row, out float fAx, out float fSlope);
            float[][] prof = afProfile(p.FirTree.WTopMm, p.FirTree.DepthMm,
                                       p.FirTree.NTeeth, p.FirTree.TaperDeg,
                                       p.FirTree.LobeFrac,
                                       p.FirTree.ClearanceMm);
            var msh = new Mesh();
            for (int k = 0; k < row.BladeCount; k++)
                AddSolid(msh, prof, fRPlatform, row.ZCenterMm,
                         fAx * fSlotAxialFactor,
                         row.Sections[0].StaggerDeg, fSlope,
                         2f * MathF.PI * k / row.BladeCount);
            return new Voxels(msh);
        }
    }
}
