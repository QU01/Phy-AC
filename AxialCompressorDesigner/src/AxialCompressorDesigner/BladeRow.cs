//
// BladeRow.cs — axial blade rows, two recipes:
//
// 1. SOLID PROFILE LOFT (default since 2026-07-16). The contract carries
//    the full closed profile per section (`points`: NACA-65 / DCA
//    thickness distribution, 60 CCW points, analytically guaranteed free
//    of self-intersection by blade_profiles.py). Each blade is a
//    watertight mesh: ruled side walls between densified ribs + hub/tip
//    caps triangulated LE→TE (zipper between pressure and suction
//    chains), voxelized directly. This preserves the real thickness
//    distribution — LE radius, max-thickness position and asymmetry —
//    which the camber-shell recipe destroyed (uniform thickness).
//
//    Note on the Phy-CC v3 warning ("never stitch pressure+suction"):
//    that failure mode is meshes built by OFFSETTING a camber sheet at
//    mesh level, which fold over where the camber radius < ½·thickness.
//    Here the closed contour comes analytically from Python with the
//    no-self-intersection contract invariant, so the solid loft is safe.
//
// 2. CAMBER SHEET + voxMeshShell (fallback — Phy-CC v3 recipe). Used for
//    legacy contracts without `points` and for rows thinner than two
//    voxels, where the solid path cannot guarantee printable walls; the
//    shell clamps the thickness (and warns — geometry thicker than
//    designed).
//
// Geometry: each section arrives in its chord frame (x chordwise, y
// normal, centroid at the origin, mm) plus a SIGNED stagger (rotor +,
// stator −). For each span rib the polyline is rotated by the stagger,
// then wrapped onto the cylinder of its stacking radius:
//
//   phi = phi0 + y'/r      (arc length preserved)
//   P   = (r·cos phi, r·sin phi, zCenter + x')
//
// Rotor rows root at the hub (root rib sunk into the drum, free end
// pulled in by the running clearance); stator rows hang from the casing
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

        /// <summary>Radial over-extension of the raw loft before it is
        /// trimmed against the flow path. Must comfortably exceed the
        /// clearances and the annulus contraction along the chord; the
        /// CadQuery path uses the same 4 mm.</summary>
        const float fTrimExtendMm = 4f;

        /// <summary>
        /// Builds all blades of one row as a single voxel field, TRIMMED
        /// against the annulus.
        ///
        /// Until phase 10 each row was simply laid on the radii its own
        /// sections carried. That is wrong twice over: the annulus also
        /// contracts ALONG THE CHORD, so a tip that fits at the row centre
        /// pokes out of the casing at the trailing edge; and it left the
        /// running clearances to be faked by shifting the end sections.
        /// Now the blade is lofted with radial margin and cut against the
        /// real hub/tip polylines — the same recipe as the CadQuery route
        /// (geometry_generator._cq_blade_trimmed), so both describe the
        /// same blade.
        ///
        /// Clearances follow the contract convention: the ROTOR is born at
        /// the hub and ends at the tip line (the running gap lives between
        /// it and the casing bore, which sits at tip + eps); the STATOR
        /// hangs from the casing bore and leaves eps over the drum.
        /// </summary>
        /// <param name="p">Machine parameters (annulus, clearance, fir-tree).</param>
        /// <param name="row">Row parameters (sections hub→tip).</param>
        /// <param name="fMinThicknessMm">Lower bound on printable wall
        /// thickness (2 voxels).</param>
        public static Voxels voxBuildRow(AxialCompressorParameters p,
                                         RowParams row,
                                         float fMinThicknessMm)
        {
            bool bSolid = true;
            foreach (SectionParams s in row.Sections)
                if (s.ProfilePoints == null || s.ProfilePoints.Length < 8)
                {
                    bSolid = false;
                    break;
                }
            if (bSolid && row.MeanThicknessMm < fMinThicknessMm)
            {
                Library.Log(
                    $"WARNING row stage {row.StageIndex}: mean thickness " +
                    $"{row.MeanThicknessMm:F2} mm < 2 voxels " +
                    $"({fMinThicknessMm:F2} mm) — solid profile would not " +
                    "survive voxelization; falling back to the camber " +
                    "shell.");
                bSolid = false;
            }

            Voxels vox = bSolid
                ? voxBuildRowSolid(row, fMinThicknessMm)
                : voxBuildRowShell(row, fMinThicknessMm);

            vox &= voxFlowPath(p, row, fMinThicknessMm);

            // Fir-tree roots hang BELOW the hub line, so they are added
            // after the trim (the trim would cut them off). The disc
            // broaches the matching slot with the same profile plus the
            // assembly fit — see FirTree.cs.
            if (row.Rotating && p.FirTree.Enabled)
                vox += FirTree.voxRootRing(
                    p, row, RotorDrum.fRHubAt(p.HubLine, row.ZCenterMm));
            return vox;
        }

        /// <summary>The solid of revolution the row is cut against: the
        /// flow path between the hub and tip polylines over the row's
        /// axial slot, with the running clearances applied.</summary>
        static Voxels voxFlowPath(AxialCompressorParameters p, RowParams row,
                                  float fMinThicknessMm)
        {
            // The voxelized surface still rounds by about half a voxel, so
            // the free end is pulled in by that much: without it a rotor
            // tip trimmed exactly at the tip line welds itself to the
            // casing bore whenever eps <= the voxel size.
            float fEdge = 0.25f * fMinThicknessMm;      // = ½ voxel
            float fEps = p.TipClearanceMm;
            if (row.Rotating && fEps < 2f * fEdge)
                Library.Log(
                    $"WARNING row stage {row.StageIndex}: tip clearance " +
                    $"{fEps:F2} mm is below one voxel " +
                    $"({2f * fEdge:F2} mm) — the printed gap is set by the " +
                    "voxel, not by the design. Use a finer voxel for a " +
                    "faithful running clearance.");

            float fRHubOff = row.Rotating ? 0f : fEps + fEdge;
            float fRTipOff = row.Rotating ? -fEdge : fEps;
            float fZ0 = row.ZLeMm - 2f;
            float fZ1 = row.ZTeMm + 2f;

            Voxels vox = LatticeUtils.voxRevolveZ(t =>
            {
                float fZ = fZ0 + t * (fZ1 - fZ0);
                return new Vector2(
                    Casing.fRTipAt(p.TipLine, fZ) + fRTipOff, fZ);
            }, 0f, 64);
            return vox - LatticeUtils.voxRevolveZ(t =>
            {
                float fZ = (fZ0 - 2f) + t * ((fZ1 + 2f) - (fZ0 - 2f));
                return new Vector2(
                    RotorDrum.fRHubAt(p.HubLine, fZ) + fRHubOff, fZ);
            }, 0f, 64);
        }

        // ══════════════════════════════════════════════════════════════
        // Recipe 1 — solid profile loft (real thickness distribution)
        // ══════════════════════════════════════════════════════════════
        static Voxels voxBuildRowSolid(RowParams row, float fMinThicknessMm)
        {
            int nSecs = row.Sections.Count;
            var afR = new float[nSecs];
            for (int i = 0; i < nSecs; i++)
                afR[i] = row.Sections[i].RMm;
            // Radial margin at BOTH ends; the flow-path trim decides where
            // the blade actually stops (voxBuildRow).
            afR[0] -= fTrimExtendMm;
            afR[nSecs - 1] += fTrimExtendMm;

            Rib[] aRibs = DensifyRibs(row, afR, bProfile: true);
            Mesh msh = new Mesh();
            for (int nBlade = 0; nBlade < row.BladeCount; nBlade++)
            {
                float fPhi0 = 2f * MathF.PI * nBlade / row.BladeCount;
                AddBladeSolid(msh, aRibs, row.ZCenterMm, fPhi0);
            }
            return new Voxels(msh);
        }

        /// <summary>Closed watertight blade: ruled side walls around the
        /// closed contour + hub/tip caps (LE→TE zipper triangulation,
        /// independent of where the polygon starts or its winding).</summary>
        static void AddBladeSolid(Mesh msh, Rib[] aRibs, float fZCenter,
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

            // side walls: closed ruled tube between consecutive ribs
            for (int i = 0; i < nRibs - 1; i++)
                for (int j = 0; j < nPts; j++)
                {
                    int j1 = (j + 1) % nPts;
                    msh.AddQuad(an[i, j], an[i + 1, j],
                                an[i + 1, j1], an[i, j1]);
                }

            // caps on the root and free-end ribs close the solid
            AddCap(msh, aRibs[0], an, 0, nPts);
            AddCap(msh, aRibs[nRibs - 1], an, nRibs - 1, nPts);
        }

        /// <summary>Triangulates the airfoil polygon of one rib. The two
        /// surfaces are x-monotone by construction (blade_profiles.py), so
        /// a zipper between the LE→TE chains is always a valid, crossing-
        /// free triangulation — no assumption on start index or winding
        /// (the stator sections arrive mirrored and re-ordered).</summary>
        static void AddCap(Mesh msh, Rib rib, int[,] an, int iRib, int nPts)
        {
            // LE/TE = extreme points along the chord (profile frame x)
            int iLe = 0, iTe = 0;
            for (int j = 1; j < nPts; j++)
            {
                if (rib.afPts[j][0] < rib.afPts[iLe][0]) iLe = j;
                if (rib.afPts[j][0] > rib.afPts[iTe][0]) iTe = j;
            }

            // chain A: iLe → iTe forward (cyclic); chain B: the other way
            int nA = (iTe - iLe + nPts) % nPts + 1;
            int nB = (iLe - iTe + nPts) % nPts + 1;
            var aA = new int[nA];
            var aB = new int[nB];
            for (int k = 0; k < nA; k++)
                aA[k] = (iLe + k) % nPts;
            for (int k = 0; k < nB; k++)
                aB[k] = (iLe - k + nPts) % nPts;   // reversed → also LE→TE

            int a = 0, b = 0;
            while (a < nA - 1 || b < nB - 1)
            {
                bool bAdvanceA =
                    b == nB - 1
                    || (a < nA - 1
                        && rib.afPts[aA[a + 1]][0] <= rib.afPts[aB[b + 1]][0]);
                if (bAdvanceA)
                {
                    AddTri(msh, an[iRib, aA[a]], an[iRib, aA[a + 1]],
                           an[iRib, aB[b]]);
                    a++;
                }
                else
                {
                    AddTri(msh, an[iRib, aA[a]], an[iRib, aB[b + 1]],
                           an[iRib, aB[b]]);
                    b++;
                }
            }
        }

        static void AddTri(Mesh msh, int n0, int n1, int n2)
        {
            if (n0 == n1 || n1 == n2 || n0 == n2)
                return;   // degenerate at the shared LE/TE vertices
            msh.nAddTriangle(n0, n1, n2);
        }

        // ══════════════════════════════════════════════════════════════
        // Recipe 2 — camber sheet + voxMeshShell (Phy-CC v3 fallback)
        // ══════════════════════════════════════════════════════════════
        static Voxels voxBuildRowShell(RowParams row, float fMinThicknessMm)
        {
            float fThick = MathF.Max(row.MeanThicknessMm, fMinThicknessMm);
            if (row.MeanThicknessMm < fMinThicknessMm)
                Library.Log(
                    $"WARNING row stage {row.StageIndex}: shell thickness " +
                    $"clamped {row.MeanThicknessMm:F2} → {fThick:F2} mm " +
                    "(2 voxels) — the printed blade is THICKER than the " +
                    "verified design; use a finer voxel for faithful parts.");
            float fHalf = 0.5f * fThick;

            int nSecs = row.Sections.Count;

            // Radial margin at both ends; the flow-path trim in
            // voxBuildRow decides where the blade stops. The shell also
            // inflates by fHalf, which the trim removes as well.
            var afR = new float[nSecs];
            for (int i = 0; i < nSecs; i++)
                afR[i] = row.Sections[i].RMm;
            afR[0] -= fTrimExtendMm + fHalf;
            afR[nSecs - 1] += fTrimExtendMm + fHalf;

            // densify the span: interpolated ribs (r, stagger, camber)
            var aRibs = DensifyRibs(row, afR, bProfile: false);

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

        static Rib[] DensifyRibs(RowParams row, float[] afR, bool bProfile)
        {
            int nSecs = row.Sections.Count;
            float[][] afSrc(SectionParams s)
                => bProfile ? s.ProfilePoints : s.CamberPoints;
            int nPts = afSrc(row.Sections[0]).Length;
            int nRibs = (nSecs - 1) * nSpanSubdiv + 1;
            var aRibs = new Rib[nRibs];
            int k = 0;
            for (int i = 0; i < nSecs - 1; i++)
            {
                SectionParams a = row.Sections[i];
                SectionParams b = row.Sections[i + 1];
                float[][] afA = afSrc(a), afB = afSrc(b);
                for (int s = 0; s < nSpanSubdiv; s++)
                {
                    float t = s / (float)nSpanSubdiv;
                    var pts = new float[nPts][];
                    for (int j = 0; j < nPts; j++)
                        pts[j] = new float[]
                        {
                            afA[j][0] * (1 - t) + afB[j][0] * t,
                            afA[j][1] * (1 - t) + afB[j][1] * t,
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
                afPts = afSrc(last),
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
