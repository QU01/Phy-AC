//
// Parameters.cs — data-driven parameters for the axial compressor builder.
//
// Unlike Phy-CC's CentrifugalCompressorParameters (a handful of scalars),
// the axial builder is DATA-DRIVEN: the Python layer 5a emits fully
// resolved blade-row sections (camber polylines per span station) inside
// axial_compressor.json, and this record simply carries them. All lengths
// are in millimeters; angles in degrees.
//
// Sign convention of the contract: rotor staggers/angles are positive,
// stator ones negative (they turn the flow in opposite senses around the
// axis); the builder applies the SIGNED rotation uniformly.
//

using System.Collections.Generic;

namespace AxialCompressorDesigner
{
    /// <summary>One blade section at a span station.</summary>
    public sealed class SectionParams
    {
        public float RMm;                 // stacking radius of this section
        public float ChordMm;
        public float StaggerDeg;          // signed (rotor +, stator −)
        public float ThicknessMm;         // t_max — voxMeshShell thickness
        public float[][] CamberPoints;    // chord frame, mm, centroid at 0,0
        public float[][] ProfilePoints;   // closed CCW profile (contract
                                          // `points`) — enables the solid
                                          // loft with the REAL NACA-65/DCA
                                          // thickness distribution; null in
                                          // legacy contracts (shell fallback)
    }

    /// <summary>Fir-tree blade retention, mirroring the contract's
    /// `assembly.firtree` block one to one. Defaults match
    /// geometry_generator.build_contract so a legacy contract without the
    /// block still builds something sensible — but `Enabled` is false
    /// there, and the builder then falls back to the old sunk root.</summary>
    public sealed class FirTreeParams
    {
        public bool Enabled = false;
        public int NTeeth = 3;
        public float WTopMm = 6f;         // total width at the platform
        public float DepthMm = 10f;       // radial depth of the root
        public float TaperDeg = 8f;       // half-angle of the tree
        public float LobeFrac = 0.22f;    // neck depth / local half-width
        public float ClearanceMm = 0.08f; // broaching/assembly fit per side
        public float AxialFrac = 0.85f;   // root length / row axial chord
        public float MinWidthMm = 3f;
    }

    /// <summary>One blade row (rotor or stator) with its span sections.</summary>
    public sealed class RowParams
    {
        public int StageIndex;
        public bool Rotating;
        public int BladeCount;
        public float ZCenterMm;           // axial position of the stacking axis
        public float ZLeMm;               // axial slot of the row
        public float ZTeMm;
        public List<SectionParams> Sections = new List<SectionParams>();

        /// <summary>Mean blade thickness over the span [mm].</summary>
        public float MeanThicknessMm
        {
            get
            {
                float f = 0f;
                foreach (SectionParams s in Sections)
                    f += s.ThicknessMm;
                return Sections.Count > 0 ? f / Sections.Count : 1f;
            }
        }
    }

    /// <summary>Full machine parameters consumed by AxialCompressor.</summary>
    public sealed class AxialCompressorParameters
    {
        public string OutputName = "PhyAC";
        public float VoxelSizeMm = 0.5f;

        // Annulus polylines [z_mm, r_mm], inlet → exit
        public float[][] HubLine;
        public float[][] TipLine;
        public float TipClearanceMm = 0.5f;

        public float DrumInnerRadiusMm = 20f;   // shaft radius under the disks
        public float CasingWallMm = 5f;         // casing shell thickness
        public float RootSinkMm = 1.5f;         // blade root sink into body.
                                                // Legacy: only used when the
                                                // contract carries no fir-tree
                                                // (FirTree.Enabled == false).
                                                // With retention the root is a
                                                // real fir-tree in a broached
                                                // slot, not a sunk rib.
        public float BladeFilletRMm = 2f;       // root fillet radius (0 = off;
                                                // same value structures_core
                                                // uses for the Peterson K_t)

        // Assembly detail (turbodesigner-style; contract `assembly` block)
        public float ShaftStubMm = 30f;         // shaft stub past each end
        public float ShaftBoreFrac = 0.35f;     // center bore / shaft radius
        public float HubWallMm = 5f;            // hub flow-path shell wall
        public float DiskWebFrac = 0.5f;        // disk web / row axial chord
        public float DiskWebMinMm = 4f;
        public float DiskWebMaxMm = 14f;
        public int FlangeBoltCount = 12;
        public float FlangeBoltDMm = 6f;
        public float FlangeWMm = 12f;           // flange radial width
        public float FlangeTMm = 8f;            // flange axial thickness
        public int BleedStage = 1;              // 1-based stage of bleed port
        public float BleedHoleDMm = 18f;
        public float BleedBossDMm = 32f;
        public float BleedBossHMm = 14f;
        public int TieBoltCount = 12;           // tie bolts through the discs
        public float TieBoltDMm = 10f;
        public float RimReliefMm = 0.10f;       // disc rim sits this much
                                                // below the blade platform —
                                                // the platform forms the gas
                                                // path, the rim does not

        /// <summary>Fir-tree blade retention (contract `assembly.firtree`).
        /// The SAME parametric profile the CadQuery path uses; see
        /// FirTree.cs and geometry_generator.firtree_profile.</summary>
        public FirTreeParams FirTree = new FirTreeParams();

        public List<RowParams> Rows = new List<RowParams>();

        // Optional guide-vane rows (casing-mounted, non-rotating). The
        // physics assumes stage-1 pre-swirl (implicit IGV) and an OGV that
        // returns the exit flow to axial; contracts >= 2026-07-16 carry
        // both rows so the printed machine matches the verified triangles.
        public RowParams IgvRow;
        public RowParams OgvRow;

        /// <summary>Returns an error message, or null when valid.</summary>
        public string Validate()
        {
            if (HubLine == null || HubLine.Length < 2)
                return "HubLine must have at least 2 points";
            if (TipLine == null || TipLine.Length < 2)
                return "TipLine must have at least 2 points";
            if (Rows.Count == 0)
                return "no blade rows";
            if (VoxelSizeMm <= 0.05f)
                return "VoxelSizeMm too small";
            foreach (RowParams row in Rows)
            {
                string strErr = strValidateRow(row);
                if (strErr != null)
                    return strErr;
            }
            if (IgvRow != null)
            {
                string strErr = strValidateRow(IgvRow);
                if (strErr != null)
                    return "IGV: " + strErr;
            }
            if (OgvRow != null)
            {
                string strErr = strValidateRow(OgvRow);
                if (strErr != null)
                    return "OGV: " + strErr;
            }
            return null;
        }

        static string strValidateRow(RowParams row)
        {
            if (row.BladeCount < 3)
                return $"row stage {row.StageIndex}: BladeCount < 3";
            if (row.Sections.Count < 2)
                return $"row stage {row.StageIndex}: needs >= 2 sections";
            int n = row.Sections[0].CamberPoints.Length;
            foreach (SectionParams s in row.Sections)
                if (s.CamberPoints == null || s.CamberPoints.Length != n)
                    return $"row stage {row.StageIndex}: camber point " +
                           "count differs between sections (loft needs " +
                           "equal counts)";
            int nProf = -1;
            foreach (SectionParams s in row.Sections)
            {
                if (s.ProfilePoints == null)
                    continue;
                if (nProf < 0)
                    nProf = s.ProfilePoints.Length;
                else if (s.ProfilePoints.Length != nProf)
                    return $"row stage {row.StageIndex}: profile point " +
                           "count differs between sections (solid loft " +
                           "needs equal counts)";
            }
            return null;
        }
    }
}
