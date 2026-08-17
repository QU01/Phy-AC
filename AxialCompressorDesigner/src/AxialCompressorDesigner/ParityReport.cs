//
// ParityReport.cs — machine-readable dump of what layer 5c actually
// builds, so the STL route and the CadQuery/STEP route can be compared
// by a test instead of by eye.
//
// Two levels, because they cost very different amounts:
//
//   ProfileJson  pure math (fir-tree profile, root length and platform
//                slope per row). No PicoGK, no GPU, milliseconds. This is
//                what catches the two ports drifting apart.
//   Write        builds every part and records volume and bounding box.
//                Needs the PicoGK runtime and takes minutes at a fine
//                voxel; it is what catches the two routes describing
//                different machines.
//
// Consumed by validation/parity_stl_step.py.
//

using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;
using PicoGK;

namespace AxialCompressorDesigner
{
    public static class ParityReport
    {
        static string F(float f)
            => f.ToString("R", CultureInfo.InvariantCulture);

        /// <summary>Fir-tree profile and per-row root geometry as JSON.
        /// Pure math — safe to call without Library.Go.</summary>
        public static string strProfileJson(AxialCompressorParameters p)
        {
            var sb = new StringBuilder();
            sb.Append("{\n");
            sb.Append("  \"schema_major\": ")
              .Append(PhyACImport.nSCHEMA_MAJOR).Append(",\n");
            sb.Append("  \"firtree_enabled\": ")
              .Append(p.FirTree.Enabled ? "true" : "false").Append(",\n");
            sb.Append("  \"n_sections\": ")
              .Append(FirTree.nSections).Append(",\n");
            sb.Append("  \"slot_axial_factor\": ")
              .Append(F(FirTree.fSlotAxialFactor)).Append(",\n");
            sb.Append("  \"tie_bolt_count\": ")
              .Append(p.TieBoltCount).Append(",\n");
            sb.Append("  \"tie_bolt_d_mm\": ")
              .Append(F(p.TieBoltDMm)).Append(",\n");
            sb.Append("  \"rim_relief_mm\": ")
              .Append(F(p.RimReliefMm)).Append(",\n");

            AppendProfile(sb, "root_profile",
                          FirTree.afProfile(p.FirTree.WTopMm,
                                            p.FirTree.DepthMm,
                                            p.FirTree.NTeeth,
                                            p.FirTree.TaperDeg,
                                            p.FirTree.LobeFrac, 0f));
            sb.Append(",\n");
            AppendProfile(sb, "slot_profile",
                          FirTree.afProfile(p.FirTree.WTopMm,
                                            p.FirTree.DepthMm,
                                            p.FirTree.NTeeth,
                                            p.FirTree.TaperDeg,
                                            p.FirTree.LobeFrac,
                                            p.FirTree.ClearanceMm));
            sb.Append(",\n  \"rows\": [\n");
            bool bFirst = true;
            foreach (RowParams row in p.Rows)
            {
                if (!row.Rotating)
                    continue;
                FirTree.RootAxialAndSlope(p, row, out float fAx,
                                          out float fSlope);
                if (!bFirst)
                    sb.Append(",\n");
                bFirst = false;
                sb.Append("    {\"stage\": ").Append(row.StageIndex)
                  .Append(", \"z_center_mm\": ").Append(F(row.ZCenterMm))
                  .Append(", \"r_hub_c_mm\": ")
                  .Append(F(RotorDrum.fRHubAt(p.HubLine, row.ZCenterMm)))
                  .Append(", \"root_axial_mm\": ").Append(F(fAx))
                  .Append(", \"platform_slope\": ").Append(F(fSlope))
                  .Append(", \"n_blades\": ").Append(row.BladeCount)
                  .Append("}");
            }
            sb.Append("\n  ]\n}\n");
            return sb.ToString();
        }

        static void AppendProfile(StringBuilder sb, string strName,
                                  float[][] af)
        {
            sb.Append("  \"").Append(strName).Append("\": [");
            for (int i = 0; i < af.Length; i++)
            {
                if (i > 0)
                    sb.Append(", ");
                sb.Append('[').Append(F(af[i][0])).Append(", ")
                  .Append(F(af[i][1])).Append(']');
            }
            sb.Append(']');
        }

        /// <summary>Builds every part and writes volume + bounding box per
        /// part. Must run INSIDE Library.Go (it creates Voxels).</summary>
        public static void WritePartsJson(AxialCompressorParameters p,
                                          string strPath)
        {
            var oComp = new AxialCompressor(p);
            var sb = new StringBuilder();
            sb.Append("{\n  \"output_name\": \"").Append(p.OutputName)
              .Append("\",\n  \"voxel_mm\": ").Append(F(p.VoxelSizeMm))
              .Append(",\n  \"parts\": {\n");
            bool bFirst = true;
            foreach ((string strName, Func<Voxels> fnBuild)
                     in oComp.PartBuilders())
            {
                Voxels vox = fnBuild();
                vox.CalculateProperties(out float fVol, out BBox3 oBB);
                // r_max of a solid of revolution is the bbox half-extent;
                // r_min is NOT recoverable from a bbox (a bore is inside),
                // so only the outer envelope and the axial span travel.
                float fRMax = MathF.Max(
                    MathF.Max(MathF.Abs(oBB.vecMin.X), MathF.Abs(oBB.vecMax.X)),
                    MathF.Max(MathF.Abs(oBB.vecMin.Y), MathF.Abs(oBB.vecMax.Y)));
                if (!bFirst)
                    sb.Append(",\n");
                bFirst = false;
                sb.Append("    \"").Append(strName).Append("\": {")
                  .Append("\"volume_mm3\": ").Append(F(fVol))
                  .Append(", \"r_max_mm\": ").Append(F(fRMax))
                  .Append(", \"z_min_mm\": ").Append(F(oBB.vecMin.Z))
                  .Append(", \"z_max_mm\": ").Append(F(oBB.vecMax.Z))
                  .Append('}');
                Library.Log($"Parity {strName}: {fVol:F0} mm^3, "
                            + $"r_max {fRMax:F2} mm");
            }
            sb.Append("\n  }\n}\n");
            File.WriteAllText(strPath, sb.ToString());
            Library.Log($"Parity report → {strPath}");
        }
    }
}
