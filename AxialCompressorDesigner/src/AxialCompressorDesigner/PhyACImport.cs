//
// PhyACImport.cs — bridge from the Quasar Phy-AC design pipeline.
//
// Phy-AC (physics_core.py / neural_optimizer.py / geometry_generator.py)
// optimizes a 13-D design vector and writes the verified result as
// `axial_compressor.json` (schema phyac-axial-1) inside <outdir>/geometry/.
// This class maps that file onto AxialCompressorParameters so the
// optimized design can be built as printable rotor + casing STLs:
//
//   var p = PhyACImport.FromContractJson(@"runs\x\geometry\axial_compressor.json");
//   AxialCompressorBuilder.Build(p, outputDir);
//
// The STL is built on the VERIFIED geometry: annulus polylines, per-row
// camber mid-surfaces at 7 span stations (thickened by voxMeshShell —
// the Phy-CC v3 vane pattern), and the structural drum bore.
//

using System;
using System.IO;
using System.Text.Json;

namespace AxialCompressorDesigner
{
    public static class PhyACImport
    {
        public static AxialCompressorParameters FromContractJson(
            string jsonPath, string outputName = null)
        {
            using JsonDocument oDoc = JsonDocument.Parse(File.ReadAllText(jsonPath));
            JsonElement oRoot = oDoc.RootElement;

            if (!oRoot.TryGetProperty("schema", out JsonElement eSchema)
                || !eSchema.GetString().StartsWith("phyac-axial"))
                throw new ArgumentException(
                    $"'{jsonPath}' is not a Phy-AC axial_compressor.json");

            var p = new AxialCompressorParameters();

            JsonElement oAnn = oRoot.GetProperty("annulus");
            p.HubLine = afPolyline(oAnn.GetProperty("hub"));
            p.TipLine = afPolyline(oAnn.GetProperty("tip"));
            if (oAnn.TryGetProperty("tip_clearance_mm", out JsonElement eClr))
                p.TipClearanceMm = (float)eClr.GetDouble();

            if (oRoot.TryGetProperty("structural", out JsonElement oStr)
                && oStr.ValueKind == JsonValueKind.Object
                && oStr.TryGetProperty("drum_inner_r_mm", out JsonElement eDr))
                p.DrumInnerRadiusMm = (float)eDr.GetDouble();

            if (oRoot.TryGetProperty("assembly", out JsonElement oAsm)
                && oAsm.ValueKind == JsonValueKind.Object)
            {
                p.ShaftStubMm = fOpt(oAsm, "shaft_stub_mm", p.ShaftStubMm);
                p.ShaftBoreFrac = fOpt(oAsm, "shaft_bore_frac", p.ShaftBoreFrac);
                p.HubWallMm = fOpt(oAsm, "hub_wall_mm", p.HubWallMm);
                p.DiskWebFrac = fOpt(oAsm, "disk_web_frac", p.DiskWebFrac);
                p.DiskWebMinMm = fOpt(oAsm, "disk_web_min_mm", p.DiskWebMinMm);
                p.DiskWebMaxMm = fOpt(oAsm, "disk_web_max_mm", p.DiskWebMaxMm);
                p.FlangeBoltCount = (int)fOpt(oAsm, "flange_bolt_count",
                                              p.FlangeBoltCount);
                p.FlangeBoltDMm = fOpt(oAsm, "flange_bolt_d_mm", p.FlangeBoltDMm);
                p.FlangeWMm = fOpt(oAsm, "flange_w_mm", p.FlangeWMm);
                p.FlangeTMm = fOpt(oAsm, "flange_t_mm", p.FlangeTMm);
                p.BleedStage = (int)fOpt(oAsm, "bleed_stage", p.BleedStage);
                p.BleedHoleDMm = fOpt(oAsm, "bleed_hole_d_mm", p.BleedHoleDMm);
                p.BleedBossDMm = fOpt(oAsm, "bleed_boss_d_mm", p.BleedBossDMm);
                p.BleedBossHMm = fOpt(oAsm, "bleed_boss_h_mm", p.BleedBossHMm);
            }

            foreach (JsonElement oStage in oRoot.GetProperty("stages").EnumerateArray())
            {
                int iStage = oStage.GetProperty("index").GetInt32();
                foreach (string strKind in new[] { "rotor", "stator" })
                    p.Rows.Add(rowFromJson(oStage.GetProperty(strKind),
                                           iStage));
            }

            // Guide-vane rows (contracts >= 2026-07-16; optional for
            // backward compatibility with older axial_compressor.json)
            if (oRoot.TryGetProperty("igv", out JsonElement oIgv)
                && oIgv.ValueKind == JsonValueKind.Object)
                p.IgvRow = rowFromJson(oIgv, -1);
            if (oRoot.TryGetProperty("ogv", out JsonElement oOgv)
                && oOgv.ValueKind == JsonValueKind.Object)
                p.OgvRow = rowFromJson(oOgv, p.Rows.Count / 2);

            string strName = outputName;
            if (strName == null)
            {
                int nStages = oRoot.GetProperty("derived")
                                   .GetProperty("n_stages").GetInt32();
                strName = $"PhyAC_{nStages}stg";
            }
            p.OutputName = strName;
            return p;
        }

        static RowParams rowFromJson(JsonElement oRow, int iStage)
        {
            var row = new RowParams
            {
                StageIndex = iStage,
                Rotating = oRow.GetProperty("rotating").GetBoolean(),
                BladeCount = oRow.GetProperty("n_blades").GetInt32(),
                ZCenterMm = (float)oRow.GetProperty("z_center_mm").GetDouble(),
                ZLeMm = (float)oRow.GetProperty("z_le_mm").GetDouble(),
                ZTeMm = (float)oRow.GetProperty("z_te_mm").GetDouble(),
            };
            foreach (JsonElement oSec in oRow.GetProperty("sections").EnumerateArray())
            {
                var sec = new SectionParams
                {
                    RMm = (float)oSec.GetProperty("r_mm").GetDouble(),
                    ChordMm = (float)oSec.GetProperty("chord_mm").GetDouble(),
                    StaggerDeg = (float)oSec.GetProperty("stagger_deg").GetDouble(),
                    ThicknessMm = (float)oSec.GetProperty("thickness_mm").GetDouble(),
                    CamberPoints = afPolyline(oSec.GetProperty("camber_points")),
                };
                // full closed profile (real thickness distribution) —
                // enables the solid loft; optional for legacy contracts
                if (oSec.TryGetProperty("points", out JsonElement ePts)
                    && ePts.ValueKind == JsonValueKind.Array)
                    sec.ProfilePoints = afPolyline(ePts);
                row.Sections.Add(sec);
            }
            return row;
        }

        static float fOpt(JsonElement o, string strName, float fDefault)
        {
            return o.TryGetProperty(strName, out JsonElement e)
                ? (float)e.GetDouble() : fDefault;
        }

        static float[][] afPolyline(JsonElement e)
        {
            var af = new float[e.GetArrayLength()][];
            int i = 0;
            foreach (JsonElement pt in e.EnumerateArray())
            {
                af[i] = new float[2];
                int j = 0;
                foreach (JsonElement v in pt.EnumerateArray())
                {
                    if (j < 2) af[i][j] = (float)v.GetDouble();
                    j++;
                }
                i++;
            }
            return af;
        }
    }
}
