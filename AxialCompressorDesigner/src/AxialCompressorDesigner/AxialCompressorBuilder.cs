//
// AxialCompressorBuilder.cs — public entry point (Phy-CC Build pattern).
//
//   1. Validates the parameters
//   2. Calls PicoGK.Library.Go(voxelSize, task, logFile, ...)
//   3. The task builds the geometry and exports the STLs
//   4. Returns the absolute path of the assembly STL
//

using System;
using System.IO;
using PicoGK;

namespace AxialCompressorDesigner
{
    public static class AxialCompressorBuilder
    {
        /// <summary>
        /// Build every part STL (shaft, per-stage bladed discs, casing
        /// rings) plus the union views into <paramref name="outputDir"/>.
        /// Returns the assembly STL path. <paramref name="partFilter"/>
        /// (optional, case-insensitive substring, e.g. "RotorStage3" or
        /// "Stator") builds only the matching parts and skips the unions —
        /// the fast path for inspecting a single piece.
        /// </summary>
        public static string Build(AxialCompressorParameters p,
                                   string outputDir,
                                   bool showViewer = false,
                                   string partFilter = null)
        {
            string strErr = p.Validate();
            if (strErr != null)
                throw new ArgumentException("Invalid parameters: " + strErr);

            Directory.CreateDirectory(outputDir);
            string strOutputStl = Path.GetFullPath(
                Path.Combine(outputDir, p.OutputName + ".stl"));

            string strLogFile = Path.Combine(outputDir, "PicoGK.log");
            PicoGK.Library.Go(
                p.VoxelSizeMm,
                () => BuildTask(p, strOutputStl, showViewer, partFilter),
                strLogFile,
                !showViewer,
                $"AxialCompressorDesigner — {p.OutputName}");

            if (partFilter == null && !File.Exists(strOutputStl))
                throw new Exception(
                    "Build failed: STL was not written. Check the log at " +
                    strLogFile);
            return strOutputStl;
        }

        /// <summary>Builds every part and writes a machine-readable
        /// parity report (volume + bounding box per part) instead of
        /// STLs. Used by validation/parity_stl_step.py to check that the
        /// PicoGK route and the CadQuery route describe the SAME machine.
        /// </summary>
        public static void BuildParityReport(AxialCompressorParameters p,
                                             string outputDir,
                                             string strReportPath)
        {
            string strErr = p.Validate();
            if (strErr != null)
                throw new ArgumentException("Invalid parameters: " + strErr);
            Directory.CreateDirectory(outputDir);
            Directory.CreateDirectory(
                Path.GetDirectoryName(Path.GetFullPath(strReportPath)));
            PicoGK.Library.Go(
                p.VoxelSizeMm,
                () => ParityReport.WritePartsJson(p, strReportPath),
                Path.Combine(outputDir, "PicoGK.log"),
                true,
                $"AxialCompressorDesigner — parity {p.OutputName}");
        }

        static void BuildTask(AxialCompressorParameters p,
                              string strOutputStl, bool showViewer,
                              string strFilter)
        {
            Library.Log($"Building axial compressor '{p.OutputName}' at " +
                        $"voxel {p.VoxelSizeMm} mm ({p.Rows.Count} rows)" +
                        (strFilter != null ? $" — filter '{strFilter}'" : ""));

            var oComp = new AxialCompressor(p);
            string strDir = Path.GetDirectoryName(strOutputStl) ?? ".";
            string strBase = Path.GetFileNameWithoutExtension(strOutputStl);

            // ── streaming: una pieza a la vez (construir → exportar →
            //    acumular unión → liberar). Mantiene el pico de memoria
            //    en UNA pieza y reporta progreso/tiempos por pieza. ────
            Voxels voxRotor = null;    // unión de piezas rotativas
            Voxels voxCasing = null;   // unión de piezas estáticas
            int nPart = 0;
            var oWatch = System.Diagnostics.Stopwatch.StartNew();
            foreach ((string strName, Func<Voxels> fnBuild)
                     in oComp.PartBuilders())
            {
                if (strFilter != null && !strName.Contains(
                        strFilter, StringComparison.OrdinalIgnoreCase))
                    continue;

                long lMs0 = oWatch.ElapsedMilliseconds;
                Voxels vox = fnBuild();
                if (vox.bIsEmpty())
                    throw new Exception($"part '{strName}' is empty — " +
                                        "check parameters.");
                vox.CalculateProperties(out float fVol, out BBox3 oBB);
                string strPartStl = Path.Combine(
                    strDir, $"{strBase}_{strName}.stl");
                new Mesh(vox).SaveToStlFile(strPartStl);
                Library.Log($"Part {strName}: {fVol:F0} mm^3 in " +
                            $"{(oWatch.ElapsedMilliseconds - lMs0) / 1000f:F1} s " +
                            $"→ {strPartStl}");

                if (strFilter == null)
                {
                    bool bRotating = strName == "Shaft"
                                  || strName.StartsWith("RotorStage");
                    if (bRotating)
                        voxRotor = voxRotor == null ? vox : voxRotor + vox;
                    else
                        voxCasing = voxCasing == null ? vox : voxCasing + vox;
                }
                if (showViewer)
                    Library.oViewer().Add(vox, nPart++ % 4);
            }

            if (strFilter != null)
            {
                Library.Log("Part filter active — union views skipped.");
                return;
            }

            // ── uniones en posición de marcha + ensamble ──────────────
            string strRotorStl = Path.Combine(strDir, strBase + "_Rotor.stl");
            string strCasingStl = Path.Combine(strDir, strBase + "_Casing.stl");
            new Mesh(voxRotor).SaveToStlFile(strRotorStl);
            Library.Log($"Exported STL: {strRotorStl}");
            new Mesh(voxCasing).SaveToStlFile(strCasingStl);
            Library.Log($"Exported STL: {strCasingStl}");

            Voxels voxAssembly = voxRotor + voxCasing;
            new Mesh(voxAssembly).SaveToStlFile(strOutputStl);
            Library.Log($"Exported STL: {strOutputStl}");
        }
    }
}
