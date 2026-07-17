//
// Program.cs — Phy-AC bridge: axial_compressor.json → rotor + casing STLs.
//
//   Example.dll <ruta\a\axial_compressor.json> [outDir] [voxelMm]
//
// Byte-compatible with the subprocess call in phyac_cli.generate_stls.
//

using System;
using System.Globalization;
using System.IO;
using AxialCompressorDesigner;

namespace AxialCompressorDesigner.Example
{
    internal static class Program
    {
        public static int Main(string[] args)
        {
            if (args.Length == 0)
            {
                Console.Error.WriteLine(
                    "Uso: Example.dll <axial_compressor.json> [outDir] " +
                    "[voxelMm] [filtroPieza]");
                Console.Error.WriteLine(
                    "  filtroPieza: subcadena, p.ej. 'RotorStage3', 'Stator',");
                Console.Error.WriteLine(
                    "  'Shaft' — construye solo esas piezas (sin uniones).");
                return 1;
            }

            string strJson = Path.GetFullPath(args[0]);
            if (!File.Exists(strJson))
            {
                Console.Error.WriteLine($"No existe: {strJson}");
                return 1;
            }
            string strOutDir = args.Length > 1
                ? Path.GetFullPath(args[1])
                : Path.Combine(Environment.CurrentDirectory, "output");

            var p = PhyACImport.FromContractJson(strJson);

            if (args.Length > 2 && float.TryParse(args[2], NumberStyles.Float,
                    CultureInfo.InvariantCulture, out float fVoxel))
                p.VoxelSizeMm = fVoxel;

            string strFilter = args.Length > 3 ? args[3] : null;

            Console.WriteLine($"[PhyAC] {strJson}");
            int nExtra = (p.IgvRow != null ? 1 : 0) + (p.OgvRow != null ? 1 : 0);
            Console.WriteLine($"        {p.Rows.Count} filas" +
                (nExtra > 0 ? $" + {nExtra} guide vanes (IGV/OGV)" : "") +
                ", hub " +
                $"{p.HubLine[0][1]:0.#}→{p.HubLine[^1][1]:0.#} mm, tip " +
                $"{p.TipLine[0][1]:0.#}→{p.TipLine[^1][1]:0.#} mm, " +
                $"holgura {p.TipClearanceMm:0.##} mm, voxel {p.VoxelSizeMm} mm" +
                (strFilter != null ? $", filtro '{strFilter}'" : ""));

            string strStl = AxialCompressorBuilder.Build(
                p, strOutDir, showViewer: false, partFilter: strFilter);
            Console.WriteLine($"        → {strStl}");
            return 0;
        }
    }
}
