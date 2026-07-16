//
// AxialCompressor.cs — multi-part assembly (real turbojet construction):
//
//   Shaft                    drive shaft with stubs and center bore
//   RotorStage{i}            bladed disc (blisk): disc web + hub-shell
//                            segment (incl. spacer under its stator) +
//                            rotor blades of stage i
//   StatorRing{i}            casing ring: shell segment + stator vanes;
//                            first/last rings carry the bolted flanges;
//                            the BleedStage ring carries the bleed port
//   Rotor / Casing           unions of the above (running position)
//   Assembly                 everything
//
// Split planes sit mid-gap between one stage's stator TE and the next
// stage's rotor LE — where real bolted discs and casing rings join.
//

using System;
using System.Collections.Generic;
using PicoGK;

namespace AxialCompressorDesigner
{
    public sealed class AxialCompressor
    {
        readonly AxialCompressorParameters m_p;
        readonly float m_fMinThick;
        readonly List<RowParams> m_aRotors = new List<RowParams>();
        readonly List<RowParams> m_aStators = new List<RowParams>();

        public AxialCompressor(AxialCompressorParameters p)
        {
            m_p = p;
            // keep blade walls >= 2 voxels so they do not vanish
            m_fMinThick = 2f * p.VoxelSizeMm;
            foreach (RowParams row in p.Rows)
                (row.Rotating ? m_aRotors : m_aStators).Add(row);
            m_aRotors.Sort((a, b) => a.ZCenterMm.CompareTo(b.ZCenterMm));
            m_aStators.Sort((a, b) => a.ZCenterMm.CompareTo(b.ZCenterMm));
        }

        public int nStages => m_aRotors.Count;

        /// <summary>Split plane between stage i and i+1 (0-based i).</summary>
        float fSplitAfter(int i)
        {
            float fTe = m_aStators[i].ZTeMm;
            float fNextLe = i + 1 < m_aRotors.Count
                ? m_aRotors[i + 1].ZLeMm : fTe + 4f;
            return 0.5f * (fTe + fNextLe);
        }

        // ── rotating parts ────────────────────────────────────────────
        public Voxels voxShaft() => RotorDrum.voxShaft(m_p);

        public Voxels voxRotorStage(int i)
        {
            float fZ0 = i == 0 ? m_p.HubLine[0][0] : fSplitAfter(i - 1);
            float fZ1 = i == nStages - 1
                ? m_p.HubLine[m_p.HubLine.Length - 1][0] : fSplitAfter(i);
            Voxels vox = RotorDrum.voxHubShellSegment(m_p, fZ0, fZ1);
            vox += RotorDrum.voxDiscWeb(m_p, m_aRotors[i]);
            vox += BladeRow.voxBuildRow(m_aRotors[i], m_p.TipClearanceMm,
                                        m_p.RootSinkMm, m_fMinThick);
            return vox;
        }

        // ── static parts ──────────────────────────────────────────────
        public Voxels voxStatorRing(int i)
        {
            float fZC0 = m_p.TipLine[0][0] - Casing.fLipMm;
            float fZC1 = m_p.TipLine[m_p.TipLine.Length - 1][0]
                       + Casing.fLipMm;
            float fZ0 = i == 0 ? fZC0 : fSplitAfter(i - 1);
            float fZ1 = i == nStages - 1 ? fZC1 : fSplitAfter(i);

            Voxels vox = Casing.voxShellSegment(m_p, fZ0, fZ1);
            if (i == 0)
                vox += Casing.voxFlange(m_p, fZC0, bForward: true);
            if (i == nStages - 1)
                vox += Casing.voxFlange(m_p, fZC1, bForward: false);
            if (i + 1 == Math.Clamp(m_p.BleedStage, 1, nStages))
            {
                float fZB = m_aStators[i].ZCenterMm;
                vox += Casing.voxBleedBoss(m_p, fZB);
                vox -= Casing.voxBleedHole(m_p, fZB);
            }
            vox += BladeRow.voxBuildRow(m_aStators[i], m_p.TipClearanceMm,
                                        m_p.RootSinkMm, m_fMinThick);
            return vox;
        }

        // ── lazy part builders: the caller builds ONE part at a time,
        //    exports it and releases the field (streaming keeps peak
        //    memory at one part instead of eleven and gives per-part
        //    progress in the log) ─────────────────────────────────────
        public List<(string strName, Func<Voxels> fnBuild)> PartBuilders()
        {
            var parts = new List<(string, Func<Voxels>)>
            {
                ("Shaft", voxShaft),
            };
            for (int i = 0; i < nStages; i++)
            {
                int n = i;   // capture by value
                parts.Add(($"RotorStage{n + 1}", () => voxRotorStage(n)));
            }
            for (int i = 0; i < nStages; i++)
            {
                int n = i;
                parts.Add(($"StatorRing{n + 1}", () => voxStatorRing(n)));
            }
            return parts;
        }
    }
}
