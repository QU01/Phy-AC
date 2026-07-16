//
// LatticeUtils.cs — minimal helpers built on top of PicoGK primitives.
// No ShapeKernel dependency.
//

using System.Numerics;
using PicoGK;

namespace AxialCompressorDesigner
{
    public static class LatticeUtils
    {
        /// <summary>
        /// Solid axis-aligned box (centered at vecCenter) with full
        /// dimensions fX × fY × fZ. Voxelized from a watertight cube
        /// mesh — the result is a true solid, not a shell.
        /// </summary>
        public static Voxels voxBox(Vector3 vecCenter, float fX, float fY, float fZ)
        {
            Mesh msh = Utils.mshCreateCube(new Vector3(fX, fY, fZ), vecCenter);
            return new Voxels(msh);
        }

        /// <summary>
        /// Solid cylinder of length fL along +Z, radius fR, with the
        /// bottom face at vecBaseBottom (so the centre of the base is
        /// at vecBaseBottom + (0, 0, 0), top face at vecBaseBottom +
        /// (0, 0, fL)). A single flat-capped lattice beam along the
        /// axis with the full radius — a true solid, not a shell.
        /// </summary>
        public static Voxels voxCylinderZ(Vector3 vecBaseBottom, float fR, float fL)
        {
            Lattice lat = new Lattice();
            lat.AddBeam(vecBaseBottom, fR,
                        vecBaseBottom + new Vector3(0, 0, fL), fR,
                        false);
            return new Voxels(lat);
        }

        /// <summary>
        /// Solid cylinder along the local +Z of a LocalFrame, length
        /// fL, radius fR, base at the frame origin.
        /// </summary>
        public static Voxels voxCylinderAlongAxis(LocalFrame oFrame, float fR, float fL)
        {
            Vector3 vecBase = oFrame.vecGetPosition();
            Vector3 vecZ    = oFrame.vecGetLocalZ();

            Lattice lat = new Lattice();
            lat.AddBeam(vecBase, fR, vecBase + fL * vecZ, fR, false);
            return new Voxels(lat);
        }

        /// <summary>
        /// Solid of revolution about the Z axis. fnContour maps
        /// t in [0,1] to a (radius, z) point; the solid fills everything
        /// inside the revolved contour (r &lt;= contour radius at each z),
        /// built as a stack of <b>tapered (conical) frustums</b> — each band
        /// uses the true contour radius at each of its two ends, so the
        /// surface follows the contour smoothly instead of terracing into
        /// constant-radius cylinders. fROffset is added to the contour
        /// radius (use it for clearance / wall offsets). The contour must
        /// be monotonic in z.
        /// </summary>
        public static Voxels voxRevolveZ(System.Func<float, Vector2> fnContour,
                                         float fROffset = 0f, int nBands = 64)
        {
            Lattice lat = new Lattice();
            for (int i = 0; i < nBands; i++)
            {
                Vector2 a = fnContour(i / (float)nBands);
                Vector2 b = fnContour((i + 1) / (float)nBands);
                float fRa = a.X + fROffset;
                float fRb = b.X + fROffset;
                if (fRa <= 0.1f && fRb <= 0.1f)
                    continue;
                // Conical band from (a.Y, fRa) to (b.Y, fRb). Consecutive
                // bands share an endpoint (same z, same radius), so the
                // frustums weld into a continuous, terrace-free surface.
                lat.AddBeam(new Vector3(0f, 0f, a.Y), System.MathF.Max(fRa, 0.05f),
                            new Vector3(0f, 0f, b.Y), System.MathF.Max(fRb, 0.05f),
                            false);
            }
            return new Voxels(lat);
        }
    }
}
