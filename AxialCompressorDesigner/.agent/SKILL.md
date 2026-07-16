---
name: picogk-geometry-authoring
description: |
  Author voxel-based 3D geometries in C# using the PicoGK voxel geometry kernel
  (and its companion LEAP 71 ShapeKernel). Use this skill when the user asks to
  build, generate, or author an engineering geometry (heat exchanger, compressor
  impeller, rocket nozzle, manifold, lattice structure, fin pack, screw thread
  support, etc.) as C# code targeting the PicoGK `Library.Go()` task. Covers the
  full RISC pattern: `Library.Go` → `Task()` → `Lattice`/`BaseShape`/implicit →
  `Voxels` boolean pipeline → `Mesh` STL export. Triggers on phrases like "build
  a [part] with PicoGK", "voxel geometry for...", "shape kernel [part]",
  "compresor centrífugo en PicoGK", "genera la geometría de un intercambiador",
  "computational engineering model". Do NOT use this skill for: FreeCAD, OpenSCAD,
  OpenVDB raw C++ authoring, traditional B-Rep CAD, or generic CAD questions
  unrelated to the PicoGK/ShapeKernel stack.
---

# PicoGK Geometry Authoring

A Computational Engineering Model (CEM) in PicoGK is a small C# program that
combines a handful of primitives, runs them through boolean voxel operations,
and exports an STL. The whole architecture is deliberately RISC: you have very
few primitive instructions, and you compose everything from them.

## Inputs to collect

Before writing code, confirm:

1. **Target geometry** — what physical part (e.g. centrifugal compressor, heat
   exchanger, impeller, nozzle, manifold). If unclear, ask the user.
2. **Bounding box / envelope** — overall dimensions in mm. If the user gives
   performance specs (flow rate, pressure ratio, RPM), translate them into
   geometric envelopes yourself and confirm.
3. **Voxel size** — default 0.5 mm for parts in the 10–100 cm range. Lower to
   0.3 mm or 0.2 mm for fine internal features (fins < 0.5 mm thick); bump up
   to 1 mm for outer-shell iteration drafts.
4. **Output folder** — an existing absolute path where the viewer can write
   screenshots and the STL.
5. **ShapeKernel available?** — PicoGK alone is enough, but `BaseBox` /
   `BaseCylinder` / `BasePipe` / `BaseSphere` / `BaseRing` / `BaseLens` /
   `LatticeManifold` / `ScrewHole` / `ThreadReinforcement` / `ThreadCutter` from
   the LEAP 71 ShapeKernel cut hundreds of lines. Assume ShapeKernel is present
   unless the user says otherwise.

## Procedure

### 1. Skeleton

Every PicoGK program is a `static void Task()` invoked by `Library.Go`:

```csharp
using System.Numerics;
using Leap71.ShapeKernel;     // optional but recommended
using PicoGK;

public static class MyPart
{
    public static void Task()
    {
        PicoGK.Library.Log("Starting Task.");
        var oPart = new MyPart();
        Voxels voxResult = oPart.voxConstruct();
        Sh.ExportVoxelsToSTLFile(voxResult, Sh.strGetExportPath(Sh.EExport.STL, "MyPart"));
        PicoGK.Library.Log("Finished Task.");
    }
}
```

`Program.cs` calls it:

```csharp
try
{
    PicoGK.Library.Go(0.5f, MyPart.Task, strOutputFolder);
}
catch (Exception e) { Console.WriteLine(e); }
```

Always wrap in `try/catch` so library load errors surface in the log.

### 2. The RISC instruction set (what to reach for)

You only need ~10 primitive calls. Internalize this table:

| Intent                              | PicoGK call                                     |
|-------------------------------------|--------------------------------------------------|
| Spawn empty voxel field             | `new Voxels(...)` (via `Library`)                |
| Sphere / box / beam voxel field     | `Voxels.voxSphere(lib, center, r)`, `Utils.mshCreateCube(lib, bbox)` |
| Boolean union / sub / intersect     | `voxA + voxB`, `voxA - voxB`, `voxA & voxB`      |
| In-place boolean                    | `voxA.BoolAdd(voxB)`, `voxA.BoolSubtract(voxB)` |
| Offset / shell / fillet / smoothen  | `vox.Offset(d)`, `vox.voxShell(neg, pos)`, `vox.Fillet(r)`, `vox.Smoothen(d)` |
| Slice / project                     | `vox.ProjectZSlice(zStart, zEnd)`                |
| Implicit field (gyroid, signed fn)  | `vox.RenderImplicit(IImplicit, bbox)`           |
| Lattice (beams + spheres)           | `new Lattice(); lat.AddBeam(a, ra, b, rb)`       |
| Export STL                          | `msh = new Mesh(vox); msh.SaveToStlFile(path)`   |
| Preview / screenshot                | `Sh.PreviewVoxels(vox, color)`, `Library.oViewer().RequestScreenShot(...)` |
| Frame / coord system                | `new LocalFrame(pos, xDir, zDir)`                |

For everything else (screw threads, fins, support lattice, base primitives),
use the ShapeKernel helpers. See `references/api-quickref.md`.

### 3. The canonical pattern (use this skeleton for any complex part)

This is exactly the HelixHeatX pattern. Adopt it for compressors, pumps, etc.:

```csharp
public partial class MyPart
{
    // 1. INPUTS (LocalFrames, dimensions, material thicknesses)
    LocalFrame m_oInlet, m_oOutlet;
    float m_fWallThickness;

    // 2. CONSTRUCTOR — define inputs + bounding box
    public MyPart() { /* set fields, build m_voxBounding */ }

    // 3. ONE METHOD PER SUB-COMPONENT (each returns Voxels)
    Voxels voxGetFlowPath()        { /* lattice of beam spirals */ }
    Voxels voxGetFins()            { /* rooftop fins for printability */ }
    Voxels voxGetShell()           { /* outer shell with offset + fillet */ }
    Voxels voxGetInletOutlet()     { /* pipe transitions */ }
    Voxels voxGetSupports()        { /* lattice supports for overhangs */ }

    // 4. ASSEMBLY — combine via boolean ops
    Voxels voxConstruct()
    {
        Voxels voxFluid  = voxGetFlowPath();
        Voxels voxFins   = voxGetFins();
        Voxels voxShell  = voxGetShell();
        voxShell         += voxGetInletOutlet();
        voxShell         += voxGetSupports();
        Voxels voxResult = voxShell - voxFluid;   // shell MINUS fluid void
        voxResult        += voxFins;              // add internal fins back
        voxResult        &= m_voxBounding;        // trim to envelope
        return voxResult;
    }
}
```

**Why this pattern works**: it matches the "primitive operations, free
combination" idea of RISC. Each sub-component is a small function returning
`Voxels`, so you can comment out a line during development and re-run quickly
at coarse voxel size to iterate on geometry. Final assembly is a few lines
that read like a sentence.

### 4. The "shell minus fluid" inversion (critical mental model)

In traditional CAD you build walls. In PicoGK you build **the void first**, then
derive the shell:

```csharp
Voxels voxInnerVolume  = voxGetFlowPath();        // the fluid void
Voxels voxOuterVolume  = voxInnerVolume.voxOffset(0.9f);  // shell
Voxels voxResult       = voxOuterVolume - voxInnerVolume;   // wall only
```

This is counterintuitive but it's the heat-exchanger pattern verbatim, and
it scales to compressors (impeller + volute), pumps, nozzles — anywhere you
have a flowing fluid inside a solid part. See
`references/heat-exchanger-pattern.md` for the full annotated analysis.

### 5. Lattice for sweeps, splines, and supports

`Lattice` is the workhorse for anything curved, swept, or fine-structured
(fins, support webs, screw threads, manifolds):

```csharp
Lattice lat = new Lattice();
for (float t = 0; t < 1f; t += 0.001f)
{
    Vector3 a = SampleCurve(t);
    Vector3 b = a + 0.4f * NormalAt(t);
    lat.AddBeam(a, fBeamA, b, fBeamB);
}
Voxels vox = new Voxels(lat);   // voxelize the lattice
```

`AddBeam(a, rA, b, rB)` puts a tapered capsule from `a` to `b`. Add
`AddSphere` for end caps. Voxelize with `new Voxels(lat)`.

**Tip**: For pipes, manifolds, and splitters, use `Sh.latFromBeam(...)` or
build a `Lattice` and sweep `AddBeam` along sample points of a curve. See
`references/api-quickref.md` for the full lattice API.

### 6. Implicits for lattices and TPMS

For gyroids, Schwarz-P, and other triply-periodic minimal surfaces, use
`IImplicit`:

```csharp
public class Gyroid : IBoundedImplicit
{
    public BBox3 oBounds => new(...);
    public float fSignedDistance(in Vector3 vec)
        => 0.5f * MathF.Sin(vec.X) * MathF.Cos(vec.Y)
         + 0.5f * MathF.Sin(vec.Y) * MathF.Cos(vec.Z)
         + 0.5f * MathF.Sin(vec.Z) * MathF.Cos(vec.X);
}
// usage:
Voxels vox = new Voxels(lib, new Gyroid(), bbox);
```

Use `IntersectImplicit` to mask a field with an implicit (gyroid inside a box).

### 7. Voxel size rules (fast iteration + final print)

The HelixHeatX tutorial gives a clean rule: **smallest feature ≥ 1 voxel**.

| Part size          | Start voxel | Final voxel   |
|--------------------|-------------|----------------|
| 1 m tall           | 1.0 mm      | 0.5 mm         |
| 10–100 cm          | 0.5 mm      | 0.3 mm         |
| < 10 cm (small)    | 0.2 mm      | 0.1 mm         |

Iteration loop: build outer shell at 1 mm → swap to 0.3 mm only when
ready to print. Hot path: write a flag (`bool bRunOuterOnly`) so you can
comment out interior components during dev.

### 8. Printability guards

If the part will be 3D-printed, every PicoGK part should encode:

- **Rooftop fins**: for any fin thinner than ~0.4 mm, taper the height
  (e.g. use `fBeam` that drops to 0.2 mm at the tip). See
  `InternalFins.cs` in HelixHeatX for the canonical pattern.
- **Bottom web / powder removal**: at the build-plate contact surface, add
  a lattice web (see `PrintWeb.cs`) and subtract it for grooves.
- **Supports for overhangs**: use `Lattice` beams with kink angle > 45°
  to support pipe stubs and IO threads (see `IOSupports.cs`).
- **Fillet + smoothen** the outer shell at the end:
  `voxShell.Fillet(5f); voxShell.Smoothen(0.5f);` (the HelixHeatX final
  step, gives a polished look).

## Output contract

A working PicoGK project is a Visual Studio solution with:

```
MyPart.sln
Program.cs                       # calls Library.Go(voxelSize, Task, logPath)
PicoGK/                          # library source
LEAP71_ShapeKernel/              # optional but recommended
src/MyPart/
    MyPart.cs                    # main class with Task() + voxConstruct()
    FlowPath.cs                  # one sub-component per file (partial class)
    Fins.cs
    Shell.cs
    IO.cs
    Supports.cs
```

Code style follows the LEAP 71 conventions:

- `m_` prefix for member fields, `f` for floats, `o` for objects, `a` for arrays.
- `m_voxBounding` for the envelope, `m_fWallThickness` for global parameters.
- `LocalFrame` for all coordinate systems.
- `Voxels` is the only return type of sub-component methods.
- `Sh.PreviewVoxels(...)` + `Library.oViewer().RequestScreenShot(...)` between
  sub-components so the user sees the build sequence.

## Failure handling

Common failure modes and how to fix them:

- **`PicoGKLights.zip` not found** — viewer comes up dark, but geometry still
  builds. Add `assets/ViewerEnvironment.zip` next to the executable or pass
  `strLightsFile` to `Library.Go`.
- **Native lib not found** — exception on startup. Make sure
  `PicoGKLib.dll` / `PicoGKLib.dylib` is on the executable path. On Windows
  it lives under `PicoGK/native/win-x64/`.
- **Features vanish at fine voxel** — you're under the smallest-feature
  threshold. Either thicken the feature (`m_fWallThickness += 0.2f`) or
  drop voxel size to 0.2 mm.
- **Fins/thread walls print as solid blocks** — voxel size is too coarse
  for the wall thickness. Lower voxel size or thicken the wall.
- **Stl export is empty** — the final `Voxels` field is empty. Check that
  `voxResult.bIsEmpty()` is false before exporting. Likely cause: subtracted
  the void AND the shell in the wrong order.
- **Slice count / memory blow-up** — voxel field too large at 0.1 mm. Back
  off to 0.3 mm or trim the bounding box.

## Examples

### Canonical: a centrifugal compressor

The user asks: *"build the geometry of a small centrifugal compressor
(60 mm inlet OD, 80 mm impeller OD, 50 mm axial length)"*. Full annotated
recipe in `references/centrifugal-compressor-pattern.md`. Minimal skeleton:

```csharp
public partial class CentrifugalCompressor
{
    LocalFrame m_oInletFrame, m_oImpellerFrame, m_oVoluteFrame;
    float m_fImpellerOD, m_fInletOD, m_fAxialLength, m_fWallThickness;
    Voxels m_voxBounding;

    public CentrifugalCompressor()
    {
        m_fInletOD       = 60f;
        m_fImpellerOD    = 80f;
        m_fAxialLength   = 50f;
        m_fWallThickness = 0.8f;
        m_oInletFrame    = new LocalFrame(new Vector3(-25, 0, 0), Vector3.UnitX);
        m_oImpellerFrame = new LocalFrame(new Vector3(0, 0, 0), Vector3.UnitZ);
        m_oVoluteFrame   = new LocalFrame(new Vector3(0, 0, 0), Vector3.UnitZ);
        m_voxBounding    = new BaseBox(new LocalFrame(), m_fImpellerOD + 20, m_fImpellerOD + 20, m_fAxialLength + 10)
                              .voxConstruct();
    }

    public Voxels voxConstruct()
    {
        Voxels voxFlowPath = voxGetImpellerFlowPath();      // curved blades
        Voxels voxShroud   = voxGetShroud();                // top + bottom cover
        Voxels voxInlet    = voxGetInletPipe();             // axial pipe
        Voxels voxVolute   = voxGetVolute();                // spiral collector
        Voxels voxHub      = voxGetHub();                   // central hub

        Voxels voxOuter    = voxShroud + voxInlet + voxVolute + voxHub;
        voxOuter           += voxGetOuterSupports();

        Voxels voxResult   = voxOuter - voxFlowPath;        // shell minus void
        voxResult          += voxGetInternalBlades();       // add blades back
        voxResult          &= m_voxBounding;

        voxResult.Fillet(3f);
        voxResult.Smoothen(0.5f);
        return voxResult;
    }
}
```

This is the heat-exchanger pattern inverted (volute instead of helix) plus
a shroud/hub/volute assembly. See the reference for the full impeller-blade
lattice, volute spiral, and printability details.

### Second example: minimal Hello World (PicoGK_Examples style)

```csharp
public class HelloWorld
{
    public static void Task()
    {
        Mesh msh = Utils.mshCreateCube();      // 1×1×1 mm cube at origin
        Library.oViewer().Add(msh);
    }
}
// Program.cs:
Library.Go(0.5f, HelloWorld.Task);
```

### Anti-example: do not do this

```csharp
// BAD: building walls directly, then adding features one-by-one
Voxels vox = voxOuterBox;
vox -= voxInnerCylinder;       // void 1
vox -= voxSecondVoid;          // void 2 — already lost track of what's solid
vox.Fillet(5f);                // fills holes
// 200 lines later: no idea what's solid anymore
```

Always design the **void first**, then take `outer.Offset(wall) - inner`.

## Companion references

- `references/api-quickref.md` — PicoGK + ShapeKernel API cheat-sheet
- `references/heat-exchanger-pattern.md` — annotated walk-through of the
  HelixHeatX example (the canonical "shell minus fluid" pattern)
- `references/centrifugal-compressor-pattern.md` — full recipe for a
  centrifugal compressor (impeller blades, volute, shroud, supports)
