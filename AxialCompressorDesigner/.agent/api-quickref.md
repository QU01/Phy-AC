# PicoGK + ShapeKernel API Quick Reference

A condensed lookup of the calls you actually use. Method signatures abbreviated;
consult the source for full docs.

## PicoGK.Library (the global façade)

Created once via `Library.Go(voxelSizeMM, TaskDelegate, [logPath], [title])`.

```csharp
// entry point
Library.Go(float fVoxelSizeMM, ThreadStart fnTask,
           string strLogFilePath = "",
           bool   bEndAppWithTask = false,
           string strWindowTitle  = "PicoGK");

// runtime accessors (only valid after Go())
Library.oViewer()         // Viewer instance
Library.oLibrary()        // Library instance
Library.Log(string, args) // log to file + console
Library.strLogFolder      // current output dir
Library.fVoxelSizeMM      // active voxel size
Library.bContinueTask()   // returns false if user closed the viewer
Library.EndTask()         // request task end
```

## Viewer

```csharp
Library.oViewer().Add(Voxels vox)               // default material
Library.oViewer().Add(Voxels vox, int group)    // group 0..N
Library.oViewer().Add(Mesh msh)
Library.oViewer().RemoveAllObjects()
Library.oViewer().SetGroupMaterial(int group, string hexColor, float roughness, float metalness)
Library.oViewer().RequestScreenShot(string pathWithoutExt)
```

## Voxels (the core boolean field)

Construction:

```csharp
new Voxels(Library lib)                          // empty
new Voxels(Voxels source)                        // copy
new Voxels(ScalarField field)                    // from signed distance
new Voxels(Library lib, IImplicit x, BBox3 box)  // render implicit in bbox
new Voxels(Library lib, IBoundedImplicit x)      // implicit with own bounds
new Voxels(Mesh msh)                             // voxelize a mesh
new Voxels(Lattice lat)                          // voxelize a lattice

// static factories
Voxels.voxSphere(Library lib, Vector3 center, float radius)
Voxels.voxLatticeBeam(Library lib, Vector3 a, float ra, Vector3 b, float rb)
Voxels.voxMeshShell(Library lib, Mesh msh, float radius)
Voxels.voxCombine(Voxels a, Voxels b)
Voxels.voxCombineAll(Library lib, IEnumerable<Voxels> list)
```

Boolean (operators and in-place):

```csharp
voxA + voxB          // union
voxA - voxB          // difference (A minus B)
voxA & voxB          // intersect
voxA.BoolAdd(voxB)
voxA.BoolSubtract(voxB)
voxA.BoolIntersect(voxB)
voxA.BoolAddAll(IEnumerable<Voxels>)
voxA.BoolSubtractAll(IEnumerable<Voxels>)
```

Offset / shell / fillet / smooth:

```csharp
vox.Offset(float dMM)        // in-place; positive=expand, negative=shrink
vox.voxOffset(float dMM)     // returns new copy
vox.DoubleOffset(d1, d2)     // expand d1 then contract d2
vox.TripleOffset(d)          // smooth filter: shrink, grow 2x, shrink back
vox.Smoothen(d)              // alias for TripleOffset
vox.OverOffset(dFirst, dFinalSurface=0)  // round corners / fillet
vox.Fillet(d)                // alias: OverOffset(d, 0)
vox.voxShell(fOffset)                       // hollw from outside
vox.voxShell(fNegOffset, fPosOffset, smoothInner=0)  // double-walled
```

Trim, project, mask:

```csharp
vox.Trim(BBox3 box)              // intersect with box
vox.ProjectZSlice(float zStart, float zEnd)  // project all slices from zStart..zEnd
vox.RenderImplicit(IImplicit x, BBox3 box)   // overwrite with implicit field
vox.IntersectImplicit(IImplicit x)           // mask: intersect with implicit
vox.RenderMesh(Mesh msh)
vox.RenderLattice(Lattice lat)
```

Queries:

```csharp
vox.bIsEmpty()
vox.bIsInside(Vector3 pt)        // true if inside or on surface
vox.bIsEqual(Voxels other)
vox.oBoundingBox()               // returns BBox3 (via mesh roundtrip)
vox.GetVoxelDimensions(out int x, y, z)
vox.fSurfaceNormal(Vector3 pt)
vox.bClosestPointOnSurface(in Vector3 search, out Vector3 surfacePt)
vox.vecSurfaceNormal(in Vector3 surfacePt)
vox.CalculateProperties(out float volumeCubicMM, out BBox3 bbox)
```

Slicing (for inspection, 2D image export):

```csharp
vox.imgAllocateSlice(out int nSlices, ESliceAxis.X|Y|Z)
vox.GetVoxelSlice(int sliceIdx, ref ImageGrayScale img, ESliceMode mode, ESliceAxis)
```

## Lattice (beams + spheres)

```csharp
new Lattice(Library lib)

lat.AddSphere(Vector3 center, float radius)
lat.AddBeam(Vector3 a, float rA, Vector3 b, float rB, bool bRoundCap = true)
lat.AddBeam(Vector3 a, Vector3 b, float rA, float rB, bool bRoundCap = true)  // args swapped
```

Then convert to voxels with `new Voxels(lat)`. Lattices are efficient for
sweeps, splines, threads, and any feature where you'd otherwise build a
mesh by hand.

## Mesh

```csharp
new Mesh(Library lib)                            // empty
new Mesh(Voxels vox)                             // mesh from voxels
msh.AddVertices(IEnumerable<Vector3>, out int[] indices)
msh.nAddVertex(Vector3 v) -> int idx
msh.nAddTriangle(Vector3 A, B, C) -> int triIdx
msh.nAddTriangle(Triangle t) -> int triIdx
msh.nAddTriangle(int idxA, idxB, idxC)
msh.AddQuad(int n0, n1, n2, n3, bool bFlipped=false)
msh.AddQuad(Vector3 v0, v1, v2, v3, bool bFlipped=false)
msh.nTriangleCount(), msh.nVertexCount()
msh.GetTriangle(int idx, out Vector3 A, B, C)
msh.vecVertexAt(int idx)
msh.oBoundingBox() -> BBox3
msh.mshCreateTransformed(Vector3 scale, Vector3 offset)
msh.mshCreateTransformed(Matrix4x4 mat)
msh.mshCreateMirrored(Vector3 planePoint, Vector3 planeNormal)
msh.Append(Mesh other)
msh.SaveToStlFile(string path, EStlUnit unit=AUTO, Vector3? offset=null, float scale=1.0f)
```

## IO

```csharp
MeshIo.mshFromStlFile(string path, out EStlUnit detectedUnit)
MeshIo.mshFromStlFile(FileStream)
msh.SaveToStlFile(string path)

ImageIo  // 2D image <-> grayscale slice conversions
VdbFile  // load/save OpenVDB files for interop with Houdini etc.
VoxelsIo // save/load the voxel field directly (VDB-backed)
```

## Utils

```csharp
Utils.mshCreateCube(Library lib, BBox3 bbox)
Utils.mshCreateCube(Library lib, Vector3? scale, Vector3? offsetMM)
Utils.strDocumentsFolder()                // cross-platform
Utils.strProjectRootFolder()
Utils.strPicoGKSourceCodeFolder()
Utils.strExecutableFolder()
Utils.strDateTimeFilename(prefix, postfix)
Utils.strShorten(string, int maxChars)
```

## Types

```csharp
Vector3, Vector2, Matrix4x4, BBox3, BBox2, Triangle, Coord, ColorFloat
LocalFrame  // (see ShapeKernel below)
```

## ShapeKernel (LEAP 71, optional but strongly recommended)

### `BaseShape` derivatives — the parametric primitives

Each takes a `LocalFrame` and dimensions, has `voxConstruct() -> Voxels`.

```csharp
new BaseBox(LocalFrame,        float x, y, z)             // mm-aligned box
new BaseCylinder(LocalFrame,   float length, float radius)
new BaseSphere(LocalFrame,     float radius)
new BaseRing(LocalFrame,       float length, float radius)
new BaseLens(LocalFrame,       float radius, float height)
new BasePipe(LocalFrame,       float length)               // pipe with mod radius
new BasePipeSegment(LocalFrame, float radius, float arcDeg, float length)
```

### Pipe modulation

```csharp
oPipe.SetRadius(SurfaceModulation inner, SurfaceModulation outer)
// SurfaceModulation: delegate float (float fPhi, float fLengthRatio)
```

A `SurfaceModulation` is a `(phi, lengthRatio) -> radius` function. Use this
to taper pipe ends, add bulbous sections, etc.

### Construction modules

```csharp
new ScrewHole(LocalFrame, float threadLen, float threadRadius,
              float headLen, float headRadius).voxConstruct()
new ThreadReinforcement(LocalFrame, float length, float innerR, float outerR).voxConstruct()
new ThreadCutter(LocalFrame, float length, float maxR, float coreR, float pitch).voxConstruct()
new LatticeManifold(LocalFrame, float length, float radius).voxConstruct()
```

### Static helpers (`Sh`)

```csharp
Sh.PreviewVoxels(Voxels vox, Color c, float alpha = 1f)
Sh.PreviewLattice(Lattice lat, Color c, float alpha = 1f)
Sh.PreviewBoxWireframe(BaseBox box, Color c)
Sh.PreviewCylinderWireframe(BaseCylinder cyl, Color c)
Sh.ExportVoxelsToSTLFile(Voxels vox, string path)
Sh.strGetExportPath(Sh.EExport.STL, "filename")       // returns full path under log folder
Sh.strGetExportPath(Sh.EExport.TGA, "screenshot")     // screenshot path
Sh.latFromBeam(Vector3 a, Vector3 b, float ra, float rb, bool bRoundCap)
Sh.oGetBoundingBox(Voxels vox) -> BBox3
Sh.ExportMeshToSTLFile(Mesh msh, string path)
```

### Vector / frame ops (`VecOperations`, `LocalFrame`)

```csharp
// LocalFrame: position + local X / Y / Z axes
LocalFrame oFrame = new LocalFrame(Vector3 pos)        // default axes
LocalFrame oFrame = new LocalFrame(Vector3 pos, Vector3 xDir, Vector3 zDir)
LocalFrame.vecGetPosition() -> Vector3
LocalFrame.vecGetLocalX() / Y() / Z() -> Vector3

LocalFrame.oGetTranslatedFrame(LocalFrame source, Vector3 translation)
LocalFrame.oGetInvertFrame(LocalFrame source, bool bFlipX, bool bFlipY)

// VecOperations
VecOperations.vecGetCylPoint(float r, float phi, float z) -> Vector3
VecOperations.vecTranslatePointOntoFrame(LocalFrame frame, Vector3 pt) -> Vector3
VecOperations.vecRotateAroundZ(Vector3 pt, float angle, Vector3 pivot) -> Vector3
VecOperations.vecRotateAroundAxis(Vector3 v, float angle, Vector3 axis) -> Vector3
```

### Uf (utility functions)

```csharp
Uf.fTransFixed(float a, float b, float ratio)            // smooth a→b
Uf.fGetSuperShapeRadius(float phi, Uf.ESuperShape.ROUND|QUAD|...)
Uf.Wait(float seconds)
```

### Cp (color palette)

```csharp
Cp.clrRed, Cp.clrBlue, Cp.clrBlack, Cp.clrGray, Cp.clrRock,
Cp.clrPitaya, Cp.clrFrozen, Cp.clrWarning, Cp.clrToothpaste,
Cp.clrRacingGreen, Cp.clrYellow
```

## Color hex format

`"RRGGBB"` or `"AARRGGBB"`. Examples: `"FF000033"` (red, 20% alpha),
`"555577"` (dark blue-gray).

## Common gotchas

- **Library mismatch**: all `Voxels`/`Lattice`/`Mesh` must come from the same
  `Library` instance. Inside `Library.Go`, that's the global one; use
  `Library.oLibrary()` if you need it explicitly.
- **Voxel size vs feature size**: if a wall is 0.4 mm, voxel must be ≤ 0.2 mm
  or the wall will partially vanish.
- **Mesh → Voxels → Mesh roundtrip loses zero-thickness surfaces**. Always
  use Offset to make a wall have a positive thickness.
- **bIsEmpty check**: call `vox.bIsEmpty()` before exporting to avoid
  empty STLs.
- **BoolAdd vs +**: both work, but `+` returns a new copy; `BoolAdd` is
  in-place. For pipelines of >5 unions, prefer `BoolAdd` to save memory.
