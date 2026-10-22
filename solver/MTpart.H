// ************************************************************************* //
// 9. alpha1, epsilon1 -> find keyhole depth / melt depth
//     + deepest keyhole position (z)
// ************************************************************************* //

const vectorField& cellC = mesh.C();

label nIntF = mesh.nInternalFaces();

depth_y *= 0.0;
melt_depth_y *= 0.0;


// ------------------------------------------------------------------------- //
// Variables for deepest keyhole interface position
// ------------------------------------------------------------------------- //

scalar localMaxKeyholeY = -GREAT;
scalar localZatMaxKeyholeY = GREAT;
label localKeyholeInterfaceCount = 0;


// ------------------------------------------------------------------------- //
// Find alpha1 = 0.5 interface (keyhole)
// and epsilon1 = 0.5 interface (melt depth)
// ------------------------------------------------------------------------- //

for (label faceI = 0; faceI < nIntF; ++faceI)
{
    label own = mesh.owner()[faceI];
    label nei = mesh.neighbour()[faceI];

    scalar cOwn = alpha1[own];
    scalar cNei = alpha1[nei];

    scalar eOwn = epsilon1[own];
    scalar eNei = epsilon1[nei];


    // =====================================================================
    // Keyhole interface : alpha1 = 0.5
    // =====================================================================
    if ((cOwn - 0.5) * (cNei - 0.5) < 0)
    {
        scalar yOwn = cellC[own].y();
        scalar yNei = cellC[nei].y();

        scalar zOwn = cellC[own].z();
        scalar zNei = cellC[nei].z();

        scalar alpha_con =
            (0.5 - cOwn) / (cNei - cOwn);

        scalar yInt =
            yOwn + alpha_con * (yNei - yOwn);

        scalar zInt =
            zOwn + alpha_con * (zNei - zOwn);


        // Existing keyhole depth calculation
        if (depth_y[own] < yInt)
        {
            depth_y[own] = yInt;
        }


        // Deepest keyhole interface point on this processor
        if (yInt > localMaxKeyholeY)
        {
            localMaxKeyholeY = yInt;
            localZatMaxKeyholeY = zInt;
        }

        localKeyholeInterfaceCount++;
    }


    // =====================================================================
    // Melt interface : epsilon1 = 0.5
    // =====================================================================
    if ((eOwn - 0.5) * (eNei - 0.5) < 0)
    {
        scalar yOwn = cellC[own].y();
        scalar yNei = cellC[nei].y();

        scalar alpha_con =
            (0.5 - eOwn) / (eNei - eOwn);

        scalar yInt =
            yOwn + alpha_con * (yNei - yOwn);

        if (melt_depth_y[own] < yInt)
        {
            melt_depth_y[own] = yInt;
        }
    }
}


// ------------------------------------------------------------------------- //
// Global depth
// ------------------------------------------------------------------------- //

scalar globalMaxY =
    gMax(depth_y.internalField());

scalar globalMaxMeltY =
    gMax(melt_depth_y.internalField());

scalar y_ini = 1e-4;

Info << "keyhole depth: "
     << globalMaxY - y_ini
     << endl;

Info << "melt depth: "
     << globalMaxMeltY - y_ini
     << endl;


// ------------------------------------------------------------------------- //
// Find global deepest keyhole position
// ------------------------------------------------------------------------- //

label globalKeyholeInterfaceCount =
    localKeyholeInterfaceCount;

reduce
(
    globalKeyholeInterfaceCount,
    sumOp<label>()
);


scalar globalMaxKeyholeY =
    localMaxKeyholeY;

reduce
(
    globalMaxKeyholeY,
    maxOp<scalar>()
);


scalar globalZatMaxKeyholeY = GREAT;

scalar tolerance = 1e-12;


// Only the processor containing the global deepest point
// supplies its z position
if
(
    localKeyholeInterfaceCount > 0
 && mag(localMaxKeyholeY - globalMaxKeyholeY) < tolerance
)
{
    globalZatMaxKeyholeY =
        localZatMaxKeyholeY;
}


// If multiple processors share the same maximum,
// take the minimum z among them.
// Normally this should not matter except for ties.
reduce
(
    globalZatMaxKeyholeY,
    minOp<scalar>()
);


if (globalKeyholeInterfaceCount > 0)
{
    Info << "deepest keyhole y: "
         << globalMaxKeyholeY
         << endl;

    Info << "deepest keyhole z: "
         << globalZatMaxKeyholeY
         << endl;
}
else
{
    Info << "deepest keyhole y: 0"
         << endl;

    Info << "deepest keyhole z: 0"
         << endl;
}



// ************************************************************************* //
// 10. liquid volume fraction -> find width / length / head / tail
//
// x direction : transverse direction -> melt pool width
// z direction : scan direction       -> melt pool length
//
// Melt interface criterion:
//     epsilon1 * alpha1 = 0.5
//
// Assumption:
//     laser scanning direction = +z
//
// Therefore:
//     zMin = rear end (tail end)
//     zMax = front end (head end)
//
// ************************************************************************* //


// ------------------------------------------------------------------------- //
// Local extrema
// ------------------------------------------------------------------------- //

scalar localMinX = GREAT;
scalar localMaxX = -GREAT;

scalar localMinZ = GREAT;
scalar localMaxZ = -GREAT;


// Coordinates at width extrema
scalar localYatMinX = GREAT;
scalar localYatMaxX = GREAT;


// Optional coordinates at longitudinal extrema
scalar localYatMinZ = GREAT;
scalar localYatMaxZ = GREAT;

scalar localXatMinZ = GREAT;
scalar localXatMaxZ = GREAT;


label localMeltInterfaceCount = 0;


// ------------------------------------------------------------------------- //
// Find liquid melt interface
// ------------------------------------------------------------------------- //

for (label faceI = 0; faceI < nIntF; ++faceI)
{
    label own = mesh.owner()[faceI];
    label nei = mesh.neighbour()[faceI];

    scalar lOwn =
        epsilon1[own] * alpha1[own];

    scalar lNei =
        epsilon1[nei] * alpha1[nei];


    if ((lOwn - 0.5) * (lNei - 0.5) < 0)
    {
        scalar xOwn = cellC[own].x();
        scalar xNei = cellC[nei].x();

        scalar yOwn = cellC[own].y();
        scalar yNei = cellC[nei].y();

        scalar zOwn = cellC[own].z();
        scalar zNei = cellC[nei].z();


        scalar alpha_con =
            (0.5 - lOwn) / (lNei - lOwn);


        scalar xInt =
            xOwn + alpha_con * (xNei - xOwn);

        scalar yInt =
            yOwn + alpha_con * (yNei - yOwn);

        scalar zInt =
            zOwn + alpha_con * (zNei - zOwn);


        // ================================================================
        // Width direction : X
        // ================================================================

        if (xInt < localMinX)
        {
            localMinX = xInt;
            localYatMinX = yInt;
        }

        if (xInt > localMaxX)
        {
            localMaxX = xInt;
            localYatMaxX = yInt;
        }


        // ================================================================
        // Length direction : Z
        // ================================================================

        if (zInt < localMinZ)
        {
            localMinZ = zInt;

            localYatMinZ = yInt;
            localXatMinZ = xInt;
        }

        if (zInt > localMaxZ)
        {
            localMaxZ = zInt;

            localYatMaxZ = yInt;
            localXatMaxZ = xInt;
        }


        localMeltInterfaceCount++;
    }
}



// ************************************************************************* //
// MPI reduction
// ************************************************************************* //


// ------------------------------------------------------------------------- //
// Global extrema
// ------------------------------------------------------------------------- //

scalar globalMinX = localMinX;
scalar globalMaxX = localMaxX;

scalar globalMinZ = localMinZ;
scalar globalMaxZ = localMaxZ;


label globalMeltInterfaceCount =
    localMeltInterfaceCount;


reduce
(
    globalMinX,
    minOp<scalar>()
);

reduce
(
    globalMaxX,
    maxOp<scalar>()
);

reduce
(
    globalMinZ,
    minOp<scalar>()
);

reduce
(
    globalMaxZ,
    maxOp<scalar>()
);

reduce
(
    globalMeltInterfaceCount,
    sumOp<label>()
);



// ************************************************************************* //
// Get coordinates corresponding to global extrema
// ************************************************************************* //


// ------------------------------------------------------------------------- //
// Width extrema coordinates
// ------------------------------------------------------------------------- //

scalar globalYatMinX = GREAT;
scalar globalYatMaxX = GREAT;


if
(
    localMeltInterfaceCount > 0
 && mag(localMinX - globalMinX) < tolerance
)
{
    globalYatMinX =
        localYatMinX;
}


if
(
    localMeltInterfaceCount > 0
 && mag(localMaxX - globalMaxX) < tolerance
)
{
    globalYatMaxX =
        localYatMaxX;
}


reduce
(
    globalYatMinX,
    minOp<scalar>()
);

reduce
(
    globalYatMaxX,
    minOp<scalar>()
);



// ------------------------------------------------------------------------- //
// Length extrema coordinates
// ------------------------------------------------------------------------- //

scalar globalYatMinZ = GREAT;
scalar globalYatMaxZ = GREAT;

scalar globalXatMinZ = GREAT;
scalar globalXatMaxZ = GREAT;


if
(
    localMeltInterfaceCount > 0
 && mag(localMinZ - globalMinZ) < tolerance
)
{
    globalYatMinZ =
        localYatMinZ;

    globalXatMinZ =
        localXatMinZ;
}


if
(
    localMeltInterfaceCount > 0
 && mag(localMaxZ - globalMaxZ) < tolerance
)
{
    globalYatMaxZ =
        localYatMaxZ;

    globalXatMaxZ =
        localXatMaxZ;
}


reduce
(
    globalYatMinZ,
    minOp<scalar>()
);

reduce
(
    globalYatMaxZ,
    minOp<scalar>()
);

reduce
(
    globalXatMinZ,
    minOp<scalar>()
);

reduce
(
    globalXatMaxZ,
    minOp<scalar>()
);



// ************************************************************************* //
// Final melt pool quantities
// ************************************************************************* //

if (globalMeltInterfaceCount > 0)
{
    // ---------------------------------------------------------------------
    // Width
    // ---------------------------------------------------------------------

    scalar meltWidth =
        globalMaxX - globalMinX;


    // ---------------------------------------------------------------------
    // Length
    // ---------------------------------------------------------------------

    scalar meltLength =
        globalMaxZ - globalMinZ;


    // ---------------------------------------------------------------------
    // Head / tail
    //
    // Assumption:
    // laser moves toward +z
    //
    // globalMinZ -------- deepest KH -------- globalMaxZ
    //     tail                  |                 head
    //
    // ---------------------------------------------------------------------

    scalar meltTailLength = 0.0;
    scalar meltHeadLength = 0.0;
    scalar tailHeadRatio = 0.0;

    scalar keyholeRelativePosition = 0.0;


    if (globalKeyholeInterfaceCount > 0)
    {
        meltTailLength =
            globalZatMaxKeyholeY - globalMinZ;

        meltHeadLength =
            globalMaxZ - globalZatMaxKeyholeY;


        if (meltHeadLength > SMALL)
        {
            tailHeadRatio =
                meltTailLength / meltHeadLength;
        }


        if (meltLength > SMALL)
        {
            keyholeRelativePosition =
                (globalZatMaxKeyholeY - globalMinZ)
              / meltLength;
        }
    }



    // =====================================================================
    // Output
    // =====================================================================

    Info << "melt left x: "
         << globalMinX
         << endl;

    Info << "melt left y: "
         << globalYatMinX
         << endl;


    Info << "melt right x: "
         << globalMaxX
         << endl;

    Info << "melt right y: "
         << globalYatMaxX
         << endl;


    Info << "melt width: "
         << meltWidth
         << endl;



    Info << "melt tail z: "
         << globalMinZ
         << endl;

    Info << "melt tail x: "
         << globalXatMinZ
         << endl;

    Info << "melt tail y: "
         << globalYatMinZ
         << endl;



    Info << "melt head z: "
         << globalMaxZ
         << endl;

    Info << "melt head x: "
         << globalXatMaxZ
         << endl;

    Info << "melt head y: "
         << globalYatMaxZ
         << endl;



    Info << "melt length: "
         << meltLength
         << endl;



    if (globalKeyholeInterfaceCount > 0)
    {
        Info << "deepest keyhole z: "
             << globalZatMaxKeyholeY
             << endl;

        Info << "melt tail length: "
             << meltTailLength
             << endl;

        Info << "melt head length: "
             << meltHeadLength
             << endl;

        Info << "tail/head ratio: "
             << tailHeadRatio
             << endl;

        Info << "keyhole relative position: "
             << keyholeRelativePosition
             << endl;
    }
    else
    {
        Info << "deepest keyhole z: 0"
             << endl;

        Info << "melt tail length: 0"
             << endl;

        Info << "melt head length: 0"
             << endl;

        Info << "tail/head ratio: 0"
             << endl;

        Info << "keyhole relative position: 0"
             << endl;
    }
}
else
{
    Info << "melt left x: 0"
         << endl;

    Info << "melt left y: 0"
         << endl;

    Info << "melt right x: 0"
         << endl;

    Info << "melt right y: 0"
         << endl;

    Info << "melt width: 0"
         << endl;


    Info << "melt tail z: 0"
         << endl;

    Info << "melt tail x: 0"
         << endl;

    Info << "melt tail y: 0"
         << endl;


    Info << "melt head z: 0"
         << endl;

    Info << "melt head x: 0"
         << endl;

    Info << "melt head y: 0"
         << endl;


    Info << "melt length: 0"
         << endl;

    Info << "deepest keyhole z: 0"
         << endl;

    Info << "melt tail length: 0"
         << endl;

    Info << "melt head length: 0"
         << endl;

    Info << "tail/head ratio: 0"
         << endl;

    Info << "keyhole relative position: 0"
         << endl;
}
