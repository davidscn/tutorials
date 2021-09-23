// if (runTime.writeTime())
{
    volVectorField gradT(fvc::grad(T));

    volScalarField gradTx(
        IOobject(
            "gradTx",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        gradT.component(vector::X));

    volScalarField gradTy(
        IOobject(
            "gradTy",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        gradT.component(vector::Y));

    volScalarField gradTz(
        IOobject(
            "gradTz",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        gradT.component(vector::Z));

    volVectorField DTgradT(
        IOobject(
            "flux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        DT * gradT);

    {
        double error = 0;
        const Foam::volScalarField *T_(&mesh.lookupObject<volScalarField>("T"));
        const double time = runTime.value();

        // Get the locations of the volume centered mesh vertices
        const vectorField &CellCenters = mesh.C();
        unsigned int numDataLocations = CellCenters.size();

        for (int i = 0; i < CellCenters.size(); i++)
        {
            const double coord_x = CellCenters[i].x();
            const double coord_y = CellCenters[i].y();

            const double exact_solution = 1 + (coord_x * coord_x) +
                                          (alpha * coord_y * coord_y) + beta * (time);

            //  mesh.magSf()[i] * 10 * 10 * mesh.magSf()[i]
            error += (exact_solution - T_->internalField()[i]) *
                     (exact_solution - T_->internalField()[i]);

            if (i < 10)
                Info << (exact_solution - T_->internalField()[i])
                     << endl;
        }

        // Get the locations of the mesh vertices (here: face centers)
        // for all the patches
        // We need to exclude the default (empty faces)
        // for (int j = 0; j < mesh.boundaryMesh().size() - 1; ++j)
        // {
        //     // Get the face centers of the current patch
        //     const vectorField faceCenters =
        //         mesh.boundaryMesh()[j].faceCentres();

        //     numDataLocations += faceCenters.size();

        //     // Assign the (x,y,z) locations to the vertices
        //     for (int i = 0; i < faceCenters.size(); i++)
        //     {
        //         const double coord_x = faceCenters[i].x();
        //         const double coord_y = faceCenters[i].y();

        //         const double exact_solution = 1 + (coord_x * coord_x) +
        //                                       (alpha * coord_y * coord_y) + beta * time;

        // error += (exact_solution - T_->boundaryField()[j][i]) *
        //      (exact_solution - T_->boundaryField()[j][i]);
        //         }
        // }
        Info << "\n\nLocal square root error: " << std::sqrt(error / numDataLocations) << " at t = " << time << endl;
        Info << "Global absolute error: " << error << " at t = " << time << "\n"
             << endl;
    }

    runTime.write();
}
