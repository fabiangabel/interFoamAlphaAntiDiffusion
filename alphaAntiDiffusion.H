//#################################################################################
// This header performs a interface steepening method according to So et al.
// Comments contaning "(XX)" refer to the respective equations in So et al. paper
// diffusion equation is solved explicitely, boundary values remain unconsidered
//#################################################################################

//#################### Variable Definition/Initialization #########################
{
    label own,neigh,faceI;
    scalar mirror;

    // save values of field alpha1 before postprocessing
    volScalarField alpha1Old(alpha1);

    // other Fields
    volVectorField      gradAlphaNCd(fvc::grad(alpha1));
    surfaceScalarField  alphaNf(fvc::interpolate(alpha1));
    surfaceScalarField  gradAlphaNProjectedf(alphaNf);
    surfaceScalarField  dAlphaNDcf(fvc::snGrad(alpha1));
    surfaceScalarField  alphaApostNf(alphaNf);
    volVectorField      gradAlphaNLim(fvc::grad(alphaApostNf));
    surfaceVectorField  gradAlphaNf(fvc::interpolate(gradAlphaNCd));
    //surfaceScalarField  gradAlphaNfProjectedf(gradAlphaNf & mesh.Sf());
    volScalarField      diffGradAlpha(fvc::div(-mag(phi/(mesh.magSf())) * gradAlphaNf & mesh.Sf()));

    /*volScalarField alpha1OldOut
    (
        IOobject
        (
            "alpha1Old",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1Old
    );
    alpha1OldOut=alpha1Old;
    //alpha1OldOut.write();

    volVectorField UOldOut
    (
        IOobject
        (
            "UOld",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U
    );

    UOldOut=U;
    UOldOut.write();
    */

    // initialize control variables (geplant: maxIter in fVSolution angeben)
    int maxIter=10;
    scalar TOL1=0.0;
    scalar TOL2=0.0;
    scalar numerator1 = 0.0;
    scalar denominator1 = 0.0;
    scalar numerator2 = 0.0;
    scalar denominator2 = 0.0;
    volScalarField numerator(mag(fvc::div(gradAlphaNLim)));
    volScalarField denominator(Foam::pow(mag(gradAlphaNCd),2));

//##################### Find uMax and deltaXMin and set pseudo time step #########

    // Two options either use U or phi/mesh.magSf() to determine max-velocity
    // No difference for inflow outflow, since inflow velocity equals max-velocity
    dimensionedScalar uMax("uMax",dimensionSet(0,1,-1,0,0,0,0),
                //Foam::max(mag(U)).value());
                Foam::max(phi/mesh.magSf()).value());
    dimensionedScalar deltaXMin("deltaXMin",dimensionSet(0,1,0,0,0,0,0),
                1.0/Foam::max(mesh.deltaCoeffs()).value());

    Foam::reduce(uMax, maxOp<dimensionedScalar>());
    Foam::reduce(deltaXMin, minOp<dimensionedScalar>());

    //##### (32) set pseudo time step ############################################
    dimensionedScalar deltaTau("deltaTau",dimensionSet(0,1,1,0,0,0,0),
                0.25*Foam::pow(deltaXMin.value(),2)/(uMax.value()+VSMALL));

    //Info << "Solving Anti-diffusion eqation: " << endl;
    Info << "Solving Anti-diffusion eqation: " ;
    //Sout << "Solving Anti-diffusion eqation: " << endl;

//##################### Perform interface steepening method ########################

    for (int i=0; i < maxIter; i++)
    {
        // save last iteration value
        alpha1Old = alpha1;

        // (24) Gauss approximated gradient at cell centre (CDS) ###################
        // adapt interpolate if working on non-orthogonal grids ####################
        //##########################################################################
        gradAlphaNCd = fvc::grad(alpha1);
        alphaNf = fvc::interpolate(alpha1);

        // (25) CDS approximated gradient at cell face #############################
        //##########################################################################
        dAlphaNDcf = fvc::snGrad(alpha1);

        // (26) perform regulation of pde ##########################################
        //##########################################################################
        //alphaApostNf = alphaNf;
        forAll (mesh.cells(),cellI)
        {
            cell faces = mesh.cells()[cellI];
            gradAlphaNLim[cellI] *= 0.0;
            forAll (faces, i)
            {
                faceI = faces[i];
                label patchI = mesh.boundaryMesh().whichPatch(faceI);
                if (patchI == -1)
                {
                    // Project cell centered gradient onto cell-face-normal direction
                    gradAlphaNProjectedf[faceI] =
                        gradAlphaNCd[cellI] & mesh.Sf()[faceI]/(mesh.magSf()[faceI]);
                    //Info << cellI << " " << gradAlphaNCd[cellI] << " " << mesh.Sf()[faceI] << " " << dAlphaNDcf[faceI] << endl;
                    own = mesh.owner()[faceI];
                    neigh = mesh.neighbour()[faceI];
                    // switch owner and neighbour if necessary ######################
                    if (neigh == cellI)
                    {
                        neigh = own;
                        own = cellI;
                        // mirror takes into account, that the surface vector always should point outwards
                        mirror = -1.0;
                    }
                    else mirror = 1.0;
                    if (gradAlphaNProjectedf[faceI]*dAlphaNDcf[faceI] > 0.0 &&
                        mag(dAlphaNDcf[faceI]) < mag(gradAlphaNProjectedf[faceI]))
                    {
                        alphaApostNf[faceI] = alpha1[neigh];
                        //Info << "case1, alphaApostNf[faceI] = " << alphaApostNf[faceI] << endl;
                    }
                    else if (gradAlphaNProjectedf[faceI]*dAlphaNDcf[faceI] > 0.0 &&
                             mag(gradAlphaNProjectedf[faceI]) == mag(dAlphaNDcf[faceI]))
                    {
                        alphaApostNf[faceI] = alphaNf[faceI];
                        //Info << "case2, alphaApostNf[faceI] = " << alphaApostNf[faceI] << endl;
                    }
                    else
                    {
                        alphaApostNf[faceI] = alpha1[own];
                        //Info << "case3, alphaApostNf[faceI] = " << alphaApostNf[faceI] << endl;
                    }
                    /*Sout << "INT:  " << mesh.C()[0].component(0)
                        << " " << own
                        << " " << alpha1[own]
                        << " " << alpha1[neigh]
                        << " " << mirror*mesh.Sf()[faceI]
                        //<< " " << alphaApostNf.boundaryField()[patchI][faceI]
                        << " " << mirror*gradAlphaNProjectedf[faceI]
                        << " " << mirror*dAlphaNDcf[faceI]
                        << " " << gradAlphaNLim[own]
                        //<< " " << endl;
                        << " " ;
                        */
                    gradAlphaNLim[cellI] += alphaApostNf[faceI]*(mirror*mesh.Sf()[faceI]);
                    //Sout << gradAlphaNLim[cellI] << endl;
                    //gradAlphaNLim[cellI] += alphaApostNf[faceI]*mesh.Sf()[faceI];
                    //Info << cellI << " " << faceI << " " << gradAlphaNLim[cellI] << endl;
                }
            }
        }
        //Sout << " " << endl;

        // Boundary Treatment (implemented: fixedValue, zeroGradient)
        forAll (mesh.boundaryMesh(), patchI)
        {
            vectorField nf = mesh.boundary()[patchI].nf();
            const labelList & faceCells = mesh.boundary()[patchI].patch().faceCells();
	    if (alpha1.boundaryField()[patchI].type() == "fixedValue")
            // if (alpha1.boundaryField()[patchI].fixesValue())
            {
                //forAll(alpha1.boundaryField()[patchI], faceI)
                forAll (faceCells, faceI)
                {
                    //label faceIGlobal = faceI+mesh.boundaryMesh()[patchI].start();
                    //label own = mesh.owner()[faceIGlobal];
                    //gradAlphaNProjectedf.boundaryField()[patchI][faceI] = gradAlphaNCd[own] & nf[faceI];
                    gradAlphaNProjectedf.boundaryField()[patchI][faceI] = gradAlphaNCd[faceCells[faceI]] & nf[faceI];
                    if (gradAlphaNProjectedf.boundaryField()[patchI][faceI]*dAlphaNDcf.boundaryField()[patchI][faceI] > 0.0 &&
                        mag(dAlphaNDcf.boundaryField()[patchI][faceI]) < mag(gradAlphaNProjectedf.boundaryField()[patchI][faceI]))
                    {
                        alphaApostNf.boundaryField()[patchI][faceI] = alpha1.boundaryField()[patchI][faceI];
                        //Info << "case1, alphaApostNf[faceI] = " << alphaApostNf.boundaryField()[patchI][faceI] << endl;
                    }
                    else if ((gradAlphaNProjectedf.boundaryField()[patchI][faceI]*dAlphaNDcf.boundaryField()[patchI][faceI] > 0.0 &&
                        mag(dAlphaNDcf.boundaryField()[patchI][faceI]) == mag(gradAlphaNProjectedf.boundaryField()[patchI][faceI])))
                    {
                        alphaApostNf.boundaryField()[patchI][faceI] = alphaNf.boundaryField()[patchI][faceI];
                        //Info << "case2, alphaApostNf[faceI] = " << alphaApostNf.boundaryField()[patchI][faceI] << endl;
                    }
                    else
                    {
                        //alphaApostNf.boundaryField()[patchI][faceI] = alpha1[own];
                        alphaApostNf.boundaryField()[patchI][faceI] = alpha1[faceCells[faceI]];
                        //Info << "case3, alphaApostNf[faceI] = " << alphaApostNf.boundaryField()[patchI][faceI] << endl;
                    }
                    // no mirror needed since faces of patches always store an outward pointing normal vector
                    gradAlphaNLim[faceCells[faceI]] += alphaApostNf.boundaryField()[patchI][faceI]*mesh.Sf().boundaryField()[patchI][faceI];
                    //Sout << gradAlphaNLim[own] << endl;
                }
            }
            //same as fixedValue?
            else if (alpha1.boundaryField()[patchI].type() == "processor")
            {
                //forAll(alpha1.boundaryField()[patchI], faceI)
                forAll (faceCells, faceI)
                {
                    //Info << nf[faceI] << endl;
                    //label faceIGlobal = faceI+mesh.boundaryMesh()[patchI].start();
                    //label own = mesh.owner()[faceIGlobal];
                    gradAlphaNProjectedf.boundaryField()[patchI][faceI] = gradAlphaNCd[faceCells[faceI]] & nf[faceI];
                    if (gradAlphaNProjectedf.boundaryField()[patchI][faceI]*dAlphaNDcf.boundaryField()[patchI][faceI] > 0.0 &&
                        mag(dAlphaNDcf.boundaryField()[patchI][faceI]) < mag(gradAlphaNProjectedf.boundaryField()[patchI][faceI]))
                    {
                        alphaApostNf.boundaryField()[patchI][faceI] = alpha1.boundaryField()[patchI][faceI];
                        //Info << "case1, alphaApostNf[faceI] = " << alphaApostNf[faceI] << endl;
                    }
                    else if (gradAlphaNProjectedf.boundaryField()[patchI][faceI]*dAlphaNDcf.boundaryField()[patchI][faceI] > 0.0 &&
                             mag(gradAlphaNProjectedf.boundaryField()[patchI][faceI]) == mag(dAlphaNDcf.boundaryField()[patchI][faceI]))
                    {
                        alphaApostNf.boundaryField()[patchI][faceI] = alphaNf.boundaryField()[patchI][faceI];
                        //Info << "case2, alphaApostNf[faceI] = " << alphaApostNf[faceI] << endl;
                    }
                    else
                    {
                        alphaApostNf.boundaryField()[patchI][faceI] = alpha1[faceCells[faceI]];
                        //Info << "case3, alphaApostNf[faceI] = " << alphaApostNf[faceI] << endl;
                    }
                    gradAlphaNLim[faceCells[faceI]] += alphaApostNf.boundaryField()[patchI][faceI]*mesh.Sf().boundaryField()[patchI][faceI];
                    //Sout << gradAlphaNLim[own] << endl;
                }
            }
            //else if (alpha1.boundaryField()[patchI].type() == "zeroGradient")
            else if (alpha1.boundaryField()[patchI].type() != "empty")
            {
                //forAll(alpha1.boundaryField()[patchI], faceI)
                forAll (faceCells, faceI)
        	{
                    //Info << faceI << " " << faceCells[faceI] << endl;
                    gradAlphaNLim[faceCells[faceI]] += alpha1[faceCells[faceI]]*mesh.Sf().boundaryField()[patchI][faceI];
                }
            }
        }

        // (27) Gauss approximated gradient at cell centre (CDS) ###################
        //##########################################################################
        forAll (mesh.cells(),cellI)
        {
            gradAlphaNLim[cellI] = gradAlphaNLim[cellI]/mesh.V()[cellI];
        }

        // update cell centre gradient
        gradAlphaNLim.correctBoundaryConditions();

        // (28) flux limiting due to minimization of gradient al cell face #########
        //##########################################################################
        forAll (mesh.faces(),faceI)
        {
            label patchI = mesh.boundaryMesh().whichPatch(faceI);
            if (patchI == -1)
            {
                own = mesh.owner()[faceI];
                label neigh = mesh.neighbour()[faceI];
                if (mag(gradAlphaNLim[own]) <= mag(gradAlphaNLim[neigh]))
                {
                    gradAlphaNf[faceI] = gradAlphaNLim[own];
                }
                else
                {
                    gradAlphaNf[faceI] = gradAlphaNLim[neigh];
                }
            }
        }

        // Necessary to ensure boundedness and global continuity
        forAll (gradAlphaNf.boundaryField(), patchI)
        {
            const labelList & faceCells = mesh.boundary()[patchI].patch().faceCells();
            if (alpha1.boundaryField()[patchI].type() == "processor")
            {
                //forAll(alpha1.boundaryField()[patchI], faceI)
                forAll (faceCells, faceI)
                {
                    //label faceIGlobal = faceI+mesh.boundaryMesh()[patchI].start();
                    //label own = mesh.owner()[faceIGlobal];
                    if (mag(gradAlphaNLim[faceCells[faceI]]) <= mag(gradAlphaNLim.boundaryField()[patchI][faceI]))
                    {
                        //gradAlphaNf.boundaryField()[patchI][faceI] = gradAlphaNLim[own];
                        gradAlphaNf.boundaryField()[patchI][faceI] = gradAlphaNLim[faceCells[faceI]];
                    }
                    else
                    {
                        gradAlphaNf.boundaryField()[patchI][faceI] = gradAlphaNLim.boundaryField()[patchI][faceI];
                    }
                }
            }
            else gradAlphaNf.boundaryField()[patchI] *= 0.0;
        }

        // (33) Calculate and check tolerances #####################################
        //##########################################################################

        //gradAlphaNfProjectedf = gradAlphaNf & mesh.Sf();
        //not needed, since gradAlphaNLim.correctBoundaryConditions() already updates the necessary values
        //forAll (gradAlphaNfProjectedf.boundaryField(), patchI)
        //{
        //    syncTools::syncBoundaryFaceList(gradAlphaNfProjectedf.boundaryField()[patchI], mesh.boundary()[patchI].patch().faceCells(), minEqOp<scalar>());
        //}

        // using limitted divergence as numerator
	//numerator = mag(fvc::div(gradAlphaNfProjectedf));
	numerator = mag(fvc::div(gradAlphaNf & mesh.Sf()));
	denominator = Foam::pow(mag(gradAlphaNCd),2);
	numerator1 = denominator1 = numerator2 = denominator2 = 0.0;
	forAll (mesh.cells(),cellI)
	{
            numerator1 += numerator[cellI]*mesh.V()[cellI];
            denominator1 += denominator[cellI]*mesh.V()[cellI];
            if (numerator2 < numerator[cellI]*mesh.V()[cellI])
            {
                numerator2 = numerator[cellI]*mesh.V()[cellI];
            }
            if (denominator2 < denominator[cellI]*mesh.V()[cellI])
            {
                denominator2 = denominator[cellI]*mesh.V()[cellI];
            }
	}
	Foam::reduce(numerator1, sumOp<scalar>());
	Foam::reduce(denominator1, sumOp<scalar>());
	Foam::reduce(numerator2, maxOp<scalar>());
	Foam::reduce(denominator2, maxOp<scalar>());
	TOL1 = numerator1/denominator1;
	TOL2 = numerator2/denominator2;
        if ((TOL1 <= 1.0 && TOL2 <= 1.0) || i==maxIter-1)
        //if (TOL2 <= 1.0)
        {
            Info //<< "Time = " << runTime.timeName()
                << " TOL1 = " << TOL1
                << ", TOL2 = " << TOL2
                << ", No Iterations " << i << endl;
            break;
        }

        // (30) calculate diffusionTerm ###########################################
        //##########################################################################
        diffGradAlpha =
            fvc::div(-mag(phi/mesh.magSf()) * gradAlphaNf & mesh.Sf());
            //fvc::div(-mag(phi/mesh.magSf()) * gradAlphaNfProjectedf);

        //forAll (mesh.cells(),cellI)
        //{
        //    Sout << mesh.C()[0].component(0) << " " << diffGradAlpha[cellI] << endl;
        //}

        // (31) solve discretized equation explicitely ############################
        //##########################################################################
        alpha1=alpha1Old+diffGradAlpha*deltaTau;

        // update values for next iteration
        alpha1.correctBoundaryConditions();
        gradAlphaNCd.correctBoundaryConditions();
        alphaNf = fvc::interpolate(alpha1);
        dAlphaNDcf = fvc::snGrad(alpha1);
        //apparently it is not possible to update surface fields using .correctBoundaryConditions()
        //alphaNf.correctBoundaryConditions();
        //dAlphaNDcf.correctBoundaryConditions();
    }
    //alpha1.write();
    //U.write();
    //std::cout << 1/0 << endl;
}
//#################################################################################
