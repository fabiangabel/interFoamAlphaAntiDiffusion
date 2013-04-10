//#################################################################################
// This header performs a interface steepening method according to So et al.
// Comments contaning "(XX)" refer to the respective equations in So et al. paper
// diffusion equation is solved explicitely, boundary values remain unconsidered
//#################################################################################

//#################### Variable Definition/Initialization #########################

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
    surfaceScalarField  gradAlphaNfProjectedf(gradAlphaNf & mesh.Sf());
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
    int maxIter=5;
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
    Info << "Solving Anti-diffusion eqation: ";

    /*
    uMax = " << uMax << endl;
    Info << "deltaXMin = " << deltaXMin << endl;
    Info << "deltaTau = " << deltaTau << endl;
    */

//##################### Perform interface steepening method ########################

    for (int i=0; i < maxIter; i++)
    {
        //Info << "Iteration " << i << endl;

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
                    gradAlphaNLim[cellI] += alphaApostNf[faceI]*(mirror*mesh.Sf()[faceI]);
                    //gradAlphaNLim[cellI] += alphaApostNf[faceI]*mesh.Sf()[faceI];
                    //Info << cellI << " " << faceI << " " << gradAlphaNLim[cellI] << endl;
                }
            }
        }

        // Boundary Treatment (implemented: fixedValue, zeroGradient)
        forAll(mesh.boundaryMesh(), patchI)
        {   
            vectorField nf = mesh.boundary()[patchI].nf();
	        if (alpha1.boundaryField()[patchI].type() == "fixedValue")
	        // if (alpha1.boundaryField()[patchI].fixesValue())
        	{
        	    Info << "fixedValue" << endl;
        		forAll(alpha1.boundaryField()[patchI], faceI)
        		{
        			label faceIGlobal = faceI+mesh.boundaryMesh()[patchI].start();
        		    label own = mesh.owner()[faceIGlobal];
                    gradAlphaNProjectedf.boundaryField()[patchI][faceI] = gradAlphaNCd[own] & nf[faceI];
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
                        alphaApostNf.boundaryField()[patchI][faceI] = alpha1[own];
                        //Info << "case3, alphaApostNf[faceI] = " << alphaApostNf.boundaryField()[patchI][faceI] << endl;
                    }
                    // no mirror needed since faces of patches always store an outward pointing normal vector
                    gradAlphaNLim[own] += alphaApostNf.boundaryField()[patchI][faceI]*mesh.Sf()[faceIGlobal];
        	    }
        	}
            //same as fixedValue?
            else if (alpha1.boundaryField()[patchI].type() == "processor")
            {
                forAll(alpha1.boundaryField()[patchI], faceI)
                {
                    //Info << nf[faceI] << endl;
                    label faceIGlobal = faceI+mesh.boundaryMesh()[patchI].start();
                    label own = mesh.owner()[faceIGlobal];                    
                    gradAlphaNProjectedf.boundaryField()[patchI][faceI] = gradAlphaNCd[own] & nf[faceI];
                    if (gradAlphaNProjectedf.boundaryField()[patchI][faceI]*dAlphaNDcf.boundaryField()[patchI][faceI] > 0.0 &&
                        mag(dAlphaNDcf.boundaryField()[patchI][faceI]) < mag(gradAlphaNProjectedf.boundaryField()[patchI][faceI]))
                    {
                        alphaApostNf.boundaryField()[patchI][faceI] = alpha1.boundaryField()[patchI][faceI];
                        //Info << "case1, alphaApostNf[faceI] = " << alphaApostNf[faceI] << endl;
                    }
                    else if (gradAlphaNProjectedf[faceI]*dAlphaNDcf[faceI] > 0.0 &&
                             mag(gradAlphaNProjectedf[faceI]) == mag(dAlphaNDcf[faceI]))
                    {
                        alphaApostNf.boundaryField()[patchI][faceI] = alphaNf.boundaryField()[patchI][faceI];
                        //Info << "case2, alphaApostNf[faceI] = " << alphaApostNf[faceI] << endl;
                    }
                    else
                    {
                        alphaApostNf.boundaryField()[patchI][faceI] = alpha1[own];
                        //Info << "case3, alphaApostNf[faceI] = " << alphaApostNf[faceI] << endl;
                    }
                    Info << "\t" << own << " " << alpha1[own] 
                                << " " << alpha1.boundaryField()[patchI][faceI] 
                                //<< " " << alphaApostNf.boundaryField()[patchI][faceI]
                                << " " << gradAlphaNProjectedf.boundaryField()[patchI][faceI]
                                << " " << dAlphaNDcf.boundaryField()[patchI][faceI]
                                << " " << gradAlphaNLim[own]
                                //<< " " << endl;
                                << " " ;
                    gradAlphaNLim[own] += alphaApostNf.boundaryField()[patchI][faceI]*mesh.Sf()[faceIGlobal];
                    Info << gradAlphaNLim[own] << endl;
        	    }
            }
        	//else if (alpha1.boundaryField()[patchI].type() == "zeroGradient")
            else
            {
                forAll(alpha1.boundaryField()[patchI], faceI)
        	    {
        	        label faceIGlobal = faceI+mesh.boundaryMesh()[patchI].start();
        	        label own = mesh.owner()[faceIGlobal];
                    gradAlphaNLim[own] += alpha1[own]*mesh.Sf()[faceIGlobal];
                    //Info << "zeroGradient " << own << " " << mesh.Sf()[faceI] << endl;
                }
            }
        }

        // (27) Gauss approximated gradient at cell centre (CDS) ###################
        //##########################################################################
        forAll (mesh.cells(),cellI)
        {
            gradAlphaNLim[cellI] = gradAlphaNLim[cellI]/mesh.V()[cellI];
            //Info << cellI << " gradAlphaNLim[cellI]= " << gradAlphaNLim[cellI] << " " << mesh.V()[cellI] << endl;
            //Info << endl;
        }       

        // (28) flux limiting due to minimization of gradient al cell face #########
        //##########################################################################   
        forAll (mesh.faces(),faceI)
        {
            own = mesh.owner()[faceI];
            label patchI = mesh.boundaryMesh().whichPatch(faceI);
            if (patchI == -1)
            {
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
            //else if PROC minEqOp
        }
        //Info << gradAlphaNf << endl;
        // Necessary to ensure boundedness and global continuity
        forAll (gradAlphaNf.boundaryField(), patchI)
        {
            gradAlphaNf.boundaryField()[patchI] *= 0.0;
        }
        //Info << gradAlphaNf << endl;

        // (33) Calculate and check tolerances #####################################
        //##########################################################################

        gradAlphaNfProjectedf = gradAlphaNf & mesh.Sf();
        
        //Info << gradAlphaNfProjectedf << endl;
        syncTools::syncFaceList(mesh, gradAlphaNfProjectedf, minEqOp<scalar>());

        // using limitted divergence as numerator
	    numerator = mag(fvc::div(gradAlphaNfProjectedf));
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
	    Foam::reduce(numerator2, sumOp<scalar>());
	    Foam::reduce(denominator2, sumOp<scalar>());
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
            //fvc::div(-mag(phi/mesh.magSf()) * gradAlphaNf & mesh.Sf());
            fvc::div(-mag(phi/mesh.magSf()) * gradAlphaNfProjectedf);

        // (31) solve discretized equation explicitely ############################
        //##########################################################################
        alpha1=alpha1Old+diffGradAlpha*deltaTau;

        // output results #########################################################
        //##########################################################################
        /*
        Info << "32 - after regulation" << endl;
        forAll (mesh.cells(),cellI)
        {
            Info << cellI << " "
                 << alpha1OOld[cellI] << " "
                 << alpha1Old[cellI] << " "
                 << alpha1[cellI] << " "
                 << gradAlphaNCd[cellI] << " "
                 << gradAlphaNLim[cellI] <<  " "
                 << diffGradAlpha[cellI] << endl;
        }
        */
    }

    //alpha1.write();
    //U.write();
    //std::cout << 1/0 << endl;

//#################################################################################
