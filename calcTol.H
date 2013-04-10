//######## using limitted divergence as numerator############

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
	TOL1 = numerator1/denominator1;
	TOL2 = numerator2/denominator2;
/*
//######## using divergence of limitted gradient as numerator ############

	numerator = mag(fvc::div(gradAlphaNLim));
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
	TOL1 = numerator1/denominator1;
	TOL2 = numerator2/denominator2;
	Info << "divlimgrad" << endl;
	Info << "numerator1 = " << numerator1 << " denominator1 = " << denominator1 << endl;
	Info << "numerator2 = " << numerator2 << " denominator2 = " << denominator2 << endl;
	Info << "TOL1 " << TOL1 << " TOL2 " << TOL2 << endl;

//######## using divergence of limitted gradient as numerator ############
//######## using limitted gradient as denominator ############

	numerator = mag(fvc::div(gradAlphaNLim));
	denominator = Foam::pow(mag(gradAlphaNLim),2);
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
	TOL1 = numerator1/denominator1;
	TOL2 = numerator2/denominator2;
	Info << "limlim" << endl;
	Info << "numerator1 = " << numerator1 << " denominator1 = " << denominator1 << endl;
	Info << "numerator2 = " << numerator2 << " denominator2 = " << denominator2 << endl;
	Info << "TOL1 " << TOL1 << " TOL2 " << TOL2 << endl;

*/
