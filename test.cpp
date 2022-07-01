// A numerical validation for the paper "Microsurface Transformations"
// by Asen Atanasov, Vladimir Koylazov, Rossen Dimov and Alexander Wilkie, EGSR'22

// A number of microsurfaces (GTR, GGX, Beckmann, Phong, Sheen, STD, Discrete) with random parameters are created.
// Their normalization and shadowing constraints are integrated numerically, then these microsurfaces
// are transformed by a random matrix M and the constraints are computed again. Note that in some cases like 
// Phong, Sheen, STD, Discrete the Smith shadowing approximations are not always accurate. When a shadowing
// constraint is not accurate for the original surface it is also not accurate for the transformed surface.

// Author: Asen Atanasov

#include <iostream>
#include <math.h>
#include <string.h>
#include "microsurface_transformations.h"

using namespace std;

int main() {

	const int numTests=10;
	const int numSamples=1000000;

	for(int i=1; i<=numTests; i++) {
		cout << "**********************************************************************************************" << endl;
		cout << "Test " << i << ":" << endl;

		const float alphaX=rnd();
		const float alphaY=rnd();
		const float gamma=rnd(0.0f, 4.0f);

		cout << "Random roughness and tail parameters:" << endl;
		cout << "AlphaX: " << alphaX << endl;
		cout << "AlphaY: " << alphaY << endl;
		cout << "Gamma: " << gamma << endl;
		cout << endl;

		const Vector outgoing=getSphereDir(rnd(), rnd());
		cout << "Random shadowing direction = (" << outgoing.x << ", " << outgoing.y << ", " << outgoing.z << ")" << endl << endl;

		Matrix m;
		m.setCol(0, Vector(rnd(-10.0f, 10.0f), rnd(-10.0f, 10.0f), 0.0f));
		m.setCol(1, Vector(rnd(-10.0f, 10.0f), rnd(-10.0f, 10.0f), 0.0f));
		m.setCol(2, Vector(0.0f, 0.0f, 1.0f));

		cout << "Generate matrix M with random entries in (-10, 10):" << endl;
		cout << "M = (" << m[0][0] << ", " << m[0][1] << ") (" << m[1][0] << ", " << m[1][1] << ")" << endl << endl;

		GTR gtr(alphaX, gamma);
		GGX ggx(alphaX, alphaY);
		Beckmann beckmann(alphaX, alphaY);
		Phong phong(alphaX);
		Sheen sheen(alphaX);
		STD std(alphaX, alphaY, gamma+1.5f);
		Discrete discrete(alphaX, gamma, 100);

		int numSurfaces=0;
		Microsurface *microsurfaces[8];
		microsurfaces[numSurfaces++]=&gtr;
		microsurfaces[numSurfaces++]=&ggx;
		microsurfaces[numSurfaces++]=&beckmann;
		microsurfaces[numSurfaces++]=&phong;
		microsurfaces[numSurfaces++]=&sheen;
		microsurfaces[numSurfaces++]=&std;
		microsurfaces[numSurfaces++]=&discrete;
		
		cout << "Test normalization and shadowing constraints for the original and transformed by M surfaces:" << endl;
		for(int j=0; j<numSurfaces; j++) {
			Microsurface *surface=microsurfaces[j];
			const char *name=surface->getName();
			const float norm=surface->integrateDistribution(numSamples);
			const float shadow=surface->shadowingConstraint(outgoing, numSamples);
			surface->initTransform(m);
			const float normM=surface->integrateDistribution(numSamples);
			const float shadowM=surface->shadowingConstraint(outgoing, numSamples);

			cout << name << " Normalization constraint:\toriginal=" << norm << ",\ttransformed=" << normM << endl;
			cout << name << "     Shadowing constraint:\toriginal=" << shadow << ",\ttransformed=" << shadowM << endl;
		}
	}

	return 0;
}