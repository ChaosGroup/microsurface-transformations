#pragma once

// A numerical validation for the paper "Microsurface Transformations"
// by Asen Atanasov, Vladimir Koylazov, Rossen Dimov and Alexander Wilkie, EGSR'22

// The basic technique to transform microsurfaces is implemented in the base class Microsurface.
// A list of implementations for specific distributions are provided afterwards:
// 1. GTR - "Physically Based Shading at Disney" by Burley (2012) and "Deriving the Smith shadowing function G_1 for \gamma \in (0,4]" by Dimov (2015)
// 2. GGX - Anisotropic GGX, "Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs" by Heitz (2014)
// 3. Beckmann - Anisotropic Beckmann, "Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs" by Heitz (2014)
// 4. Phong - "Microfacet Models for Refraction through Rough Surfaces" by Walter et al. (2007)
// 5. Sheen - "Production Friendly Microfacet Sheen BRDF", Estevez and Kulla (2017)
// 6. STD - "STD: Student's t-Distribution of Slopes for Microfacet Based BSDFs" by Ribardiere et al. (2017)
// 7. Discrete - discrete GGX with classical Smith shadowing, "Discrete Stochastic Microfacet Models", Jakob et al. (2014)

// Author: Asen Atanasov
// The code for GTR D and G1 is copied from V-Ray SDK

#include "math_utils.h"

class Microsurface {
	Matrix m, im;
	float det;
	virtual float getD(Vector h) const = 0;
	virtual float getG1(Vector dir) const = 0;

public:

	virtual const char* getName() const = 0;

	// Inilialize a matrix of the form given in Equation (9)
	void initTransform(Matrix m) {
		this->m=m;
		this->im=m;
		im.makeInverse();
		det=fabs(mixed(m[0], m[1], m[2]));
	}

	// The transformed distribution, Equation (22)
	float getMicrofacetDistribution(Vector h) const {
		Vector hTransformed=h*m; // multiply with the transpose of m
		const float normalization=1.0f/hTransformed.lengthSqr();
		hTransformed*=sqrtf(normalization);
		return det*sqr(normalization)*getD(hTransformed);
	}

	// The transformed shadowing function, Equation (12)
	float getSmithG1(Vector dir) const {
		const Vector dirTransformed=normalize(im*dir);
		return getG1(dirTransformed);
	}

	// Normalization constraint, Equation (3)
	float integrateDistribution(int numSamples) const {
		double sum=0.0;
		for(int i=0; i<numSamples; i++) {
			const float u = float(i)/float(numSamples);
			const float v = radicalInverse2(i);
			const Vector h=getSphereDir(u, v);
			const float hz=h.z;
			sum += double(getMicrofacetDistribution(h)*hz);
		}
		const float res = 2.0f*PI*float(sum)/float(numSamples);
		return res;
	}

	// Shadowing function constraint, Equation (2)
	float shadowingConstraint(const Vector &outgoing, int numSamples) const {
		double sum=0.0f;
		for(int i=0; i<numSamples; i++) {
			const float u = float(i)/float(numSamples);
			const float v = radicalInverse2(i);
			const Vector h=getSphereDir(u, v);
			sum += double(std::max(0.0f, dotf(h, outgoing)*getMicrofacetDistribution(h)*getSmithG1(outgoing)));
		}
		const float res = 2.0f*PI*float(sum)/float(numSamples);
		return res/outgoing.z;
	}
};

template <int KNOTS>
static float evalEquidistantNaturalSpline(float knots[KNOTS], float t) {
	float f[KNOTS-2];
	float m[KNOTS];
	m[0] = 0.0f;
	m[KNOTS-1] = 0.0f;

	// Initialize moments m and result vector f.
	for(int i=0; i<KNOTS-2; i++) {
		m[i+1] = 4.0f;
		f[i] = knots[i+2]-2.0f*knots[i+1]+knots[i];
	}

	// Solve tridiagonal system. The solution is in m.
	for(int i=0; i<KNOTS-3; i++) {
		const float k = 1.0f/m[i+1];
		m[i+2] -= k;
		f[i+1] -= k*f[i];
	}
	for(int i=KNOTS-3; i>=0; i--) {
		m[i+1] = (f[i]-m[i+2])/m[i+1];
	}

	const int i=clamp(int(floorf(t)), 0, KNOTS-2);

	// Calculate polynomial coefficients.
	const float a = m[i+1]-m[i];
	const float b = 3.0f*m[i];
	const float c = knots[i+1]-knots[i]-2.0f*m[i]-m[i+1];
	const float d = knots[i];

	const float x=t-float(i);
	const float res=((a*x+b)*x+c)*x+d;
	return res;
}

class GTR: public Microsurface {
	float alpha;
	float gamma;

	float getD(Vector h) const {
		float res=0.0f;
		const float cosTheta=h.z;
		if(cosTheta<=1e-3f)
			return res;

		const float cosTheta2=sqr(cosTheta);
		const float tanTheta2=(1.0f/cosTheta2)-1.0f;
		const float roughness2=sqr(alpha);

		if(fabsf(gamma-1.0f)<1e-6f) {
			const float div=PI*logf(roughness2)*cosTheta2*(roughness2+tanTheta2);
			res=(div<(roughness2-1.0f)*1e-6f)? (roughness2-1.0f)/div : 0.0f;
		}
		else {
			const float divisor=PI*(1.0f-powf(roughness2, 1.0f-gamma))*powf(cosTheta2*(roughness2+tanTheta2), gamma);
			const float dividend=(gamma-1.0f)*(roughness2-1.0f);
			res=(fabsf(divisor)>fabsf(dividend)*1e-6f)? dividend/divisor : 0.0f;
		}

		return res;
	}

	float getG1Gamma1(float alpha, float cotTheta) const {
		const float cotTheta2=sqr(cotTheta);
		const float alpha2=sqr(alpha);
		const float a=sqrtf(cotTheta2+alpha2);
		const float b=sqrtf(cotTheta2+1.0f);
		return cotTheta*logf(alpha2)/(a-b+cotTheta*logf(alpha2*(cotTheta+b)/(cotTheta+a)));
	}

	float getG1Gamma2(float alpha, float cotTheta) const {
		return 2.f/(1.f+sqrtf(1.f+sqr(alpha/cotTheta)));
	}

	float getG1Gamma3(float alpha, float cotTheta) const {
		const float cotTheta2=sqr(cotTheta);
		const float alpha2=sqr(alpha);
		const float a=sqrtf(cotTheta2+alpha2);
		const float b=alpha2+1.0f;
		return 4.0f*cotTheta*a*b/(2.0f*cotTheta*b*(cotTheta+a)+alpha2*(3.0f*alpha2+1.0f));
	}

	float getG1Gamma4(float alpha, float cotTheta) const {
		const float cotTheta2=sqr(cotTheta);
		const float alpha2=sqr(alpha);
		const float alpha4=sqr(alpha2);
		const float a=8.f*(alpha4+alpha2+1.f);
		const float b=sqrtf(cotTheta2+alpha2);
		const float b3=b*(cotTheta2+alpha2);
		return 2.f*cotTheta*a*b3/(a*cotTheta*(b3+cotTheta*cotTheta2)+3.f*alpha2*(4.f*cotTheta2*(2.f*alpha4+alpha2+1.f)+alpha2*(5.f*alpha4+2.f*alpha2+1.f)));
	}

	float getG1(Vector dir) const {
		const float cosTheta = dir.z;
		if(cosTheta<=1e-3f) return 0.0f;
		if(cosTheta>=1.0f-1e-6f) return 1.0f;
		const float cotTheta = cosTheta/sqrtf(1.f-sqr(cosTheta));

		float res=0.0f;
		// when gamma is any of the integer values 0, 1, 2, 3, 4 apply analytical solution
		if(gamma<=0.01f)
			res=getG1Gamma2(1.0f, cotTheta);
		else if(fabsf(gamma-1.0f)<=1e-2f)
			res=getG1Gamma1(alpha, cotTheta);
		else if(fabsf(gamma-2.0f)<=1e-2f)
			res=getG1Gamma2(alpha, cotTheta);
		else if(fabsf(gamma-3.0f)<=1e-2f)
			res=getG1Gamma3(alpha, cotTheta);
		else if(gamma>=4.0f-1e-2f)
			res=getG1Gamma4(alpha, cotTheta);
		else {
			float knots[5];
			knots[0]=getG1Gamma2(1.0f, cotTheta);
			knots[1]=getG1Gamma1(alpha, cotTheta);
			knots[2]=getG1Gamma2(alpha, cotTheta);
			knots[3]=getG1Gamma3(alpha, cotTheta);
			knots[4]=getG1Gamma4(alpha, cotTheta);
			res=evalEquidistantNaturalSpline<5>(knots, gamma);
		}
		return clamp(res, 0.0f, 1.0f);
	}

public:

	GTR(float alpha=0.1f, float gamma=2.0f, Matrix m=Matrix(1.0f)): alpha(alpha), gamma(gamma) {
		initTransform(m);
	}

	const char* getName() const  {
		return "GTR     ";
	}
};

inline float getAlpha(Vector dir, float alphaX, float alphaY) {
	const float cosTheta2    = sqr(dir.z);
	const float invSinTheta2 = 1.0f/(1.0f-cosTheta2);

	if(alphaX==alphaY || invSinTheta2<=0) {
		return alphaX;
	}

	const float cosPhi2 = sqr(dir.x)*invSinTheta2;
	const float sinPhi2 = sqr(dir.y)*invSinTheta2;

	return sqrt(cosPhi2*sqr(alphaX) + sinPhi2*sqr(alphaY));
}

class GGX: public Microsurface {
	float alphaX;
	float alphaY;

	float getD(Vector h) const {
		if(h.z<=1e-3f) return 0.0f;
		const float denominator=PI*alphaX*alphaY*sqr(sqr(h.x/alphaX)+sqr(h.y/alphaY)+sqr(h.z));
		return 1.0f/denominator;
	}

	float getG1(Vector dir) const {
		const float cosTheta = dir.z;
		if(cosTheta<=1e-3f) return 0.0f;
		if(cosTheta>=1.0f-1e-6f) return 1.0f;
		const float alpha=getAlpha(dir, alphaX, alphaY);
		const float tanTheta2=(1.0f/sqr(dir.z))-1.0f;
		const float denominator=1.0f+sqrtf(1.0f+tanTheta2*sqr(alpha));
		return 2.0f/denominator;
	}

public:

	GGX(float alphaX, float alphaY, Matrix m=Matrix(1.0f)): alphaX(alphaX), alphaY(alphaY) {
		initTransform(m);
	}

	const char* getName() const  {
		return "GGX     ";
	}
};

class Beckmann: public Microsurface {
	float alphaX;
	float alphaY;

	float getD(Vector h) const {
		const float cosTheta=h.z;
		if(cosTheta<=1e-3f) return 0.0f;
		const float cosTheta2=sqr(cosTheta);
		const float sinTheta2=1.0f-cosTheta2;
		float numerator=1.0f;
		if(sinTheta2>1e-6f) {
			const float tanTheta2=sinTheta2/cosTheta2;
			const float cosPhi2=sqr(h.x)/sinTheta2;
			const float sinPhi2=sqr(h.y)/sinTheta2;
			numerator=expf(-tanTheta2*(cosPhi2/sqr(alphaX)+sinPhi2/sqr(alphaY)));
		}
		const float denominator=PI*alphaX*alphaY*sqr(cosTheta2);
		return numerator/denominator;
	}

	float getG1(Vector dir) const {
		const float cosTheta = dir.z;
		if(cosTheta<=1e-3f) return 0.0f;
		if(cosTheta>=1.0f-1e-6f) return 1.0f;
		const float alpha=getAlpha(dir, alphaX, alphaY);
		const float tanTheta=sqrtf((1.0f/sqr(cosTheta))-1.0f);
		const float a=1.0f/(alpha*tanTheta);
		if(a<1.6f) {
			return (3.535f*a+2.181f*sqr(a))/(1.0f+2.276f*a+2.577f*sqr(a));
		}
		return 1.0f;
	}

public:

	Beckmann(float alphaX, float alphaY, Matrix m=Matrix(1.0f)): alphaX(alphaX), alphaY(alphaY) {
		initTransform(m);
	}

	const char* getName() const  {
		return "Beckmann";
	}
};

class Phong: public Microsurface {
	float alpha;

	float getD(Vector h) const {
		const float cosTheta=h.z;
		if(cosTheta<=1e-3f) return 0.0f;
		return (alpha+2.0f)*powf(cosTheta, alpha)/(2.0f*PI);
	}

	float getG1(Vector dir) const {
		const float cosTheta = dir.z;
		if(cosTheta<=1e-3f) return 0.0f;
		if(cosTheta>=1.0f-1e-6f) return 1.0f;
		const float tanTheta=sqrtf((1.0f/sqr(cosTheta))-1.0f);
		const float a=sqrtf(0.5f*alpha+1.0f)/tanTheta;
		if(a<1.6f) {
			return (3.535f*a+2.181f*sqr(a))/(1.0f+2.276f*a+2.577f*sqr(a));
		}
		return 1.0f;
	}

public:

	Phong(float alpha, Matrix m=Matrix(1.0f)): alpha(alpha) {
		initTransform(m);
	}

	const char* getName() const  {
		return "Phong   ";
	}
};

class Sheen: public Microsurface {
	float alpha;
	float a, b, c, d, e;

	float getD(Vector h) const {
		const float cosTheta=h.z;
		if(cosTheta<=1e-3f) return 0.0f;
		const float sinTheta=sqrtf(1.0f-sqr(cosTheta));
		const float power=alpha>1e-6f?1.0f/alpha:1e18f;
		return (2.0f+power)*powf(sinTheta, power)/(2.0f*PI);
	}

	void initG1() {
		const float t=sqr(1.0f-alpha);
		a=t*25.3245f+(1-t)*21.5473f;
		b=t*3.32435f+(1-t)*3.82987f;
		c=t*0.16801f+(1-t)*0.19823f;
		d=t*-1.27393f+(1-t)*-1.97760f;
		e=t*-4.85967f+(1-t)*-4.32054f;
	}

	float getL(float x) const {
		return a/(1.0f+b*powf(x, c))+d*x+e;
	}

	float getG1(Vector dir) const {
		const float cosTheta = dir.z;
		if(cosTheta<=1e-3f) return 0.0f;
		if(cosTheta>=1.0f-1e-6f) return 1.0f;
		if(dir.z<0.5f) {
			const float lambda=expf(getL(cosTheta));
			return 1.0f/(1.0f+lambda);
		}
		const float lambda=expf(2.0f*getL(0.5)-getL(1.0f-cosTheta));
		return 1.0f/(1.0f+lambda);
	}

public:

	Sheen(float alpha, Matrix m=Matrix(1.0f)): alpha(alpha) {
		initTransform(m);
		this->alpha=clamp(alpha, 1e-6f, 1.0f);
		initG1();
	}

	const char* getName() const  {
		return "Sheen   ";
	}
};

class STD: public Microsurface {
	float alphaX;
	float alphaY;
	float gamma;

	float getD(Vector h) const {
		const float cosTheta=h.z;
		if(cosTheta<=1e-3f) return 0.0f;
		const float cosTheta2=sqr(cosTheta);
		const float sinTheta2=1.0f-cosTheta2;
		float denominator=PI*alphaX*alphaY;
		if(sinTheta2>1e-6f) {
			const float tanTheta2=sinTheta2/cosTheta2;
			const float cosPhi2=sqr(h.x)/sinTheta2;
			const float sinPhi2=sqr(h.y)/sinTheta2;
			denominator*=sqr(cosTheta2)*powf(1.0f+tanTheta2*(cosPhi2/sqr(alphaX)+sinPhi2/sqr(alphaY))/(gamma-1.0f), gamma);
		}
		return 1.0f/denominator;
	}

	float getS1(float alpha, float mu) const {
		return alpha*powf(gamma-1.0f+sqr(mu/alpha), 1.5f-gamma)/mu;
	}

	float getF21(float z) const {
		return z*(1.066f+z*(2.655f+4.892f*z))/(1.038f+z*(2.969f+z*(4.305f+4.418f*z)));
	}

	float getF22() const {
		return (-14.402f+gamma*(-27.145f+gamma*(20.574f-2.745f*gamma)))/(-30.612f+gamma*(86.567f+gamma*(-84.341f+29.938f*gamma)));
	}

	float getF23() const {
		return (-129.404f+gamma*(324.987f+gamma*(-299.305f+93.268f*gamma)))/(-92.609f+gamma*(256.006f+gamma*(-245.663f+86.064f*gamma)));
	}

	float getF24(float z) const {
		return (6.537f+z*(6.074f+z*(-0.623f+5.223f*z)))/(6.538f+z*(6.103f+z*(-3.218f+6.347f*z)));
	}

	float getS2(float alpha, float mu) const {
		const float z=mu/alpha;
		return getF21(z)*(getF22()+getF23()*getF24(z));
	}

	float getLambda(float alpha, float mu) const {
		const float invSqrtPi=0.56418958354775628694807945156077f;
		return invSqrtPi*tgammaf(gamma-0.5f)*(powf(gamma-1.0f, gamma)*getS1(alpha, mu)/(2.0f*gamma-3.0f)+sqrtf(gamma-1.0f)*getS2(alpha, mu))/tgammaf(gamma)-0.5f;
	}

	float getG1(Vector dir) const {
		const float cosTheta = dir.z;
		if(cosTheta<=1e-3f) return 0.0f;
		if(cosTheta>=1.0f-1e-6f) return 1.0f;
		const float alpha=getAlpha(dir, alphaX, alphaY);
		const float cosTheta2=sqr(cosTheta);
		const float cotTheta=sqrtf(cosTheta2/(1.0f-cosTheta2));
		const float lambda=getLambda(alpha, cotTheta);
		return 1.0f/(1.0f+lambda);
	}

public:

	STD(float alphaX, float alphaY, float gamma, Matrix m=Matrix(1.0f)): alphaX(alphaX), alphaY(alphaY), gamma(gamma) {
		initTransform(m);
	}

	const char* getName() const  {
		return "STD     ";
	}
};

class Discrete: public Microsurface {
	float alpha;
	float gamma;
	float cosGamma;
	int count;
	Vector *normals;

	float getD(Vector h) const {
		if (h.z<=1e-3f) return 0.0f;
		int inside=0;
		for (int i=0; i<count; i++) {
			if (dotf(normals[i], h)>cosGamma) {
				inside++;
			}
		}
		const float condeSolidAngle=2.0f*PI*(1.0f-cosGamma);
		return float(inside)/(float(count)*condeSolidAngle*h.z);
	}

	float getG1(Vector dir) const {
		const float cosTheta = dir.z;
		if(cosTheta<=1e-3f) return 0.0f;
		if(cosTheta>=1.0f-1e-6f) return 1.0f;
		const float tanTheta2=(1.0f/sqr(cosTheta))-1.0f;
		const float denominator=1.0f+sqrtf(1.0f+tanTheta2*sqr(alpha));
		return 2.0f/denominator;
	}

	Vector sampleGGXNormal(float alpha, float u, float v) {
		float cosTheta2 = (1.0f-u)/(1.0f+(sqr(alpha)-1.0f)*u);
		float cosTheta = sqrtf(std::max(cosTheta2, 0.0f));
		float sinTheta = sqrtf(std::max(1.0f - cosTheta2, 0.0f));

		float phi = 2.f*PI*v;
		return Vector(cosf(phi)*sinTheta, sinf(phi)*sinTheta, cosTheta);
	}

public:

	Discrete(float alpha, float gamma, int count, Matrix m=Matrix(1.0f)): alpha(alpha), gamma(gamma), count(count) {
		initTransform(m);

		cosGamma=cosf(PI*gamma/180.0f);

		normals=new Vector[count];
		for (int i=0; i<count; i++) {
			const float u=rnd();
			const float v=rnd();
			normals[i]=sampleGGXNormal(alpha, u, v);
		}
	}

	~Discrete() {
		delete[] normals;
	}

	const char* getName() const  {
		return "Discrete";
	}
};