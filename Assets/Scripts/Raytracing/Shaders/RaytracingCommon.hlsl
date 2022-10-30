
/////////////////////////
/// Precompiler Stuff ///
/////////////////////////

#include "Noise.hlsl"
#include "HelperFunctions.hlsl"

#define RAYMARCH_STEPS 500
#define EPSILON 0.0001
#define MIN_STOP_DISTANCE 400.0
#define MAX_RAY_DISTANCE 1000.0

#define PI 3.1415926

////////////////////////////
/// Black hole variables ///
////////////////////////////

Texture2D _Skybox;
SamplerState sampler_Skybox;

Texture2D _AccretionDiskTexture;
SamplerState sampler_AccretionDiskTexture;

#define BLACK_HOLE_MASS 0.33
#define EVENT_HORIZON_RADIUS 3 * BLACK_HOLE_MASS

////////////////////////
/// Required Structs ///
////////////////////////

class Ray {
	float3 origin;
	float3 direction;
};

class IntegratedRay {
	float3 samplePosition;
	float3 sampleVelocity;
};

///////////////////////////////
/// Distance Estimator Base ///
///////////////////////////////

interface IDistanceEstimator {
	float GetDistance(float3 p);

	IDistanceEstimator GetFirst();
	IDistanceEstimator GetSecond();
	IDistanceEstimator GetThird();
};

class DistanceEstimatorBase : IDistanceEstimator {
	IDistanceEstimator GetFirst() {
		DistanceEstimatorBase base;
		return base;
	}

	IDistanceEstimator GetSecond() {
		DistanceEstimatorBase base;
		return base;
	}

	IDistanceEstimator GetThird() {
		DistanceEstimatorBase base;
		return base;
	}

	float GetDistance(float3 p) {
		return length(p) - 1;
	}
};

///////////////////////
/// Base primitives ///
///////////////////////

class BoxDE : DistanceEstimatorBase {
	float3 dimensions;

	float GetDistance(float3 p) {
		float3 q = abs(p) - dimensions;
		return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0);
	}

	static IDistanceEstimator New(float3 dimensions) {
		BoxDE box;
		box.dimensions = dimensions;
		return box;
	}
};

class BoxFrameDE : DistanceEstimatorBase {
	float3 dimensions;
	float r;

	float GetDistance(float3 p) {
		p = abs(p) - dimensions;
		float3 q = abs(p + r) - r;
		return min(min(
			length(max(float3(p.x, q.y, q.z), 0.0)) + min(max(p.x, max(q.y, q.z)), 0.0),
			length(max(float3(q.x, p.y, q.z), 0.0)) + min(max(q.x, max(p.y, q.z)), 0.0)),
			length(max(float3(q.x, q.y, p.z), 0.0)) + min(max(q.x, max(q.y, p.z)), 0.0));
	}

	static IDistanceEstimator New(float3 dimensions, float r) {
		BoxFrameDE box;
		box.dimensions = dimensions;
		box.r = r;
		return box;
	}
};

class SphereDE : DistanceEstimatorBase {
	float radius;

	float GetDistance(float3 p) {
		return length(p) - radius;
	}

	static IDistanceEstimator New(float radius) {
		SphereDE sphere;
		sphere.radius = radius;
		return sphere;
	}
};

class CapsuleDE : DistanceEstimatorBase {
	float3 a;
	float3 b;
	float r;

	float GetDistance(float3 p) {
		float3 pa = p - a, ba = b - a;
		float h = clamp(dot(pa, ba) / dot(ba, ba), 0.0, 1.0);
		return length(pa - ba * h) - r;
	}

	static IDistanceEstimator New(float3 a, float3 b, float radius) {
		CapsuleDE capsule;
		capsule.a = a;
		capsule.b = b;
		capsule.r = radius;
		return capsule;
	}
};

class PlaneDE : DistanceEstimatorBase {

	float3 normal;

	float GetDistance(float3 p) {
		return dot(p, normal);
	}

	static IDistanceEstimator New(float3 normal) {
		PlaneDE plane;
		plane.normal = normal;
		return plane;
	}

};

class TorusDE : DistanceEstimatorBase {

	float innerRadius;
	float outerRadius;

	float GetDistance(float3 p) {
		float2 q = float2(length(p.xz) - outerRadius, p.y);
		return length(q) - innerRadius;
	}

	static IDistanceEstimator New(float outerRadius, float innerRadius) {
		TorusDE torus;
		torus.innerRadius = innerRadius;
		torus.outerRadius = outerRadius;
		return torus;
	}

};

class RingDE : DistanceEstimatorBase {
	
	float innerRadius;
	float outerRadius;

	float angle;
	float ringSize;

	// Made more complicated to get a "real" SDF
	float GetDistance(float3 p) {
		float3 closestPoint = float3(p.x, 0, p.z);

		closestPoint.xz = normalize(closestPoint.xz);

		float2 a = float2(cos(angle), sin(angle));
		float threshold = cos(ringSize);

		float dp = dot(closestPoint.xz, a);
		dp = max(dp, threshold);

		float b = dot(closestPoint.xz, float2(-a.y, a.x)) > 0 ? 1 : -1;

		float na = angle + acos(dp) * b;

		closestPoint.xz = float2(cos(na), sin(na)) * outerRadius;

		return length(p - closestPoint) - innerRadius;
	}

	static IDistanceEstimator New(float outerRadius, float innerRadius, float angle, float ringSize) {
		RingDE ring;
		ring.outerRadius = outerRadius;
		ring.innerRadius = innerRadius;
		ring.angle = angle;
		ring.ringSize = ringSize;
		return ring;
	}

};

///////////////////////////////////
/// Fractal Distance Estimators ///
///////////////////////////////////

class JuliaDE : DistanceEstimatorBase {

	float4 offset;

	float GetDistance(float3 p) {

		float4 z = float4(p, 0.0);
		float md2 = 1.0;
		float mz2 = dot(z, z);

		for (int i = 0; i < 15; i++) {
			md2 *= mz2 * 4;
			z = qsqr(z) + offset;
			mz2 = dot(z, z);

			if (mz2 > 1000.0) break;
		}

		return 0.25 * sqrt(mz2 / md2) * log(mz2) - 0.0001;
	}

	static IDistanceEstimator New(float4 offset) {
		JuliaDE julia;
		julia.offset = offset;
		return julia;
	}

};

class CubicJuliaDE : DistanceEstimatorBase {

	float4 offset;

	float GetDistance(float3 p) {

		float4 z = float4(p, 0.0);
		float md2 = 1.0;
		float mz2 = dot(z, z);

		for (int i = 0; i < 15; i++) {
			md2 *= dot2(qsqr(z)) * 9;
			z = qcube(z) + offset;
			mz2 = dot(z, z);

			if (mz2 > 256.0) break;
		}

		return 0.25 * sqrt(mz2 / md2) * log(mz2) - 0.0001;
	}

	static IDistanceEstimator New(float4 offset) {
		CubicJuliaDE julia;
		julia.offset = offset;
		return julia;
	}

};

class MandelbulbDE : DistanceEstimatorBase {
	
	float GetDistance(float3 p) {
		float power = 3;

		float3 z = p;
		float dr = 1;
		float r;

		for (int i = 0; i < 15; i++) {
			r = length(z);
			if (r > 2) break;

			float theta = acos(z.z / r) * power;
			float phi = atan2(z.y, z.x) * power;
			float zr = pow(r, power);
			dr = pow(r, power - 1) * power * dr + 1;

			z = zr * float3(sin(theta) * cos(phi), sin(phi) * sin(theta), cos(theta));
			z += p;
		}

		return 0.5 * log(r) * r / dr;
	}

	static IDistanceEstimator New() {
		MandelbulbDE bulb;
		return bulb;
	}

};

////////////////////////////
/// Primitive Operations ///
////////////////////////////

class TransformOP : DistanceEstimatorBase {
	float3x3 m;

	float GetDistance(float3 p) {
		float leastScale = GetLeastScale(m);
		return GetFirst().GetDistance(mul(m, p)) * leastScale;
	}

	static IDistanceEstimator New(IDistanceEstimator primitive, float3x3 m) {
		class Local : TransformOP {
			IDistanceEstimator GetFirst() {
				return primitive;
			}
		} local;

		local.m = m;
		return local;
	}
};

class TwistOP : DistanceEstimatorBase {
	float k;

	float GetDistance(float3 p) {
		float c = cos(k * p.y);
		float s = sin(k * p.y);

		float2x2 m = float2x2(c, -s, s, c);

		float3 q = float3(mul(m, p.xz), p.y);

		return GetFirst().GetDistance(q) * 0.5f;
	}

	static IDistanceEstimator New(IDistanceEstimator primitive, float k) {
		class Local : TwistOP {
			IDistanceEstimator GetFirst() { 
				return primitive; 
			} 
		} local;
		local.k = k;
		return local;
	}
};

class BendOP : DistanceEstimatorBase {
	float k;

	float GetDistance(float3 p) {
		float c = cos(k * p.x);
		float s = sin(k * p.x);

		float2x2 m = { c, -s, s, c };

		float3 q = float3(mul(m, p.xy), p.z);

		return GetFirst().GetDistance(q) * 0.5f;
	}

	static IDistanceEstimator New(IDistanceEstimator primitive, float k) {
		class Local : BendOP {
			IDistanceEstimator GetFirst() {
				return primitive;
			}
		} local;
		local.k = k;
		return local;
	}
};

class UnionOP : DistanceEstimatorBase {

	float smoothingConstant;

	float GetDistance(float3 p) {
		float d1 = GetFirst().GetDistance(p);
		float d2 = GetSecond().GetDistance(p);

		float k = smoothingConstant;

		float h = clamp(0.5 + 0.5 * (d2 - d1) / k, 0.0, 1.0);
		return mix(d2, d1, h) - k * h * (1.0 - h);
	}

	static IDistanceEstimator New(IDistanceEstimator primitiveA, IDistanceEstimator primitiveB, float smoothingConstant) {		
		class Local : UnionOP {
			IDistanceEstimator GetFirst() {
				return primitiveA;
			}
			IDistanceEstimator GetSecond() {
				return primitiveB;
			}
		} local;
		local.smoothingConstant = smoothingConstant;
		return local;
	}

};

class SubtractionOP : DistanceEstimatorBase {

	float smoothingConstant;

	float GetDistance(float3 p) {
		float d1 = GetFirst().GetDistance(p);
		float d2 = GetSecond().GetDistance(p);

		float k = smoothingConstant;

		float h = clamp(0.5 - 0.5 * (d2 + d1) / k, 0.0, 1.0);
		return mix(d2, -d1, h) + k * h * (1.0 - h);
	}

	static IDistanceEstimator New(IDistanceEstimator primitiveA, IDistanceEstimator primitiveB, float smoothingConstant) {		
		class Local : SubtractionOP {
			IDistanceEstimator GetFirst() {
				return primitiveA;
			}
			IDistanceEstimator GetSecond() {
				return primitiveB;
			}
		} local;
		local.smoothingConstant = smoothingConstant;
		return local;
	}

};

class SubtractionRigidOP : DistanceEstimatorBase {
	float GetDistance(float3 p) {
		float d1 = GetFirst().GetDistance(p);
		float d2 = GetSecond().GetDistance(p);
		return max(d2, -d1);
	}

	static IDistanceEstimator New(IDistanceEstimator primitiveA, IDistanceEstimator primitiveB) {
		class Local : SubtractionRigidOP {
			IDistanceEstimator GetFirst() {
				return primitiveA;
			}
			IDistanceEstimator GetSecond() {
				return primitiveB;
			}
		} local;
		return local;
	}
};

class IntersectionOP : DistanceEstimatorBase {

	float smoothingConstant;

	float GetDistance(float3 p) {
		float d1 = GetFirst().GetDistance(p);
		float d2 = GetSecond().GetDistance(p);

		float k = smoothingConstant;

		float h = clamp(0.5 - 0.5 * (d2 - d1) / k, 0.0, 1.0);
		return mix(d2, d1, h) + k * h * (1.0 - h);
	}

	static IDistanceEstimator New(IDistanceEstimator primitiveA, IDistanceEstimator primitiveB, float smoothingConstant) {
		class Local : IntersectionOP {
			IDistanceEstimator GetFirst() {
				return primitiveA;
			}
			IDistanceEstimator GetSecond() {
				return primitiveB;
			}
		} local;
		local.smoothingConstant = smoothingConstant;
		return local;
	}

};

class IntersectionRigidOP : DistanceEstimatorBase {
	float GetDistance(float3 p) {
		float d1 = GetFirst().GetDistance(p);
		float d2 = GetSecond().GetDistance(p);
		return max(d1, d2);
	}

	static IDistanceEstimator New(IDistanceEstimator primitiveA, IDistanceEstimator primitiveB) {
		class Local : IntersectionRigidOP {
			IDistanceEstimator GetFirst() {
				return primitiveA;
			}
			IDistanceEstimator GetSecond() {
				return primitiveB;
			}
		} local;
		return local;
	}
};

class RoundOP : DistanceEstimatorBase {
	float r;

	float GetDistance(float3 p) {
		return GetFirst().GetDistance(p) - r;
	}

	static IDistanceEstimator New(IDistanceEstimator primitive, float r) {
		class Local : RoundOP {
			IDistanceEstimator GetFirst() {
				return primitive;
			}
		} local;

		local.r = r;

		return local;
	}
};

class DisplaceOP : DistanceEstimatorBase {
	float3 displacement;

	float GetDistance(float3 p) {
		return GetFirst().GetDistance(p - displacement);
	}

	static IDistanceEstimator New(IDistanceEstimator primitive, float3 displacement) {
		class Local : DisplaceOP {
			IDistanceEstimator GetFirst() {
				return primitive;
			}
		} local;

		local.displacement = displacement;

		return local;
	}
};

class ModOP : DistanceEstimatorBase {
	float3 modulus;

	float GetDistance(float3 p) {
		return GetFirst().GetDistance(mod(p, modulus));
	}

	static IDistanceEstimator New(IDistanceEstimator primitive, float3 modulus) {
		class Local : ModOP {
			IDistanceEstimator GetFirst() {
				return primitive;
			}
		} local;

		local.modulus = modulus;

		return local;
	}
};

class RepeatOP : DistanceEstimatorBase {

    float3 modulus;

    float GetDistance(float3 p) {
        // Check a bunch
        float d = 10000.0;

        float3 lp = mod(p - modulus * 0.5, modulus);

        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                for (int dz = -1; dz <= 1; dz++) {
                    // Todo: replace with smoothmin
                    d = min(d, GetFirst().GetDistance(lp + modulus * float3(dx, dy, dz)));
                }
            }
        }

        return d;
    }

    static IDistanceEstimator New(IDistanceEstimator primitive, float3 modulus) {
        class Local : ModOP {
            IDistanceEstimator GetFirst() {
                return primitive;
            }
        } local;

        local.modulus = modulus;

        return local;
    }

};

class HorizontalRepeatOP : DistanceEstimatorBase {
	float2 modulus;

	float GetDistance(float3 p) {

		float2 h = mod(p.xz, modulus);

		return GetFirst().GetDistance(float3(h.x, p.y, h.y));
	}

	static IDistanceEstimator New(IDistanceEstimator primitive, float2 modulus) {
		class Local : HorizontalRepeatOP {
			IDistanceEstimator GetFirst() {
				return primitive;
			}
		} local;

		local.modulus = modulus;

		return local;
	}
};

class SolidifyOP : DistanceEstimatorBase {
	float GetDistance(float3 p) {
		return GetFirst().GetDistance(max(p, 0));
	}

	static IDistanceEstimator New(IDistanceEstimator primitive) {
		class Local : SolidifyOP {
			IDistanceEstimator GetFirst() {
				return primitive;
			}
		} local;
		return local;
	}
};

class InvertOP : DistanceEstimatorBase {
	float GetDistance(float3 p) {
		return -GetFirst().GetDistance(p);
	}

	static IDistanceEstimator New(IDistanceEstimator primitive) {
		class Local : InvertOP {
			IDistanceEstimator GetFirst() {
				return primitive;
			}
		} local;
		return local;
	}
};

//////////////////////////////////////
/// Black hole specific Primitives ///
//////////////////////////////////////

class AccretionDiskDE : DistanceEstimatorBase {
	float GetDistance(float3 p) {
		
		float3 closestPoint = p;

		float m = length(p);

		m = min(m, 4.0f);
		m = max(m, 1.4f);

		closestPoint = m * normalize(closestPoint);

		closestPoint.y = min(closestPoint.y, 0.01f);
		closestPoint.y = max(closestPoint.y, -0.01f);

		return 0.5 * length(p - closestPoint);
	}

	static IDistanceEstimator New() {
		AccretionDiskDE disk;
		return disk;
	}
};

class AccretionDiskParticlesDE : DistanceEstimatorBase {


	float GetDistance(float3 p) {		
		IDistanceEstimator ring = RingDE::New(2.5f, 0.1f, 0.0f, PI / 2);

		float d = ring.GetDistance(p);

        for (int i = 0; i < 50; i++) {

            float outer = 2.0f + cos(float(i)) * 0.5f;
            float inner = 0.05f;

            float ang = i / 10.0 * PI * 2;

            float wid = PI / 8;

            IDistanceEstimator newRing = RingDE::New(outer, inner, ang, wid);

            float h = cos(float(i) * 2) * 0.1f;

            d = min(d, newRing.GetDistance(p + float3(0, h, 0)));

        }

		return d;
	}

	static IDistanceEstimator New() {
		AccretionDiskParticlesDE disk;
		return disk;
	}

};

////////////////////////////////////////////////
/// Different Kinds of Rays and Integrations ///
////////////////////////////////////////////////

interface Integration {
	void Integrate(inout IntegratedRay ray, float timestep);
};

class SimpleIntegration : Integration {
	void Integrate(inout IntegratedRay ray, float timestep) {
		ray.samplePosition += ray.sampleVelocity * timestep;
	}
};

float3 BlackHoleAcceleration(float3 p) {
	// Black hole constant
	const float k = 1.0f;

	return -1.5 * k * k * p * pow(length(p), -6);
}

class BlackHoleEuler : Integration {
	void Step(inout IntegratedRay ray, float dt) {
		dt /= length(ray.sampleVelocity);

		float3 a = BlackHoleAcceleration(ray.samplePosition);
		ray.samplePosition += ray.sampleVelocity * dt;
		ray.sampleVelocity += a * dt;
	}
	void Integrate(inout IntegratedRay ray, float dt) {
		for (int i = 0; i < 10; i ++) Step(ray, dt / 10);
	}
};

class BlackHoleEulerMidpoint : Integration {
	void Step(inout IntegratedRay ray, float dt) {
		dt /= length(ray.sampleVelocity);
		float3 a = BlackHoleAcceleration(ray.samplePosition);
		ray.samplePosition += ray.sampleVelocity * dt;
		a = a / 2 + BlackHoleAcceleration(ray.samplePosition) / 2;
		ray.sampleVelocity += a * dt;
    }

	void Integrate(inout IntegratedRay ray, float dt) {
		const int numSteps = 10;
		for (int i = 0; i < numSteps; i++) Step(ray, dt / numSteps);
	}
};

class BlackHoleImplicitEuler : Integration {
	void Step(inout IntegratedRay ray, float dt) {
		float3 a = BlackHoleAcceleration(ray.samplePosition);
		float3 predictedVelocity = a * dt;
		dt /= max(1, length(predictedVelocity));
		ray.sampleVelocity += a * dt;
		ray.samplePosition += ray.sampleVelocity * dt;
	}
	void Integrate(inout IntegratedRay ray, float dt) {
		for (int i = 0; i < 20; i ++) Step(ray, dt / 20);
	}
};

class BlackHoleRK4 : Integration {
	void Evaluate(float3 ix, float3 iv, float dt, float3 dx, float3 dv, out float3 odx, out float3 odv) {
		float3 sx = ix + dx * dt;
		float3 sv = iv + dv * dt;
		odx = sv;
		odv = BlackHoleAcceleration(sx);
	}

	void Integrate(inout IntegratedRay ray, float dt) {
		float3 sx = ray.samplePosition;
		float3 sv = ray.sampleVelocity;

		dt /= length(sv);

		float3 ax, av, bx, bv, cx, cv, dx, dv;

		Evaluate(sx, sv, 0, 
			float3(0, 0, 0), float3(0, 0, 0), 
			ax, av);
		
		Evaluate(sx, sv, dt * 0.5f, 
			ax, av, 
			bx, bv);

		Evaluate(sx, sv, dt * 0.5f, 
			bx, bv, 
			cx, cv);

		Evaluate(sx, sv, dt, 
			cx, cv, 
			dx, dv);

		float3 dxdt = (ax + 2.0f * (bx + cx) + dx) / 6;
		float3 dvdt = (av + 2.0f * (bv + cv) + dv) / 6;

		ray.samplePosition += dxdt * dt;
		ray.sampleVelocity += dvdt * dt;
	}

};

//////////////////////////////////////////
/// Different kinds of skybox sampling ///
//////////////////////////////////////////

interface SkyboxSampler {
    float3 GetSkybox(IntegratedRay ray);
};

class BasicSkybox : SkyboxSampler {

    float rampBrightness(float b) {
        return b * b * (3 - 2 * b);
    }

    float3 GetSkybox(IntegratedRay ray) {
        float3 endDirection = normalize(ray.sampleVelocity);

        float x = atan2(endDirection.x, endDirection.z) / PI * 0.5 + 0.5;
        float y = acos(endDirection.y) / PI;

        float3 c = _Skybox.SampleLevel(sampler_Skybox, float2(x, y), 0).rgb;

        // float b = max(c.r, max(c.g, c.b));
        // b = sqrt(b);

        return c;
    }
};

class AdvancedSkybox : SkyboxSampler {
    float3 GetSkybox(IntegratedRay ray) {
        float3 ro = ray.samplePosition;
        float3 rd = normalize(ray.sampleVelocity);

        float x = atan2(rd.x, rd.z) / PI * 0.5 + 0.5;
        float y = acos(rd.y) / PI;

        float3 base = _Skybox.SampleLevel(sampler_Skybox, float2(x, y), 0).rgb;

        // float3 base2 = _Skybox.SampleLevel(sampler_Skybox, float2(x * 5, clamp(y * 5, -1, 1)), 0).rgb;

        float n = snoise(float3(x, y, 0) * 2000.0);
        n = n * n * (3 - 2 * n);
        n *= exp(-y * y);
        n = 0.7 + 0.3 * n;
        return n * base;
    }
};

class DanielsSkybox : SkyboxSampler {
    float3 GetSkybox(IntegratedRay ray) {

        float3 ro = ray.samplePosition;
        float3 rd = normalize(ray.sampleVelocity);

        float3 absRd = abs(rd);
        float maxComp = max(max(absRd.x, absRd.y), absRd.z);

        float2 uv = maxComp == absRd.x ? rd.yz / absRd.x :
            (maxComp == absRd.y ? rd.xz / absRd.y : rd.xy / absRd.z);

        uv = uv * 0.5 + 0.5;

        float3 col = float3(0, 0, 0);

        // Main Background
        col += fnoise(rd * 0.25) * float3(0.7, 0.0, 0.4);
        col += fnoise(rd * 0.25 + float3(0.1, 0.2, -0.3)) * float3(0.0, 0.2, 0.4);
        col += fnoise(rd * 0.25 + float3(-0.4, 0.5, 0.6)) * float3(0.1, 0.2, 0.1);

        col *= 0.8;

        // Highlights
        col += 3.0 * pow(fnoise(rd * 0.5 + float3(0.7, -0.8, 0.9)), 3.0) * float3(0.3, 0.0, 0.5);

        // Stars
        float stars = fnoise(rd * 3.0);
        stars = smoothstep(0.6, 0.8, stars);

        float starMask = pow(dot(1.0 - col, float3(0.333, 0.333, 0.333)), 8.0);
        col += starMask * stars * 1.1;

        return col;
    }
};

//////////////
/// Scenes ///
//////////////

interface IScene {
	float GetDistance(float3 p);
	float3 GetColor(float3 p);
	float3 GetNormal(float3 p);
};


class BlackHoleScene : IScene {

	float GetDistance(float3 p) {
		IDistanceEstimator eventHorizon = SphereDE::New(0.8f);

		return eventHorizon.GetDistance(p);
	}

	float3 GetColor(float3 p) {
		// Event horizon
		return float3(0.0, 0.0, 0.0);
	}

	float3 GetNormal(float3 p) {
		const float h = EPSILON;
		const float2 k = float2(1, -1);
		return normalize(k.xyy * GetDistance(p + k.xyy * h) +
			k.yyx * GetDistance(p + k.yyx * h) +
			k.yxy * GetDistance(p + k.yxy * h) +
			k.xxx * GetDistance(p + k.xxx * h));
	}

};

class FractalScene : IScene {

    float GetDistance(float3 p) {

        float4 c = float4(-2, 6, 15, -6) / 22.0;

        float4 z = float4(p, 0.0);
        float md2 = 1.0;
        float mz2 = dot(z, z);

        float o = 1e10;

        for (int i = 0; i < 15; i++) {
            md2 *= dot2(qsqr(z)) * 9;
            z = qcube(z) + c;
            mz2 = dot(z, z);

            o = min(o, length(z.xz - float2(0.45, 0.55)) - 0.1);

            if (mz2 > 256.0) break;
        }

        float d = 0.25 * sqrt(mz2 / md2) * log(mz2) - 0.0001;

        d = min(o, d);

        d = max(d, p.y);

        return min(d, 0.5);

    }

    float3 GetColor(float3 p) {
        float4 c = float4(-2, 6, 15, -6) / 22.0;

        float4 z = float4(p, 0.0);
        float md2 = 1.0;
        float mz2 = dot(z, z);

        float md = 1.0;

        for (int i = 0; i < 15; i++) {
            md2 *= dot2(qsqr(z)) * 9;
            z = qcube(z) + c;
            mz2 = dot(z, z);

            md = min(md, length(z.xz - float2(0.45, 0.55)) - 0.1);

            if (mz2 > 256.0) break;
        }

        float3 col = 0.5 + 0.5 * cos(log2(md) * 0.9 + float3(0.0, 0.6, 1.0) + 3.5);

        if (p.y > 0.0) col = mix(col, float3(1, 1, 1), 0.2);

        float inside = smoothstep(14.0, 15.0, md);

        col *= float3(0.45, 0.42, 0.40) + float3(0.55, 0.58, 0.60) * inside;

        col = mix(col * col * (3.0 - 2.0 * col), col, inside);

        /*
        float a = md;
        return 0.5 + 0.5 * float3(
            cos(a + 3.1415 / 3 * 0),
            cos(a + 3.1415 / 3 * 2),
            cos(a + 3.1415 / 3 * 4)
        );
        */
        return col;
    }

    float3 GetNormal(float3 p) {

        const float h = EPSILON;
        const float2 k = float2(1, -1);
        return normalize(k.xyy * GetDistance(p + k.xyy * h) +
            k.yyx * GetDistance(p + k.yyx * h) +
            k.yxy * GetDistance(p + k.yxy * h) +
            k.xxx * GetDistance(p + k.xxx * h));

    }

};

class TestScene : IScene {

	float GetDistance(float3 p) {
		float3x3 transformation = {
			1, 1, 0,
			0, -1, 0,
			0, 1, 1
		};

		IDistanceEstimator shape1 = BoxDE::New(float3(4, 4, 4));

		IDistanceEstimator shape2 = SphereDE::New(5.0f);

		IDistanceEstimator shape3 = BoxFrameDE::New(float3(5, 2, 2), 0.3f);
		shape3 = RoundOP::New(shape3, 0.1);

		IDistanceEstimator combined = IntersectionOP::New(shape2, shape1, 0.1f);

		IDistanceEstimator twisted = TwistOP::New(combined, 0.1);

		IDistanceEstimator bent = BendOP::New(twisted, 0.1);

		IDistanceEstimator transformed = TransformOP::New(bent, transformation);

		IDistanceEstimator final = UnionOP::New(transformed, shape3, 0.1f);

		return final.GetDistance(p);
	}

	float3 GetColor(float3 p) {
		return float3(1, 1, 1);
	}

	float3 GetNormal(float3 p) {
		const float h = 0.0001f;
		const float2 k = float2(1, -1);
		return normalize(k.xyy * GetDistance(p + k.xyy * h) +
			k.yyx * GetDistance(p + k.yyx * h) +
			k.yxy * GetDistance(p + k.yxy * h) +
			k.xxx * GetDistance(p + k.xxx * h));
	}
};