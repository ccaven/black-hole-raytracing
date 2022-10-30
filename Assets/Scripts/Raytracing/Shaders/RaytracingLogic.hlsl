
#define ANTI_ALIASING 0

///////////////////////
/// Input variables ///
///////////////////////

RWTexture2D<float4> _Result;

float4x4 _CameraToWorld;
float4x4 _CameraInverseProjection;

//////////////////////////////
/// Global scene variables ///
//////////////////////////////

BlackHoleScene scene;
BlackHoleImplicitEuler integration;
BasicSkybox skybox;

float MAX_STEP_SIZE = 100.0f;
float CLIP_DISTANCE = 500.0f;

/////////////////////
/// Ray functions ///
/////////////////////

static inline Ray CreateRay(float3 origin, float3 direction) {
	Ray ray;
	ray.origin = origin;
	ray.direction = direction;
	return ray;
}

static inline Ray CreateCameraRay(float2 uv) {
	
	float3 origin = mul(_CameraToWorld, float4(0, 0, 0, 1)).xyz;
	float3 direction = mul(_CameraInverseProjection, float4(uv, 0, 1)).xyz;
	direction = mul(_CameraToWorld, float4(direction, 0)).xyz;
	direction = normalize(direction);
	return CreateRay(origin, direction);
}

static inline IntegratedRay CreateIntegratedRay(Ray ray) {
	IntegratedRay integratedRay;

	integratedRay.samplePosition = ray.origin;
	integratedRay.sampleVelocity = ray.direction;

	return integratedRay;
}

static inline float3 GetSkybox(float3 endDirection) {
	float x = atan2(endDirection.x, endDirection.z) / 3.1415926 * 0.5 + 0.5;
	float y = acos(endDirection.y) / 3.1415926;

	return _Skybox.SampleLevel(sampler_Skybox, float2(x, y), 0).rgb;
}


//////////////////////////////
/// Main raymarch function ///
//////////////////////////////

static inline float4 Raymarch(Ray ray, IScene scene, Integration integration, SkyboxSampler skybox) {
	
	IntegratedRay integratedRay = CreateIntegratedRay(ray);

	int i = 0;
	
	float distanceToScene = 10;

	for (; i < RAYMARCH_STEPS && distanceToScene > EPSILON && distanceToScene < MIN_STOP_DISTANCE; i++) {
		distanceToScene = scene.GetDistance(integratedRay.samplePosition);
		integration.Integrate(integratedRay, distanceToScene);
	}

    float4 col = float4(0, 0, 0, 0);

	if (i > RAYMARCH_STEPS - 2 || distanceToScene >= MIN_STOP_DISTANCE)
		col = float4(skybox.GetSkybox(integratedRay), 1);
    else
	    col.rgb = scene.GetColor(integratedRay.samplePosition);

	return col;

}

/////////////////////////////////////
/// Compute kernel for raytracing ///
/////////////////////////////////////

[numthreads(8, 8, 1)]
void CSMain(uint3 id : SV_DispatchThreadID) {
	uint width, height;
	_Result.GetDimensions(width, height);

#if ANTI_ALIASING == 0
	float2 uv0 = float2(id.xy + 0.25f) / float2(width, height) * 2.0f - 1.0f;
	float2 uv1 = float2(id.xy + float2(0.25f, 0.75f)) / float2(width, height) * 2.0f - 1.0f;
	float2 uv2 = float2(id.xy + float2(0.75f, 0.25f)) / float2(width, height) * 2.0f - 1.0f;
	float2 uv3 = float2(id.xy + 0.75f) / float2(width, height) * 2.0f - 1.0f;

	Ray ray0 = CreateCameraRay(uv0);
	Ray ray1 = CreateCameraRay(uv1);
	Ray ray2 = CreateCameraRay(uv2);
	Ray ray3 = CreateCameraRay(uv3);

	float4 r0 = Raymarch(ray0, scene, integration, skybox);
	float4 r1 = Raymarch(ray1, scene, integration, skybox);
	float4 r2 = Raymarch(ray2, scene, integration, skybox);
	float4 r3 = Raymarch(ray3, scene, integration, skybox);

    _Result[id.xy] = (r0 + r1 + r2 + r3) / 4;
#elif ANTI_ALIASING == 1
    // Antialiasing
    float2 uv0 = float2(id.xy + 0.25f) / float2(width, height) * 2.0f - 1.0f;
    float2 uv3 = float2(id.xy + 0.75f) / float2(width, height) * 2.0f - 1.0f;

    Ray ray0 = CreateCameraRay(uv0);
    Ray ray3 = CreateCameraRay(uv3);

    float4 r0 = Raymarch(ray0, scene, integration, skybox);
    float4 r3 = Raymarch(ray3, scene, integration, skybox);

    _Result[id.xy] = (r0 + r3) / 2;

#else    
    // No antialiasing
    _Result[id.xy] = Raymarch(CreateCameraRay(float2(id.xy + 0.5f) / float2(width, height) * 2.0f - 1.0f), scene, integration, skybox);
#endif
}

