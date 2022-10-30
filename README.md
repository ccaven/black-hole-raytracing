# black-hole-raytracing
Geodesic raytracing of Schwarzschild black holes.

Built with [Unity](https://unity.com/download)\
Version 2020.3.11f1

## Scene Abstraction
```hlsl
class TestScene : IScene {

  float GetDistance(float3 p) {
    float3x3 transformation = {
      1, 1, 0,
      0, -1, 0,
      0, 1, 1
    };

    IDistanceEstimator shape1 = BoxDE::New(float3(4, 4, 4));

    IDistanceEstimator shape2 = SphereDE::New(5.0f);

    IDistanceEstimator shape3 = RoundOP::New(BoxFrameDE::New(float3(5, 2, 2), 0.3f), 0.1);

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
```

## Cartesian Schwarzschild approximation
```hlsl
float3 BlackHoleAcceleration(float3 p, float k) {
  return -1.5 * k * k * p * pow(length(p), -6);
}
```
