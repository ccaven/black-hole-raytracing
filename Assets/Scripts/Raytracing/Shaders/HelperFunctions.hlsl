////////////////////////
/// Helper Functions ///
////////////////////////

float dot2(in float2 v) {
    return dot(v, v);
}
float dot2(in float3 v) {
    return dot(v, v);
}
float dot2(in float4 v) {
    return dot(v, v);
}

float mod(in float a, in float b) {
    return fmod(fmod(a, b) + b, b);
}
float2 mod(in float2 a, in float2 b) {
    return fmod(fmod(a, b) + b, b);
}
float3 mod(in float3 a, in float3 b) {
    return fmod(fmod(a, b) + b, b);
}

float ndot(float2 a, float2 b) {
    return a.x * b.x - a.y * b.y;
}

float mix(float a, float b, float k) {
    return a + (b - a) * k;
}
float3 mix(float3 a, float3 b, float k) {
    return a + (b - a) * k;
}

float GetLeastScale(float3x3 m) {

    const float SQRT2 = sqrt(2);
    const float SQRT3 = sqrt(3);

    float3 dx = m._11_12_13;
    float3 dy = m._21_22_23;
    float3 dz = m._31_32_33;

    float d = dot2(dx);

    d = min(d, dot2(dy));
    d = min(d, dot2(dz));

    d = min(d, dot2(dx + dy) / SQRT2);
    d = min(d, dot2(dy + dz) / SQRT2);
    d = min(d, dot2(dz + dx) / SQRT2);

    d = min(d, dot2(dx + dy + dz) / SQRT3);

    return sqrt(d);
}

/////////////////////////////////
/// Complex Number Operations ///
/////////////////////////////////

float4 qsqr(float4 a) {
    return float4(a.x * a.x - a.y * a.y - a.z * a.z - a.w * a.w, 2.0 * a.x * a.y, 2.0 * a.x * a.z, 2.0 * a.x * a.w);
}
float4 qcube(in float4 q) {
    float4  q2 = q * q;
    return float4(q.x * (q2.x - 3.0 * q2.y - 3.0 * q2.z - 3.0 * q2.w),
        q.yzw * (3.0 * q2.x - q2.y - q2.z - q2.w));
}
float4 qmul(float4 a, float4 b) {
    return float4(b.xyz * a.w + a.xyz * b.w + cross(a.xyz, b.xyz), a.w * b.w - dot(a.xyz, b.xyz));
}
float4 qinv(float4 a) {
    return float4(a.x, -a.yzw) / dot2(a);
}
float4 qdiv(float4 a, float4 b) {
    return qmul(a, qinv(b));
}
float4 qexp(float4 a) {
    float t = exp(a.x);
    float l = length(a.yzw);
    return float4(t * cos(l), t * a.yzw * sin(l) / l);
}
float4 qln(float4 a) {
    float l = length(a);
    float t = log(l);
    return float4(t, a.yzw * acos(a.x / l) / length(a.yzw));
}
float4 qpow(float4 a, float n) {
    return qexp(n * qln(a));
}
