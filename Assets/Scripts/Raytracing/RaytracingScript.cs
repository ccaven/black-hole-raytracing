using UnityEngine;

[RequireComponent(typeof(Camera))]
public class RaytracingScript : MonoBehaviour {

    [Header("References")]
    [SerializeField] private ComputeShader raytracingShader;

    [Header("Raymarch Settings")]
    [SerializeField] private Texture2D skybox;
    [SerializeField] private Texture2D accretionDiskTexture;

    [Header("Bloom Settings")]
    [SerializeField] private bool bloomEnabled;
    [SerializeField] [Range(0, 2)] private float bloomThreshold;
    [SerializeField] [Range(0, 2)] private float bloomIntensity;

    private RenderTexture _temp;

    private Camera _camera;

    private void Awake () {
        _camera = GetComponent<Camera>();

        /*
        _temp = new RenderTexture(Screen.width, Screen.height, 0, RenderTextureFormat.ARGBFloat, RenderTextureReadWrite.Linear) {
            enableRandomWrite = true
        };
        _temp.Create();

        raytracingShader.SetTexture(0, "_Result", _temp);
        */

        raytracingShader.SetTexture(0, "_AccretionDiskTexture", accretionDiskTexture);
    }

    private void InitializeTexture () {
        _temp = new RenderTexture(Screen.width, Screen.height, 0, RenderTextureFormat.ARGBFloat, RenderTextureReadWrite.Linear) {
            enableRandomWrite = true
        };
        _temp.Create();

        raytracingShader.SetTexture(0, "_Result", _temp);
    }

    private void Start () {
        raytracingShader.SetTexture(0, "_Skybox", skybox);
    }

    private void OnRenderImage ( RenderTexture source, RenderTexture destination ) {

        if ( _temp == null ) InitializeTexture();

        raytracingShader.SetMatrix("_CameraToWorld", _camera.cameraToWorldMatrix);
        raytracingShader.SetMatrix("_CameraInverseProjection", _camera.projectionMatrix.inverse);

        raytracingShader.Dispatch(0,
            Mathf.CeilToInt(Screen.width / 8f),
            Mathf.CeilToInt(Screen.height / 8f), 1);

        if ( bloomEnabled ) {
            BloomMaster.threshold = bloomThreshold;
            BloomMaster.intensity = bloomIntensity;
            BloomMaster.Bloom(_temp, destination);
        }
        else Graphics.Blit(_temp, destination);
    }
}
