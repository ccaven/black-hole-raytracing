using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class DiskGenerator : MonoBehaviour {

    public int angleDensity = 256;
    public int radiusDensity = 64;
    public int verticalDensity = 32;

    public RenderTexture diskTexture;

    private void Awake () {

        // Initialize texture
        diskTexture = new RenderTexture(angleDensity, radiusDensity, 0, RenderTextureFormat.ARGBFloat) {
            volumeDepth = verticalDensity,
            dimension = UnityEngine.Rendering.TextureDimension.Tex3D,
            filterMode = FilterMode.Bilinear,
            enableRandomWrite = true
        };

        diskTexture.Create();

        print(diskTexture);
    }

    void Update () {

    }
}
