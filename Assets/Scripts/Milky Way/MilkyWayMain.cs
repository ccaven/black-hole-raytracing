using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class MilkyWayMain : MonoBehaviour {

    [Header("References")]
    [SerializeField] private ComputeShader noiseShader;
    [SerializeField] private Transform cameraTransform;

    public RenderTexture milkyWayTexture;

    private Camera _camera;

    private void Awake () {

        _camera = cameraTransform.GetComponent<Camera>();

    }


}