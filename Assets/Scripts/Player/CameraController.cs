using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CameraController : MonoBehaviour {

    public enum CameraPath {
        Manual,
        Introduction,
        Second
    }

    [Header("References")]
    [SerializeField] private Transform playerTransform;

    [Header("Camera Settings")]
    [SerializeField] private float sensitivity;
    [SerializeField] private int cameraPathID;
    [SerializeField] private float cameraSeconds;

    private float pitch;
    private float yaw;
    private Camera _camera;

    private float pathTime = 0;

    // Camera controlled by player
    private void PlayerControls (float time) {
        pitch -= Input.GetAxisRaw("Mouse Y") * sensitivity;
        yaw += Input.GetAxisRaw("Mouse X") * sensitivity;

        pitch = Mathf.Clamp(pitch, -90f, 90f);

        transform.position = playerTransform.position;

        transform.rotation = Quaternion.Euler(pitch, yaw, 0);
        playerTransform.rotation = Quaternion.Euler(0, yaw, 0);
    }

    readonly private Action<float, Transform>[] paths = {
        // First path
        // Start close and zoom out from black hole
        (t, transform) => {
            float r = 0.75f + t * 5.0f / 18.0f;
            transform.position = new Vector3(0, 0, -r);
        },

        // Second path
        // Start far and orbit closer
        (t, transform) => {
            float r = 10.0f * Mathf.Exp(-t * 0.1f);
            float a = t * 10.0f;
            transform.position = Quaternion.Euler(0, a, 0) * new Vector3(0, 0, -r);

            transform.rotation = Quaternion.Euler(0, a, 0);
        },

        // Third path
        // Move far away
        (t, transform) => {

            float d = 10.0f + 10.0f * t;

            d = Mathf.Min(d, 400.0f);

            float b = 5.0f;

            float FOV = Mathf.Atan(b / d);

            float a = t * 10.0f;

            transform.gameObject.GetComponent<Camera>().fieldOfView = FOV * Mathf.Rad2Deg;
            transform.position = Quaternion.Euler(0, a, 0) *  new Vector3(0, 0, -d);
            transform.rotation = Quaternion.Euler(0, a, 0);

        },

        // Fourth path
        // flyby
        (t, transform) => {

            Vector3 p1 = new Vector3(1, 0, -10f);
            Vector3 p2 = new Vector3(1, 0, 10f);

            transform.position = Vector3.Lerp(p1, p2, t / 10f);

        },

        // Final path
        // Falling into black hole
        // Start far then get close and turn backwards as you fall in
        (t, transform) => {

            // r 10 -> 0.75
            // a 0 -> 180

            float k = t / 10.0f;

            k = Mathf.Clamp(k, 0.0f, 1.0f);

            // Smoothstep
            k = k * k * (3 - 2 * k);
            k = k * k * (3 - 2 * k);

            float r = Mathf.Lerp(10.0f, 0.94f, k);

            float a = Mathf.Lerp(0.0f, 180.0f, k);

            float angleOffset = Mathf.Sin(k * Mathf.PI) * 90.0f;

            transform.position = Quaternion.Euler(0, a, 0) * new Vector3(0, 0, -r);
            transform.rotation = Quaternion.Euler(0, angleOffset, 0);

        },
    };

    private void Update () {
        paths[cameraPathID](pathTime, transform);
        pathTime += Time.deltaTime;
    }

    public void LockMouse () {
        Cursor.lockState = CursorLockMode.Locked;
    }

    private void Start () {
        pathTime = 0;
    }

    private void Awake () {
        _camera = GetComponent<Camera>();
    }

    
}