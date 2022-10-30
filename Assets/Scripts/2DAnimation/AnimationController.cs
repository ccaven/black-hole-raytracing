using System;
using System.Linq;
using System.Collections.Generic;
using UnityEngine;

[Serializable]
public struct ColorSettings {
    public Color backgroundColor, gridColor, sphereColor, v_color, e0_color, e1_color;
};

[Serializable]
public struct GridSettings {
    public int gridWidth, gridHeight, gridPrecision;
    public float gridThickness, gridSpacing;
};

public class AnimationController : MonoBehaviour {

    // Camera reference
    [Header("References")]
    public new Camera camera;
    public GameObject background;
    public GameObject sphere;

    // Play percentage
    [Header("Play options")]
    [Range(0, 1)]
    public float time;
    public bool isPlaying;

    [Header("Settings")]
    public GridSettings gridSettings;
    public ColorSettings colorSettings;

    // Materials of each line
    private Material gridLineMaterial, backgroundMaterial, vMaterial, e0Material, e1Material, sphereMaterial;

    // The lines we are going to change
    private LineRenderer vLine, e0Line, e1Line;

    // Initialize objects
    private void Awake () {
        // Initialize materials
        gridLineMaterial = new Material(Shader.Find("Unlit/Color")) { color = colorSettings.gridColor };
        backgroundMaterial = new Material(Shader.Find("Unlit/Color")) { color = colorSettings.backgroundColor };
        vMaterial = new Material(Shader.Find("Unlit/Color")) { color = colorSettings.v_color };
        e0Material = new Material(Shader.Find("Unlit/Color")) { color = colorSettings.e0_color };
        e1Material = new Material(Shader.Find("Unlit/Color")) { color = colorSettings.e1_color };
        sphereMaterial = new Material(Shader.Find("Unlit/Color")) { color = colorSettings.sphereColor };

        // Correct background
        background.GetComponent<MeshRenderer>().material = backgroundMaterial;
        background.transform.localScale = new Vector3(16f / 9f, 1, 1.0f);

        sphere.GetComponent<MeshRenderer>().material = sphereMaterial;
    }

    // Add a line with multiple positions
    private LineRenderer AddLine ( string name, Vector3[] positions, Material material, float thickness, float z) {
        GameObject obj = new GameObject(name);

        obj.transform.parent = transform;

        LineRenderer renderer = obj.AddComponent<LineRenderer>();

        renderer.numCornerVertices = 6;
        renderer.numCapVertices = 6;
        renderer.widthMultiplier = thickness;

        renderer.positionCount = positions.Length;

        obj.transform.position = Vector3.back * z;
        renderer.SetPositions(positions);
        //renderer.SetPositions(positions.Select(v => new Vector3(v.x, v.y, z)).ToArray());

        renderer.material = material;

        return renderer;
    }

    // Add a line with only two positions
    private LineRenderer AddLine ( string name, Vector3 start, Vector3 end, Material material, float thickness, float z) {
        return AddLine(name, new Vector3[] { start, end }, material, thickness, z);
    }

    // Calculate the displacement of any point on the grid
    private Vector2 Displacement ( Vector2 original ) {
        const float gs = 0.5f;
        const float wf = 0.5f;

        return original + gs * new Vector2(
            Mathf.Sin(original.y * wf),
            Mathf.Sin(original.x * wf));
    }

    private Vector2 Derivative ( Func<Vector2, Vector2> f, Vector2 o, Vector2 d ) {
        const float h = 0.001f;
        return (f(o + d * h) - f(o)) / h;
    }

    private Func<Vector2, Vector2> Derivative ( Func<Vector2, Vector2> f, Vector2 d ) {
        const float h = 0.001f;
        return ( Vector2 i ) => {
            return (f(i + d * h) - f(i)) / h;
        };
    }

    private Func<Vector2, Vector2> NDerivative ( Func<Vector2, Vector2> f, Vector2[] ds ) {
        Func<Vector2, Vector2>[] fs = new Func<Vector2, Vector2>[ds.Length];

        for ( int i = 0; i < ds.Length; i++ ) {
            if ( i == 0 ) fs[0] = Derivative(f, ds[0]);
            else fs[i] = Derivative(fs[i - 1], ds[i]);
        }

        return fs[fs.Length - 1];
    }

    private Vector2 Inverse ( Func<Vector2, Vector2> f, Vector2 o ) {

        // expand f(x) as a taylor series

        return o;
    }

    private Vector2 InverseSpecific ( Vector2 o ) {

        // Use newtons method

        return o;

    }


    // Generate the positions of each gridline
    private void GenerateGridLines () {
        // Generate Y gridlines
        for ( int x = 0; x <= gridSettings.gridWidth; x++ ) {
            Vector3[] p = new Vector3[gridSettings.gridPrecision + 1];

            // Compute each position
            for ( int iy = 0; iy <= gridSettings.gridPrecision; iy++ ) {
                float y = (float)iy / gridSettings.gridPrecision * gridSettings.gridHeight;

                p[iy] = Displacement(gridSettings.gridSpacing * new Vector2(x, y)) -
                    gridSettings.gridSpacing * new Vector2(gridSettings.gridWidth, gridSettings.gridHeight) / 2;
            }
            AddLine("Gridline Y", p, gridLineMaterial, gridSettings.gridThickness, 0);
        }

        // Generate X gridlines
        for ( int y = 0; y <= gridSettings.gridHeight; y++ ) {
            Vector3[] p = new Vector3[gridSettings.gridPrecision + 1];

            // Compute each position
            for ( int ix = 0; ix <= gridSettings.gridPrecision; ix++ ) {
                float x = (float)ix / gridSettings.gridPrecision * gridSettings.gridWidth;

                p[ix] = Displacement(gridSettings.gridSpacing * new Vector2(x, y)) -
                    gridSettings.gridSpacing * new Vector2(gridSettings.gridWidth, gridSettings.gridHeight) / 2; ;
            }
            AddLine("Gridline X", p, gridLineMaterial, gridSettings.gridThickness, 0);
        }
    }

    // Initialize lines
    void Start () {
        // Initialize at final position
        GenerateGridLines();

        vLine = AddLine("Velocity Vector", new Vector3(-10, -5), new Vector3(10, 5), vMaterial, 0.2f, 0.01f);

        // Initialize at zero
        e0Line = AddLine("e0 Vector", Vector3.zero, Vector3.zero, e0Material, 0.1f, 0.02f);
        e1Line = AddLine("e1 Vector", Vector3.zero, Vector3.zero, e1Material, 0.1f, 0.02f);
    }

    // Update the lines that move
    void Update () {
        // Update sphere
        Vector2 spherePos = Vector3.Lerp(vLine.GetPosition(0), vLine.GetPosition(1), time);


        // Vector2 inverseSpherePos = InverseGridDisplacement(spherePos);

        sphere.transform.position = new Vector3(spherePos.x, spherePos.y, sphere.transform.position.z);

        // Update lines
        Vector2 v = Vector3.Normalize(vLine.GetPosition(1) - vLine.GetPosition(0));

        Vector2 e0 = Derivative(Displacement, spherePos, Vector2.right);
        Vector2 e1 = Derivative(Displacement, spherePos, Vector2.up);

        float v0 = Vector2.Dot(v, e0);
        float v1 = Vector2.Dot(v, e1);

        e0Line.SetPosition(0, spherePos);
        e0Line.SetPosition(1, spherePos + e0 * v0);

        e1Line.SetPosition(0, spherePos);
        e1Line.SetPosition(1, spherePos + e1 * v1);


        backgroundMaterial.color = colorSettings.backgroundColor;
        gridLineMaterial.color = colorSettings.gridColor;
        vMaterial.color = colorSettings.v_color;
    }

}
