using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class BloomMaster {

    public static float threshold;
    public static float intensity;

    static readonly RenderTexture[] textures = new RenderTexture[3];

    static Shader bloomShader;

    static Material bloomMaterial;

    public static void Bloom ( RenderTexture source, RenderTexture destination ) {

        if ( bloomShader == null ) {
            bloomShader = Shader.Find("Hidden/BloomShader");
            bloomMaterial = new Material(bloomShader);
        }

        bloomMaterial.SetFloat("_Intensity", intensity);
        bloomMaterial.SetFloat("_Threshold", threshold);

        textures[0] = RenderTexture.GetTemporary(source.width, source.height, 0, source.format);

        Graphics.Blit(source, textures[0], bloomMaterial, 0);

        int i;

        for ( i = 1; i < textures.Length; i++ ) {
            textures[i] = RenderTexture.GetTemporary(textures[i - 1].width / 2, textures[i - 1].height / 2, 0, source.format);
            Graphics.Blit(textures[i - 1], textures[i], bloomMaterial, 1);            
        }

        for ( i = textures.Length - 1 ; i >= 1; i -- ) {
            Graphics.Blit(textures[i], textures[i - 1], bloomMaterial, 2);
            RenderTexture.ReleaseTemporary(textures[i]);
            textures[i] = null;
        }

        bloomMaterial.SetTexture("_SourceTex", source);
        Graphics.Blit(textures[0], destination, bloomMaterial, 3);

        RenderTexture.ReleaseTemporary(textures[0]);
    }
}
