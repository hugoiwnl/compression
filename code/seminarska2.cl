float edist(float prvi, float drugi, float treci){
    return prvi*prvi + drugi*drugi + treci*treci;
}

int random_number(uint x, int y, int gid, int width, int height){
    uint seed = x - gid;
    uint t = seed ^ (seed << 11);
    uint result = y ^ (y >> 19) ^ (t ^ (t >> 8));

    return result % (height*width);
}


__kernel void kompresijaKernel(__global  unsigned char *imageIn,
                            __global  int *indeksi_centroida,
                            __global  int *centroidi_b,
                            __global  int *centroidi_g,
                            __global  int *centroidi_r,
                            __global  int *sums_r,
                            __global  int *sums_g,
                            __global  int *sums_b,
                            __global  int *counters,
                            int width,
                            int height,
                            int st_gruc
                            )
{
    __local int blue_local[16];
    __local int green_local[16];
    __local int red_local[16];
    int gid = get_global_id(0);
    int lid = get_local_id(0);
    //lokalni centroidi
    if(lid < st_gruc){
        //prekopiraj elemente
        blue_local[lid] = centroidi_b[lid];
        green_local[lid] = centroidi_g[lid];
        red_local[lid] = centroidi_r[lid];

        barrier(CLK_LOCAL_MEM_FENCE);
    }

    

    if(gid < width * height){

        

        //uzmi values od tog piksela u slici
        int blue = imageIn[gid*4];
        int green = imageIn[gid*4 + 1];
        int red = imageIn[gid*4 + 2];
        //int index_gruce = najbliziCentroid(centroidi_r, centroidi_g, centroidi_b, red, green, blue, sums_r, sums_g, sums_b,counters);
        
        //nadji najblizi centroid
        int centroid_index = 0;
        //float minDistance = edist((centroidi_r[0] - red), (centroidi_g[0] - green), (centroidi_b[0] - blue));
        float minDistance = edist((red_local[0] - red), (green_local[0] - green), (blue_local[0] - blue));
        for (int i = 1; i < st_gruc; i++)
        {
            //float tmp_distanca = edist((centroidi_r[i] - red), (centroidi_g[i] - green), (centroidi_b[i] - blue));
            float tmp_distanca = edist((red_local[i] - red), (green_local[i] - green), (blue_local[i] - blue));
            if(tmp_distanca < minDistance){
                minDistance = tmp_distanca;
                centroid_index = i;

            }
        }
        int index_gruce = centroid_index;
        atomic_add(&sums_r[index_gruce], red);
        atomic_add(&sums_g[index_gruce], green);
        atomic_add(&sums_b[index_gruce], blue);
        atomic_add(&counters[index_gruce], 1);
        indeksi_centroida[gid] = index_gruce;
    }
    
}

__kernel void refresh(__global  unsigned char *imageIn,
                    __global  int *indeksi_centroida,
                    __global  int *centroidi_b,
                    __global  int *centroidi_g,
                    __global  int *centroidi_r,
                    __global  int *sums_r,
                    __global  int *sums_g,
                    __global  int *sums_b,
                    __global  int *counters,
                    int width,
                    int height,
                    int seed,
                    int st_gruc
                    )
{
    int gid = get_global_id(0);
    
    if(gid < st_gruc){
        if(counters[gid] != 0){
            centroidi_r[gid] = sums_r[gid] / counters[gid];
            centroidi_g[gid] = sums_g[gid] / counters[gid];
            centroidi_b[gid] = sums_b[gid] / counters[gid];
        }else{
            int randomPixel = random_number(seed, 2, gid, width, height);
            indeksi_centroida[randomPixel] = gid;
            centroidi_b[gid] = imageIn[gid * 4];
            centroidi_g[gid] = imageIn[gid * 4 + 1];
            centroidi_r[gid] = imageIn[gid * 4 + 2];
        }
    }
    
    
}