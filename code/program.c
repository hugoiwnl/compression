#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "FreeImage.h"
#include <math.h>
#include <time.h>
#include <omp.h>
#include <CL/cl.h>


#define MAX_SOURCE_SIZE	16384
#define ITERACIJE 50
#define ST_GRUC 64
#define NTHREADS 16

void prvi_centroidi(int *centroidi_r, int *centroidi_g, int *centroidi_b, unsigned char *imageIn,int height, int width,int prvi, int drugi, int treci, int cetvrti, int deo){
    
    //kopiraj(imageIn, image_out,width,height);
    //shrani centroide, razdelim sliko na 16 vrstic po 4 redov (64 gruci)
    int popunjeni = 0;
    int delovi = 0;
    for (int i = 0; i < height; i++)
    {
        if(delovi<ST_GRUC/4 & (i==0 || i%deo==0 || i==height-1)){
            //shrani celicu
            
            //prvi piksel
            centroidi_b[popunjeni] = imageIn[(i * width + prvi) * 4];
            centroidi_g[popunjeni] = imageIn[(i * width + prvi) * 4 + 1];
            centroidi_r[popunjeni] = imageIn[(i * width + prvi) * 4 + 2];
            popunjeni++;
            //drugi
            centroidi_b[popunjeni] = imageIn[(i * width + drugi) * 4];
            centroidi_g[popunjeni] = imageIn[(i * width + drugi) * 4 + 1];
            centroidi_r[popunjeni] = imageIn[(i * width + drugi) * 4 + 2];
            popunjeni++;
            //treci
            centroidi_b[popunjeni] = imageIn[(i * width + treci) * 4];
            centroidi_g[popunjeni] = imageIn[(i * width + treci) * 4 + 1];
            centroidi_r[popunjeni] = imageIn[(i * width + treci) * 4 + 2];
            popunjeni++;
            //cetvrti
            centroidi_b[popunjeni] = imageIn[(i * width + cetvrti) * 4];
            centroidi_g[popunjeni] = imageIn[(i * width + cetvrti) * 4 + 1];
            centroidi_r[popunjeni] = imageIn[(i * width + cetvrti) * 4 + 2];
            popunjeni++;

            //dodam broj delova horizontalno
            delovi++;

        }
         
    }
}



float distance(float x1, float x2, float x3)
{
    //return sqrt(pow(x1, 2.0) + pow(x2, 2.0) + pow(x3, 2.0));
    return pow(x1, 2.0) + pow(x2, 2.0) + pow(x3, 2.0);
}



void izprazni(int* arr){
    for (int i = 0; i < ST_GRUC; i++)
    {
        arr[i] = 0;
    }
    
}

int najbliziCentroid(int *centroid_r,int *centroid_g,int *centroid_b, int red, int green, int blue, int* sum_r, int* sum_g, int* sum_b, int* counters){
    int centroid_index = 0;
    float minDistance = distance((centroid_r[0] - red), (centroid_g[0] - green), (centroid_b[0] - blue));
    for (int i = 1; i < ST_GRUC; i++)
    {
        float tmp_distanca = distance((centroid_r[i] - red), (centroid_g[i] - green), (centroid_b[i] - blue));
        if(tmp_distanca < minDistance){
            minDistance = tmp_distanca;
            centroid_index = i;

        }
    }
    sum_r[centroid_index] += red;
    sum_g[centroid_index] += green;
    sum_b[centroid_index] += blue;
    counters[centroid_index] += 1;
    return centroid_index;

}

void kompresija(int height, int width, unsigned char *imageIn, int *indeksi_centroida, int *centroidi_r, int *centroidi_g, int *centroidi_b){
    //razvrstaj piksele v gruce
    
    int *sums_r = (int*)malloc(ST_GRUC * sizeof(int)); //vsote r vrednosti
    int *sums_g = (int*)malloc(ST_GRUC * sizeof(int)); //vsote g vrednosti
    int *sums_b = (int*)malloc(ST_GRUC * sizeof(int)); //vsote b vrednosti
    int *counters = (int*)malloc(ST_GRUC * sizeof(int)); //counteri piksel v gruci
    double start=omp_get_wtime();
    for (int iteracija = 0; iteracija < ITERACIJE; iteracija++)
    {
        izprazni(sums_r); //empty the array
        izprazni(sums_g); //empty the array
        izprazni(sums_b); //empty the array
        izprazni(counters);
        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < width; j++)
            {
                //za svaki piksel izracunaj najblizu grucu
                int blue = imageIn[(i * width + j) * 4];
                int green = imageIn[(i * width + j) * 4 + 1];
                int red = imageIn[(i * width + j) * 4 + 2];
                int index_gruce = najbliziCentroid(centroidi_r, centroidi_g, centroidi_b, red, green, blue, sums_r, sums_g, sums_b, counters);
                indeksi_centroida[i*width+j] = index_gruce;
            }        
        }

        //izracunaj avg vrednosti za svaku gruc
        if(iteracija < ITERACIJE-1){
            for (int i = 0; i < ST_GRUC; i++)
            {
                if(counters[i] != 0){
                    //stavi nove average vrednosti gruc
                    centroidi_r[i] = sums_r[i] / counters[i];
                    centroidi_g[i] = sums_g[i] / counters[i];
                    centroidi_b[i] = sums_b[i] / counters[i];
                }else{
                    int randomPixel = rand()%((height*width+1)-0) + 1;
                    indeksi_centroida[randomPixel] = i;
                    centroidi_r[i] = imageIn[randomPixel * 4];
                    centroidi_g[i] = imageIn[randomPixel * 4 + 1];
                    centroidi_b[i] = imageIn[randomPixel * 4 + 2];
                }
            }
        }
        
        
    }
    double stop=omp_get_wtime();
    double ptime=stop-start;
    printf("Pretecen sekvencni cas: %f s\n",ptime); 
}



int najbliziCentroidCPU(int *centroid_r,int *centroid_g,int *centroid_b, int red, int green, int blue, int* sum_r, int* sum_g, int* sum_b, int* counters){
    int centroid_index = 0;
    float minDistance = distance((centroid_r[0] - red), (centroid_g[0] - green), (centroid_b[0] - blue));
    for (int i = 1; i < ST_GRUC; i++)
    {
        float tmp_distanca = distance((centroid_r[i] - red), (centroid_g[i] - green), (centroid_b[i] - blue));
        if(tmp_distanca < minDistance){
            minDistance = tmp_distanca;
            centroid_index = i;

        }
    }
    sum_r[centroid_index] += red;
    sum_g[centroid_index] += green;
    sum_b[centroid_index] += blue;
    #pragma omp atomic update
        counters[centroid_index] += 1;
    return centroid_index;

}

void kompresijaCPU(int height, int width, unsigned char *imageIn, int *indeksi_centroida, int *centroidi_r, int *centroidi_g, int *centroidi_b){
    //razvrstaj piksele v gruce
    int *sums_r = (int*)malloc(ST_GRUC * sizeof(int)); //vsote r vrednosti
    int *sums_g = (int*)malloc(ST_GRUC * sizeof(int)); //vsote g vrednosti
    int *sums_b = (int*)malloc(ST_GRUC * sizeof(int)); //vsote b vrednosti
    int *counters = (int*)malloc(ST_GRUC * sizeof(int)); //counteri piksel v gruci
    int index_gruce;
    double start=omp_get_wtime();

    for (int iteracija = 0; iteracija < ITERACIJE; iteracija++)
    {
        izprazni(sums_r); //empty the array
        izprazni(sums_g); //empty the array
        izprazni(sums_b); //empty the array
        izprazni(counters);

        #pragma omp parallel for private(index_gruce) shared(sums_r, sums_g, sums_b, counters)
        for (int i = 0; i < width * height; i++)
        {
            int blue = imageIn[i*4];
            int green = imageIn[i*4 + 1];
            int red = imageIn[i*4 + 2];
            index_gruce = najbliziCentroidCPU(centroidi_r, centroidi_g, centroidi_b, red, green, blue, sums_r, sums_g, sums_b, counters);
            indeksi_centroida[i] = index_gruce;
        }

        //sum
        

        //izracunaj avg vrednosti za svaku gruc
    
        if(iteracija < ITERACIJE-1){
            for (int i = 0; i < ST_GRUC; i++)
            {
                if(counters[i] != 0){
                    //stavi nove average vrednosti gruc
                    centroidi_r[i] = sums_r[i] / counters[i];
                    centroidi_g[i] = sums_g[i] / counters[i];
                    centroidi_b[i] = sums_b[i] / counters[i];
                }else{
                    int randomPixel = rand()%((height*width+1)-0) + 1;
                    indeksi_centroida[randomPixel] = i;
                    centroidi_r[i] = imageIn[randomPixel * 4];
                    centroidi_g[i] = imageIn[randomPixel * 4 + 1];
                    centroidi_b[i] = imageIn[randomPixel * 4 + 2];
                }
            }
        }
    
        
    }
    double stop=omp_get_wtime();
    double ptime=stop-start;
    printf("Pretecen CPU2 cas: %f s\n",ptime); 
}



void kompresijaGPU(int height, int width, int pitch, FIBITMAP *imageBitmap32, int *indeksi_centroida, int *centroidi_r, int *centroidi_g, int *centroidi_b){
    char ch;
    int i;
	cl_int ret;
	// Branje datoteke
    FILE *fp;
    char *source_str;
    size_t source_size;

    fp = fopen("seminarska.cl", "r");
	
    if (!fp) 
	{
		fprintf(stderr, ":-(#\n");
        exit(1);
    }
    source_str = (char*)malloc(MAX_SOURCE_SIZE);
    source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
	source_str[source_size] = '\0';
    fclose( fp );

    //alociram memoriju za input picture array
    unsigned char *image_in_gpu = (unsigned char *)malloc(height * pitch * sizeof(unsigned char));
    FreeImage_ConvertToRawBits(image_in_gpu, imageBitmap32, pitch, 32, 0xFF, 0xFF, 0xFF, TRUE);

    //alociram memoriju za output array
    //unsigned char *image_out_gpu = (unsigned char *)malloc(height * pitch * sizeof(unsigned char));

    int randomPixel = rand()%((height*width+1)-0) + 1; //random pixel koji passujem u drugi kernel ako je stevilo gruc 0
    int gruci = ST_GRUC;

    int *sums_r = (int*)malloc(ST_GRUC * sizeof(int)); //vsote r vrednosti
    int *sums_g = (int*)malloc(ST_GRUC * sizeof(int)); //vsote g vrednosti
    int *sums_b = (int*)malloc(ST_GRUC * sizeof(int)); //vsote b vrednosti
    int *counters = (int*)malloc(ST_GRUC * sizeof(int)); //counteri piksel v gruci

    izprazni(sums_r);
    izprazni(sums_g);
    izprazni(sums_b);
    izprazni(counters);

    //platforma
    cl_platform_id	platform_id[10];
    cl_uint			ret_num_platforms;
	char			*buf;
	size_t			buf_len;
	ret = clGetPlatformIDs(10, platform_id, &ret_num_platforms);
			// max. "stevilo platform, kazalec na platforme, dejansko "stevilo platform

	// Podatki o napravi
	cl_device_id	device_id[10];
	cl_uint			ret_num_devices;
	// Delali bomo s platform_id[0] na GPU
	ret = clGetDeviceIDs(platform_id[0], CL_DEVICE_TYPE_GPU, 10,	
						 device_id, &ret_num_devices);				
			// izbrana platforma, tip naprave, koliko naprav nas zanima
			// kazalec na naprave, dejansko "stevilo naprav
    // Kontekst
	cl_context context = clCreateContext(NULL,
     1,
      &device_id[0],
       NULL,
        NULL,
         &ret);
			// kontekst: vklju"cene platforme - NULL je privzeta, "stevilo naprav, 
			// kazalci na naprave, kazalec na call-back funkcijo v primeru napake
			// dodatni parametri funkcije, "stevilka napake
    // Ukazna vrsta
    cl_command_queue command_queue = clCreateCommandQueue(context, device_id[0], 0, &ret);

    size_t local_item_size=256;
    size_t num_groups = (((width * height) - 1) / local_item_size + 1);
    size_t global_item_size = num_groups * local_item_size;

    int gruc_st = ST_GRUC;

    //mem obj za input image
    cl_mem img_in_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
									  height * pitch * sizeof(unsigned char), image_in_gpu, &ret); //mem obj za input image
    //mem obj za indekse gruci od piksela
    cl_mem indeksi_mem_obj = clCreateBuffer(context, CL_MEM_WRITE_ONLY, 
                                      height * width * sizeof(int), NULL, &ret); 
    //mem objekti za centroide gruc
    cl_mem b_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
                                      ST_GRUC * sizeof(int), centroidi_b, &ret); //mem obj za b gruce
    cl_mem g_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
                                      ST_GRUC * sizeof(int), centroidi_g, &ret); //mem obj za g gruce
    cl_mem r_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
                                      ST_GRUC * sizeof(int), centroidi_r, &ret); //mem obj za r gruce
    cl_mem sum_r_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
                                      ST_GRUC * sizeof(int), sums_r, &ret); //mem obj za r sums
    cl_mem sum_g_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
                                      ST_GRUC * sizeof(int), sums_g, &ret); //mem obj za g sums
    cl_mem sum_b_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
                                      ST_GRUC * sizeof(int), sums_b, &ret); //mem obj za b sums
    cl_mem counters_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
                                      ST_GRUC * sizeof(int), counters, &ret); //mem obj za broj piksela u gruci
    
    
    cl_program program = clCreateProgramWithSource(context,	1, (const char **)&source_str,  
												   NULL, &ret);
			// kontekst, "stevilo kazalcev na kodo, kazalci na kodo,		
			// stringi so NULL terminated, napaka													
 
    // Prevajanje
    ret = clBuildProgram(program, 1, &device_id[0], NULL, NULL, NULL);
			// program, "stevilo naprav, lista naprav, opcije pri prevajanju,
			// kazalec na funkcijo, uporabni"ski argumenti
    // Log
	size_t build_log_len;
	char *build_log;
	ret = clGetProgramBuildInfo(program, device_id[0], CL_PROGRAM_BUILD_LOG, 
								0, NULL, &build_log_len);
			// program, "naprava, tip izpisa, 
			// maksimalna dol"zina niza, kazalec na niz, dejanska dol"zina niza
	build_log = (char *)malloc(sizeof(char)*(build_log_len+1));
	ret = clGetProgramBuildInfo(program, device_id[0], CL_PROGRAM_BUILD_LOG, 
							    build_log_len, build_log, NULL);
    
	free(build_log);
    double start=omp_get_wtime();
    //ustvarim kernel
    cl_kernel kernel = clCreateKernel(program, "kompresijaKernel", &ret); //prvi kernel za razvrscanje piksela u gruce

    cl_kernel kernel2 = clCreateKernel(program, "refresh", &ret); //drugi kernel za racunanje average v gruc


    ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&img_in_mem_obj);
    ret |= clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&indeksi_mem_obj);
    ret |= clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&b_mem_obj);
    ret |= clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&g_mem_obj);
    ret |= clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&r_mem_obj);
    ret |= clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&sum_r_mem_obj);
    ret |= clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&sum_g_mem_obj);
    ret |= clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *)&sum_b_mem_obj);
    ret |= clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *)&counters_mem_obj);
    ret |= clSetKernelArg(kernel, 9, sizeof(cl_int), (void *)&width);
    ret |= clSetKernelArg(kernel, 10, sizeof(cl_int), (void *)&height);
    ret |= clSetKernelArg(kernel, 11, sizeof(cl_int), (void *)&gruci);


    ret = clSetKernelArg(kernel2, 0, sizeof(cl_mem), (void *)&img_in_mem_obj);
    ret |= clSetKernelArg(kernel2, 1, sizeof(cl_mem), (void *)&indeksi_mem_obj);
    ret |= clSetKernelArg(kernel2, 2, sizeof(cl_mem), (void *)&b_mem_obj);
    ret |= clSetKernelArg(kernel2, 3, sizeof(cl_mem), (void *)&g_mem_obj);
    ret |= clSetKernelArg(kernel2, 4, sizeof(cl_mem), (void *)&r_mem_obj);
    ret |= clSetKernelArg(kernel2, 5, sizeof(cl_mem), (void *)&sum_r_mem_obj);
    ret |= clSetKernelArg(kernel2, 6, sizeof(cl_mem), (void *)&sum_g_mem_obj);
    ret |= clSetKernelArg(kernel2, 7, sizeof(cl_mem), (void *)&sum_b_mem_obj);
    ret |= clSetKernelArg(kernel2, 8, sizeof(cl_mem), (void *)&counters_mem_obj);
    ret |= clSetKernelArg(kernel2, 9, sizeof(cl_int), (void *)&width);
    ret |= clSetKernelArg(kernel2, 10, sizeof(cl_int), (void *)&height);
    ret |= clSetKernelArg(kernel2, 11, sizeof(cl_int), (void *)&randomPixel);
    ret |= clSetKernelArg(kernel2, 12, sizeof(cl_int), (void *)&gruci);



    for (int i = 0; i < ITERACIJE; i++)
    {
        ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL,						
								 &global_item_size, &local_item_size, 0, NULL, NULL);
    

    
    


        ret = clEnqueueNDRangeKernel(command_queue, kernel2, 1, NULL,						
	   							 &global_item_size, &local_item_size, 0, NULL, NULL);
    }
    
	

    ret = clEnqueueReadBuffer(command_queue, indeksi_mem_obj, CL_TRUE, 0,						
							  height * width * sizeof(int), indeksi_centroida, 0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, b_mem_obj, CL_TRUE, 0,						
							 ST_GRUC * sizeof(int), centroidi_b, 0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, g_mem_obj, CL_TRUE, 0,						
							  ST_GRUC * sizeof(int), centroidi_g, 0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, r_mem_obj, CL_TRUE, 0,						
							  ST_GRUC * sizeof(int), centroidi_r, 0, NULL, NULL);

    double stop=omp_get_wtime();
    double ptime=stop-start;
    printf("Pretecen GPU cas: %f s\n",ptime);

    ret = clFlush(command_queue);
    ret = clFinish(command_queue);
    ret = clReleaseKernel(kernel);
    ret = clReleaseKernel(kernel2);
    ret = clReleaseProgram(program);
    ret = clReleaseMemObject(img_in_mem_obj);
    ret = clReleaseMemObject(indeksi_mem_obj);
    ret = clReleaseMemObject(b_mem_obj);
    ret = clReleaseMemObject(g_mem_obj);
    ret = clReleaseMemObject(r_mem_obj);
    ret = clReleaseMemObject(sum_b_mem_obj);
    ret = clReleaseMemObject(sum_g_mem_obj);
    ret = clReleaseMemObject(sum_r_mem_obj);
    ret = clReleaseMemObject(counters_mem_obj);
    ret = clReleaseCommandQueue(command_queue);
    ret = clReleaseContext(context);

    
    
    
}

void kompresijaGPU2(int height, int width, int pitch, FIBITMAP *imageBitmap32, int *indeksi_centroida, int *centroidi_r, int *centroidi_g, int *centroidi_b){
    char ch;
    int i;
	cl_int ret;
	// Branje datoteke
    FILE *fp;
    char *source_str;
    size_t source_size;

    fp = fopen("seminarska2.cl", "r");
	
    if (!fp) 
	{
		fprintf(stderr, ":-(#\n");
        exit(1);
    }
    source_str = (char*)malloc(MAX_SOURCE_SIZE);
    source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
	source_str[source_size] = '\0';
    fclose( fp );

    //alociram memoriju za input picture array
    unsigned char *image_in_gpu = (unsigned char *)malloc(height * pitch * sizeof(unsigned char));
    FreeImage_ConvertToRawBits(image_in_gpu, imageBitmap32, pitch, 32, 0xFF, 0xFF, 0xFF, TRUE);

    //alociram memoriju za output array
    //unsigned char *image_out_gpu = (unsigned char *)malloc(height * pitch * sizeof(unsigned char));

    int randomPixel = rand()%((height*width+1)-0) + 1; //random pixel koji passujem u drugi kernel ako je stevilo gruc 0
    int gruci = ST_GRUC;

    int *sums_r = (int*)malloc(ST_GRUC * sizeof(int)); //vsote r vrednosti
    int *sums_g = (int*)malloc(ST_GRUC * sizeof(int)); //vsote g vrednosti
    int *sums_b = (int*)malloc(ST_GRUC * sizeof(int)); //vsote b vrednosti
    int *counters = (int*)malloc(ST_GRUC * sizeof(int)); //counteri piksel v gruci

    izprazni(sums_r);
    izprazni(sums_g);
    izprazni(sums_b);
    izprazni(counters);

    //platforma
    cl_platform_id	platform_id[10];
    cl_uint			ret_num_platforms;
	char			*buf;
	size_t			buf_len;
	ret = clGetPlatformIDs(10, platform_id, &ret_num_platforms);
			// max. "stevilo platform, kazalec na platforme, dejansko "stevilo platform

	// Podatki o napravi
	cl_device_id	device_id[10];
	cl_uint			ret_num_devices;
	// Delali bomo s platform_id[0] na GPU
	ret = clGetDeviceIDs(platform_id[0], CL_DEVICE_TYPE_GPU, 10,	
						 device_id, &ret_num_devices);				
			// izbrana platforma, tip naprave, koliko naprav nas zanima
			// kazalec na naprave, dejansko "stevilo naprav
    // Kontekst
	cl_context context = clCreateContext(NULL,
     1,
      &device_id[0],
       NULL,
        NULL,
         &ret);
			// kontekst: vklju"cene platforme - NULL je privzeta, "stevilo naprav, 
			// kazalci na naprave, kazalec na call-back funkcijo v primeru napake
			// dodatni parametri funkcije, "stevilka napake
    // Ukazna vrsta
    cl_command_queue command_queue = clCreateCommandQueue(context, device_id[0], 0, &ret);

    size_t local_item_size=256;
    size_t num_groups = (((width * height) - 1) / local_item_size + 1);
    size_t global_item_size = num_groups * local_item_size;

    int gruc_st = ST_GRUC;

    //mem obj za input image
    cl_mem img_in_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
									  height * pitch * sizeof(unsigned char), image_in_gpu, &ret); //mem obj za input image
    //mem obj za indekse gruci od piksela
    cl_mem indeksi_mem_obj = clCreateBuffer(context, CL_MEM_WRITE_ONLY, 
                                      height * width * sizeof(int), NULL, &ret); 
    //mem objekti za centroide gruc
    cl_mem b_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
                                      ST_GRUC * sizeof(int), centroidi_b, &ret); //mem obj za b gruce
    cl_mem g_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
                                      ST_GRUC * sizeof(int), centroidi_g, &ret); //mem obj za g gruce
    cl_mem r_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
                                      ST_GRUC * sizeof(int), centroidi_r, &ret); //mem obj za r gruce
    cl_mem sum_r_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
                                      ST_GRUC * sizeof(int), sums_r, &ret); //mem obj za r sums
    cl_mem sum_g_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
                                      ST_GRUC * sizeof(int), sums_g, &ret); //mem obj za g sums
    cl_mem sum_b_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
                                      ST_GRUC * sizeof(int), sums_b, &ret); //mem obj za b sums
    cl_mem counters_mem_obj = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
                                      ST_GRUC * sizeof(int), counters, &ret); //mem obj za broj piksela u gruci
    
    
    cl_program program = clCreateProgramWithSource(context,	1, (const char **)&source_str,  
												   NULL, &ret);
			// kontekst, "stevilo kazalcev na kodo, kazalci na kodo,		
			// stringi so NULL terminated, napaka													
 
    // Prevajanje
    ret = clBuildProgram(program, 1, &device_id[0], NULL, NULL, NULL);
			// program, "stevilo naprav, lista naprav, opcije pri prevajanju,
			// kazalec na funkcijo, uporabni"ski argumenti
    // Log
	size_t build_log_len;
	char *build_log;
	ret = clGetProgramBuildInfo(program, device_id[0], CL_PROGRAM_BUILD_LOG, 
								0, NULL, &build_log_len);
			// program, "naprava, tip izpisa, 
			// maksimalna dol"zina niza, kazalec na niz, dejanska dol"zina niza
	build_log = (char *)malloc(sizeof(char)*(build_log_len+1));
	ret = clGetProgramBuildInfo(program, device_id[0], CL_PROGRAM_BUILD_LOG, 
							    build_log_len, build_log, NULL);
    
	free(build_log);
    double start=omp_get_wtime();
    //ustvarim kernel
    cl_kernel kernel = clCreateKernel(program, "kompresijaKernel", &ret); //prvi kernel za razvrscanje piksela u gruce

    cl_kernel kernel2 = clCreateKernel(program, "refresh", &ret); //drugi kernel za racunanje average v gruc


    ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&img_in_mem_obj);
    ret |= clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&indeksi_mem_obj);
    ret |= clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&b_mem_obj);
    ret |= clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&g_mem_obj);
    ret |= clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&r_mem_obj);
    ret |= clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&sum_r_mem_obj);
    ret |= clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&sum_g_mem_obj);
    ret |= clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *)&sum_b_mem_obj);
    ret |= clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *)&counters_mem_obj);
    ret |= clSetKernelArg(kernel, 9, sizeof(cl_int), (void *)&width);
    ret |= clSetKernelArg(kernel, 10, sizeof(cl_int), (void *)&height);
    ret |= clSetKernelArg(kernel, 11, sizeof(cl_int), (void *)&gruci);


    ret = clSetKernelArg(kernel2, 0, sizeof(cl_mem), (void *)&img_in_mem_obj);
    ret |= clSetKernelArg(kernel2, 1, sizeof(cl_mem), (void *)&indeksi_mem_obj);
    ret |= clSetKernelArg(kernel2, 2, sizeof(cl_mem), (void *)&b_mem_obj);
    ret |= clSetKernelArg(kernel2, 3, sizeof(cl_mem), (void *)&g_mem_obj);
    ret |= clSetKernelArg(kernel2, 4, sizeof(cl_mem), (void *)&r_mem_obj);
    ret |= clSetKernelArg(kernel2, 5, sizeof(cl_mem), (void *)&sum_r_mem_obj);
    ret |= clSetKernelArg(kernel2, 6, sizeof(cl_mem), (void *)&sum_g_mem_obj);
    ret |= clSetKernelArg(kernel2, 7, sizeof(cl_mem), (void *)&sum_b_mem_obj);
    ret |= clSetKernelArg(kernel2, 8, sizeof(cl_mem), (void *)&counters_mem_obj);
    ret |= clSetKernelArg(kernel2, 9, sizeof(cl_int), (void *)&width);
    ret |= clSetKernelArg(kernel2, 10, sizeof(cl_int), (void *)&height);
    ret |= clSetKernelArg(kernel2, 11, sizeof(cl_int), (void *)&randomPixel);
    ret |= clSetKernelArg(kernel2, 12, sizeof(cl_int), (void *)&gruci);



    for (int i = 0; i < ITERACIJE; i++)
    {
        ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL,						
								 &global_item_size, &local_item_size, 0, NULL, NULL);
    

    
    


        ret = clEnqueueNDRangeKernel(command_queue, kernel2, 1, NULL,						
	   							 &global_item_size, &local_item_size, 0, NULL, NULL);
    }
    
	

    ret = clEnqueueReadBuffer(command_queue, indeksi_mem_obj, CL_TRUE, 0,						
							  height * width * sizeof(int), indeksi_centroida, 0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, b_mem_obj, CL_TRUE, 0,						
							 ST_GRUC * sizeof(int), centroidi_b, 0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, g_mem_obj, CL_TRUE, 0,						
							  ST_GRUC * sizeof(int), centroidi_g, 0, NULL, NULL);
    ret = clEnqueueReadBuffer(command_queue, r_mem_obj, CL_TRUE, 0,						
							  ST_GRUC * sizeof(int), centroidi_r, 0, NULL, NULL);

    double stop=omp_get_wtime();
    double ptime=stop-start;
    printf("Pretecen GPU cas: %f s\n",ptime);

    ret = clFlush(command_queue);
    ret = clFinish(command_queue);
    ret = clReleaseKernel(kernel);
    ret = clReleaseKernel(kernel2);
    ret = clReleaseProgram(program);
    ret = clReleaseMemObject(img_in_mem_obj);
    ret = clReleaseMemObject(indeksi_mem_obj);
    ret = clReleaseMemObject(b_mem_obj);
    ret = clReleaseMemObject(g_mem_obj);
    ret = clReleaseMemObject(r_mem_obj);
    ret = clReleaseMemObject(sum_b_mem_obj);
    ret = clReleaseMemObject(sum_g_mem_obj);
    ret = clReleaseMemObject(sum_r_mem_obj);
    ret = clReleaseMemObject(counters_mem_obj);
    ret = clReleaseCommandQueue(command_queue);
    ret = clReleaseContext(context);

    
    
    
}



int main(int argc, char *argv[]){
    //Load image from file
	FIBITMAP *imageBitmap = FreeImage_Load(FIF_PNG, "3840.png", 0);
	
    FIBITMAP *imageBitmap32 = FreeImage_ConvertTo32Bits(imageBitmap);
    FIBITMAP *imageBitmap8 = FreeImage_ConvertTo8Bits(imageBitmap);
	
    //Get image dimensions
    int width = FreeImage_GetWidth(imageBitmap32);
	int height = FreeImage_GetHeight(imageBitmap32);
	int pitch = FreeImage_GetPitch(imageBitmap32);

    
    unsigned char *imageIn = (unsigned char *)malloc(height*pitch * sizeof(unsigned char));
    FreeImage_ConvertToRawBits(imageIn, imageBitmap32, pitch, 32, FI_RGBA_RED_MASK, FI_RGBA_GREEN_MASK, FI_RGBA_BLUE_MASK, TRUE);
    //unsigned char *image_out = imageIn;
    unsigned char *image_out = (unsigned char *)malloc(height * pitch * sizeof(unsigned char));
    FreeImage_ConvertToRawBits(image_out, imageBitmap32, pitch, 32, FI_RGBA_RED_MASK, FI_RGBA_GREEN_MASK, FI_RGBA_BLUE_MASK, TRUE);

    //Free source image data
	
    int *centroidi_r = (int*)malloc(ST_GRUC * sizeof(int));
    int *centroidi_g = (int*)malloc(ST_GRUC * sizeof(int));
    int *centroidi_b = (int*)malloc(ST_GRUC * sizeof(int));

    double deljen = height/(ST_GRUC/4);
    //int deo = ceil(deljen);
    int deo = (int)deljen + 1;


    //inicijalizacija centroida
    int prvi = 0;
    int drugi = width/4;
    int treci = drugi*3;
    int cetvrti = width-1;
    prvi_centroidi(centroidi_r,centroidi_g,centroidi_b,imageIn,height,width,prvi,drugi,treci,cetvrti,deo);
    int *indeksi_centroida = (int*)malloc(height * width * sizeof(int));

    //ovde krece program
    
    //kompresija(height,width,imageIn,indeksi_centroida,centroidi_r,centroidi_g,centroidi_b); //sekvencni algoritem

    kompresijaCPU(height,width,imageIn,indeksi_centroida,centroidi_r,centroidi_g,centroidi_b);
 
    //kompresijaGPU(height,width,pitch,imageBitmap32,indeksi_centroida,centroidi_r,centroidi_g,centroidi_b);

    //kompresijaGPU2(height,width,pitch,imageBitmap32,indeksi_centroida,centroidi_r,centroidi_g,centroidi_b);


    //apliciraj nove boje na novu sliku  
    
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            //vidim kojoj gruci pripada taj piksel
            int gruca = indeksi_centroida[i*width+j];
            image_out[(i * width + j) * 4] = centroidi_b[gruca];
            image_out[(i * width + j) * 4 + 1] = centroidi_g[gruca];
            image_out[(i * width + j) * 4 + 2] = centroidi_r[gruca];
        }
        
    }

    
    
    


    
    //shranim sliku
    FIBITMAP *dst = FreeImage_ConvertFromRawBits(image_out, width, height, pitch,
		32, 0xFF, 0xFF, 0xFF, TRUE);
	FreeImage_Save(FIF_PNG, dst, "kompresija.png", 0);  
    
    
}

/*					Meritev casa(64 threads 64 gruc)

Rezolucija			Serijski algoritem			OpenMP			OpenCL

640x480				6.128390 s					1.769160 s		0.033635 s
800x600				9.605911 s					2.748614 s		0.075463 s
1600x900			28.698578 s					6.685092 s		0.130239 s
1920x1080			41.330940 s					11.911868 s		0.178933 s
3840x2160			165.541810 s				39.195659 s     0.749682 s
*/

