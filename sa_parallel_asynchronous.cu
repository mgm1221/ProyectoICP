#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <tuple>
#include <math.h>
#include <string.h>
#include <curand_kernel.h>

using namespace std;
using namespace std::chrono;


__device__ __constant__ char d_sigma[4];
__device__ int *d_T_table;

int sigma_size = 4;

// Calcula valor posible de la función T y los guarda en T_table_flat
void compute_T_table_flat(int* T_table_flat, int HEIGHT, int WIDTH, int MAX_L, int MAX_K) {

    for (int i = 0; i < HEIGHT * WIDTH; ++i)
        T_table_flat[i] = 0;

    T_table_flat[0 * WIDTH + MAX_K] = 1;

    for (int L = 1; L <= MAX_L; ++L) {
        for (int k = -MAX_K; k <= MAX_K; ++k) {
            int offset_k = k + MAX_K;
            int val = 0;

            if (offset_k - 1 >= 0)
                val += T_table_flat[(L - 1) * WIDTH + offset_k - 1];

            val += (sigma_size - 2) * T_table_flat[(L - 1) * WIDTH + offset_k];

            if (offset_k + 1 < WIDTH)
                val += T_table_flat[(L - 1) * WIDTH + offset_k + 1];

            T_table_flat[L * WIDTH + offset_k] = val;
        }
    }
}


// Función para recuperar valores de T_table_flat calculados previmente
__device__ int T_lookup(const int* T_table, int L, int k, int MAX_K, int WIDTH, int HEIGHT) {
    int offset_k = k + MAX_K;
    if (L < 0 || L >= HEIGHT || offset_k < 0 || offset_k >= WIDTH)
        return 0;
    return T_table[L * WIDTH + offset_k];
}


__device__ int rand_int(int max) {
     int tid = threadIdx.x + blockIdx.x * blockDim.x;

        curandState state;
        unsigned long seed = clock64() + tid;
        curand_init(seed, tid, 0, &state);

        int result = curand(&state) % max;

        return result;
}


__device__ void random_sequence(char* seq, int m) {

    for (int i = 0; i < m; i++) {
        seq[i] = d_sigma[rand_int(4)];
    }
}


__device__ void random_neighbor(char* neighbor, char* current, int m) {

    int pos = rand_int(4);

    for (int i = 0; i < m; ++i) {
        neighbor[i] = current[i];
    }

    char old_char = neighbor[pos];
    char new_char;
    do {
        new_char = d_sigma[rand_int(4)];
    } while (new_char == old_char);

    neighbor[pos] = new_char;

}


__device__ int hamming_distance(const char* a, const char* b, int m) {

    int dist = 0;
    for (int i = 0; i < m; i++) {
        if (a[i] != b[i]){
            dist++;
        }
    }
    return dist;
}


// Función de evaluación
__device__ int h(char* s, char* S, int MAX_K, int WIDTH, int HEIGHT, int n,
                 int m, float th, int* d_array, int* c_array, int tid){


    int near = 0;
    char* seq;
    int d = m * th;


    memset(&d_array[tid * n], 0, sizeof(int) * n);
    memset(&c_array[tid * n], 0, sizeof(int) * n);

    for (int i = 0; i < n; i++) {
        seq = &S[i * m];
        d_array[(tid * n) + i] = hamming_distance(s, seq, m);
        c_array[(tid * n) + i] = m - d_array[(tid * n) + i];

        if (d_array[(tid * n) + i] < d) {
            near++;
        }
    }

    int f = n - near;
    int GpC;
    int gi = 0;
    double sumGpC;

    double sumP;

    if (near == 0){
        GpC = 0;
    }
    else{
        sumGpC = 0;

        for (int i = 0; i < n; i++){
            if (d_array[(tid * n) + i] < d){

                gi = 1;

                for (int j = 0; j < n; j++){
                    if (i != j){
                        sumP = 0;

                        for (int c = c_array[(tid * n) + j]; c < c_array[(tid * n) + i]; c++){
                            sumP += T_lookup(d_T_table, c_array[(tid * n) + i], c, MAX_K, WIDTH, HEIGHT) / pow(4, c_array[(tid * n) + i]);
                        }
                        gi += sumP;
                    }
                }
                sumGpC += static_cast<double>(gi) / c_array[(tid * n) + i];
            }
        }
        GpC = sumGpC / near;
    }

    return (n + 1) * f + GpC;
}


__global__ void simulated_annealing(char* flat_sequences, int MAX_K, int WIDTH, int HEIGHT, int n, int m, float th, char* results, float* qualities,
                                    int* d_array,int* c_array, char* neighbors, char* best_s, 
                                    char* current_s, int num_threads, float T_init = 1000, float alpha = 0.95, int max_iter = 1000) {

    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid >= num_threads) return;

    
    curandState state;
    curand_init(1234, tid, 0, &state);

    random_sequence(&current_s[tid * m], m);

    float current_cost = h(&current_s[tid * m], flat_sequences, MAX_K, WIDTH, HEIGHT , n, m, th, d_array, c_array, tid);


    for (int i = 0; i < m; ++i) {
        best_s[(tid * m) + i] = current_s[(tid * m) + i];
    }

    float best_cost = current_cost;
    float T = T_init;

    for (int iter = 0; iter < max_iter && T > 1e-3; iter++) {

        random_neighbor(&neighbors[tid * m] , current_s, m);

        float neighbor_cost = h(&neighbors[tid * m], flat_sequences, MAX_K, WIDTH, HEIGHT, n, m, th, d_array, c_array, tid);
        float delta = neighbor_cost - current_cost;

        float r = curand_uniform(&state);

        if (delta > 0 || (r < expf(-delta / T))) {

            for (int i = 0; i < m; ++i) {
                current_s[(tid * m) + i] = neighbors[(tid * m) + i];
            }

            current_cost = neighbor_cost;
            if (best_cost < current_cost) {

                //best_s = corrent_s
                for (int i = 0; i < m; ++i) {
                    best_s[(tid * m) + i] = current_s[(tid * m) + i];
                }

                best_cost = current_cost;
            }
        }
        
        T *= alpha;
    }

    char* seq;
    int count = 0;
    for (int i = 0; i < n; i++) {
        seq = &flat_sequences[i * m];

        if (hamming_distance(&best_s[tid * m], seq, m) >= th * m) count++;
    }

    float quality = (count / (float)n) * 100;

    for (int i = 0; i < m; ++i) {
        results[(tid * m) + i] = best_s[(tid * m) + i];

    }

    qualities[tid] = quality;

}


int main(int argc, char *argv[])
{

    char sigma[] = {'A', 'C', 'G', 'T'};

    srand(time(NULL));

    if(argc < 5){
        cout << "faltan argumentos" << endl;
        return -1;
    }

    string i_arg = argv[1];
    string str_file = argv[2];

    string th_arg = argv[3];
    float th = stof(argv[4]);

    string T_arg = argv[5];
    float T = stof(argv[6]);

    string alpha_arg = argv[7];
    float alpha = stof(argv[8]);

    string num_threads_arg = argv[9];
    int num_threads = stoi(argv[10]);

    string threads_per_block_arg = argv[11];
    int threads_per_block = stoi(argv[12]);

    auto start = high_resolution_clock::now();

    ifstream myfile(str_file);

    size_t lastSlash = str_file.find_last_of('/');
    string filename = (lastSlash == string::npos) ? str_file : str_file.substr(lastSlash + 1);

    size_t dot = filename.find('.');
    if (dot != string::npos) {
        filename = filename.substr(0, dot);
    }

    int dash1 = filename.find('-');
    int dash2 = filename.find('-', dash1 + 1);

    int n = stoi(filename.substr(0, dash1));
    int m = stoi(filename.substr(dash1 + 1, dash2 - dash1 - 1));
    int MAX_L = n;
    int MAX_K = n;
    int HEIGHT = MAX_L + 1;
    int WIDTH = 2 * MAX_K + 1;

    int* T_table_flat = new int[WIDTH * HEIGHT];
    
    char* host_flat_sequences = new char[n * m];
    string line;
    int counter = 0 ;

    while (getline(myfile, line) && counter < n) {
        if (line.size() != (m+1)) {
            cerr << "Sequence length mismatch on line " << counter << endl;
            return -1;
        }
        memcpy(&host_flat_sequences[counter * m], line.c_str(), m);
        counter++;
    }
    myfile.close();

    compute_T_table_flat(T_table_flat, HEIGHT, WIDTH, MAX_L, MAX_K);

    /*Paralelizacion*/
    int blocks_per_grid = (num_threads + threads_per_block - 1) / threads_per_block;

    // Vectores de resultados y calidad
    char* device_flat_sequences;
    char* d_results;
    float* d_quality;
    int* device_d_array ;
    int* device_c_array;
    char* d_neighbors_flat;
    char* d_best_s_flat;
    char* d_current_s_flat;

    cudaMalloc(&d_T_table, sizeof(int) * HEIGHT * WIDTH);

    //reservar memoria para los resultados
    cudaMalloc(&device_flat_sequences, n * m * sizeof(char));
    cudaMalloc(&d_results, num_threads * sizeof(char) * m);
    cudaMalloc(&d_quality, num_threads * sizeof(float));
    cudaMalloc(&device_d_array, num_threads * n * sizeof(int));
    cudaMalloc(&device_c_array, num_threads * n * sizeof(int));
    cudaMalloc(&d_neighbors_flat, num_threads * sizeof(char) * m);
    cudaMalloc(&d_best_s_flat, num_threads* m * sizeof(char));
    cudaMalloc(&d_current_s_flat, num_threads* m * sizeof(char));

    cudaMemcpyToSymbol(d_sigma, sigma, sizeof(char) * 4);
    
    int* temp_device_T_table;
    cudaMalloc(&temp_device_T_table, sizeof(int) * HEIGHT * WIDTH);
    cudaMemcpy(device_flat_sequences, host_flat_sequences, n * m * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(temp_device_T_table, T_table_flat, sizeof(int) * HEIGHT * WIDTH, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(d_T_table, &temp_device_T_table, sizeof(int*));

    simulated_annealing<<<blocks_per_grid, threads_per_block>>>(device_flat_sequences,
                                                                MAX_K, WIDTH, HEIGHT, 
                                                                n, m, th, d_results, d_quality, 
                                                                device_d_array, device_c_array, d_neighbors_flat, 
                                                                d_best_s_flat, d_current_s_flat, num_threads,
                                                                T, alpha, (int)T);

    char* h_results = (char*)malloc(num_threads * m * sizeof(char));
    float* h_quality = (float*)malloc(num_threads * sizeof(float));
    
    cudaDeviceSynchronize();
    
    // Copia de resultados
    cudaMemcpy(h_results, d_results, num_threads * m * sizeof(char), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_quality, d_quality, num_threads * sizeof(float), cudaMemcpyDeviceToHost);

    int max_quality = 0;
    int best_index = 0;
    for (int i = 0; i < num_threads; i++){  

        if (max_quality < h_quality[i]){
            max_quality = h_quality[i];
            best_index = i;

        }
    }

    auto end = high_resolution_clock::now();

    cout << "Result: ";
    for (int j = 0; j < m;j++){
        printf("%c", h_results[(best_index * m) + j]);
    }

    std::chrono::duration<double> duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);

    cout << " Cardinality: " << h_quality[best_index] << "%" <<
    " Time taken: " << duration.count() << "s" << endl;
    
    cudaFree(d_results);
    cudaFree(d_quality);
    cudaFree(device_flat_sequences);
    cudaFree(temp_device_T_table);
    delete[] T_table_flat;
    delete[] host_flat_sequences;
    free(h_results);
    free(h_quality);

    return 0;
}
