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
    #include <random>

    using namespace std;
    using namespace std::chrono;

    vector<vector<int>> memo;
    vector<vector<bool>> computed;

    vector<char> sigma = {'A', 'C', 'G', 'T'};
    int sigma_size = 4;

    string random_sequence(int m){

        string str;
        str.resize(m);

        for(int i = 0; i < m; i++){
            str[i] = sigma[rand() % sigma_size];
        }

        return str;
    }

    string random_neighbor(string current_solution) {
        
        int n = current_solution.length();
        
        int pos = rand()% n;

        string neighbor_solution = current_solution;
        
        char old_char = neighbor_solution[pos];
        char new_char;
        do {
            new_char = sigma[rand() % sigma_size];
        } while (new_char == old_char);

        neighbor_solution[pos] = new_char;

        return neighbor_solution;

    }

    int hammingDist(string str1, string str2) 
    { 
        int i = 0, count = 0; 
        while (str1[i] != '\0') { 
            if (str1[i] != str2[i]) 
                count++; 
            i++; 
        } 
        return count; 
    }

    int T(int L, int k, int MAX_K, int WIDTH, int HEIGHT) {
        int offset_k = k + MAX_K;

        if (offset_k < 0 || offset_k > 2 * MAX_K)
            return 0;

        if (L == 0 && k == 0)
            return 1;
        if (L == 0)
            return 0;

        if (computed[L][offset_k])
            return memo[L][offset_k];

        int val = T(L - 1, k - 1, MAX_K, WIDTH, HEIGHT) +
                (sigma_size - 2) * T(L - 1, k, MAX_K, WIDTH, HEIGHT) +
                T(L - 1, k + 1, MAX_K, WIDTH, HEIGHT);

        memo[L][offset_k] = val;
        computed[L][offset_k] = true;
        return val;
    }

    int h(string s,vector<string> S, int MAX_K, int WIDTH, int HEIGHT,int n, float th){

        int near = 0;
        string seq;
        int m = s.length();
        float d = m * th;


        int d_array[n];
        int c_array[n];
        int g_array[n];

        memset(d_array, 0, sizeof(int) * n);
        memset(c_array, 0, sizeof(int) * n);
        memset(g_array, 0, sizeof(int) * n);

        for (int i = 0; i < n; i++) {
            seq = S[i];
            d_array[i] = hammingDist(s, seq);
            c_array[i] = m - d_array[i];

            if (d_array[i] < d) {
                near++;
            }
        }
        
        int f = n - near;
        int GpC;
        double sumGpC;

        double sumP;

        if (near == 0){
            GpC = 0;
        }
        else{
            sumGpC = 0;

            for (int i = 0; i < n; i++){
                if (d_array[i] < d){

                    g_array[i] = 1;

                    for (int j = 0; j < n; j++){
                        if (i != j){
                            sumP = 0;
                        
                            for (int c = c_array[j]; c < c_array[i]; c++){

                                sumP += T(c_array[i], c, MAX_K, WIDTH, HEIGHT) / pow(sigma_size, c_array[i]);
                            }
                            g_array[i] += sumP;
                        }
                    }
                    
                    sumGpC += static_cast<double>(g_array[i]) / c_array[i];
                }
            }
            GpC = sumGpC / near;
        }
        
        return (n + 1) * f + GpC;
    }


    tuple<string, int, float> simulated_annealing(vector<string> S, string initial_s, int MAX_K, int WIDTH, int HEIGHT, int n, float th,mt19937& gen, uniform_real_distribution<float>& dist, float T_init = 5000, float alpha = 0.95, int max_iter = 5000) {
        vector<char> sigma = {'A', 'C', 'G', 'T'};
        string current_s = initial_s;
        float current_cost = h(current_s, S, MAX_K, WIDTH, HEIGHT, n, th);
        string best_s = current_s;
        float best_cost = current_cost;
        float T = T_init;

        for (int iter = 0; iter < max_iter && T > 1e-3; iter++) {
            string neighbor = random_neighbor(current_s);

            float neighbor_cost = h(neighbor, S, MAX_K, WIDTH, HEIGHT, n, th);
            float delta = neighbor_cost - current_cost;
            float r = dist(gen);
            
            if (delta > 0 || r < exp(-delta / T)) {
                
                current_s = neighbor;
                current_cost = neighbor_cost;

                if (best_cost < current_cost) {
                    best_s = current_s;
                    best_cost = current_cost;
                }
                
            }
            T *= alpha;
        }


        string seq;
        int count = 0;
        for (int i = 0; i < n; i++) {
            seq = S[i];
            if (hammingDist(best_s, seq) >= th * best_s.length()) count++;
        }

        float quality = ((float)count / (float)n) * 100;
        return make_tuple(best_s, hammingDist(initial_s, best_s), quality);
    }


    int main(int argc, char *argv[])
    {

        srand(time(NULL));
        random_device rd;                         // Non-deterministic seed
        mt19937 gen(rd());  

        uniform_real_distribution<float> dist(0.0f, 1.0f);
        if(argc < 5){
            cout << "faltan argumentos" << endl;
            return -1;
        }
        string initial_s, s;


        string i_arg = argv[1];
        string str_file = argv[2];

        string th_arg = argv[3];
        float th = stof(argv[4]);

        string T_arg = argv[5];
        float T = stof(argv[6]);

        string alpha_arg = argv[7];
        float alpha = stof(argv[8]);

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

        int MAX_L = m;
        int MAX_K = m;
        int HEIGHT = MAX_L + 1;
        int WIDTH = 2 * MAX_K + 1;

        memo = vector<vector<int>>(HEIGHT, vector<int>(WIDTH, 0));
        computed = vector<vector<bool>>(HEIGHT, vector<bool>(WIDTH, false));

        vector<string> sequences;
        string line;
        while (myfile >> line) {
            sequences.push_back(line);
        }
        
        string init_s = sequences[0];
        
        
        auto [final_s, changes, quality] = simulated_annealing(sequences, init_s, MAX_K, WIDTH, HEIGHT, n, th, gen, dist, T, alpha, (int) T);
        auto end = high_resolution_clock::now();

        std::chrono::duration<double> duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);

        cout << "Result: " << final_s;
    

        cout << " Cardinality: " << quality << "%" <<
        " Time taken: " << duration.count() << "s" << endl;

        myfile.close();
        return 0;
    }
